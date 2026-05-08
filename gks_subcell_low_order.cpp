#include "gks_subcell_low_order.h"

#include "gks_kfvs_adapter.h"

#include <algorithm>
#include <cmath>
#include <limits>

GKSSubcellLowOrderType low_order_type = KFVS1;

GKSSubcellMUSCLHancockStats1D::GKSSubcellMUSCLHancockStats1D()
	: reconstructed_cells(0),
	  zero_slope_fallback_cells(0),
	  predicted_state_fallback_cells(0),
	  face_flux_fallback_faces(0),
	  final_update_fallback_cells(0),
	  min_rho_update(1.0e300),
	  min_p_update(1.0e300)
{
}

namespace
{
	const double kMUSCLDensityFloor = 1.0e-12;
	const double kMUSCLPressureFloor = 1.0e-12;

	void ReflectNormalVelocityState(const double in_Q[3], double out_Q[3])
	{
		out_Q[0] = in_Q[0];
		out_Q[1] = -in_Q[1];
		out_Q[2] = in_Q[2];
	}

	void ZeroState(double out[3])
	{
		for (int m = 0; m < 3; ++m)
		{
			out[m] = 0.0;
		}
	}

	bool StateFinite1D(const double Q[3])
	{
		return std::isfinite(Q[0]) && std::isfinite(Q[1]) && std::isfinite(Q[2]);
	}

	double PressureFromConservative1D(const double Q[3])
	{
		if (Q[0] <= 0.0)
		{
			return -std::numeric_limits<double>::infinity();
		}
		return Pressure(Q[0], Q[1], Q[2]);
	}

	bool ConservativeAdmissible1D(const double Q[3])
	{
		return StateFinite1D(Q) &&
			Q[0] > kMUSCLDensityFloor &&
			PressureFromConservative1D(Q) > kMUSCLPressureFloor;
	}

	bool PrimitiveAdmissible1D(const double W[3])
	{
		return std::isfinite(W[0]) && std::isfinite(W[1]) && std::isfinite(W[2]) &&
			W[0] > kMUSCLDensityFloor && W[2] > kMUSCLPressureFloor;
	}

	void ConservativeToPrimitiveSafe1D(const double Q[3], double W[3], bool& ok)
	{
		double local_Q[3] = { Q[0], Q[1], Q[2] };
		Convar_to_primvar_1D(W, local_Q);
		ok = PrimitiveAdmissible1D(W);
	}

	void PrimitiveToConservativeSafe1D(const double W[3], double Q[3])
	{
		double local_W[3] = { W[0], W[1], W[2] };
		Primvar_to_convar_1D(Q, local_W);
	}

	void PhysicalFluxFromConservative1D(const double Q[3], double flux[3])
	{
		const double rho = Q[0];
		const double momentum = Q[1];
		const double energy = Q[2];
		const double u = momentum / rho;
		const double p = PressureFromConservative1D(Q);
		flux[0] = momentum;
		flux[1] = momentum * u + p;
		flux[2] = u * (energy + p);
	}

	double Minmod(double a, double b)
	{
		if (a * b <= 0.0)
		{
			return 0.0;
		}
		const double sign = (a > 0.0) ? 1.0 : -1.0;
		return sign * std::min(std::fabs(a), std::fabs(b));
	}

	bool FluxFinite1D(const double flux[3])
	{
		return std::isfinite(flux[0]) && std::isfinite(flux[1]) && std::isfinite(flux[2]);
	}

	void FirstOrderKFVSFaceFlux1D(
		const double left_Q[3],
		const double right_Q[3],
		double dt,
		double flux_avg[3])
	{
		KFVS1_TimeAveragedFlux1D(left_Q, right_Q, dt, flux_avg);
	}

	void ComputeFirstOrderInternalFluxes1D(
		GKSSubcellBranch1D& branch,
		double dt)
	{
		for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
		{
			for (int f = 0; f < 2; ++f)
			{
				FirstOrderKFVSFaceFlux1D(
					branch.cell[e].low_dof[f],
					branch.cell[e].low_dof[f + 1],
					dt,
					branch.cell[e].internal_flux[f]);
			}
		}
	}

	void ComputeMUSCLHancockInternalFluxes1D(
		GKSSubcellBranch1D& branch,
		double dt)
	{
		branch.muscl_stats = GKSSubcellMUSCLHancockStats1D();

		const double big_dx =
			branch.geom.subcell_length_x[0] +
			branch.geom.subcell_length_x[1] +
			branch.geom.subcell_length_x[2];
		double center_x[3];
		double face_x[4];
		for (int j = 0; j < 3; ++j)
		{
			center_x[j] = 0.5 * big_dx * branch.geom.subcell_center_s[j];
		}
		for (int f = 0; f < 4; ++f)
		{
			face_x[f] = 0.5 * big_dx * branch.geom.subcell_face_s[f];
		}

		for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
		{
			double W[3][3];
			bool center_ok[3];
			for (int j = 0; j < 3; ++j)
			{
				ConservativeToPrimitiveSafe1D(branch.cell[e].low_dof[j], W[j], center_ok[j]);
			}

			double Q_left_star[3][3];
			double Q_right_star[3][3];
			for (int j = 0; j < 3; ++j)
			{
				branch.muscl_stats.reconstructed_cells++;

				double WL[3] = { W[j][0], W[j][1], W[j][2] };
				double WR[3] = { W[j][0], W[j][1], W[j][2] };
				bool use_zero_slope = !center_ok[j];

				if (!use_zero_slope && j > 0 && j < 2 && center_ok[j - 1] && center_ok[j + 1])
				{
					for (int m = 0; m < 3; ++m)
					{
						const double slope_left = (W[j][m] - W[j - 1][m]) / (center_x[j] - center_x[j - 1]);
						const double slope_right = (W[j + 1][m] - W[j][m]) / (center_x[j + 1] - center_x[j]);
						const double slope = Minmod(slope_left, slope_right);
						WL[m] = W[j][m] - slope * (center_x[j] - face_x[j]);
						WR[m] = W[j][m] + slope * (face_x[j + 1] - center_x[j]);
					}
					use_zero_slope = !PrimitiveAdmissible1D(WL) || !PrimitiveAdmissible1D(WR);
				}

				if (use_zero_slope)
				{
					branch.muscl_stats.zero_slope_fallback_cells++;
					for (int m = 0; m < 3; ++m)
					{
						Q_left_star[j][m] = branch.cell[e].low_dof[j][m];
						Q_right_star[j][m] = branch.cell[e].low_dof[j][m];
					}
					continue;
				}

				double QL[3], QR[3];
				PrimitiveToConservativeSafe1D(WL, QL);
				PrimitiveToConservativeSafe1D(WR, QR);

				double FL[3], FR[3];
				PhysicalFluxFromConservative1D(QL, FL);
				PhysicalFluxFromConservative1D(QR, FR);

				for (int m = 0; m < 3; ++m)
				{
					const double half_update =
						0.5 * dt / branch.geom.subcell_length_x[j] * (FR[m] - FL[m]);
					Q_left_star[j][m] = QL[m] - half_update;
					Q_right_star[j][m] = QR[m] - half_update;
				}

				if (!ConservativeAdmissible1D(Q_left_star[j]) ||
					!ConservativeAdmissible1D(Q_right_star[j]))
				{
					branch.muscl_stats.predicted_state_fallback_cells++;
					for (int m = 0; m < 3; ++m)
					{
						Q_left_star[j][m] = branch.cell[e].low_dof[j][m];
						Q_right_star[j][m] = branch.cell[e].low_dof[j][m];
					}
				}
			}

			for (int f = 0; f < 2; ++f)
			{
				bool use_first_order_flux =
					!ConservativeAdmissible1D(Q_right_star[f]) ||
					!ConservativeAdmissible1D(Q_left_star[f + 1]);
				if (!use_first_order_flux)
				{
					FirstOrderKFVSFaceFlux1D(
						Q_right_star[f],
						Q_left_star[f + 1],
						dt,
						branch.cell[e].internal_flux[f]);
					use_first_order_flux = !FluxFinite1D(branch.cell[e].internal_flux[f]);
				}

				if (use_first_order_flux)
				{
					branch.muscl_stats.face_flux_fallback_faces++;
					FirstOrderKFVSFaceFlux1D(
						branch.cell[e].low_dof[f],
						branch.cell[e].low_dof[f + 1],
						dt,
						branch.cell[e].internal_flux[f]);
				}
			}
		}
	}

	void FaceFluxByMode2D_X(
		const double left_Q[4],
		const double right_Q[4],
		double dt,
		double flux_avg[4])
	{
		KFVS1_TimeAveragedFlux2D_X(left_Q, right_Q, dt, flux_avg);
	}

	void FaceFluxByMode2D_Y(
		const double bottom_Q[4],
		const double top_Q[4],
		double dt,
		double flux_avg[4])
	{
		KFVS1_TimeAveragedFlux2D_Y(bottom_Q, top_Q, dt, flux_avg);
	}

	int XFaceIndex2D(int cells_y, int face_i, int cell_j)
	{
		return face_i * cells_y + cell_j;
	}

	int YFaceIndex2D(int cells_y, int cell_i, int face_j)
	{
		return cell_i * (cells_y + 1) + face_j;
	}
}

void GKSSubcellBranchResize1D(
	GKSSubcellBranch1D& branch,
	const GKSFRMesh1D& mesh)
{
	GKSSubcellBuildGeometry1D(branch.geom, mesh.dx);
		branch.cell.resize(mesh.cells);
	branch.muscl_stats = GKSSubcellMUSCLHancockStats1D();
}

void GKSSubcellInitializeFromCurrentDofs1D(
	const GKSFRMesh1D& mesh,
	GKSSubcellBranch1D& branch)
{
	GKSSubcellBranchResize1D(branch, mesh);
	for (int e = 0; e < mesh.cells; ++e)
	{
		// LWFR first-order branch:
		// the current solution DOFs u_j^e themselves are the low-order subcell
		// unknowns, one per weighted subcell, without any extra projection.
		for (int j = 0; j < 3; ++j)
		{
			for (int m = 0; m < 3; ++m)
			{
				branch.cell[e].low_dof[j][m] = mesh.cell[e].Q[j][m];
			}
		}
	}
}

void GKSSubcellComputeInternalFluxes1D(
	GKSSubcellBranch1D& branch,
	double dt)
{
	GKSSubcellComputeInternalFluxes1D(branch, dt, low_order_type);
}

void GKSSubcellComputeInternalFluxes1D(
	GKSSubcellBranch1D& branch,
	double dt,
	GKSSubcellLowOrderType mode)
{
	if (mode == MUSCL_HANCOCK)
	{
		ComputeMUSCLHancockInternalFluxes1D(branch, dt);
		return;
	}
	branch.muscl_stats = GKSSubcellMUSCLHancockStats1D();
	ComputeFirstOrderInternalFluxes1D(branch, dt);
}

void GKSSubcellComputeLowFaceFluxes1D(
	const GKSSubcellBranch1D& branch,
	double dt,
	GKSFRBoundary1D boundary,
	std::vector<GKSFRFaceFlux1D>& face_fluxes)
{
	const int cells = static_cast<int>(branch.cell.size());
	face_fluxes.assign(cells + 1, GKSFRFaceFlux1D());

	for (int e = 1; e < cells; ++e)
	{
		FirstOrderKFVSFaceFlux1D(
			branch.cell[e - 1].low_dof[2],
			branch.cell[e].low_dof[0],
			dt,
			face_fluxes[e].F);
	}

	if (boundary == gksfr_periodic)
	{
		FirstOrderKFVSFaceFlux1D(
			branch.cell[cells - 1].low_dof[2],
			branch.cell[0].low_dof[0],
			dt,
			face_fluxes[0].F);
		for (int m = 0; m < 3; ++m)
		{
			face_fluxes[cells].F[m] = face_fluxes[0].F[m];
		}
	}
	else if (boundary == gksfr_reflective)
	{
		double mirrored_Q[3];
		ReflectNormalVelocityState(branch.cell[0].low_dof[0], mirrored_Q);
		FirstOrderKFVSFaceFlux1D(
			mirrored_Q,
			branch.cell[0].low_dof[0],
			dt,
			face_fluxes[0].F);

		ReflectNormalVelocityState(branch.cell[cells - 1].low_dof[2], mirrored_Q);
		FirstOrderKFVSFaceFlux1D(
			branch.cell[cells - 1].low_dof[2],
			mirrored_Q,
			dt,
			face_fluxes[cells].F);
	}
	else
	{
		// Free boundary: use transmissive states at the physical boundaries.
		FirstOrderKFVSFaceFlux1D(
			branch.cell[0].low_dof[0],
			branch.cell[0].low_dof[0],
			dt,
			face_fluxes[0].F);
		FirstOrderKFVSFaceFlux1D(
			branch.cell[cells - 1].low_dof[2],
			branch.cell[cells - 1].low_dof[2],
			dt,
			face_fluxes[cells].F);
	}
}

void GKSSubcellComputeElementResiduals1D(
	GKSSubcellBranch1D& branch)
{
	const double* dx = branch.geom.subcell_length_x;
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		// First-order FV residual on the weighted subcell geometry:
		// j=0: -(f_{1/2}^e - F_{e-1/2}) / (w_0 dx_e)
		// j=1: -(f_{3/2}^e - f_{1/2}^e) / (w_1 dx_e)
		// j=2: -(F_{e+1/2} - f_{3/2}^e) / (w_2 dx_e)
		// The boundary-face contributions are added separately in
		// GKSSubcellAddFaceResiduals1D.
		for (int m = 0; m < 3; ++m)
		{
			branch.cell[e].residual_subcell[0][m] = -branch.cell[e].internal_flux[0][m] / dx[0];
			branch.cell[e].residual_subcell[1][m] =
				-(branch.cell[e].internal_flux[1][m] - branch.cell[e].internal_flux[0][m]) / dx[1];
			branch.cell[e].residual_subcell[2][m] = branch.cell[e].internal_flux[1][m] / dx[2];
		}
	}
}

void GKSSubcellAddFaceResiduals1D(
	GKSSubcellBranch1D& branch,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes)
{
	const double* dx = branch.geom.subcell_length_x;
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		for (int m = 0; m < 3; ++m)
		{
			branch.cell[e].residual_subcell[0][m] += face_fluxes[e].F[m] / dx[0];
			branch.cell[e].residual_subcell[2][m] -= face_fluxes[e + 1].F[m] / dx[2];
		}
	}
}

void GKSSubcellScaleResiduals1D(
	GKSSubcellBranch1D& branch,
	const std::vector<double>& alpha_cell)
{
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int m = 0; m < 3; ++m)
			{
				branch.cell[e].residual_subcell[i][m] *= alpha_cell[e];
			}
		}
	}
}

void GKSSubcellUpdateLowOrderDofs1D(
	GKSSubcellBranch1D& branch,
	double dt)
{
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int m = 0; m < 3; ++m)
			{
				branch.cell[e].low_dof_new[i][m] =
					branch.cell[e].low_dof[i][m] + dt * branch.cell[e].residual_subcell[i][m];
			}
		}
		GKSSubcellBigCellAverageFromLowOrderDofs1D(
			branch.geom,
			branch.cell[e].low_dof_new,
			branch.cell[e].cell_avg_new);
	}
}

void GKSSubcellFallbackBadMUSCLUpdates1D(
	GKSSubcellBranch1D& branch,
	const GKSSubcellBranch1D& kfvs_fallback_branch,
	double dt)
{
	if (low_order_type != MUSCL_HANCOCK)
	{
		return;
	}
	const int cells = static_cast<int>(branch.cell.size());
	for (int e = 0; e < cells; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			double candidate[3];
			for (int m = 0; m < 3; ++m)
			{
				candidate[m] = branch.cell[e].low_dof[i][m]
					+ dt * branch.cell[e].residual_subcell[i][m];
			}

			if (StateFinite1D(candidate) && candidate[0] > 0.0)
			{
				branch.muscl_stats.min_rho_update =
					std::min(branch.muscl_stats.min_rho_update, candidate[0]);
				branch.muscl_stats.min_p_update =
					std::min(branch.muscl_stats.min_p_update, PressureFromConservative1D(candidate));
			}

			if (!ConservativeAdmissible1D(candidate))
			{
				branch.muscl_stats.final_update_fallback_cells++;
				for (int m = 0; m < 3; ++m)
				{
					branch.cell[e].residual_subcell[i][m] =
						kfvs_fallback_branch.cell[e].residual_subcell[i][m];
				}
			}
		}
	}
}

void GKSSubcellAdvanceWithFaceFluxes1D(
	GKSSubcellBranch1D& branch,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes,
	double dt)
{
	GKSSubcellComputeElementResiduals1D(branch);
	GKSSubcellAddFaceResiduals1D(branch, face_fluxes);
	GKSSubcellUpdateLowOrderDofs1D(branch, dt);
}

void GKSSubcellBranchResize2D(
	GKSSubcellBranch2D& branch,
	const GKSFRMesh2D& mesh)
{
	GKSSubcellBuildGeometry2D(branch.geom, mesh.dx, mesh.dy);
	branch.cells_x = mesh.cells_x;
	branch.cells_y = mesh.cells_y;
	branch.cell.resize(mesh.cells_x * mesh.cells_y);
}

int GKSSubcellCellIndex2D(const GKSSubcellBranch2D& branch, int i, int j)
{
	return i * branch.cells_y + j;
}

void GKSSubcellInitializeFromCurrentDofs2D(
	const GKSFRMesh2D& mesh,
	GKSSubcellBranch2D& branch)
{
	GKSSubcellBranchResize2D(branch, mesh);
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			const GKSFRCell2D& src = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
			GKSSubcellCell2D& dst = branch.cell[GKSSubcellCellIndex2D(branch, i, j)];
			for (int p = 0; p < 3; ++p)
			{
				for (int q = 0; q < 3; ++q)
				{
					for (int m = 0; m < 4; ++m)
					{
						dst.low_dof[p][q][m] = src.Q[p][q][m];
						dst.residual_subcell[p][q][m] = 0.0;
						dst.low_dof_new[p][q][m] = src.Q[p][q][m];
					}
				}
			}
		}
	}
}

void GKSSubcellComputeInternalFluxes2D(
	GKSSubcellBranch2D& branch,
	double dt)
{
	for (int i = 0; i < branch.cells_x; ++i)
	{
		for (int j = 0; j < branch.cells_y; ++j)
		{
			GKSSubcellCell2D& cell = branch.cell[GKSSubcellCellIndex2D(branch, i, j)];
			for (int f = 0; f < 2; ++f)
			{
				for (int q = 0; q < 3; ++q)
				{
					FaceFluxByMode2D_X(
						cell.low_dof[f][q],
						cell.low_dof[f + 1][q],
						dt,
						cell.internal_flux_x[f][q]);
				}
			}
			for (int p = 0; p < 3; ++p)
			{
				for (int f = 0; f < 2; ++f)
				{
					FaceFluxByMode2D_Y(
						cell.low_dof[p][f],
						cell.low_dof[p][f + 1],
						dt,
						cell.internal_flux_y[p][f]);
				}
			}
		}
	}
}

void GKSSubcellComputeLowFaceFluxes2D(
	const GKSSubcellBranch2D& branch,
	double dt,
	GKSFRBoundary2D boundary,
	std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& y_face_fluxes)
{
	x_face_fluxes.assign((branch.cells_x + 1) * branch.cells_y, GKSFRFaceFlux2D());
	y_face_fluxes.assign(branch.cells_x * (branch.cells_y + 1), GKSFRFaceFlux2D());
	if (boundary != gksfr2d_periodic)
	{
		return;
	}
	for (int face_i = 0; face_i <= branch.cells_x; ++face_i)
	{
		const int left_i = (face_i + branch.cells_x - 1) % branch.cells_x;
		const int right_i = face_i % branch.cells_x;
		for (int j = 0; j < branch.cells_y; ++j)
		{
			const GKSSubcellCell2D& left = branch.cell[GKSSubcellCellIndex2D(branch, left_i, j)];
			const GKSSubcellCell2D& right = branch.cell[GKSSubcellCellIndex2D(branch, right_i, j)];
			for (int q = 0; q < 3; ++q)
			{
				FaceFluxByMode2D_X(
					left.low_dof[2][q],
					right.low_dof[0][q],
					dt,
					x_face_fluxes[XFaceIndex2D(branch.cells_y, face_i, j)].F[q]);
			}
		}
	}
	for (int i = 0; i < branch.cells_x; ++i)
	{
		for (int face_j = 0; face_j <= branch.cells_y; ++face_j)
		{
			const int bottom_j = (face_j + branch.cells_y - 1) % branch.cells_y;
			const int top_j = face_j % branch.cells_y;
			const GKSSubcellCell2D& bottom = branch.cell[GKSSubcellCellIndex2D(branch, i, bottom_j)];
			const GKSSubcellCell2D& top = branch.cell[GKSSubcellCellIndex2D(branch, i, top_j)];
			for (int p = 0; p < 3; ++p)
			{
				FaceFluxByMode2D_Y(
					bottom.low_dof[p][2],
					top.low_dof[p][0],
					dt,
					y_face_fluxes[YFaceIndex2D(branch.cells_y, i, face_j)].F[p]);
			}
		}
	}
}

void GKSSubcellComputeElementResiduals2D(
	GKSSubcellBranch2D& branch)
{
	const double* dx = branch.geom.subcell_length_x;
	const double* dy = branch.geom.subcell_length_y;
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		GKSSubcellCell2D& cell = branch.cell[e];
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					double rx = 0.0;
					if (i == 0)
					{
						rx = -cell.internal_flux_x[0][j][m] / dx[0];
					}
					else if (i == 1)
					{
						rx = -(cell.internal_flux_x[1][j][m] - cell.internal_flux_x[0][j][m]) / dx[1];
					}
					else
					{
						rx = cell.internal_flux_x[1][j][m] / dx[2];
					}
					double ry = 0.0;
					if (j == 0)
					{
						ry = -cell.internal_flux_y[i][0][m] / dy[0];
					}
					else if (j == 1)
					{
						ry = -(cell.internal_flux_y[i][1][m] - cell.internal_flux_y[i][0][m]) / dy[1];
					}
					else
					{
						ry = cell.internal_flux_y[i][1][m] / dy[2];
					}
					cell.residual_subcell[i][j][m] = rx + ry;
				}
			}
		}
	}
}

void GKSSubcellAddFaceResiduals2D(
	GKSSubcellBranch2D& branch,
	const std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& y_face_fluxes)
{
	const double* dx = branch.geom.subcell_length_x;
	const double* dy = branch.geom.subcell_length_y;
	for (int i = 0; i < branch.cells_x; ++i)
	{
		for (int j = 0; j < branch.cells_y; ++j)
		{
			GKSSubcellCell2D& cell = branch.cell[GKSSubcellCellIndex2D(branch, i, j)];
			for (int q = 0; q < 3; ++q)
			{
				const double* left_flux = x_face_fluxes[XFaceIndex2D(branch.cells_y, i, j)].F[q];
				const double* right_flux = x_face_fluxes[XFaceIndex2D(branch.cells_y, i + 1, j)].F[q];
				for (int m = 0; m < 4; ++m)
				{
					cell.residual_subcell[0][q][m] += left_flux[m] / dx[0];
					cell.residual_subcell[2][q][m] -= right_flux[m] / dx[2];
				}
			}
			for (int p = 0; p < 3; ++p)
			{
				const double* bottom_flux = y_face_fluxes[YFaceIndex2D(branch.cells_y, i, j)].F[p];
				const double* top_flux = y_face_fluxes[YFaceIndex2D(branch.cells_y, i, j + 1)].F[p];
				for (int m = 0; m < 4; ++m)
				{
					cell.residual_subcell[p][0][m] += bottom_flux[m] / dy[0];
					cell.residual_subcell[p][2][m] -= top_flux[m] / dy[2];
				}
			}
		}
	}
}

void GKSSubcellUpdateLowOrderDofs2D(
	GKSSubcellBranch2D& branch,
	double dt)
{
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		GKSSubcellCell2D& cell = branch.cell[e];
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					cell.low_dof_new[i][j][m] = cell.low_dof[i][j][m]
						+ dt * cell.residual_subcell[i][j][m];
				}
			}
		}
		GKSSubcellBigCellAverageFromLowOrderDofs2D(branch.geom, cell.low_dof_new, cell.cell_avg_new);
	}
}

void GKSSubcellAdvanceWithFaceFluxes2D(
	GKSSubcellBranch2D& branch,
	const std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& y_face_fluxes,
	double dt)
{
	GKSSubcellComputeElementResiduals2D(branch);
	GKSSubcellAddFaceResiduals2D(branch, x_face_fluxes, y_face_fluxes);
	GKSSubcellUpdateLowOrderDofs2D(branch, dt);
}
