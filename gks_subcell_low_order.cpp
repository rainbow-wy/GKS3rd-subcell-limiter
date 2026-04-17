#include "gks_subcell_low_order.h"

#include "gks_kfvs_adapter.h"

GKSSubcellLowOrderType low_order_type = KFVS1;

namespace
{
	void ZeroState(double out[3])
	{
		for (int m = 0; m < 3; ++m)
		{
			out[m] = 0.0;
		}
	}

	void FaceFluxByMode(
		const double left_Q[3],
		const double right_Q[3],
		double dt,
		double flux_avg[3])
	{
		// The framework is already organized by low_order_type. In the current
		// stage, KFVS1 is the fully implemented robust branch. MUSCL-Hancock will
		// reuse the same face-flux kernel later after reconstructed face states are
		// added, so for now both modes share this KFVS1 core.
		KFVS1_TimeAveragedFlux1D(left_Q, right_Q, dt, flux_avg);
	}
}

void GKSSubcellBranchResize1D(
	GKSSubcellBranch1D& branch,
	const GKSFRMesh1D& mesh)
{
	GKSSubcellBuildGeometry1D(branch.geom, mesh.dx);
		branch.cell.resize(mesh.cells);
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
	for (int e = 0; e < static_cast<int>(branch.cell.size()); ++e)
	{
		for (int f = 0; f < 2; ++f)
		{
			FaceFluxByMode(
				branch.cell[e].low_dof[f],
				branch.cell[e].low_dof[f + 1],
				dt,
				branch.cell[e].internal_flux[f]);
		}
	}
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
		FaceFluxByMode(
			branch.cell[e - 1].low_dof[2],
			branch.cell[e].low_dof[0],
			dt,
			face_fluxes[e].F);
	}

	if (boundary == gksfr_periodic)
	{
		FaceFluxByMode(
			branch.cell[cells - 1].low_dof[2],
			branch.cell[0].low_dof[0],
			dt,
			face_fluxes[0].F);
		for (int m = 0; m < 3; ++m)
		{
			face_fluxes[cells].F[m] = face_fluxes[0].F[m];
		}
	}
	else
	{
		// Free boundary: use transmissive states at the physical boundaries.
		FaceFluxByMode(
			branch.cell[0].low_dof[0],
			branch.cell[0].low_dof[0],
			dt,
			face_fluxes[0].F);
		FaceFluxByMode(
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

void GKSSubcellAdvanceWithFaceFluxes1D(
	GKSSubcellBranch1D& branch,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes,
	double dt)
{
	GKSSubcellComputeElementResiduals1D(branch);
	GKSSubcellAddFaceResiduals1D(branch, face_fluxes);
	GKSSubcellUpdateLowOrderDofs1D(branch, dt);
}
