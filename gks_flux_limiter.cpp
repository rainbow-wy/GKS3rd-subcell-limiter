#include "gks_flux_limiter.h"

#include "gks_basic.h"

#include <algorithm>
#include <cmath>
#include <limits>

GKSFluxLimiterParam1D::GKSFluxLimiterParam1D()
	: rho_floor(1.0e-12),
	  p_floor(1.0e-12)
{
}

GKSFluxLimiterParam2D::GKSFluxLimiterParam2D()
	: rho_floor(1.0e-12),
	  p_floor(1.0e-12),
	  kx(0.5),
	  ky(0.5)
{
}

GKSFluxLimiterDiag2D::GKSFluxLimiterDiag2D()
	: limited_faces_x(0),
	  limited_faces_y(0),
	  min_theta_x(1.0),
	  min_theta_y(1.0)
{
}

namespace
{
	void ConservativeToDiagnostic(const double U[3], double& rho, double& u, double& p)
	{
		rho = U[0];
		if (!std::isfinite(U[0]) || std::fabs(U[0]) <= 1.0e-300)
		{
			u = std::numeric_limits<double>::quiet_NaN();
			p = std::numeric_limits<double>::quiet_NaN();
			return;
		}
		u = U[1] / U[0];
		p = Pressure(U[0], U[1], U[2]);
	}

	bool StateFinite(const double U[3])
	{
		return std::isfinite(U[0]) && std::isfinite(U[1]) && std::isfinite(U[2]);
	}

	bool StateAdmissible(const double U[3], double rho_floor, double p_floor)
	{
		if (!StateFinite(U))
		{
			return false;
		}
		if (U[0] <= rho_floor)
		{
			return false;
		}
		return Pressure(U[0], U[1], U[2]) > p_floor;
	}

	double DensityTheta(double rho_safe, double rho_candidate, double rho_floor)
	{
		if (!std::isfinite(rho_safe) || !std::isfinite(rho_candidate))
		{
			return 0.0;
		}
		if (rho_safe <= rho_floor)
		{
			return 0.0;
		}
		if (rho_candidate >= rho_floor)
		{
			return 1.0;
		}
		const double denom = rho_safe - rho_candidate;
		if (denom <= 1.0e-14)
		{
			return 0.0;
		}
		double theta = (rho_safe - rho_floor) / denom;
		theta = std::max(0.0, std::min(1.0, theta));
		return theta;
	}

	double PressureFromConservative(const double U[3])
	{
		return Pressure(U[0], U[1], U[2]);
	}

	bool StateFinite2D(const double U[4])
	{
		return std::isfinite(U[0]) && std::isfinite(U[1]) &&
			std::isfinite(U[2]) && std::isfinite(U[3]);
	}

	double Pressure2DFromConservative(const double U[4])
	{
		if (U[0] <= 0.0)
		{
			return -std::numeric_limits<double>::infinity();
		}
		return (Gamma - 1.0) *
			(U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
	}

	bool StateAdmissible2D(const double U[4], double rho_floor, double p_floor)
	{
		return StateFinite2D(U) && U[0] > rho_floor &&
			Pressure2DFromConservative(U) > p_floor;
	}

	double PressureTheta(const double U_safe[3], const double U_candidate[3], double p_floor)
	{
		if (!StateFinite(U_safe) || !StateFinite(U_candidate))
		{
			return 0.0;
		}
		if (U_safe[0] <= 0.0 || PressureFromConservative(U_safe) <= p_floor)
		{
			return 0.0;
		}
		if (PressureFromConservative(U_candidate) >= p_floor)
		{
			return 1.0;
		}

		const double eps_bar = p_floor / (Gamma - 1.0);
		const double dr = U_candidate[0] - U_safe[0];
		const double dm = U_candidate[1] - U_safe[1];
		const double dE = U_candidate[2] - U_safe[2];

		const double A = 2.0 * dr * dE - dm * dm;
		const double B = 2.0 * U_safe[0] * dE
			+ 2.0 * dr * (U_safe[2] - eps_bar)
			- 2.0 * U_safe[1] * dm;
		const double C = 2.0 * U_safe[0] * (U_safe[2] - eps_bar)
			- U_safe[1] * U_safe[1];

		if (std::fabs(A) < 1.0e-14)
		{
			if (std::fabs(B) < 1.0e-14)
			{
				return 0.0;
			}
			const double root = -C / B;
			return std::max(0.0, std::min(1.0, root));
		}

		double disc = B * B - 4.0 * A * C;
		if (disc < 0.0)
		{
			disc = 0.0;
		}
		const double sqrt_disc = std::sqrt(disc);
		const double r1 = (-B - sqrt_disc) / (2.0 * A);
		const double r2 = (-B + sqrt_disc) / (2.0 * A);

		double theta = 0.0;
		bool found = false;
		if (r1 >= 0.0 && r1 <= 1.0)
		{
			theta = r1;
			found = true;
		}
		if (r2 >= 0.0 && r2 <= 1.0)
		{
			theta = found ? std::min(theta, r2) : r2;
			found = true;
		}
		if (!found)
		{
			return 0.0;
		}
		return std::max(0.0, std::min(1.0, theta));
	}

	double PressureTheta2D(const double U_safe[4], const double U_candidate[4], double p_floor)
	{
		if (!StateFinite2D(U_safe) || !StateFinite2D(U_candidate))
		{
			return 0.0;
		}
		if (U_safe[0] <= 0.0 || Pressure2DFromConservative(U_safe) <= p_floor)
		{
			return 0.0;
		}
		if (Pressure2DFromConservative(U_candidate) >= p_floor)
		{
			return 1.0;
		}
		const double eta_p = p_floor / (Gamma - 1.0);
		const double dr = U_candidate[0] - U_safe[0];
		const double dmx = U_candidate[1] - U_safe[1];
		const double dmy = U_candidate[2] - U_safe[2];
		const double dE = U_candidate[3] - U_safe[3];
		const double A = 2.0 * dr * dE - dmx * dmx - dmy * dmy;
		const double B = 2.0 * U_safe[0] * dE
			+ 2.0 * dr * (U_safe[3] - eta_p)
			- 2.0 * U_safe[1] * dmx
			- 2.0 * U_safe[2] * dmy;
		const double C = 2.0 * U_safe[0] * (U_safe[3] - eta_p)
			- U_safe[1] * U_safe[1]
			- U_safe[2] * U_safe[2];

		if (std::fabs(A) < 1.0e-14)
		{
			if (std::fabs(B) < 1.0e-14)
			{
				return 0.0;
			}
			return std::max(0.0, std::min(1.0, -C / B));
		}
		double disc = B * B - 4.0 * A * C;
		if (disc < 0.0)
		{
			disc = 0.0;
		}
		const double sqrt_disc = std::sqrt(disc);
		const double r1 = (-B - sqrt_disc) / (2.0 * A);
		const double r2 = (-B + sqrt_disc) / (2.0 * A);
		double theta = 0.0;
		bool found = false;
		if (r1 >= 0.0 && r1 <= 1.0)
		{
			theta = r1;
			found = true;
		}
		if (r2 >= 0.0 && r2 <= 1.0)
		{
			theta = found ? std::min(theta, r2) : r2;
			found = true;
		}
		return found ? std::max(0.0, std::min(1.0, theta)) : 0.0;
	}

	void BlendFlux4(const double high[4], const double low[4], double alpha, double out[4])
	{
		bool high_finite = true;
		bool low_finite = true;
		for (int m = 0; m < 4; ++m)
		{
			high_finite = high_finite && std::isfinite(high[m]);
			low_finite = low_finite && std::isfinite(low[m]);
		}
		if (!high_finite && low_finite)
		{
			for (int m = 0; m < 4; ++m)
			{
				out[m] = low[m];
			}
			return;
		}
		if (alpha <= 1.0e-14)
		{
			for (int m = 0; m < 4; ++m)
			{
				out[m] = high[m];
			}
			return;
		}
		if (alpha >= 1.0 - 1.0e-14)
		{
			for (int m = 0; m < 4; ++m)
			{
				out[m] = low[m];
			}
			return;
		}
		for (int m = 0; m < 4; ++m)
		{
			out[m] = (1.0 - alpha) * high[m] + alpha * low[m];
		}
	}

	int XFaceIndex2D(int cells_y, int face_i, int cell_j)
	{
		return face_i * cells_y + cell_j;
	}

	int YFaceIndex2D(int cells_y, int cell_i, int face_j)
	{
		return cell_i * (cells_y + 1) + face_j;
	}

	void UpdateXLeftState2D(
		const GKSSubcellBranch2D& branch,
		int e_left,
		int q,
		const double face_flux[4],
		double dt,
		double kx,
		double U_out[4])
	{
		const GKSSubcellCell2D& cell = branch.cell[e_left];
		const double dx = branch.geom.subcell_length_x[2];
		for (int m = 0; m < 4; ++m)
		{
			U_out[m] = cell.low_dof[2][q][m]
				- dt / (kx * dx) * (face_flux[m] - cell.internal_flux_x[1][q][m]);
		}
	}

	void UpdateXRightState2D(
		const GKSSubcellBranch2D& branch,
		int e_right,
		int q,
		const double face_flux[4],
		double dt,
		double kx,
		double U_out[4])
	{
		const GKSSubcellCell2D& cell = branch.cell[e_right];
		const double dx = branch.geom.subcell_length_x[0];
		for (int m = 0; m < 4; ++m)
		{
			U_out[m] = cell.low_dof[0][q][m]
				- dt / (kx * dx) * (cell.internal_flux_x[0][q][m] - face_flux[m]);
		}
	}

	void UpdateYBottomState2D(
		const GKSSubcellBranch2D& branch,
		int e_bottom,
		int p,
		const double face_flux[4],
		double dt,
		double ky,
		double U_out[4])
	{
		const GKSSubcellCell2D& cell = branch.cell[e_bottom];
		const double dy = branch.geom.subcell_length_y[2];
		for (int m = 0; m < 4; ++m)
		{
			U_out[m] = cell.low_dof[p][2][m]
				- dt / (ky * dy) * (face_flux[m] - cell.internal_flux_y[p][1][m]);
		}
	}

	void UpdateYTopState2D(
		const GKSSubcellBranch2D& branch,
		int e_top,
		int p,
		const double face_flux[4],
		double dt,
		double ky,
		double U_out[4])
	{
		const GKSSubcellCell2D& cell = branch.cell[e_top];
		const double dy = branch.geom.subcell_length_y[0];
		for (int m = 0; m < 4; ++m)
		{
			U_out[m] = cell.low_dof[p][0][m]
				- dt / (ky * dy) * (cell.internal_flux_y[p][0][m] - face_flux[m]);
		}
	}

	void BlendFlux(
		const double F_high[3],
		const double F_low[3],
		double alpha,
		double F_out[3])
	{
		if (alpha <= 1.0e-14)
		{
			for (int m = 0; m < 3; ++m)
			{
				F_out[m] = F_high[m];
			}
			return;
		}
		if (alpha >= 1.0 - 1.0e-14)
		{
			for (int m = 0; m < 3; ++m)
			{
				F_out[m] = F_low[m];
			}
			return;
		}
		for (int m = 0; m < 3; ++m)
		{
			F_out[m] = (1.0 - alpha) * F_high[m] + alpha * F_low[m];
		}
	}

	void UpdateLeftAdjacentState(
		const GKSSubcellBranch1D& branch,
		int e_left,
		const double face_flux[3],
		double dt,
		double U_out[3])
	{
		const int j = 2;
		const double dx = branch.geom.subcell_length_x[j];
		for (int m = 0; m < 3; ++m)
		{
			U_out[m] = branch.cell[e_left].low_dof[j][m]
				- dt / dx * (face_flux[m] - branch.cell[e_left].internal_flux[1][m]);
		}
	}

	void UpdateRightAdjacentState(
		const GKSSubcellBranch1D& branch,
		int e_right,
		const double face_flux[3],
		double dt,
		double U_out[3])
	{
		const int j = 0;
		const double dx = branch.geom.subcell_length_x[j];
		for (int m = 0; m < 3; ++m)
		{
			U_out[m] = branch.cell[e_right].low_dof[j][m]
				- dt / dx * (branch.cell[e_right].internal_flux[0][m] - face_flux[m]);
		}
	}

	void StoreFaceStateDiagnostics(
		GKSFluxLimiterDiag1D& diag,
		int f,
		const GKSFluxLimiterParam1D& param,
		const double UL_safe[3],
		const double UR_safe[3],
		const double UL_candidate[3],
		const double UR_candidate[3])
	{
		ConservativeToDiagnostic(UL_safe, diag.left_safe_rho[f], diag.left_safe_u[f], diag.left_safe_p[f]);
		ConservativeToDiagnostic(UR_safe, diag.right_safe_rho[f], diag.right_safe_u[f], diag.right_safe_p[f]);
		ConservativeToDiagnostic(UL_candidate, diag.left_candidate_rho[f], diag.left_candidate_u[f], diag.left_candidate_p[f]);
		ConservativeToDiagnostic(UR_candidate, diag.right_candidate_rho[f], diag.right_candidate_u[f], diag.right_candidate_p[f]);

		diag.left_safe_rho_bad[f] = !StateFinite(UL_safe) || UL_safe[0] <= param.rho_floor;
		diag.left_safe_p_bad[f] = !StateFinite(UL_safe) || PressureFromConservative(UL_safe) <= param.p_floor;
		diag.right_safe_rho_bad[f] = !StateFinite(UR_safe) || UR_safe[0] <= param.rho_floor;
		diag.right_safe_p_bad[f] = !StateFinite(UR_safe) || PressureFromConservative(UR_safe) <= param.p_floor;

		diag.left_candidate_rho_bad[f] = !StateFinite(UL_candidate) || UL_candidate[0] <= param.rho_floor;
		diag.left_candidate_p_bad[f] = !StateFinite(UL_candidate) || PressureFromConservative(UL_candidate) <= param.p_floor;
		diag.right_candidate_rho_bad[f] = !StateFinite(UR_candidate) || UR_candidate[0] <= param.rho_floor;
		diag.right_candidate_p_bad[f] = !StateFinite(UR_candidate) || PressureFromConservative(UR_candidate) <= param.p_floor;
	}
}

void GKSFluxLimiterApply1D(
	const GKSSubcellBranch1D& branch,
	const std::vector<double>& alpha_cell,
	const std::vector<GKSFRFaceFlux1D>& high_face_fluxes,
	const std::vector<GKSFRFaceFlux1D>& low_face_fluxes,
	double dt,
	GKSFRBoundary1D boundary,
	const GKSFluxLimiterParam1D& param,
	std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
	GKSFluxLimiterDiag1D& diag)
{
	const int cells = static_cast<int>(branch.cell.size());
	final_face_fluxes.assign(cells + 1, GKSFRFaceFlux1D());
	diag.alpha_face.assign(cells + 1, 0.0);
	diag.theta_rho.assign(cells + 1, 1.0);
	diag.theta_p.assign(cells + 1, 1.0);
	diag.theta.assign(cells + 1, 1.0);
	diag.left_safe_rho.assign(cells + 1, 0.0);
	diag.left_safe_u.assign(cells + 1, 0.0);
	diag.left_safe_p.assign(cells + 1, 0.0);
	diag.right_safe_rho.assign(cells + 1, 0.0);
	diag.right_safe_u.assign(cells + 1, 0.0);
	diag.right_safe_p.assign(cells + 1, 0.0);
	diag.left_candidate_rho.assign(cells + 1, 0.0);
	diag.left_candidate_u.assign(cells + 1, 0.0);
	diag.left_candidate_p.assign(cells + 1, 0.0);
	diag.right_candidate_rho.assign(cells + 1, 0.0);
	diag.right_candidate_u.assign(cells + 1, 0.0);
	diag.right_candidate_p.assign(cells + 1, 0.0);
	diag.left_candidate_rho_bad.assign(cells + 1, 0);
	diag.left_candidate_p_bad.assign(cells + 1, 0);
	diag.right_candidate_rho_bad.assign(cells + 1, 0);
	diag.right_candidate_p_bad.assign(cells + 1, 0);
	diag.left_safe_rho_bad.assign(cells + 1, 0);
	diag.left_safe_p_bad.assign(cells + 1, 0);
	diag.right_safe_rho_bad.assign(cells + 1, 0);
	diag.right_safe_p_bad.assign(cells + 1, 0);

	for (int f = 0; f <= cells; ++f)
	{
		double alpha_face = 0.0;
		if (f == 0)
		{
			alpha_face = (boundary == gksfr_periodic)
				? 0.5 * (alpha_cell[cells - 1] + alpha_cell[0])
				: alpha_cell[0];
		}
		else if (f == cells)
		{
			alpha_face = (boundary == gksfr_periodic)
				? diag.alpha_face[0]
				: alpha_cell[cells - 1];
		}
		else
		{
			alpha_face = 0.5 * (alpha_cell[f - 1] + alpha_cell[f]);
		}
		diag.alpha_face[f] = alpha_face;

		double F_ig[3];
		BlendFlux(high_face_fluxes[f].F, low_face_fluxes[f].F, alpha_face, F_ig);

		if (boundary != gksfr_periodic &&
			(f == 0 || f == cells))
		{
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[f].F[m] = F_ig[m];
			}
			continue;
		}

		const int e_left = (f == 0) ? cells - 1 : f - 1;
		const int e_right = (f == cells) ? 0 : f;

		double UL_safe[3], UR_safe[3];
		double UL_ig[3], UR_ig[3];
		UpdateLeftAdjacentState(branch, e_left, low_face_fluxes[f].F, dt, UL_safe);
		UpdateRightAdjacentState(branch, e_right, low_face_fluxes[f].F, dt, UR_safe);
		UpdateLeftAdjacentState(branch, e_left, F_ig, dt, UL_ig);
		UpdateRightAdjacentState(branch, e_right, F_ig, dt, UR_ig);
		StoreFaceStateDiagnostics(diag, f, param, UL_safe, UR_safe, UL_ig, UR_ig);

		const bool safe_ok =
			StateAdmissible(UL_safe, param.rho_floor, param.p_floor) &&
			StateAdmissible(UR_safe, param.rho_floor, param.p_floor);
		if (!safe_ok)
		{
			diag.theta_rho[f] = 0.0;
			diag.theta_p[f] = 0.0;
			diag.theta[f] = 0.0;
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[f].F[m] = low_face_fluxes[f].F[m];
			}
			continue;
		}

		const bool candidate_ok =
			StateAdmissible(UL_ig, param.rho_floor, param.p_floor) &&
			StateAdmissible(UR_ig, param.rho_floor, param.p_floor);
		if (candidate_ok)
		{
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[f].F[m] = F_ig[m];
			}
			continue;
		}

		const double theta_rho_l = DensityTheta(UL_safe[0], UL_ig[0], param.rho_floor);
		const double theta_rho_r = DensityTheta(UR_safe[0], UR_ig[0], param.rho_floor);
		const double theta_rho = std::min(1.0, std::min(theta_rho_l, theta_rho_r));
		diag.theta_rho[f] = theta_rho;

		double F_rho[3];
		BlendFlux(F_ig, low_face_fluxes[f].F, 1.0 - theta_rho, F_rho);

		double UL_rho[3], UR_rho[3];
		UpdateLeftAdjacentState(branch, e_left, F_rho, dt, UL_rho);
		UpdateRightAdjacentState(branch, e_right, F_rho, dt, UR_rho);

		const double theta_p_l = PressureTheta(UL_safe, UL_rho, param.p_floor);
		const double theta_p_r = PressureTheta(UR_safe, UR_rho, param.p_floor);
		const double theta_p = std::min(1.0, std::min(theta_p_l, theta_p_r));
		diag.theta_p[f] = theta_p;

		const double theta = theta_rho * theta_p;
		diag.theta[f] = theta;

		if (theta <= 1.0e-14)
		{
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[f].F[m] = low_face_fluxes[f].F[m];
			}
		}
		else if (theta >= 1.0 - 1.0e-14)
		{
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[f].F[m] = F_ig[m];
			}
		}
		else
		{
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[f].F[m] = theta * F_ig[m] + (1.0 - theta) * low_face_fluxes[f].F[m];
			}
		}
	}

	if (boundary == gksfr_periodic)
	{
		for (int m = 0; m < 3; ++m)
		{
			final_face_fluxes[cells].F[m] = final_face_fluxes[0].F[m];
		}
	diag.theta_rho[cells] = diag.theta_rho[0];
	diag.theta_p[cells] = diag.theta_p[0];
	diag.theta[cells] = diag.theta[0];
	diag.alpha_face[cells] = diag.alpha_face[0];
	diag.left_safe_rho[cells] = diag.left_safe_rho[0];
	diag.left_safe_u[cells] = diag.left_safe_u[0];
	diag.left_safe_p[cells] = diag.left_safe_p[0];
	diag.right_safe_rho[cells] = diag.right_safe_rho[0];
	diag.right_safe_u[cells] = diag.right_safe_u[0];
	diag.right_safe_p[cells] = diag.right_safe_p[0];
	diag.left_candidate_rho[cells] = diag.left_candidate_rho[0];
	diag.left_candidate_u[cells] = diag.left_candidate_u[0];
	diag.left_candidate_p[cells] = diag.left_candidate_p[0];
	diag.right_candidate_rho[cells] = diag.right_candidate_rho[0];
	diag.right_candidate_u[cells] = diag.right_candidate_u[0];
	diag.right_candidate_p[cells] = diag.right_candidate_p[0];
	diag.left_candidate_rho_bad[cells] = diag.left_candidate_rho_bad[0];
	diag.left_candidate_p_bad[cells] = diag.left_candidate_p_bad[0];
	diag.right_candidate_rho_bad[cells] = diag.right_candidate_rho_bad[0];
	diag.right_candidate_p_bad[cells] = diag.right_candidate_p_bad[0];
	diag.left_safe_rho_bad[cells] = diag.left_safe_rho_bad[0];
	diag.left_safe_p_bad[cells] = diag.left_safe_p_bad[0];
	diag.right_safe_rho_bad[cells] = diag.right_safe_rho_bad[0];
	diag.right_safe_p_bad[cells] = diag.right_safe_p_bad[0];
	}
}

void GKSFluxLimiterApply2D(
	const GKSSubcellBranch2D& branch,
	const std::vector<double>& alpha_cell,
	const std::vector<GKSFRFaceFlux2D>& high_x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& high_y_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& low_x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& low_y_face_fluxes,
	double dt,
	GKSFRBoundary2D boundary,
	const GKSFluxLimiterParam2D& param,
	std::vector<GKSFRFaceFlux2D>& final_x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& final_y_face_fluxes,
	GKSFluxLimiterDiag2D& diag)
{
	const int nx = branch.cells_x;
	const int ny = branch.cells_y;
	final_x_face_fluxes.assign((nx + 1) * ny, GKSFRFaceFlux2D());
	final_y_face_fluxes.assign(nx * (ny + 1), GKSFRFaceFlux2D());
	const int x_count = (nx + 1) * ny * 3;
	const int y_count = nx * (ny + 1) * 3;
	diag.alpha_face_x.assign(x_count, 0.0);
	diag.alpha_face_y.assign(y_count, 0.0);
	diag.theta_rho_x.assign(x_count, 1.0);
	diag.theta_rho_y.assign(y_count, 1.0);
	diag.theta_p_x.assign(x_count, 1.0);
	diag.theta_p_y.assign(y_count, 1.0);
	diag.theta_x.assign(x_count, 1.0);
	diag.theta_y.assign(y_count, 1.0);
	diag.limited_faces_x = 0;
	diag.limited_faces_y = 0;
	diag.min_theta_x = 1.0;
	diag.min_theta_y = 1.0;
	if (nx <= 0 || ny <= 0)
	{
		return;
	}
	const double kx = std::max(1.0e-12, param.kx);
	const double ky = std::max(1.0e-12, param.ky);

	for (int face_i = 0; face_i <= nx; ++face_i)
	{
		const bool physical_boundary_face = boundary != gksfr2d_periodic && (face_i == 0 || face_i == nx);
		const int left_i = (boundary == gksfr2d_periodic)
			? (face_i + nx - 1) % nx
			: std::max(0, face_i - 1);
		const int right_i = (boundary == gksfr2d_periodic)
			? face_i % nx
			: std::min(nx - 1, face_i);
		for (int j = 0; j < ny; ++j)
		{
			const int e_left = GKSSubcellCellIndex2D(branch, left_i, j);
			const int e_right = GKSSubcellCellIndex2D(branch, right_i, j);
			const double alpha_face = physical_boundary_face
				? (face_i == 0 ? alpha_cell[e_right] : alpha_cell[e_left])
				: 0.5 * (alpha_cell[e_left] + alpha_cell[e_right]);
			const int face_index = XFaceIndex2D(ny, face_i, j);
			for (int q = 0; q < 3; ++q)
			{
				const int diag_index = (face_index * 3 + q);
				diag.alpha_face_x[diag_index] = alpha_face;
				double F_ig[4];
				BlendFlux4(high_x_face_fluxes[face_index].F[q], low_x_face_fluxes[face_index].F[q], alpha_face, F_ig);
				if (physical_boundary_face)
				{
					double U_safe[4], U_ig[4];
					if (face_i == 0)
					{
						UpdateXRightState2D(branch, e_right, q, low_x_face_fluxes[face_index].F[q], dt, kx, U_safe);
						UpdateXRightState2D(branch, e_right, q, F_ig, dt, kx, U_ig);
					}
					else
					{
						UpdateXLeftState2D(branch, e_left, q, low_x_face_fluxes[face_index].F[q], dt, kx, U_safe);
						UpdateXLeftState2D(branch, e_left, q, F_ig, dt, kx, U_ig);
					}
					if (!StateAdmissible2D(U_safe, param.rho_floor, param.p_floor))
					{
						for (int m = 0; m < 4; ++m)
						{
							final_x_face_fluxes[face_index].F[q][m] = low_x_face_fluxes[face_index].F[q][m];
						}
						diag.theta_x[diag_index] = 0.0;
						diag.theta_rho_x[diag_index] = 0.0;
						diag.theta_p_x[diag_index] = 0.0;
						diag.limited_faces_x++;
						diag.min_theta_x = 0.0;
						continue;
					}
					if (StateAdmissible2D(U_ig, param.rho_floor, param.p_floor))
					{
						for (int m = 0; m < 4; ++m)
						{
							final_x_face_fluxes[face_index].F[q][m] = F_ig[m];
						}
						continue;
					}
					const double theta_rho = DensityTheta(U_safe[0], U_ig[0], param.rho_floor);
					diag.theta_rho_x[diag_index] = theta_rho;
					double F_rho[4];
					for (int m = 0; m < 4; ++m)
					{
						F_rho[m] = low_x_face_fluxes[face_index].F[q][m]
							+ theta_rho * (F_ig[m] - low_x_face_fluxes[face_index].F[q][m]);
					}
					double U_rho[4];
					if (face_i == 0)
					{
						UpdateXRightState2D(branch, e_right, q, F_rho, dt, kx, U_rho);
					}
					else
					{
						UpdateXLeftState2D(branch, e_left, q, F_rho, dt, kx, U_rho);
					}
					const double theta_p = PressureTheta2D(U_safe, U_rho, param.p_floor);
					diag.theta_p_x[diag_index] = theta_p;
					const double theta = theta_rho * theta_p;
					diag.theta_x[diag_index] = theta;
					diag.min_theta_x = std::min(diag.min_theta_x, theta);
					if (theta < 1.0 - 1.0e-12)
					{
						diag.limited_faces_x++;
					}
					for (int m = 0; m < 4; ++m)
					{
						final_x_face_fluxes[face_index].F[q][m] =
							low_x_face_fluxes[face_index].F[q][m]
							+ theta * (F_ig[m] - low_x_face_fluxes[face_index].F[q][m]);
					}
					continue;
				}

				double UL_safe[4], UR_safe[4], UL_ig[4], UR_ig[4];
				UpdateXLeftState2D(branch, e_left, q, low_x_face_fluxes[face_index].F[q], dt, kx, UL_safe);
				UpdateXRightState2D(branch, e_right, q, low_x_face_fluxes[face_index].F[q], dt, kx, UR_safe);
				UpdateXLeftState2D(branch, e_left, q, F_ig, dt, kx, UL_ig);
				UpdateXRightState2D(branch, e_right, q, F_ig, dt, kx, UR_ig);
				if (!StateAdmissible2D(UL_safe, param.rho_floor, param.p_floor) ||
					!StateAdmissible2D(UR_safe, param.rho_floor, param.p_floor))
				{
					for (int m = 0; m < 4; ++m)
					{
						final_x_face_fluxes[face_index].F[q][m] = low_x_face_fluxes[face_index].F[q][m];
					}
					diag.theta_x[diag_index] = 0.0;
					diag.theta_rho_x[diag_index] = 0.0;
					diag.theta_p_x[diag_index] = 0.0;
					diag.limited_faces_x++;
					diag.min_theta_x = 0.0;
					continue;
				}
				if (StateAdmissible2D(UL_ig, param.rho_floor, param.p_floor) &&
					StateAdmissible2D(UR_ig, param.rho_floor, param.p_floor))
				{
					for (int m = 0; m < 4; ++m)
					{
						final_x_face_fluxes[face_index].F[q][m] = F_ig[m];
					}
					continue;
				}
				const double theta_rho = std::min(
					1.0,
					std::min(
						DensityTheta(UL_safe[0], UL_ig[0], param.rho_floor),
						DensityTheta(UR_safe[0], UR_ig[0], param.rho_floor)));
				diag.theta_rho_x[diag_index] = theta_rho;
				double F_rho[4];
				for (int m = 0; m < 4; ++m)
				{
					F_rho[m] = low_x_face_fluxes[face_index].F[q][m]
						+ theta_rho * (F_ig[m] - low_x_face_fluxes[face_index].F[q][m]);
				}
				double UL_rho[4], UR_rho[4];
				UpdateXLeftState2D(branch, e_left, q, F_rho, dt, kx, UL_rho);
				UpdateXRightState2D(branch, e_right, q, F_rho, dt, kx, UR_rho);
				const double theta_p = std::min(
					1.0,
					std::min(
						PressureTheta2D(UL_safe, UL_rho, param.p_floor),
						PressureTheta2D(UR_safe, UR_rho, param.p_floor)));
				diag.theta_p_x[diag_index] = theta_p;
				const double theta = theta_rho * theta_p;
				diag.theta_x[diag_index] = theta;
				diag.min_theta_x = std::min(diag.min_theta_x, theta);
				if (theta < 1.0 - 1.0e-12)
				{
					diag.limited_faces_x++;
				}
				for (int m = 0; m < 4; ++m)
				{
					final_x_face_fluxes[face_index].F[q][m] =
						low_x_face_fluxes[face_index].F[q][m]
						+ theta * (F_ig[m] - low_x_face_fluxes[face_index].F[q][m]);
				}
			}
		}
	}

	for (int i = 0; i < nx; ++i)
	{
		for (int face_j = 0; face_j <= ny; ++face_j)
		{
			const bool physical_boundary_face = boundary != gksfr2d_periodic && (face_j == 0 || face_j == ny);
			const int bottom_j = (boundary == gksfr2d_periodic)
				? (face_j + ny - 1) % ny
				: std::max(0, face_j - 1);
			const int top_j = (boundary == gksfr2d_periodic)
				? face_j % ny
				: std::min(ny - 1, face_j);
			const int e_bottom = GKSSubcellCellIndex2D(branch, i, bottom_j);
			const int e_top = GKSSubcellCellIndex2D(branch, i, top_j);
			const double alpha_face = physical_boundary_face
				? (face_j == 0 ? alpha_cell[e_top] : alpha_cell[e_bottom])
				: 0.5 * (alpha_cell[e_bottom] + alpha_cell[e_top]);
			const int face_index = YFaceIndex2D(ny, i, face_j);
			for (int p = 0; p < 3; ++p)
			{
				const int diag_index = face_index * 3 + p;
				diag.alpha_face_y[diag_index] = alpha_face;
				double G_ig[4];
				BlendFlux4(high_y_face_fluxes[face_index].F[p], low_y_face_fluxes[face_index].F[p], alpha_face, G_ig);
				if (physical_boundary_face)
				{
					double U_safe[4], U_ig[4];
					if (face_j == 0)
					{
						UpdateYTopState2D(branch, e_top, p, low_y_face_fluxes[face_index].F[p], dt, ky, U_safe);
						UpdateYTopState2D(branch, e_top, p, G_ig, dt, ky, U_ig);
					}
					else
					{
						UpdateYBottomState2D(branch, e_bottom, p, low_y_face_fluxes[face_index].F[p], dt, ky, U_safe);
						UpdateYBottomState2D(branch, e_bottom, p, G_ig, dt, ky, U_ig);
					}
					if (!StateAdmissible2D(U_safe, param.rho_floor, param.p_floor))
					{
						for (int m = 0; m < 4; ++m)
						{
							final_y_face_fluxes[face_index].F[p][m] = low_y_face_fluxes[face_index].F[p][m];
						}
						diag.theta_y[diag_index] = 0.0;
						diag.theta_rho_y[diag_index] = 0.0;
						diag.theta_p_y[diag_index] = 0.0;
						diag.limited_faces_y++;
						diag.min_theta_y = 0.0;
						continue;
					}
					if (StateAdmissible2D(U_ig, param.rho_floor, param.p_floor))
					{
						for (int m = 0; m < 4; ++m)
						{
							final_y_face_fluxes[face_index].F[p][m] = G_ig[m];
						}
						continue;
					}
					const double theta_rho = DensityTheta(U_safe[0], U_ig[0], param.rho_floor);
					diag.theta_rho_y[diag_index] = theta_rho;
					double G_rho[4];
					for (int m = 0; m < 4; ++m)
					{
						G_rho[m] = low_y_face_fluxes[face_index].F[p][m]
							+ theta_rho * (G_ig[m] - low_y_face_fluxes[face_index].F[p][m]);
					}
					double U_rho[4];
					if (face_j == 0)
					{
						UpdateYTopState2D(branch, e_top, p, G_rho, dt, ky, U_rho);
					}
					else
					{
						UpdateYBottomState2D(branch, e_bottom, p, G_rho, dt, ky, U_rho);
					}
					const double theta_p = PressureTheta2D(U_safe, U_rho, param.p_floor);
					diag.theta_p_y[diag_index] = theta_p;
					const double theta = theta_rho * theta_p;
					diag.theta_y[diag_index] = theta;
					diag.min_theta_y = std::min(diag.min_theta_y, theta);
					if (theta < 1.0 - 1.0e-12)
					{
						diag.limited_faces_y++;
					}
					for (int m = 0; m < 4; ++m)
					{
						final_y_face_fluxes[face_index].F[p][m] =
							low_y_face_fluxes[face_index].F[p][m]
							+ theta * (G_ig[m] - low_y_face_fluxes[face_index].F[p][m]);
					}
					continue;
				}

				double UB_safe[4], UT_safe[4], UB_ig[4], UT_ig[4];
				UpdateYBottomState2D(branch, e_bottom, p, low_y_face_fluxes[face_index].F[p], dt, ky, UB_safe);
				UpdateYTopState2D(branch, e_top, p, low_y_face_fluxes[face_index].F[p], dt, ky, UT_safe);
				UpdateYBottomState2D(branch, e_bottom, p, G_ig, dt, ky, UB_ig);
				UpdateYTopState2D(branch, e_top, p, G_ig, dt, ky, UT_ig);
				if (!StateAdmissible2D(UB_safe, param.rho_floor, param.p_floor) ||
					!StateAdmissible2D(UT_safe, param.rho_floor, param.p_floor))
				{
					for (int m = 0; m < 4; ++m)
					{
						final_y_face_fluxes[face_index].F[p][m] = low_y_face_fluxes[face_index].F[p][m];
					}
					diag.theta_y[diag_index] = 0.0;
					diag.theta_rho_y[diag_index] = 0.0;
					diag.theta_p_y[diag_index] = 0.0;
					diag.limited_faces_y++;
					diag.min_theta_y = 0.0;
					continue;
				}
				if (StateAdmissible2D(UB_ig, param.rho_floor, param.p_floor) &&
					StateAdmissible2D(UT_ig, param.rho_floor, param.p_floor))
				{
					for (int m = 0; m < 4; ++m)
					{
						final_y_face_fluxes[face_index].F[p][m] = G_ig[m];
					}
					continue;
				}
				const double theta_rho = std::min(
					1.0,
					std::min(
						DensityTheta(UB_safe[0], UB_ig[0], param.rho_floor),
						DensityTheta(UT_safe[0], UT_ig[0], param.rho_floor)));
				diag.theta_rho_y[diag_index] = theta_rho;
				double G_rho[4];
				for (int m = 0; m < 4; ++m)
				{
					G_rho[m] = low_y_face_fluxes[face_index].F[p][m]
						+ theta_rho * (G_ig[m] - low_y_face_fluxes[face_index].F[p][m]);
				}
				double UB_rho[4], UT_rho[4];
				UpdateYBottomState2D(branch, e_bottom, p, G_rho, dt, ky, UB_rho);
				UpdateYTopState2D(branch, e_top, p, G_rho, dt, ky, UT_rho);
				const double theta_p = std::min(
					1.0,
					std::min(
						PressureTheta2D(UB_safe, UB_rho, param.p_floor),
						PressureTheta2D(UT_safe, UT_rho, param.p_floor)));
				diag.theta_p_y[diag_index] = theta_p;
				const double theta = theta_rho * theta_p;
				diag.theta_y[diag_index] = theta;
				diag.min_theta_y = std::min(diag.min_theta_y, theta);
				if (theta < 1.0 - 1.0e-12)
				{
					diag.limited_faces_y++;
				}
				for (int m = 0; m < 4; ++m)
				{
					final_y_face_fluxes[face_index].F[p][m] =
						low_y_face_fluxes[face_index].F[p][m]
						+ theta * (G_ig[m] - low_y_face_fluxes[face_index].F[p][m]);
				}
			}
		}
	}
}
