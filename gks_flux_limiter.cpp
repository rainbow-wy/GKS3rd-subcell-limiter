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
