#include "gks_scaling_limiter.h"

#include "gks_basic.h"
#include "gks_subcell_geometry.h"

#include <algorithm>
#include <cmath>

GKSScalingLimiterParam1D::GKSScalingLimiterParam1D()
	: rho_floor(1.0e-12),
	  p_floor(1.0e-12)
{
}

GKSScalingLimiterParam2D::GKSScalingLimiterParam2D()
	: rho_floor(1.0e-12),
	  p_floor(1.0e-12)
{
}

GKSScalingLimiterDiag2D::GKSScalingLimiterDiag2D()
	: limited_cells(0)
{
}

namespace
{
	bool StateAdmissibleFromAverageBlend(
		const double Uavg[3],
		const double Upoint[3],
		double theta,
		double rho_floor,
		double p_floor)
	{
		double U[3];
		for (int m = 0; m < 3; ++m)
		{
			U[m] = Uavg[m] + theta * (Upoint[m] - Uavg[m]);
		}
		if (!std::isfinite(U[0]) || !std::isfinite(U[1]) || !std::isfinite(U[2]))
		{
			return false;
		}
		if (U[0] < rho_floor)
		{
			return false;
		}
		return Pressure(U[0], U[1], U[2]) >= p_floor;
	}

	double PressureThetaFromAverage(
		const double Uavg[3],
		const double Upoint[3],
		double p_floor)
	{
		if (Pressure(Upoint[0], Upoint[1], Upoint[2]) >= p_floor)
		{
			return 1.0;
		}

		const double eps_bar = p_floor / (Gamma - 1.0);
		const double dr = Upoint[0] - Uavg[0];
		const double dm = Upoint[1] - Uavg[1];
		const double dE = Upoint[2] - Uavg[2];

		const double A = 2.0 * dr * dE - dm * dm;
		const double B = 2.0 * Uavg[0] * dE
			+ 2.0 * dr * (Uavg[2] - eps_bar)
			- 2.0 * Uavg[1] * dm;
		const double C = 2.0 * Uavg[0] * (Uavg[2] - eps_bar) - Uavg[1] * Uavg[1];

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

	double Pressure2D(const double U[4])
	{
		if (U[0] <= 0.0)
		{
			return -1.0;
		}
		return (Gamma - 1.0) *
			(U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
	}

	bool StateAdmissibleFromAverageBlend2D(
		const double Uavg[4],
		const double Upoint[4],
		double theta,
		double rho_floor,
		double p_floor)
	{
		double U[4];
		for (int m = 0; m < 4; ++m)
		{
			U[m] = Uavg[m] + theta * (Upoint[m] - Uavg[m]);
		}
		for (int m = 0; m < 4; ++m)
		{
			if (!std::isfinite(U[m]))
			{
				return false;
			}
		}
		return U[0] >= rho_floor && Pressure2D(U) >= p_floor;
	}

	double PressureThetaFromAverage2D(
		const double Uavg[4],
		const double Upoint[4],
		double p_floor)
	{
		if (Pressure2D(Upoint) >= p_floor)
		{
			return 1.0;
		}
		const double eta_p = p_floor / (Gamma - 1.0);
		const double dr = Upoint[0] - Uavg[0];
		const double dmx = Upoint[1] - Uavg[1];
		const double dmy = Upoint[2] - Uavg[2];
		const double dE = Upoint[3] - Uavg[3];
		const double A = 2.0 * dr * dE - dmx * dmx - dmy * dmy;
		const double B = 2.0 * Uavg[0] * dE
			+ 2.0 * dr * (Uavg[3] - eta_p)
			- 2.0 * Uavg[1] * dmx
			- 2.0 * Uavg[2] * dmy;
		const double C = 2.0 * Uavg[0] * (Uavg[3] - eta_p)
			- Uavg[1] * Uavg[1]
			- Uavg[2] * Uavg[2];
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

	double EnforceAdmissibleThetaByPostcheck2D(
		const double Uavg[4],
		const double Upoint[4],
		double theta_hi,
		double rho_floor,
		double p_floor)
	{
		if (theta_hi <= 0.0)
		{
			return 0.0;
		}
		if (StateAdmissibleFromAverageBlend2D(Uavg, Upoint, theta_hi, rho_floor, p_floor))
		{
			return theta_hi;
		}
		double lo = 0.0;
		double hi = theta_hi;
		for (int iter = 0; iter < 80; ++iter)
		{
			const double mid = 0.5 * (lo + hi);
			if (StateAdmissibleFromAverageBlend2D(Uavg, Upoint, mid, rho_floor, p_floor))
			{
				lo = mid;
			}
			else
			{
				hi = mid;
			}
		}
		return lo;
	}

	double EnforceAdmissibleThetaByPostcheck(
		const double Uavg[3],
		const double Upoint[3],
		double theta_hi,
		double rho_floor,
		double p_floor)
	{
		if (theta_hi <= 0.0)
		{
			return 0.0;
		}
		if (StateAdmissibleFromAverageBlend(Uavg, Upoint, theta_hi, rho_floor, p_floor))
		{
			return theta_hi;
		}

		double lo = 0.0;
		double hi = theta_hi;
		for (int iter = 0; iter < 80; ++iter)
		{
			const double mid = 0.5 * (lo + hi);
			if (StateAdmissibleFromAverageBlend(Uavg, Upoint, mid, rho_floor, p_floor))
			{
				lo = mid;
			}
			else
			{
				hi = mid;
			}
		}
		return lo;
	}
}

void GKSScalingLimiterApply2D(
	GKSFRMesh2D& mesh,
	const GKSScalingLimiterParam2D& param,
	GKSScalingLimiterDiag2D& diag)
{
	const int cells = mesh.cells_x * mesh.cells_y;
	diag.theta_rho.assign(cells, 1.0);
	diag.theta_p.assign(cells, 1.0);
	diag.theta.assign(cells, 1.0);
	diag.limited_cells = 0;

	for (int e = 0; e < cells; ++e)
	{
		double Uavg[4];
		GKSSubcellBigCellAverageFromPoints2D(mesh.cell[e].Q, Uavg);

		double theta_rho = 1.0;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				const double rho_i = mesh.cell[e].Q[i][j][0];
				if (rho_i < param.rho_floor)
				{
					const double denom = Uavg[0] - rho_i;
					if (denom <= 1.0e-14)
					{
						theta_rho = 0.0;
					}
					else
					{
						theta_rho = std::min(theta_rho,
							std::max(0.0, std::min(1.0, (Uavg[0] - param.rho_floor) / denom)));
					}
				}
			}
		}
		diag.theta_rho[e] = theta_rho;

		double rho_limited[3][3][4];
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					rho_limited[i][j][m] = Uavg[m]
						+ theta_rho * (mesh.cell[e].Q[i][j][m] - Uavg[m]);
				}
			}
		}

		double theta_p = 1.0;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				theta_p = std::min(theta_p,
					PressureThetaFromAverage2D(Uavg, rho_limited[i][j], param.p_floor));
			}
		}
		diag.theta_p[e] = theta_p;

		double theta = theta_rho * theta_p;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				theta = std::min(theta,
					EnforceAdmissibleThetaByPostcheck2D(
						Uavg,
						mesh.cell[e].Q[i][j],
						theta,
						param.rho_floor,
						param.p_floor));
			}
		}
		diag.theta[e] = theta;
		if (theta < 1.0 - 1.0e-12)
		{
			diag.limited_cells++;
		}
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					mesh.cell[e].Q[i][j][m] = Uavg[m]
						+ theta * (mesh.cell[e].Q[i][j][m] - Uavg[m]);
				}
			}
		}
	}
}

void GKSScalingLimiterApply1D(
	GKSFRMesh1D& mesh,
	const std::vector<GKSCellAverage1D>& safe_cell_average,
	const GKSScalingLimiterParam1D& param,
	GKSScalingLimiterDiag1D& diag)
{
	diag.theta_rho.assign(mesh.cells, 1.0);
	diag.theta_p.assign(mesh.cells, 1.0);
	diag.theta.assign(mesh.cells, 1.0);

	for (int e = 0; e < mesh.cells; ++e)
	{
		double theta_rho = 1.0;
		for (int i = 0; i < 3; ++i)
		{
			const double rho_avg = safe_cell_average[e].U[0];
			const double rho_i = mesh.cell[e].Q[i][0];
			if (rho_i < param.rho_floor)
			{
				const double denom = rho_avg - rho_i;
				if (denom <= 1.0e-14)
				{
					theta_rho = 0.0;
				}
				else
				{
					const double theta_i = (rho_avg - param.rho_floor) / denom;
					theta_rho = std::min(theta_rho, std::max(0.0, std::min(1.0, theta_i)));
				}
			}
		}
		diag.theta_rho[e] = theta_rho;

		double rho_limited[3][3];
		for (int i = 0; i < 3; ++i)
		{
			for (int m = 0; m < 3; ++m)
			{
				if (theta_rho <= 1.0e-14)
				{
					rho_limited[i][m] = safe_cell_average[e].U[m];
				}
				else
				{
					rho_limited[i][m] = safe_cell_average[e].U[m]
						+ theta_rho * (mesh.cell[e].Q[i][m] - safe_cell_average[e].U[m]);
				}
			}
		}

		double theta_p = 1.0;
		for (int i = 0; i < 3; ++i)
		{
			theta_p = std::min(
				theta_p,
				PressureThetaFromAverage(
					safe_cell_average[e].U,
					rho_limited[i],
					param.p_floor));
		}
		diag.theta_p[e] = theta_p;

		double theta = theta_rho * theta_p;
		for (int i = 0; i < 3; ++i)
		{
			theta = std::min(
				theta,
				EnforceAdmissibleThetaByPostcheck(
					safe_cell_average[e].U,
					mesh.cell[e].Q[i],
					theta,
					param.rho_floor,
					param.p_floor));
		}
		diag.theta[e] = theta;
		for (int i = 0; i < 3; ++i)
		{
			for (int m = 0; m < 3; ++m)
			{
				if (theta <= 1.0e-14)
				{
					mesh.cell[e].Q[i][m] = safe_cell_average[e].U[m];
				}
				else
				{
					mesh.cell[e].Q[i][m] = safe_cell_average[e].U[m]
						+ theta * (mesh.cell[e].Q[i][m] - safe_cell_average[e].U[m]);
				}
			}
		}
	}
}
