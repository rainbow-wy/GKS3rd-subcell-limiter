#include "gks_smooth_indicator.h"

#include "gks_basic.h"

#include <algorithm>
#include <cmath>

GKSSmoothIndicatorParam1D::GKSSmoothIndicatorParam1D()
	: eps_energy(1.0e-12),
	  logistic_xi(9.21024),
	  alpha_min(1.0e-4)
{
}

GKSSmoothIndicatorParam2D::GKSSmoothIndicatorParam2D()
	: eps_energy(1.0e-12),
	  logistic_xi(9.21024),
	  alpha_min(1.0e-3)
{
}

namespace
{
	double Legendre0(double)
	{
		return 1.0;
	}

	double Legendre1(double s)
	{
		return s;
	}

	double Legendre2(double s)
	{
		return 0.5 * (3.0 * s * s - 1.0);
	}

	double SensorValueRhoP(const double Q[3])
	{
		double prim[3];
		Convar_to_primvar_1D(prim, const_cast<double*>(Q));
		return prim[0] * prim[2];
	}

	double SensorValueRhoP2D(const double Q[4])
	{
		double prim[4];
		double local_Q[4];
		for (int m = 0; m < 4; ++m)
		{
			local_Q[m] = Q[m];
		}
		Convar_to_primvar_2D(prim, local_Q);
		return prim[0] * prim[3];
	}

	double LegendreValue(int n, double s)
	{
		if (n == 0)
		{
			return Legendre0(s);
		}
		if (n == 1)
		{
			return Legendre1(s);
		}
		return Legendre2(s);
	}

	double AlphaFromModalEnergy(double E, double eps_energy, double logistic_xi, double alpha_min)
	{
		(void)eps_energy;
		const double N = 2.0;
		const double T = 0.5 * std::pow(10.0, -1.8 * std::pow(N + 1.0, 0.25));
		const double alpha_bar = 1.0 / (1.0 + std::exp(-(logistic_xi / T) * (E - T)));
		if (alpha_bar < alpha_min)
		{
			return 0.0;
		}
		if (alpha_bar > 1.0 - alpha_min)
		{
			return 1.0;
		}
		return alpha_bar;
	}
}

double GKSSmoothIndicatorCell1D(
	const double point_Q[3][3],
	const GKSSmoothIndicatorParam1D& param)
{
	const double s[3] = {
		-std::sqrt(3.0 / 5.0),
		0.0,
		std::sqrt(3.0 / 5.0)
	};
	const double w[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

	const double q0 = SensorValueRhoP(point_Q[0]);
	const double q1 = SensorValueRhoP(point_Q[1]);
	const double q2 = SensorValueRhoP(point_Q[2]);
	const double q[3] = { q0, q1, q2 };

	double qhat[3]{};
	for (int j = 0; j < 3; ++j)
	{
		qhat[0] += 0.5 * (2.0 * 0 + 1.0) * q[j] * Legendre0(s[j]) * w[j];
		qhat[1] += 0.5 * (2.0 * 1 + 1.0) * q[j] * Legendre1(s[j]) * w[j];
		qhat[2] += 0.5 * (2.0 * 2 + 1.0) * q[j] * Legendre2(s[j]) * w[j];
	}

	double energy_total = param.eps_energy;
	for (int n = 0; n < 3; ++n)
	{
		energy_total += qhat[n] * qhat[n];
	}
	const double energy_high = qhat[2] * qhat[2];
	const double E = energy_high / energy_total;

	return AlphaFromModalEnergy(E, param.eps_energy, param.logistic_xi, param.alpha_min);
}

void GKSSmoothIndicatorAllCells1D(
	const GKSFRMesh1D& mesh,
	const GKSSmoothIndicatorParam1D& param,
	std::vector<double>& alpha_raw,
	std::vector<double>& alpha_final)
{
	alpha_raw.assign(mesh.cells, 0.0);
	alpha_final.assign(mesh.cells, 0.0);

	for (int e = 0; e < mesh.cells; ++e)
	{
		alpha_raw[e] = GKSSmoothIndicatorCell1D(mesh.cell[e].Q, param);
	}

	for (int e = 0; e < mesh.cells; ++e)
	{
		double neighbor_smooth = 0.0;
		if (e > 0)
		{
			neighbor_smooth = std::max(neighbor_smooth, 0.5 * alpha_raw[e - 1]);
		}
		if (e + 1 < mesh.cells)
		{
			neighbor_smooth = std::max(neighbor_smooth, 0.5 * alpha_raw[e + 1]);
		}
		alpha_final[e] = std::max(alpha_raw[e], neighbor_smooth);
	}
}

double GKSSmoothIndicatorCell2D(
	const double point_Q[3][3][4],
	const GKSSmoothIndicatorParam2D& param)
{
	const double s[3] = {
		-std::sqrt(3.0 / 5.0),
		0.0,
		std::sqrt(3.0 / 5.0)
	};
	const double w[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

	double qhat[3][3]{};
	for (int m = 0; m < 3; ++m)
	{
		for (int n = 0; n < 3; ++n)
		{
			const double factor = 0.25 * (2.0 * m + 1.0) * (2.0 * n + 1.0);
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					qhat[m][n] += factor * w[i] * w[j]
						* SensorValueRhoP2D(point_Q[i][j])
						* LegendreValue(m, s[i])
						* LegendreValue(n, s[j]);
				}
			}
		}
	}

	double energy_total = param.eps_energy;
	for (int m = 0; m < 3; ++m)
	{
		for (int n = 0; n < 3; ++n)
		{
			energy_total += qhat[m][n] * qhat[m][n];
		}
	}
	double energy_x = 0.0;
	double energy_y = 0.0;
	for (int n = 0; n < 3; ++n)
	{
		energy_x += qhat[2][n] * qhat[2][n];
	}
	for (int m = 0; m < 3; ++m)
	{
		energy_y += qhat[m][2] * qhat[m][2];
	}
	const double E = std::max(energy_x, energy_y) / energy_total;
	return AlphaFromModalEnergy(E, param.eps_energy, param.logistic_xi, param.alpha_min);
}

void GKSSmoothIndicatorAllCells2D(
	const GKSFRMesh2D& mesh,
	const GKSSmoothIndicatorParam2D& param,
	std::vector<double>& alpha_raw,
	std::vector<double>& alpha_final)
{
	const int cells = mesh.cells_x * mesh.cells_y;
	alpha_raw.assign(cells, 0.0);
	alpha_final.assign(cells, 0.0);
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			const int e = GKSFR_CellIndex2D(mesh, i, j);
			alpha_raw[e] = GKSSmoothIndicatorCell2D(mesh.cell[e].Q, param);
		}
	}
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			double neighbor = 0.0;
			if (i > 0)
			{
				neighbor = std::max(neighbor, 0.5 * alpha_raw[GKSFR_CellIndex2D(mesh, i - 1, j)]);
			}
			if (i + 1 < mesh.cells_x)
			{
				neighbor = std::max(neighbor, 0.5 * alpha_raw[GKSFR_CellIndex2D(mesh, i + 1, j)]);
			}
			if (j > 0)
			{
				neighbor = std::max(neighbor, 0.5 * alpha_raw[GKSFR_CellIndex2D(mesh, i, j - 1)]);
			}
			if (j + 1 < mesh.cells_y)
			{
				neighbor = std::max(neighbor, 0.5 * alpha_raw[GKSFR_CellIndex2D(mesh, i, j + 1)]);
			}
			const int e = GKSFR_CellIndex2D(mesh, i, j);
			alpha_final[e] = std::max(alpha_raw[e], neighbor);
		}
	}
}
