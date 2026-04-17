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

	const double N = 2.0;
	const double T = 0.5 * std::pow(10.0, -1.8 * std::pow(N + 1.0, 0.25));
	const double alpha_bar = 1.0 / (1.0 + std::exp(-(param.logistic_xi / T) * (E - T)));

	if (alpha_bar < param.alpha_min)
	{
		return 0.0;
	}
	if (alpha_bar > 1.0 - param.alpha_min)
	{
		return 1.0;
	}
	return alpha_bar;
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
