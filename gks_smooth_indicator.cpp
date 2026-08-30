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
	: logistic_xi(9.21024),
	  alpha_min(1.0e-3)
{
}

GKSSmoothIndicatorCellDiag2D::GKSSmoothIndicatorCellDiag2D()
	: qhat{},
	  S0(0.0),
	  S1(0.0),
	  S2(0.0),
	  E1(0.0),
	  E2(0.0),
	  E(0.0),
	  alpha_raw(0.0)
{
}

GKSSmoothIndicatorFieldDiag2D::GKSSmoothIndicatorFieldDiag2D()
	: max_E1(0.0),
	  max_E2(0.0),
	  max_E(0.0),
	  max_alpha_raw(0.0),
	  max_alpha_final(0.0),
	  troubled_cells(0)
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

	double OrthonormalLegendre(int mode, double x)
	{
		if (mode == 0)
		{
			return std::sqrt(0.5);
		}
		if (mode == 1)
		{
			return std::sqrt(1.5) * x;
		}
		return std::sqrt(2.5) * Legendre2(x);
	}

	struct ModalTransform3
	{
		double coefficient[3][3];

		ModalTransform3()
		{
			const double node[3] = {
				-std::sqrt(3.0 / 5.0),
				0.0,
				std::sqrt(3.0 / 5.0)
			};
			const double weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
			for (int mode = 0; mode < 3; ++mode)
			{
				for (int nodal = 0; nodal < 3; ++nodal)
				{
					// At the three Gauss points this is A=V^{-1},
					// V(i,m)=phi_m(node_i), for orthonormal phi_m.
					coefficient[mode][nodal] = weight[nodal]
						* OrthonormalLegendre(mode, node[nodal]);
				}
			}
		}
	};

	const ModalTransform3& GetModalTransform3()
	{
		static const ModalTransform3 transform;
		return transform;
	}

	double Square(double value)
	{
		return value * value;
	}

	double AlphaFromModalEnergy2D(double E, double logistic_xi, double alpha_min)
	{
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

	bool ResolveFaceNeighbor(int index, int count, bool periodic, int& resolved)
	{
		if (periodic)
		{
			resolved = (index % count + count) % count;
			return true;
		}
		if (index < 0 || index >= count)
		{
			return false;
		}
		resolved = index;
		return true;
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

void GKSSmoothIndicatorCellDiagnostics2D(
	const double point_Q[3][3][4],
	const GKSSmoothIndicatorParam2D& param,
	GKSSmoothIndicatorCellDiag2D& diag)
{
	diag = GKSSmoothIndicatorCellDiag2D();
	double nodal[3][3];
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			nodal[i][j] = SensorValueRhoP2D(point_Q[i][j]);
		}
	}

	const ModalTransform3& transform = GetModalTransform3();
	double intermediate[3][3]{};
	for (int m = 0; m < 3; ++m)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int i = 0; i < 3; ++i)
			{
				intermediate[m][j] += transform.coefficient[m][i] * nodal[i][j];
			}
		}
	}
	for (int m = 0; m < 3; ++m)
	{
		for (int n = 0; n < 3; ++n)
		{
			for (int j = 0; j < 3; ++j)
			{
				diag.qhat[m][n] += transform.coefficient[n][j] * intermediate[m][j];
			}
		}
	}

	diag.S0 = Square(diag.qhat[0][0]);
	const double H1 = Square(diag.qhat[1][0])
		+ Square(diag.qhat[0][1])
		+ Square(diag.qhat[1][1]);
	diag.S1 = diag.S0 + H1;
	const double H2 = Square(diag.qhat[2][0])
		+ Square(diag.qhat[2][1])
		+ Square(diag.qhat[2][2])
		+ Square(diag.qhat[0][2])
		+ Square(diag.qhat[1][2]);
	diag.S2 = diag.S1 + H2;
	diag.E1 = diag.S1 > 0.0 ? H1 / diag.S1 : 0.0;
	diag.E2 = diag.S2 > 0.0 ? H2 / diag.S2 : 0.0;
	diag.E = std::max(diag.E1, diag.E2);
	diag.alpha_raw = AlphaFromModalEnergy2D(diag.E, param.logistic_xi, param.alpha_min);
}

double GKSSmoothIndicatorCell2D(
	const double point_Q[3][3][4],
	const GKSSmoothIndicatorParam2D& param)
{
	GKSSmoothIndicatorCellDiag2D diag;
	GKSSmoothIndicatorCellDiagnostics2D(point_Q, param, diag);
	return diag.alpha_raw;
}

void GKSSmoothIndicatorAllCells2D(
	const GKSFRMesh2D& mesh,
	const GKSSmoothIndicatorParam2D& param,
	GKSFRBoundary2D boundary,
	std::vector<double>& alpha_raw,
	std::vector<double>& alpha_final,
	GKSSmoothIndicatorFieldDiag2D* diagnostics)
{
	const int cells = mesh.cells_x * mesh.cells_y;
	alpha_raw.assign(cells, 0.0);
	alpha_final.assign(cells, 0.0);
	if (diagnostics != nullptr)
	{
		diagnostics->E1.assign(cells, 0.0);
		diagnostics->E2.assign(cells, 0.0);
		diagnostics->E.assign(cells, 0.0);
		diagnostics->max_E1 = 0.0;
		diagnostics->max_E2 = 0.0;
		diagnostics->max_E = 0.0;
		diagnostics->max_alpha_raw = 0.0;
		diagnostics->max_alpha_final = 0.0;
		diagnostics->troubled_cells = 0;
	}

	#pragma omp parallel for collapse(2)
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			const int e = GKSFR_CellIndex2D(mesh, i, j);
			GKSSmoothIndicatorCellDiag2D cell_diag;
			GKSSmoothIndicatorCellDiagnostics2D(mesh.cell[e].Q, param, cell_diag);
			alpha_raw[e] = cell_diag.alpha_raw;
			if (diagnostics != nullptr)
			{
				diagnostics->E1[e] = cell_diag.E1;
				diagnostics->E2[e] = cell_diag.E2;
				diagnostics->E[e] = cell_diag.E;
			}
		}
	}

	const bool periodic = boundary == gksfr2d_periodic;
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			double neighbor = 0.0;
			const int offset_i[4] = { -1, 1, 0, 0 };
			const int offset_j[4] = { 0, 0, -1, 1 };
			for (int face = 0; face < 4; ++face)
			{
				int neighbor_i = 0;
				int neighbor_j = 0;
				if (ResolveFaceNeighbor(i + offset_i[face], mesh.cells_x, periodic, neighbor_i)
					&& ResolveFaceNeighbor(j + offset_j[face], mesh.cells_y, periodic, neighbor_j))
				{
					neighbor = std::max(
						neighbor,
						0.5 * alpha_raw[GKSFR_CellIndex2D(mesh, neighbor_i, neighbor_j)]);
				}
			}
			const int e = GKSFR_CellIndex2D(mesh, i, j);
			alpha_final[e] = std::max(alpha_raw[e], neighbor);
		}
	}

	if (diagnostics != nullptr)
	{
		for (int e = 0; e < cells; ++e)
		{
			diagnostics->max_E1 = std::max(diagnostics->max_E1, diagnostics->E1[e]);
			diagnostics->max_E2 = std::max(diagnostics->max_E2, diagnostics->E2[e]);
			diagnostics->max_E = std::max(diagnostics->max_E, diagnostics->E[e]);
			diagnostics->max_alpha_raw = std::max(diagnostics->max_alpha_raw, alpha_raw[e]);
			diagnostics->max_alpha_final = std::max(diagnostics->max_alpha_final, alpha_final[e]);
			if (alpha_final[e] > 0.0)
			{
				++diagnostics->troubled_cells;
			}
		}
	}
}
