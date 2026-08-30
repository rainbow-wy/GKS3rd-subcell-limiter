#pragma once

#include "gks_fr_high_order.h"

#include <vector>

struct GKSSmoothIndicatorParam1D
{
	double eps_energy;
	double logistic_xi;
	double alpha_min;

	GKSSmoothIndicatorParam1D();
};

struct GKSSmoothIndicatorParam2D
{
	double logistic_xi;
	double alpha_min;

	GKSSmoothIndicatorParam2D();
};

struct GKSSmoothIndicatorCellDiag2D
{
	double qhat[3][3];
	double S0;
	double S1;
	double S2;
	double E1;
	double E2;
	double E;
	double alpha_raw;

	GKSSmoothIndicatorCellDiag2D();
};

struct GKSSmoothIndicatorFieldDiag2D
{
	std::vector<double> E1;
	std::vector<double> E2;
	std::vector<double> E;
	double max_E1;
	double max_E2;
	double max_E;
	double max_alpha_raw;
	double max_alpha_final;
	int troubled_cells;

	GKSSmoothIndicatorFieldDiag2D();
};

double GKSSmoothIndicatorCell1D(
	const double point_Q[3][3],
	const GKSSmoothIndicatorParam1D& param);

void GKSSmoothIndicatorAllCells1D(
	const GKSFRMesh1D& mesh,
	const GKSSmoothIndicatorParam1D& param,
	std::vector<double>& alpha_raw,
	std::vector<double>& alpha_final);

double GKSSmoothIndicatorCell2D(
	const double point_Q[3][3][4],
	const GKSSmoothIndicatorParam2D& param);

void GKSSmoothIndicatorCellDiagnostics2D(
	const double point_Q[3][3][4],
	const GKSSmoothIndicatorParam2D& param,
	GKSSmoothIndicatorCellDiag2D& diag);

void GKSSmoothIndicatorAllCells2D(
	const GKSFRMesh2D& mesh,
	const GKSSmoothIndicatorParam2D& param,
	GKSFRBoundary2D boundary,
	std::vector<double>& alpha_raw,
	std::vector<double>& alpha_final,
	GKSSmoothIndicatorFieldDiag2D* diagnostics = nullptr);
