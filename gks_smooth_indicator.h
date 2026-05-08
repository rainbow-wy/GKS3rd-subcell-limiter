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
	double eps_energy;
	double logistic_xi;
	double alpha_min;

	GKSSmoothIndicatorParam2D();
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

void GKSSmoothIndicatorAllCells2D(
	const GKSFRMesh2D& mesh,
	const GKSSmoothIndicatorParam2D& param,
	std::vector<double>& alpha_raw,
	std::vector<double>& alpha_final);
