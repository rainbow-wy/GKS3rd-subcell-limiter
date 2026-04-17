#pragma once

#include "gks_fr_high_order.h"

#include <vector>

struct GKSCellAverage1D
{
	double U[3];
};

struct GKSScalingLimiterParam1D
{
	double rho_floor;
	double p_floor;

	GKSScalingLimiterParam1D();
};

struct GKSScalingLimiterDiag1D
{
	std::vector<double> theta_rho;
	std::vector<double> theta_p;
	std::vector<double> theta;
};

void GKSScalingLimiterApply1D(
	GKSFRMesh1D& mesh,
	const std::vector<GKSCellAverage1D>& safe_cell_average,
	const GKSScalingLimiterParam1D& param,
	GKSScalingLimiterDiag1D& diag);
