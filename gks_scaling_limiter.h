#pragma once

#include "gks_fr_high_order.h"

#include <vector>

struct GKSCellAverage1D
{
	double U[3];
};

struct GKSCellAverage2D
{
	double U[4];
};

struct GKSScalingLimiterParam1D
{
	double rho_floor;
	double p_floor;

	GKSScalingLimiterParam1D();
};

struct GKSScalingLimiterParam2D
{
	double rho_floor;
	double p_floor;

	GKSScalingLimiterParam2D();
};

struct GKSScalingLimiterDiag1D
{
	std::vector<double> theta_rho;
	std::vector<double> theta_p;
	std::vector<double> theta;
};

struct GKSScalingLimiterDiag2D
{
	std::vector<double> theta_rho;
	std::vector<double> theta_p;
	std::vector<double> theta;
	int limited_cells;

	GKSScalingLimiterDiag2D();
};

void GKSScalingLimiterApply1D(
	GKSFRMesh1D& mesh,
	const std::vector<GKSCellAverage1D>& safe_cell_average,
	const GKSScalingLimiterParam1D& param,
	GKSScalingLimiterDiag1D& diag);

void GKSScalingLimiterApply2D(
	GKSFRMesh2D& mesh,
	const GKSScalingLimiterParam2D& param,
	GKSScalingLimiterDiag2D& diag);
