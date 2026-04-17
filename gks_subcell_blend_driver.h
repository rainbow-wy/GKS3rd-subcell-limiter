#pragma once

#include "gks_flux_limiter.h"
#include "gks_scaling_limiter.h"
#include "gks_smooth_indicator.h"
#include "gks_subcell_low_order.h"

enum GKSSubcellBlendMode1D
{
	gks_subcell_pure_high,
	gks_subcell_pure_low,
	gks_subcell_hybrid
};

struct GKSSubcellFrameworkConfig1D
{
	GKSSubcellBlendMode1D blend_mode;
	GKSSubcellLowOrderType low_mode;
	bool use_flux_limiter;
	bool use_scaling_limiter;
	GKSSmoothIndicatorParam1D smooth_param;
	GKSFluxLimiterParam1D flux_param;
	GKSScalingLimiterParam1D scaling_param;

	GKSSubcellFrameworkConfig1D();
};

struct GKSSubcellFrameworkDiag1D
{
	std::vector<double> alpha_raw;
	std::vector<double> alpha_final;
	GKSFluxLimiterDiag1D flux_diag;
	GKSScalingLimiterDiag1D scaling_diag;
	double min_rho;
	double min_p;
	int min_rho_cell;
	int min_rho_point;
	int min_p_cell;
	int min_p_point;

	GKSSubcellFrameworkDiag1D();
};

void GKSSubcellAdvanceOneStep1D(
	GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	const GKSSubcellFrameworkConfig1D& config,
	GKSSubcellFrameworkDiag1D& diag);

void accuracy_sinwave_1d_gks_subcell();
void riemann_problem_1d_gks_subcell();
