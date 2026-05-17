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

enum GKSSubcellBlendMode2D
{
	gks_subcell2d_pure_high,
	gks_subcell2d_pure_low,
	gks_subcell2d_hybrid
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

struct GKSSubcellFrameworkConfig2D
{
	GKSSubcellBlendMode2D blend_mode;
	GKSSubcellLowOrderType low_mode;
	bool use_flux_limiter;
	bool use_scaling_limiter;
	GKSSmoothIndicatorParam2D smooth_param;
	GKSFluxLimiterParam2D flux_param;
	GKSScalingLimiterParam2D scaling_param;

	GKSSubcellFrameworkConfig2D();
};

struct GKSSubcellFrameworkDiag1D
{
	std::vector<double> alpha_raw;
	std::vector<double> alpha_final;
	GKSFluxLimiterDiag1D flux_diag;
	GKSScalingLimiterDiag1D scaling_diag;
	GKSSubcellMUSCLHancockStats1D muscl_stats;
	double min_rho;
	double min_p;
	int min_rho_cell;
	int min_rho_point;
	int min_p_cell;
	int min_p_point;

	GKSSubcellFrameworkDiag1D();
};

struct GKSSubcellFrameworkDiag2D
{
	std::vector<double> alpha_raw;
	std::vector<double> alpha_final;
	GKSFluxLimiterDiag2D flux_diag;
	GKSScalingLimiterDiag2D scaling_diag;
	GKSSubcellMUSCLHancockStats2D muscl_stats;
	double min_rho;
	double min_p;
	double max_alpha;
	int troubled_cells;

	GKSSubcellFrameworkDiag2D();
};

struct GKSSubcellMask2D
{
	bool enabled;
	int cells_x;
	int cells_y;
	std::vector<int> active;

	GKSSubcellMask2D();
};

void GKSSubcellAdvanceOneStep1D(
	GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	const GKSSubcellFrameworkConfig1D& config,
	GKSSubcellFrameworkDiag1D& diag);

void GKSSubcellAdvanceOneStep2D(
	GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	const GKSSubcellFrameworkConfig2D& config,
	GKSSubcellFrameworkDiag2D& diag);

void GKSSubcellBuildRectangularCutoutMask2D(
	const GKSFRMesh2D& mesh,
	double x_min,
	double x_max,
	double y_min,
	double y_max,
	GKSSubcellMask2D& mask);

void GKSSubcellAdvanceOneStep2DMasked(
	GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	const GKSSubcellMask2D& mask,
	const GKSSubcellFrameworkConfig2D& config,
	GKSSubcellFrameworkDiag2D& diag);

void accuracy_sinwave_1d_gks_subcell();
void riemann_problem_1d_gks_subcell();
void accuracy_sinwave_2d_gks_subcell();
void riemann_problem_2d_gks_subcell();
void double_mach_reflection_2d_gks_subcell();
void detonation_shock_diffraction_2d_gks_subcell();
