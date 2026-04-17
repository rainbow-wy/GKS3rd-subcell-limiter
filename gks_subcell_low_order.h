#pragma once

#include "gks_fr_high_order.h"
#include "gks_subcell_geometry.h"

#include <vector>

enum GKSSubcellLowOrderType
{
	KFVS1,
	MUSCL_HANCOCK
};

extern GKSSubcellLowOrderType low_order_type;

struct GKSSubcellCell1D
{
	// In the LWFR low-order branch, the current solution DOFs themselves are the
	// three low-order subcell unknowns associated with the three weighted
	// subcells of the big cell.
	double low_dof[3][3];
	double internal_flux[2][3];
	double residual_subcell[3][3];
	double low_dof_new[3][3];
	double cell_avg_new[3];
};

struct GKSSubcellBranch1D
{
	GKSSubcellGeom1D geom;
	std::vector<GKSSubcellCell1D> cell;
};

void GKSSubcellBranchResize1D(
	GKSSubcellBranch1D& branch,
	const GKSFRMesh1D& mesh);

void GKSSubcellInitializeFromCurrentDofs1D(
	const GKSFRMesh1D& mesh,
	GKSSubcellBranch1D& branch);

void GKSSubcellComputeInternalFluxes1D(
	GKSSubcellBranch1D& branch,
	double dt);

void GKSSubcellComputeLowFaceFluxes1D(
	const GKSSubcellBranch1D& branch,
	double dt,
	GKSFRBoundary1D boundary,
	std::vector<GKSFRFaceFlux1D>& face_fluxes);

void GKSSubcellComputeElementResiduals1D(
	GKSSubcellBranch1D& branch);

void GKSSubcellAddFaceResiduals1D(
	GKSSubcellBranch1D& branch,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes);

void GKSSubcellScaleResiduals1D(
	GKSSubcellBranch1D& branch,
	const std::vector<double>& alpha_cell);

void GKSSubcellUpdateLowOrderDofs1D(
	GKSSubcellBranch1D& branch,
	double dt);

void GKSSubcellAdvanceWithFaceFluxes1D(
	GKSSubcellBranch1D& branch,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes,
	double dt);
