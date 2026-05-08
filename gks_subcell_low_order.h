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

struct GKSSubcellMUSCLHancockStats1D
{
	int reconstructed_cells;
	int zero_slope_fallback_cells;
	int predicted_state_fallback_cells;
	int face_flux_fallback_faces;
	int final_update_fallback_cells;
	double min_rho_update;
	double min_p_update;

	GKSSubcellMUSCLHancockStats1D();
};

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
	GKSSubcellMUSCLHancockStats1D muscl_stats;
};

struct GKSSubcellCell2D
{
	double low_dof[3][3][4];
	double internal_flux_x[2][3][4];
	double internal_flux_y[3][2][4];
	double residual_subcell[3][3][4];
	double low_dof_new[3][3][4];
	double cell_avg_new[4];
};

struct GKSSubcellBranch2D
{
	GKSSubcellGeom2D geom;
	int cells_x;
	int cells_y;
	std::vector<GKSSubcellCell2D> cell;
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

void GKSSubcellComputeInternalFluxes1D(
	GKSSubcellBranch1D& branch,
	double dt,
	GKSSubcellLowOrderType mode);

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

void GKSSubcellFallbackBadMUSCLUpdates1D(
	GKSSubcellBranch1D& branch,
	const GKSSubcellBranch1D& kfvs_fallback_branch,
	double dt);

void GKSSubcellAdvanceWithFaceFluxes1D(
	GKSSubcellBranch1D& branch,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes,
	double dt);

void GKSSubcellBranchResize2D(
	GKSSubcellBranch2D& branch,
	const GKSFRMesh2D& mesh);

int GKSSubcellCellIndex2D(const GKSSubcellBranch2D& branch, int i, int j);

void GKSSubcellInitializeFromCurrentDofs2D(
	const GKSFRMesh2D& mesh,
	GKSSubcellBranch2D& branch);

void GKSSubcellComputeInternalFluxes2D(
	GKSSubcellBranch2D& branch,
	double dt);

void GKSSubcellComputeLowFaceFluxes2D(
	const GKSSubcellBranch2D& branch,
	double dt,
	GKSFRBoundary2D boundary,
	std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& y_face_fluxes);

void GKSSubcellComputeElementResiduals2D(
	GKSSubcellBranch2D& branch);

void GKSSubcellAddFaceResiduals2D(
	GKSSubcellBranch2D& branch,
	const std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& y_face_fluxes);

void GKSSubcellUpdateLowOrderDofs2D(
	GKSSubcellBranch2D& branch,
	double dt);

void GKSSubcellAdvanceWithFaceFluxes2D(
	GKSSubcellBranch2D& branch,
	const std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& y_face_fluxes,
	double dt);
