#pragma once

#include "gks_fr_high_order.h"

#include <vector>

void GKSFRAdapterComputeHighFaceFluxes(
	const GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	std::vector<GKSFRFaceFlux1D>& face_fluxes);

void GKSFRAdapterAdvanceWithFaceFluxes(
	const GKSFRMesh1D& mesh_old,
	double dt,
	const std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
	GKSFRMesh1D& mesh_new);

void GKSFRAdapterAssembleBlendedHighResiduals(
	const GKSFRMesh1D& mesh_old,
	double dt,
	const std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
	const std::vector<double>& alpha_cell,
	GKSFRMesh1D& mesh_residual);
