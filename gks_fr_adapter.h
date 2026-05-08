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

void GKSFRAdapterComputeHighFaceFluxes2D(
	const GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& y_face_fluxes);

void GKSFRAdapterAdvanceWithFaceFluxes2D(
	const GKSFRMesh2D& mesh_old,
	double dt,
	const std::vector<GKSFRFaceFlux2D>& final_x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& final_y_face_fluxes,
	GKSFRMesh2D& mesh_new);
