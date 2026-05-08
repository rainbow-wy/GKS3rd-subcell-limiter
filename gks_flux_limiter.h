#pragma once

#include "gks_fr_high_order.h"
#include "gks_subcell_low_order.h"

#include <vector>

struct GKSFluxLimiterParam1D
{
	double rho_floor;
	double p_floor;

	GKSFluxLimiterParam1D();
};

struct GKSFluxLimiterParam2D
{
	double rho_floor;
	double p_floor;
	double kx;
	double ky;

	GKSFluxLimiterParam2D();
};

struct GKSFluxLimiterDiag1D
{
	std::vector<double> alpha_face;
	std::vector<double> theta_rho;
	std::vector<double> theta_p;
	std::vector<double> theta;
	std::vector<double> left_safe_rho;
	std::vector<double> left_safe_u;
	std::vector<double> left_safe_p;
	std::vector<double> right_safe_rho;
	std::vector<double> right_safe_u;
	std::vector<double> right_safe_p;
	std::vector<double> left_candidate_rho;
	std::vector<double> left_candidate_u;
	std::vector<double> left_candidate_p;
	std::vector<double> right_candidate_rho;
	std::vector<double> right_candidate_u;
	std::vector<double> right_candidate_p;
	std::vector<int> left_candidate_rho_bad;
	std::vector<int> left_candidate_p_bad;
	std::vector<int> right_candidate_rho_bad;
	std::vector<int> right_candidate_p_bad;
	std::vector<int> left_safe_rho_bad;
	std::vector<int> left_safe_p_bad;
	std::vector<int> right_safe_rho_bad;
	std::vector<int> right_safe_p_bad;
};

struct GKSFluxLimiterDiag2D
{
	std::vector<double> alpha_face_x;
	std::vector<double> alpha_face_y;
	std::vector<double> theta_rho_x;
	std::vector<double> theta_rho_y;
	std::vector<double> theta_p_x;
	std::vector<double> theta_p_y;
	std::vector<double> theta_x;
	std::vector<double> theta_y;
	int limited_faces_x;
	int limited_faces_y;
	double min_theta_x;
	double min_theta_y;

	GKSFluxLimiterDiag2D();
};

void GKSFluxLimiterApply1D(
	const GKSSubcellBranch1D& branch,
	const std::vector<double>& alpha_cell,
	const std::vector<GKSFRFaceFlux1D>& high_face_fluxes,
	const std::vector<GKSFRFaceFlux1D>& low_face_fluxes,
	double dt,
	GKSFRBoundary1D boundary,
	const GKSFluxLimiterParam1D& param,
	std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
	GKSFluxLimiterDiag1D& diag);

void GKSFluxLimiterApply2D(
	const GKSSubcellBranch2D& branch,
	const std::vector<double>& alpha_cell,
	const std::vector<GKSFRFaceFlux2D>& high_x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& high_y_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& low_x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& low_y_face_fluxes,
	double dt,
	GKSFRBoundary2D boundary,
	const GKSFluxLimiterParam2D& param,
	std::vector<GKSFRFaceFlux2D>& final_x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& final_y_face_fluxes,
	GKSFluxLimiterDiag2D& diag);
