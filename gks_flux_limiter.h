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
