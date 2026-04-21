#pragma once

#include "fvm_gks1d.h"

#include <vector>

enum GKSFRBoundary1D
{
	gksfr_periodic,
	gksfr_free,
	gksfr_transmissive_strict
};

struct GKSFRFaceFlux1D
{
	double F[3];
};

struct GKSFRCell1D
{
	// High-order branch unknowns: three GL solution points per element.
	double Q[3][3];

	// Cell-center quadratic expansion recovered from the three GL points.
	double Qc[3];
	double Qx[3];
	double Qxx[3];

	// Maxwellian state and smooth collision time used by the cell-internal GKS.
	double prim[3];
	double tau;

	// Cell-internal GKS coefficients in the notation of GKS-FR.pdf.
	double a[3];
	double A[3];
	double b[3];
	double c[3];
	double d[3];

	// Time-averaged internal flux polynomial coefficients:
	// Fbar_int = Fe0 + Fe1*x_tilde + 0.5*Fe2*x_tilde^2.
	double Fe0[3];
	double Fe1[3];
	double Fe2[3];

	// Time-averaged flux values at the three GL solution points.
	double flux_gl[3][3];

	// Local discontinuous flux extrapolated to left/right faces.
	double flux_face_local[2][3];

	// FR residuals at the three GL solution points.
	double residual[3][3];

	GKSFRCell1D();
};

struct GKSFRMesh1D
{
	int cells;
	double x_left;
	double x_right;
	double dx;
	std::vector<GKSFRCell1D> cell;

	GKSFRMesh1D();
};

double GKSFR_GL_Point(int point_id);
void GKSFR_ResizeUniformMesh(GKSFRMesh1D& mesh, int cells, double x_left, double x_right);

void GKSFR_ComputeCellCenterData(GKSFRCell1D& cell, double h);
void GKSFR_ComputeInternalGKSCoeffs(GKSFRCell1D& cell);
void GKSFR_ComputeInternalTimeAveragedFluxAtGLPoints(GKSFRCell1D& cell, double h, double dt);
void GKSFR_ExtrapolateLocalFluxToFaces(GKSFRCell1D& cell);
void GKSFR_PrepareCell(GKSFRCell1D& cell, double h, double dt);
void GKSFR_PrepareCells(GKSFRMesh1D& mesh, double dt);

void GKSFR_ComputeInterfaceCommonFlux(
	const GKSFRCell1D& left_cell,
	const GKSFRCell1D& right_cell,
	double h,
	double dt,
	double common_flux[3]);

void GKSFR_ComputeCommonInterfaceFluxes(
	const GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	std::vector<GKSFRFaceFlux1D>& face_fluxes);

void GKSFR_ComputeResiduals(
	GKSFRMesh1D& mesh,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes);

void GKSFR_UpdateSolutionPoints(GKSFRMesh1D& mesh, double dt);

void GKSFR_AdvanceOneStep(
	GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary);
