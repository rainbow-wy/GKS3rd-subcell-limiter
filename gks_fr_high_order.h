#pragma once

#include "fvm_gks1d.h"
#include "fvm_gks2d.h"

#include <vector>

enum GKSFRBoundary1D
{
	gksfr_periodic,
	gksfr_free,
	gksfr_transmissive_strict,
	gksfr_reflective
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

enum GKSFRBoundary2D
{
	gksfr2d_periodic,
	gksfr2d_transmissive,
	gksfr2d_reflective,
	gksfr2d_double_mach
};

enum GKSFRBoundarySide2D
{
	gksfr2d_left_side,
	gksfr2d_right_side,
	gksfr2d_bottom_side,
	gksfr2d_top_side
};

struct GKSFRFaceFlux2D
{
	double F[3][4];
};

struct GKSFRCell2D
{
	// Tensor-product GL solution points: Q[x_point][y_point][conservative_var].
	double Q[3][3][4];

	// Cell-center quadratic expansion recovered from the nine GL points.
	double Qc[4];
	double Qx[4];
	double Qy[4];
	double Qxx[4];
	double Qxy[4];
	double Qyy[4];

	double prim[4];
	double tau;

	// Cell-internal 2D one-stage third-order GKS coefficients.
	double ax[4];
	double ay[4];
	double A[4];
	double bxx[4];
	double bxy[4];
	double byy[4];
	double cx[4];
	double cy[4];
	double d[4];

	// Time-averaged internal flux polynomial coefficients in physical x/y.
	// coeff order: 00, 10, 01, 20, 11, 02.
	double Fcoef[6][4];
	double Gcoef[6][4];
	double F_gl[3][3][4];
	double G_gl[3][3][4];

	// Local discontinuous flux extrapolated to element faces.
	// x face: 0 left, 1 right, indexed by y GL point.
	// y face: 0 bottom, 1 top, indexed by x GL point.
	double F_face_local[2][3][4];
	double G_face_local[2][3][4];

	double residual[3][3][4];

	GKSFRCell2D();
};

struct GKSFRMesh2D
{
	int cells_x;
	int cells_y;
	double x_left;
	double x_right;
	double y_bottom;
	double y_top;
	double dx;
	double dy;
	std::vector<GKSFRCell2D> cell;

	GKSFRMesh2D();
};

void GKSFR_ResizeUniformMesh2D(
	GKSFRMesh2D& mesh,
	int cells_x,
	int cells_y,
	double x_left,
	double x_right,
	double y_bottom,
	double y_top);

int GKSFR_CellIndex2D(const GKSFRMesh2D& mesh, int i, int j);
double GKSFR_CellCenterX2D(const GKSFRMesh2D& mesh, int i);
double GKSFR_CellCenterY2D(const GKSFRMesh2D& mesh, int j);
double GKSFR_SolutionPointX2D(const GKSFRMesh2D& mesh, int i, int p);
double GKSFR_SolutionPointY2D(const GKSFRMesh2D& mesh, int j, int q);

void GKSFR_SetBoundaryTime2D(double t);
double GKSFR_GetBoundaryTime2D();
void GKSFR_DoubleMachPrimitive2D(double prim[4], double x, double y, double t);
void GKSFR_BoundaryGhostState2D(
	double ghost_Q[4],
	const double inner_Q[4],
	GKSFRBoundary2D boundary,
	GKSFRBoundarySide2D side,
	double x,
	double y,
	double t);

void GKSFR_ComputeCellCenterData2D(GKSFRCell2D& cell, double hx, double hy);
void GKSFR_ComputeInternalGKSCoeffs2D(GKSFRCell2D& cell);
void GKSFR_ComputeInternalTimeAveragedFluxAtGLPoints2D(GKSFRCell2D& cell, double hx, double hy, double dt);
void GKSFR_ExtrapolateLocalFluxToFaces2D(GKSFRCell2D& cell);
void GKSFR_PrepareCell2D(GKSFRCell2D& cell, double hx, double hy, double dt);
void GKSFR_PrepareCells2D(GKSFRMesh2D& mesh, double dt);

void GKSFR_ComputeCommonInterfaceFluxes2D(
	const GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& y_face_fluxes);

void GKSFR_ComputeResiduals2D(
	GKSFRMesh2D& mesh,
	const std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& y_face_fluxes);

void GKSFR_UpdateSolutionPoints2D(GKSFRMesh2D& mesh, double dt);

void GKSFR_AdvanceOneStep2D(
	GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary);
