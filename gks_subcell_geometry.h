#pragma once

#include "gks_fr_high_order.h"

struct GKSSubcellGeom1D
{
	// Fixed to the current N = 2 high-order branch: three GL solution points
	// and three low-order subcells per big cell.
	double gl_s[3];
	double gl_w[3];
	double subcell_face_s[4];
	double subcell_center_s[3];
	double subcell_length_s[3];
	double subcell_length_x[3];
};

struct GKSSubcellGeom2D
{
	double gl_s[3];
	double gl_w[3];
	double subcell_weight[3];
	double subcell_face_s[4];
	double subcell_center_s[3];
	double subcell_length_s[3];
	double subcell_length_x[3];
	double subcell_length_y[3];
};

void GKSSubcellBuildGeometry1D(GKSSubcellGeom1D& geom, double big_dx);

void GKSSubcellBigCellAverageFromPoints1D(
	const double point_Q[3][3],
	double cell_avg[3]);

void GKSSubcellBigCellAverageFromLowOrderDofs1D(
	const GKSSubcellGeom1D& geom,
	const double low_dof[3][3],
	double cell_avg[3]);

void GKSSubcellBuildGeometry2D(GKSSubcellGeom2D& geom, double big_dx, double big_dy);

void GKSSubcellBigCellAverageFromPoints2D(
	const double point_Q[3][3][4],
	double cell_avg[4]);

void GKSSubcellBigCellAverageFromLowOrderDofs2D(
	const GKSSubcellGeom2D& geom,
	const double low_dof[3][3][4],
	double cell_avg[4]);
