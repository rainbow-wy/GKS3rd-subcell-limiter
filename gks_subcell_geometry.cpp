#include "gks_subcell_geometry.h"

#include <cmath>

void GKSSubcellBuildGeometry1D(GKSSubcellGeom1D& geom, double big_dx)
{
	const double a = std::sqrt(3.0 / 5.0);
	geom.gl_s[0] = -a;
	geom.gl_s[1] = 0.0;
	geom.gl_s[2] = a;
	// The low-order branch interprets the three current solution DOFs as three
	// subcell values. Their control-volume lengths follow the quadrature-based
	// partition and therefore sum to the full big-cell length.
	geom.gl_w[0] = 5.0 / 18.0;
	geom.gl_w[1] = 4.0 / 9.0;
	geom.gl_w[2] = 5.0 / 18.0;

	geom.subcell_face_s[0] = -1.0;
	geom.subcell_face_s[1] = -4.0 / 9.0;
	geom.subcell_face_s[2] =  4.0 / 9.0;
	geom.subcell_face_s[3] =  1.0;

	for (int j = 0; j < 3; ++j)
	{
		geom.subcell_length_s[j] = geom.subcell_face_s[j + 1] - geom.subcell_face_s[j];
		geom.subcell_center_s[j] = 0.5 * (geom.subcell_face_s[j + 1] + geom.subcell_face_s[j]);
		geom.subcell_length_x[j] = geom.gl_w[j] * big_dx;
	}
}

void GKSSubcellBigCellAverageFromPoints1D(
	const double point_Q[3][3],
	double cell_avg[3])
{
	const double w[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
	for (int m = 0; m < 3; ++m)
	{
		cell_avg[m] = 0.5 * (
			w[0] * point_Q[0][m] +
			w[1] * point_Q[1][m] +
			w[2] * point_Q[2][m]);
	}
}

void GKSSubcellBigCellAverageFromLowOrderDofs1D(
	const GKSSubcellGeom1D& geom,
	const double low_dof[3][3],
	double cell_avg[3])
{
	for (int m = 0; m < 3; ++m)
	{
		cell_avg[m] = 0.0;
		for (int j = 0; j < 3; ++j)
		{
			cell_avg[m] += geom.subcell_length_x[j] * low_dof[j][m];
		}
		cell_avg[m] /= (geom.subcell_length_x[0] + geom.subcell_length_x[1] + geom.subcell_length_x[2]);
	}
}

void GKSSubcellBuildGeometry2D(GKSSubcellGeom2D& geom, double big_dx, double big_dy)
{
	geom.gl_s[0] = -std::sqrt(3.0 / 5.0);
	geom.gl_s[1] = 0.0;
	geom.gl_s[2] = std::sqrt(3.0 / 5.0);
	geom.gl_w[0] = 5.0 / 9.0;
	geom.gl_w[1] = 8.0 / 9.0;
	geom.gl_w[2] = 5.0 / 9.0;
	geom.subcell_face_s[0] = -1.0;
	geom.subcell_face_s[1] = -4.0 / 9.0;
	geom.subcell_face_s[2] = 4.0 / 9.0;
	geom.subcell_face_s[3] = 1.0;
	for (int i = 0; i < 3; ++i)
	{
		geom.subcell_weight[i] = 0.5 * geom.gl_w[i];
		geom.subcell_length_s[i] = geom.subcell_face_s[i + 1] - geom.subcell_face_s[i];
		geom.subcell_center_s[i] = 0.5 * (geom.subcell_face_s[i + 1] + geom.subcell_face_s[i]);
		geom.subcell_length_x[i] = geom.subcell_weight[i] * big_dx;
		geom.subcell_length_y[i] = geom.subcell_weight[i] * big_dy;
	}
}

void GKSSubcellBigCellAverageFromPoints2D(
	const double point_Q[3][3][4],
	double cell_avg[4])
{
	const double wtilde[3] = { 5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0 };
	for (int m = 0; m < 4; ++m)
	{
		cell_avg[m] = 0.0;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				cell_avg[m] += wtilde[i] * wtilde[j] * point_Q[i][j][m];
			}
		}
	}
}

void GKSSubcellBigCellAverageFromLowOrderDofs2D(
	const GKSSubcellGeom2D& geom,
	const double low_dof[3][3][4],
	double cell_avg[4])
{
	for (int m = 0; m < 4; ++m)
	{
		cell_avg[m] = 0.0;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				cell_avg[m] += geom.subcell_weight[i] * geom.subcell_weight[j] * low_dof[i][j][m];
			}
		}
	}
}
