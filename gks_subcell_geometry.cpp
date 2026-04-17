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
