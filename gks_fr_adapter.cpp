#include "gks_fr_adapter.h"

#include <cmath>

namespace
{
	const double kGLNode = std::sqrt(3.0 / 5.0);

	double LagrangeDerivAt(int basis_id, double s)
	{
		const double a2 = kGLNode * kGLNode;
		if (basis_id == 0)
		{
			return (2.0 * s - kGLNode) / (2.0 * a2);
		}
		if (basis_id == 1)
		{
			return -2.0 * s / a2;
		}
		return (2.0 * s + kGLNode) / (2.0 * a2);
	}

	double RadauLeftDeriv(double s)
	{
		return (-15.0 * s * s + 6.0 * s + 3.0) / 4.0;
	}

	double RadauRightDeriv(double s)
	{
		return (15.0 * s * s + 6.0 * s - 3.0) / 4.0;
	}
}

void GKSFRAdapterComputeHighFaceFluxes(
	const GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	std::vector<GKSFRFaceFlux1D>& face_fluxes)
{
	GKSFRMesh1D work = mesh;
	GKSFR_PrepareCells(work, dt);
	GKSFR_ComputeCommonInterfaceFluxes(work, dt, boundary, face_fluxes);
}

void GKSFRAdapterAdvanceWithFaceFluxes(
	const GKSFRMesh1D& mesh_old,
	double dt,
	const std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
	GKSFRMesh1D& mesh_new)
{
	mesh_new = mesh_old;
	GKSFR_PrepareCells(mesh_new, dt);
	GKSFR_ComputeResiduals(mesh_new, final_face_fluxes);
	GKSFR_UpdateSolutionPoints(mesh_new, dt);
}

void GKSFRAdapterAssembleBlendedHighResiduals(
	const GKSFRMesh1D& mesh_old,
	double dt,
	const std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
	const std::vector<double>& alpha_cell,
	GKSFRMesh1D& mesh_residual)
{
	mesh_residual = mesh_old;
	GKSFR_PrepareCells(mesh_residual, dt);

	for (int e = 0; e < mesh_residual.cells; ++e)
	{
		GKSFRCell1D& cell = mesh_residual.cell[e];
		for (int i = 0; i < 3; ++i)
		{
			const double s = GKSFR_GL_Point(i);
			const double corr_left = RadauLeftDeriv(s);
			const double corr_right = RadauRightDeriv(s);
			for (int m = 0; m < 3; ++m)
			{
				double dflux_ds = 0.0;
				for (int j = 0; j < 3; ++j)
				{
					dflux_ds += cell.flux_gl[j][m] * LagrangeDerivAt(j, s);
				}
				dflux_ds += (final_face_fluxes[e].F[m] - cell.flux_face_local[0][m]) * corr_left;
				dflux_ds += (final_face_fluxes[e + 1].F[m] - cell.flux_face_local[1][m]) * corr_right;
				cell.residual[i][m] = (1.0 - alpha_cell[e]) * (-2.0 / mesh_residual.dx * dflux_ds);
			}
		}
	}
}

void GKSFRAdapterComputeHighFaceFluxes2D(
	const GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& y_face_fluxes)
{
	GKSFRMesh2D work = mesh;
	GKSFR_PrepareCells2D(work, dt);
	GKSFR_ComputeCommonInterfaceFluxes2D(work, dt, boundary, x_face_fluxes, y_face_fluxes);
}

void GKSFRAdapterAdvanceWithFaceFluxes2D(
	const GKSFRMesh2D& mesh_old,
	double dt,
	const std::vector<GKSFRFaceFlux2D>& final_x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& final_y_face_fluxes,
	GKSFRMesh2D& mesh_new)
{
	mesh_new = mesh_old;
	GKSFR_PrepareCells2D(mesh_new, dt);
	GKSFR_ComputeResiduals2D(mesh_new, final_x_face_fluxes, final_y_face_fluxes);
	GKSFR_UpdateSolutionPoints2D(mesh_new, dt);
}
