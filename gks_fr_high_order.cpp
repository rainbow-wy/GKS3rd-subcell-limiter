#include "gks_fr_high_order.h"

#include <cmath>

namespace
{
	const double kPi = 3.14159265358979323846;
	const double kGLNode = std::sqrt(3.0 / 5.0);

	void Zero3(double* a)
	{
		for (int m = 0; m < 3; ++m)
		{
			a[m] = 0.0;
		}
	}

	void Copy3(double* dst, const double* src)
	{
		for (int m = 0; m < 3; ++m)
		{
			dst[m] = src[m];
		}
	}

	double LagrangeDerivAt(int basis_id, double s)
	{
		const double a_gl = kGLNode;
		const double a2 = a_gl * a_gl;
		if (basis_id == 0)
		{
			return (2.0 * s - a_gl) / (2.0 * a2);
		}
		if (basis_id == 1)
		{
			return -2.0 * s / a2;
		}
		return (2.0 * s + a_gl) / (2.0 * a2);
	}

	double RadauLeftDeriv(double s)
	{
		return (-15.0 * s * s + 6.0 * s + 3.0) / 4.0;
	}

	double RadauRightDeriv(double s)
	{
		return (15.0 * s * s + 6.0 * s - 3.0) / 4.0;
	}

	void AccumulateMoment(double* dst, const double* src, double scale)
	{
		for (int m = 0; m < 3; ++m)
		{
			dst[m] += scale * src[m];
		}
	}

	void EvaluateTimeAveragedFluxAtXtilde(double* flux, const GKSFRCell1D& cell, double x_tilde)
	{
		for (int m = 0; m < 3; ++m)
		{
			flux[m] = cell.Fe0[m]
				+ cell.Fe1[m] * x_tilde
				+ 0.5 * cell.Fe2[m] * x_tilde * x_tilde;
		}
	}

	void FillFaceTrace(Point1d& trace, const GKSFRCell1D& cell, double h, bool right_face)
	{
		const double x_tilde = right_face ? (0.5 * h) : (-0.5 * h);
		for (int m = 0; m < 3; ++m)
		{
			trace.convar[m] = cell.Qc[m]
				+ cell.Qx[m] * x_tilde
				+ 0.5 * cell.Qxx[m] * x_tilde * x_tilde;
			trace.der1[m] = cell.Qx[m] + cell.Qxx[m] * x_tilde;
			trace.der2[m] = cell.Qxx[m];
		}
	}

	void BuildCenterTraceFromFaces(Interface1d& interface)
	{
		// Wrapper for the existing one-stage 3rd-order interface GKS:
		// use face traces from the FR polynomial and reconstruct the interface
		// equilibrium state plus its first/second macro derivatives without
		// touching the tested GKS() implementation itself.
		double prim_left[3], prim_right[3];
		double left_convar[3], right_convar[3];
		Copy3(left_convar, interface.left.convar);
		Copy3(right_convar, interface.right.convar);
		Convar_to_ULambda_1d(prim_left, left_convar);
		Convar_to_ULambda_1d(prim_right, right_convar);

		MMDF1d ml(prim_left[1], prim_left[2]);
		MMDF1d mr(prim_right[1], prim_right[2]);
		double unit[3]{ 1.0, 0.0, 0.0 };

		double gl[3], gr[3];
		GL(0, 0, gl, unit, ml);
		GR(0, 0, gr, unit, mr);
		for (int m = 0; m < 3; ++m)
		{
			interface.center.convar[m] = prim_left[0] * gl[m] + prim_right[0] * gr[m];
		}

		double a_left[3], a_right[3];
		double wx_left[3], wx_right[3];
		Copy3(wx_left, interface.left.der1);
		Copy3(wx_right, interface.right.der1);
		SolveGKS_a1D(a_left, wx_left, prim_left);
		SolveGKS_a1D(a_right, wx_right, prim_right);

		double gx_left[3], gx_right[3];
		GL(0, 0, gx_left, a_left, ml);
		GR(0, 0, gx_right, a_right, mr);
		for (int m = 0; m < 3; ++m)
		{
			interface.center.der1[m] = prim_left[0] * gx_left[m] + prim_right[0] * gx_right[m];
		}

		double b_left[3], b_right[3];
		double wxx_left[3], wxx_right[3];
		Copy3(wxx_left, interface.left.der2);
		Copy3(wxx_right, interface.right.der2);
		SolveGKS_axx1D(b_left, wxx_left, prim_left);
		SolveGKS_axx1D(b_right, wxx_right, prim_right);

		double gxx_left[3], gxx_right[3];
		GL(0, 0, gxx_left, b_left, ml);
		GR(0, 0, gxx_right, b_right, mr);
		for (int m = 0; m < 3; ++m)
		{
			interface.center.der2[m] = prim_left[0] * gxx_left[m] + prim_right[0] * gxx_right[m];
		}
	}

	void ComputeBoundaryCommonFlux(
		const GKSFRCell1D& cell,
		double h,
		double dt,
		bool left_boundary,
		double common_flux[3])
	{
		Interface1d interface;
		Flux1d flux;
		Zero3(flux.f);

		Point1d trace;
		FillFaceTrace(trace, cell, h, !left_boundary ? true : false);
		if (left_boundary)
		{
			FillFaceTrace(trace, cell, h, false);
		}
		else
		{
			FillFaceTrace(trace, cell, h, true);
		}

		interface.left = trace;
		interface.right = trace;
		BuildCenterTraceFromFaces(interface);

		const GKS1d_type old_solver = gks1dsolver;
		gks1dsolver = gks3rd;
		GKS(flux, interface, dt);
		gks1dsolver = old_solver;

		for (int m = 0; m < 3; ++m)
		{
			common_flux[m] = flux.f[m] / dt;
		}
	}
}

GKSFRCell1D::GKSFRCell1D() : tau(0.0)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int m = 0; m < 3; ++m)
		{
			Q[i][m] = 0.0;
			flux_gl[i][m] = 0.0;
			residual[i][m] = 0.0;
		}
		Qc[i] = 0.0;
		Qx[i] = 0.0;
		Qxx[i] = 0.0;
		prim[i] = 0.0;
		a[i] = 0.0;
		A[i] = 0.0;
		b[i] = 0.0;
		c[i] = 0.0;
		d[i] = 0.0;
		Fe0[i] = 0.0;
		Fe1[i] = 0.0;
		Fe2[i] = 0.0;
		flux_face_local[0][i] = 0.0;
		flux_face_local[1][i] = 0.0;
	}
}

GKSFRMesh1D::GKSFRMesh1D() : cells(0), x_left(0.0), x_right(0.0), dx(0.0)
{
}

double GKSFR_GL_Point(int point_id)
{
	if (point_id == 0)
	{
		return -kGLNode;
	}
	if (point_id == 1)
	{
		return 0.0;
	}
	return kGLNode;
}

void GKSFR_ResizeUniformMesh(GKSFRMesh1D& mesh, int cells, double x_left, double x_right)
{
	mesh.cells = cells;
	mesh.x_left = x_left;
	mesh.x_right = x_right;
	mesh.dx = (x_right - x_left) / cells;
	mesh.cell.assign(cells, GKSFRCell1D());
}

void GKSFR_ComputeCellCenterData(GKSFRCell1D& cell, double h)
{
	// GKS-FR.pdf Eq. (75)-(78): recover the center expansion from the three
	// GL solution-point conservative states.
	const double a_gl = kGLNode;
	for (int m = 0; m < 3; ++m)
	{
		cell.Qc[m] = cell.Q[1][m];
		cell.Qx[m] = (cell.Q[2][m] - cell.Q[0][m]) / (a_gl * h);
		cell.Qxx[m] = 20.0 / (3.0 * h * h) * (cell.Q[0][m] - 2.0 * cell.Q[1][m] + cell.Q[2][m]);
	}
}

void GKSFR_ComputeInternalGKSCoeffs(GKSFRCell1D& cell)
{
	// GKS-FR.pdf Eq. (36)-(40): solve the cell-internal coefficients using the
	// existing Maxwellian moment and linear-system utilities.
	double center_convar[3];
	Copy3(center_convar, cell.Qc);
	Convar_to_ULambda_1d(cell.prim, center_convar);
	cell.tau = Get_Tau_NS(cell.prim[0], cell.prim[2]);

	SolveGKS_a1D(cell.a, cell.Qx, cell.prim);
	SolveGKS_A1D(cell.A, cell.a, cell.prim);
	SolveGKS_axx1D(cell.b, cell.Qxx, cell.prim);
	SolveGKS_axt1D(cell.c, cell.b, cell.prim);
	SolveGKS_att1D(cell.d, cell.c, cell.prim);
}

void GKSFR_ComputeInternalTimeAveragedFluxAtGLPoints(GKSFRCell1D& cell, double h, double dt)
{
	// GKS-FR.pdf Eq. (46)-(50): construct the quadratic time-averaged internal
	// flux polynomial and evaluate it at the three GL points.
	MMDF1d m(cell.prim[1], cell.prim[2]);
	const double rho = cell.prim[0];
	double unit[3]{ 1.0, 0.0, 0.0 };
	double moment[3];

	Zero3(cell.Fe0);
	Zero3(cell.Fe1);
	Zero3(cell.Fe2);

	// Fe0
	G(1, 0, moment, unit, m);
	AccumulateMoment(cell.Fe0, moment, rho);
	G(1, 0, moment, cell.A, m);
	AccumulateMoment(cell.Fe0, moment, rho * 0.5 * dt);
	G(1, 0, moment, cell.d, m);
	AccumulateMoment(cell.Fe0, moment, rho * dt * dt / 6.0);
	G(2, 0, moment, cell.a, m);
	AccumulateMoment(cell.Fe0, moment, -rho * cell.tau);
	G(1, 0, moment, cell.A, m);
	AccumulateMoment(cell.Fe0, moment, -rho * cell.tau);
	G(2, 0, moment, cell.c, m);
	AccumulateMoment(cell.Fe0, moment, -rho * cell.tau * 0.5 * dt);
	G(1, 0, moment, cell.d, m);
	AccumulateMoment(cell.Fe0, moment, -rho * cell.tau * 0.5 * dt);

	// Fe1
	G(1, 0, moment, cell.a, m);
	AccumulateMoment(cell.Fe1, moment, rho);
	G(1, 0, moment, cell.c, m);
	AccumulateMoment(cell.Fe1, moment, rho * 0.5 * dt);
	G(2, 0, moment, cell.b, m);
	AccumulateMoment(cell.Fe1, moment, -rho * cell.tau);
	G(1, 0, moment, cell.c, m);
	AccumulateMoment(cell.Fe1, moment, -rho * cell.tau);

	// Fe2
	G(1, 0, moment, cell.b, m);
	AccumulateMoment(cell.Fe2, moment, rho);

	for (int j = 0; j < 3; ++j)
	{
		const double x_tilde = 0.5 * h * GKSFR_GL_Point(j);
		EvaluateTimeAveragedFluxAtXtilde(cell.flux_gl[j], cell, x_tilde);
	}
}

void GKSFR_ExtrapolateLocalFluxToFaces(GKSFRCell1D& cell)
{
	// GKS-FR.pdf Eq. (58)-(59): evaluate the discontinuous local flux
	// polynomial at the left/right faces through the GL nodal flux values.
	const double a_gl = kGLNode;
	const double left_w[3] = {
		5.0 * (1.0 + a_gl) / 6.0,
		-2.0 / 3.0,
		5.0 * (1.0 - a_gl) / 6.0
	};
	const double right_w[3] = {
		5.0 * (1.0 - a_gl) / 6.0,
		-2.0 / 3.0,
		5.0 * (1.0 + a_gl) / 6.0
	};

	Zero3(cell.flux_face_local[0]);
	Zero3(cell.flux_face_local[1]);
	for (int j = 0; j < 3; ++j)
	{
		AccumulateMoment(cell.flux_face_local[0], cell.flux_gl[j], left_w[j]);
		AccumulateMoment(cell.flux_face_local[1], cell.flux_gl[j], right_w[j]);
	}
}

void GKSFR_PrepareCell(GKSFRCell1D& cell, double h, double dt)
{
	GKSFR_ComputeCellCenterData(cell, h);
	GKSFR_ComputeInternalGKSCoeffs(cell);
	GKSFR_ComputeInternalTimeAveragedFluxAtGLPoints(cell, h, dt);
	GKSFR_ExtrapolateLocalFluxToFaces(cell);
}

void GKSFR_PrepareCells(GKSFRMesh1D& mesh, double dt)
{
	for (int e = 0; e < mesh.cells; ++e)
	{
		GKSFR_PrepareCell(mesh.cell[e], mesh.dx, dt);
	}
}

void GKSFR_ComputeInterfaceCommonFlux(
	const GKSFRCell1D& left_cell,
	const GKSFRCell1D& right_cell,
	double h,
	double dt,
	double common_flux[3])
{
	Interface1d interface;
	Flux1d flux;
	Zero3(flux.f);

	FillFaceTrace(interface.left, left_cell, h, true);
	FillFaceTrace(interface.right, right_cell, h, false);
	BuildCenterTraceFromFaces(interface);

	const GKS1d_type old_solver = gks1dsolver;
	gks1dsolver = gks3rd;
	GKS(flux, interface, dt);
	gks1dsolver = old_solver;

	for (int m = 0; m < 3; ++m)
	{
		common_flux[m] = flux.f[m] / dt;
	}
}

void GKSFR_ComputeCommonInterfaceFluxes(
	const GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	std::vector<GKSFRFaceFlux1D>& face_fluxes)
{
	face_fluxes.assign(mesh.cells + 1, GKSFRFaceFlux1D());

	if (mesh.cells <= 0)
	{
		return;
	}

	if (boundary == gksfr_periodic)
	{
		GKSFR_ComputeInterfaceCommonFlux(
			mesh.cell[mesh.cells - 1], mesh.cell[0], mesh.dx, dt, face_fluxes[0].F);
		for (int e = 0; e < mesh.cells - 1; ++e)
		{
			GKSFR_ComputeInterfaceCommonFlux(
				mesh.cell[e], mesh.cell[e + 1], mesh.dx, dt, face_fluxes[e + 1].F);
		}
		Copy3(face_fluxes[mesh.cells].F, face_fluxes[0].F);
		return;
	}

	ComputeBoundaryCommonFlux(mesh.cell[0], mesh.dx, dt, true, face_fluxes[0].F);
	for (int e = 0; e < mesh.cells - 1; ++e)
	{
		GKSFR_ComputeInterfaceCommonFlux(
			mesh.cell[e], mesh.cell[e + 1], mesh.dx, dt, face_fluxes[e + 1].F);
	}
	ComputeBoundaryCommonFlux(mesh.cell[mesh.cells - 1], mesh.dx, dt, false, face_fluxes[mesh.cells].F);
}

void GKSFR_ComputeResiduals(
	GKSFRMesh1D& mesh,
	const std::vector<GKSFRFaceFlux1D>& face_fluxes)
{
	for (int e = 0; e < mesh.cells; ++e)
	{
		GKSFRCell1D& cell = mesh.cell[e];
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
				dflux_ds += (face_fluxes[e].F[m] - cell.flux_face_local[0][m]) * corr_left;
				dflux_ds += (face_fluxes[e + 1].F[m] - cell.flux_face_local[1][m]) * corr_right;

				cell.residual[i][m] = -2.0 / mesh.dx * dflux_ds;
			}
		}
	}
}

void GKSFR_UpdateSolutionPoints(GKSFRMesh1D& mesh, double dt)
{
	for (int e = 0; e < mesh.cells; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int m = 0; m < 3; ++m)
			{
				mesh.cell[e].Q[i][m] += dt * mesh.cell[e].residual[i][m];
			}
		}
	}
}

void GKSFR_AdvanceOneStep(
	GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary)
{
	std::vector<GKSFRFaceFlux1D> face_fluxes;
	GKSFR_PrepareCells(mesh, dt);
	GKSFR_ComputeCommonInterfaceFluxes(mesh, dt, boundary, face_fluxes);
	GKSFR_ComputeResiduals(mesh, face_fluxes);
	GKSFR_UpdateSolutionPoints(mesh, dt);
}
