#include "gks_fr_high_order.h"

#include <cmath>

namespace
{
	const double kPi = 3.14159265358979323846;
	const double kGLNode = std::sqrt(3.0 / 5.0);
	double g_boundary_time_2d = 0.0;

	void Zero3(double* a)
	{
		for (int m = 0; m < 3; ++m)
		{
			a[m] = 0.0;
		}
	}

	void Zero4(double* a)
	{
		for (int m = 0; m < 4; ++m)
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

	void Copy4(double* dst, const double* src)
	{
		for (int m = 0; m < 4; ++m)
		{
			dst[m] = src[m];
		}
	}

	void ZeroPointDerivatives2D(Point2d& point)
	{
		Zero4(point.der1x);
		Zero4(point.der1y);
		Zero4(point.der2xx);
		Zero4(point.der2xy);
		Zero4(point.der2yy);
	}

	void SetPointConservativeNoDerivatives2D(Point2d& point, const double* convar)
	{
		Copy4(point.convar, convar);
		ZeroPointDerivatives2D(point);
		point.x = 0.0;
		point.y = 0.0;
		point.normal[0] = 1.0;
		point.normal[1] = 0.0;
		point.is_reduce_order = false;
	}

	double LagrangeAt(int basis_id, double s)
	{
		const double a_gl = kGLNode;
		const double a2 = a_gl * a_gl;
		if (basis_id == 0)
		{
			return s * (s - a_gl) / (2.0 * a2);
		}
		if (basis_id == 1)
		{
			return 1.0 - s * s / a2;
		}
		return s * (s + a_gl) / (2.0 * a2);
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

	double LagrangeSecondDerivAt(int basis_id)
	{
		const double a2 = kGLNode * kGLNode;
		if (basis_id == 1)
		{
			return -2.0 / a2;
		}
		return 1.0 / a2;
	}

	double RadauLeft(double s)
	{
		return (-5.0 * s * s * s + 3.0 * s * s + 3.0 * s - 1.0) / 4.0;
	}

	double RadauRight(double s)
	{
		return (5.0 * s * s * s + 3.0 * s * s - 3.0 * s - 1.0) / 4.0;
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

	void AccumulateMoment4(double* dst, const double* src, double scale)
	{
		for (int m = 0; m < 4; ++m)
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

	void ZeroTraceDerivatives(Point1d& trace)
	{
		for (int m = 0; m < 3; ++m)
		{
			trace.der1[m] = 0.0;
			trace.der2[m] = 0.0;
		}
	}

	void ReflectNormalVelocity(Point1d& trace)
	{
		trace.convar[1] = -trace.convar[1];
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

	void ComputeStrictTransmissiveBoundaryCommonFlux(
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

		// Strict transmissive closure: extend the boundary face state as a
		// constant state outside the domain instead of copying the interior
		// polynomial derivatives across the physical boundary.
		ZeroTraceDerivatives(trace);
		interface.left = trace;
		interface.right = trace;
		interface.center = trace;

		const GKS1d_type old_solver = gks1dsolver;
		gks1dsolver = gks3rd;
		GKS(flux, interface, dt);
		gks1dsolver = old_solver;

		for (int m = 0; m < 3; ++m)
		{
			common_flux[m] = flux.f[m] / dt;
		}
	}

	void ComputeReflectiveBoundaryCommonFlux(
		const GKSFRCell1D& cell,
		double h,
		double dt,
		bool left_boundary,
		double common_flux[3])
	{
		Interface1d interface;
		Flux1d flux;
		Zero3(flux.f);

		Point1d interior_trace;
		FillFaceTrace(interior_trace, cell, h, !left_boundary ? true : false);
		if (left_boundary)
		{
			FillFaceTrace(interior_trace, cell, h, false);
		}
		else
		{
			FillFaceTrace(interior_trace, cell, h, true);
		}

		Point1d exterior_trace = interior_trace;
		ZeroTraceDerivatives(interior_trace);
		ZeroTraceDerivatives(exterior_trace);
		ReflectNormalVelocity(exterior_trace);

		interface.left = exterior_trace;
		interface.right = interior_trace;
		interface.center = interior_trace;

		const GKS1d_type old_solver = gks1dsolver;
		gks1dsolver = gks3rd;
		GKS(flux, interface, dt);
		gks1dsolver = old_solver;

		for (int m = 0; m < 3; ++m)
		{
			common_flux[m] = flux.f[m] / dt;
		}
	}

	void Moment2D(double* moment, int no_u, int no_v, const double* coeff, const double* prim)
	{
		double coeff_local[4];
		Copy4(coeff_local, coeff);
		MMDF m(prim[1], prim[2], prim[3]);
		G_address(no_u, no_v, 0, moment, coeff_local, m);
	}

	void AddMoment2D(double* accum, int no_u, int no_v, const double* coeff, const double* prim, double scale)
	{
		if (std::fabs(scale) <= 1.0e-30)
		{
			return;
		}
		double moment[4];
		Moment2D(moment, no_u, no_v, coeff, prim);
		AccumulateMoment4(accum, moment, prim[0] * scale);
	}

	void SolveCoeffFromMoment2D(double* coeff, const double* moment, const double* prim)
	{
		double rhs[4];
		double prim_local[4];
		Copy4(rhs, moment);
		Copy4(prim_local, prim);
		A(coeff, rhs, prim_local);
	}

	void SolveTimeCoeff2D(double* coeff, const double* x_coeff, const double* y_coeff, const double* prim)
	{
		double rhs[4]{ 0.0, 0.0, 0.0, 0.0 };
		AddMoment2D(rhs, 1, 0, x_coeff, prim, -1.0);
		AddMoment2D(rhs, 0, 1, y_coeff, prim, -1.0);
		SolveCoeffFromMoment2D(coeff, rhs, prim);
	}

	void RotateGlobalToYLocal(double* dst, const double* src)
	{
		dst[0] = src[0];
		dst[1] = src[2];
		dst[2] = -src[1];
		dst[3] = src[3];
	}

	void RotateYLocalToGlobal(double* dst, const double* src)
	{
		dst[0] = src[0];
		dst[1] = -src[2];
		dst[2] = src[1];
		dst[3] = src[3];
	}

	void RotateGlobalToYLocalScaled(double* dst, const double* src, double scale)
	{
		double tmp[4];
		RotateGlobalToYLocal(tmp, src);
		for (int m = 0; m < 4; ++m)
		{
			dst[m] = scale * tmp[m];
		}
	}

	void BuildBoundaryPointFromInner2D(
		Point2d& ghost,
		const Point2d& inner,
		GKSFRBoundary2D boundary,
		GKSFRBoundarySide2D side,
		bool x_normal,
		double x,
		double y)
	{
		double inner_global[4];
		if (x_normal)
		{
			Copy4(inner_global, inner.convar);
		}
		else
		{
			RotateYLocalToGlobal(inner_global, inner.convar);
		}

		double ghost_global[4];
		GKSFR_BoundaryGhostState2D(
			ghost_global,
			inner_global,
			boundary,
			side,
			x,
			y,
			g_boundary_time_2d);

		if (x_normal)
		{
			SetPointConservativeNoDerivatives2D(ghost, ghost_global);
		}
		else
		{
			double ghost_local[4];
			RotateGlobalToYLocal(ghost_local, ghost_global);
			SetPointConservativeNoDerivatives2D(ghost, ghost_local);
		}
	}

	void EvaluateCellStateAtXY(Point2d& trace, const GKSFRCell2D& cell, double hx, double hy, double xbar, double ybar)
	{
		const double r = 2.0 * xbar / hx;
		const double s = 2.0 * ybar / hy;
		for (int m = 0; m < 4; ++m)
		{
			trace.convar[m] = 0.0;
			trace.der1x[m] = 0.0;
			trace.der1y[m] = 0.0;
			trace.der2xx[m] = 0.0;
			trace.der2xy[m] = 0.0;
			trace.der2yy[m] = 0.0;
			for (int i = 0; i < 3; ++i)
			{
				const double li = LagrangeAt(i, r);
				const double lix = (2.0 / hx) * LagrangeDerivAt(i, r);
				const double lixx = (4.0 / (hx * hx)) * LagrangeSecondDerivAt(i);
				for (int j = 0; j < 3; ++j)
				{
					const double lj = LagrangeAt(j, s);
					const double ljy = (2.0 / hy) * LagrangeDerivAt(j, s);
					const double ljyy = (4.0 / (hy * hy)) * LagrangeSecondDerivAt(j);
					const double q = cell.Q[i][j][m];
					trace.convar[m] += q * li * lj;
					trace.der1x[m] += q * lix * lj;
					trace.der1y[m] += q * li * ljy;
					trace.der2xx[m] += q * lixx * lj;
					trace.der2xy[m] += q * lix * ljy;
					trace.der2yy[m] += q * li * ljyy;
				}
			}
		}
		trace.x = xbar;
		trace.y = ybar;
		trace.normal[0] = 1.0;
		trace.normal[1] = 0.0;
		trace.is_reduce_order = false;
	}

	int PeriodicIndex(int i, int n)
	{
		int r = i % n;
		if (r < 0)
		{
			r += n;
		}
		return r;
	}

	void FillXFaceTrace(Point2d& trace, const GKSFRCell2D& cell, double hx, double hy, bool right_face, double s)
	{
		const double xbar = right_face ? 0.5 * hx : -0.5 * hx;
		const double ybar = 0.5 * hy * s;
		EvaluateCellStateAtXY(trace, cell, hx, hy, xbar, ybar);
	}

	void FillYFaceTraceLocal(Point2d& trace, const GKSFRCell2D& cell, double hx, double hy, bool top_face, double r)
	{
		Point2d global_trace;
		const double xbar = 0.5 * hx * r;
		const double ybar = top_face ? 0.5 * hy : -0.5 * hy;
		EvaluateCellStateAtXY(global_trace, cell, hx, hy, xbar, ybar);

		RotateGlobalToYLocal(trace.convar, global_trace.convar);
		RotateGlobalToYLocal(trace.der1x, global_trace.der1y);
		RotateGlobalToYLocalScaled(trace.der1y, global_trace.der1x, -1.0);
		RotateGlobalToYLocal(trace.der2xx, global_trace.der2yy);
		RotateGlobalToYLocalScaled(trace.der2xy, global_trace.der2xy, -1.0);
		RotateGlobalToYLocal(trace.der2yy, global_trace.der2xx);
		trace.x = ybar;
		trace.y = -xbar;
		trace.normal[0] = 1.0;
		trace.normal[1] = 0.0;
		trace.is_reduce_order = false;
	}

	void BuildXFaceTraceFR(
		Recon2d& interface,
		const GKSFRMesh2D& mesh,
		int face_i,
		int cell_j,
		GKSFRBoundary2D boundary,
		double s)
	{
		const int j = PeriodicIndex(cell_j, mesh.cells_y);
		const double x = mesh.x_left + face_i * mesh.dx;
		const double y = GKSFR_CellCenterY2D(mesh, j) + 0.5 * mesh.dy * s;
		if (boundary == gksfr2d_periodic || (face_i > 0 && face_i < mesh.cells_x))
		{
			const int left_i = (boundary == gksfr2d_periodic)
				? PeriodicIndex(face_i - 1, mesh.cells_x)
				: face_i - 1;
			const int right_i = (boundary == gksfr2d_periodic)
				? PeriodicIndex(face_i, mesh.cells_x)
				: face_i;
			FillXFaceTrace(interface.left, mesh.cell[GKSFR_CellIndex2D(mesh, left_i, j)], mesh.dx, mesh.dy, true, s);
			FillXFaceTrace(interface.right, mesh.cell[GKSFR_CellIndex2D(mesh, right_i, j)], mesh.dx, mesh.dy, false, s);
			return;
		}

		if (face_i == 0)
		{
			FillXFaceTrace(interface.right, mesh.cell[GKSFR_CellIndex2D(mesh, 0, j)], mesh.dx, mesh.dy, false, s);
			BuildBoundaryPointFromInner2D(interface.left, interface.right, boundary, gksfr2d_left_side, true, x, y);
			return;
		}

		FillXFaceTrace(interface.left, mesh.cell[GKSFR_CellIndex2D(mesh, mesh.cells_x - 1, j)], mesh.dx, mesh.dy, true, s);
		BuildBoundaryPointFromInner2D(interface.right, interface.left, boundary, gksfr2d_right_side, true, x, y);
	}

	void BuildYFaceTraceFRLocal(
		Recon2d& interface,
		const GKSFRMesh2D& mesh,
		int cell_i,
		int face_j,
		GKSFRBoundary2D boundary,
		double r)
	{
		const int i = PeriodicIndex(cell_i, mesh.cells_x);
		const double x = GKSFR_CellCenterX2D(mesh, i) + 0.5 * mesh.dx * r;
		const double y = mesh.y_bottom + face_j * mesh.dy;
		if (boundary == gksfr2d_periodic || (face_j > 0 && face_j < mesh.cells_y))
		{
			const int bottom_j = (boundary == gksfr2d_periodic)
				? PeriodicIndex(face_j - 1, mesh.cells_y)
				: face_j - 1;
			const int top_j = (boundary == gksfr2d_periodic)
				? PeriodicIndex(face_j, mesh.cells_y)
				: face_j;
			FillYFaceTraceLocal(interface.left, mesh.cell[GKSFR_CellIndex2D(mesh, i, bottom_j)], mesh.dx, mesh.dy, true, r);
			FillYFaceTraceLocal(interface.right, mesh.cell[GKSFR_CellIndex2D(mesh, i, top_j)], mesh.dx, mesh.dy, false, r);
			return;
		}

		if (face_j == 0)
		{
			FillYFaceTraceLocal(interface.right, mesh.cell[GKSFR_CellIndex2D(mesh, i, 0)], mesh.dx, mesh.dy, false, r);
			BuildBoundaryPointFromInner2D(interface.left, interface.right, boundary, gksfr2d_bottom_side, false, x, y);
			return;
		}

		FillYFaceTraceLocal(interface.left, mesh.cell[GKSFR_CellIndex2D(mesh, i, mesh.cells_y - 1)], mesh.dx, mesh.dy, true, r);
		BuildBoundaryPointFromInner2D(interface.right, interface.left, boundary, gksfr2d_top_side, false, x, y);
	}

	void BuildCenterDerivativeFromFaces(
		double* center_der,
		const double* left_der,
		const double* right_der,
		const double* prim_left,
		const double* prim_right,
		MMDF& ml,
		MMDF& mr)
	{
		double left_coeff[4];
		double right_coeff[4];
		SolveCoeffFromMoment2D(left_coeff, left_der, prim_left);
		SolveCoeffFromMoment2D(right_coeff, right_der, prim_right);

		double left_moment[4];
		double right_moment[4];
		GL_address(0, 0, 0, left_moment, left_coeff, ml);
		GR_address(0, 0, 0, right_moment, right_coeff, mr);
		for (int m = 0; m < 4; ++m)
		{
			center_der[m] = prim_left[0] * left_moment[m] + prim_right[0] * right_moment[m];
		}
	}

	void BuildCenterTraceFromFaces2D(Recon2d& interface)
	{
		double prim_left[4];
		double prim_right[4];
		double convar_left[4];
		double convar_right[4];
		Copy4(convar_left, interface.left.convar);
		Copy4(convar_right, interface.right.convar);
		Convar_to_ULambda_2d(prim_left, convar_left);
		Convar_to_ULambda_2d(prim_right, convar_right);

		MMDF ml(prim_left[1], prim_left[2], prim_left[3]);
		MMDF mr(prim_right[1], prim_right[2], prim_right[3]);
		Collision(interface.center.convar, prim_left[0], prim_right[0], ml, mr);

		BuildCenterDerivativeFromFaces(
			interface.center.der1x,
			interface.left.der1x,
			interface.right.der1x,
			prim_left,
			prim_right,
			ml,
			mr);
		BuildCenterDerivativeFromFaces(
			interface.center.der1y,
			interface.left.der1y,
			interface.right.der1y,
			prim_left,
			prim_right,
			ml,
			mr);
		BuildCenterDerivativeFromFaces(
			interface.center.der2xx,
			interface.left.der2xx,
			interface.right.der2xx,
			prim_left,
			prim_right,
			ml,
			mr);
		BuildCenterDerivativeFromFaces(
			interface.center.der2xy,
			interface.left.der2xy,
			interface.right.der2xy,
			prim_left,
			prim_right,
			ml,
			mr);
		BuildCenterDerivativeFromFaces(
			interface.center.der2yy,
			interface.left.der2yy,
			interface.right.der2yy,
			prim_left,
			prim_right,
			ml,
			mr);
		interface.center.normal[0] = 1.0;
		interface.center.normal[1] = 0.0;
		interface.center.is_reduce_order = false;
	}

	void ComputeInterfaceCommonFluxPreparedCenter2D(
		Recon2d& interface,
		double dt,
		bool x_normal,
		double common_flux[4])
	{
		interface.normal[0] = 1.0;
		interface.normal[1] = 0.0;

		Flux2d flux;
		Zero4(flux.f);
		Zero4(flux.derf);
		Zero4(flux.der2f);
		flux.normal[0] = 1.0;
		flux.normal[1] = 0.0;

		const GKS2d_type old_solver = gks2dsolver;
		gks2dsolver = gks3rd_2d;
		GKS2D(flux, interface, dt);
		gks2dsolver = old_solver;

		double local_average[4];
		for (int m = 0; m < 4; ++m)
		{
			local_average[m] = flux.f[m] / dt;
		}
		if (x_normal)
		{
			Copy4(common_flux, local_average);
		}
		else
		{
			RotateYLocalToGlobal(common_flux, local_average);
		}
	}

	void ComputeInterfaceCommonFluxFromRecon2D(
		Recon2d& interface,
		double dt,
		bool x_normal,
		double common_flux[4])
	{
		BuildCenterTraceFromFaces2D(interface);
		ComputeInterfaceCommonFluxPreparedCenter2D(interface, dt, x_normal, common_flux);
	}

	void ComputeInterfaceCommonFlux2D(
		const GKSFRCell2D& left_cell,
		const GKSFRCell2D& right_cell,
		double hx,
		double hy,
		double dt,
		bool x_normal,
		double tangential_point,
		double common_flux[4])
	{
		Recon2d interface;
		if (x_normal)
		{
			FillXFaceTrace(interface.left, left_cell, hx, hy, true, tangential_point);
			FillXFaceTrace(interface.right, right_cell, hx, hy, false, tangential_point);
		}
		else
		{
			FillYFaceTraceLocal(interface.left, left_cell, hx, hy, true, tangential_point);
			FillYFaceTraceLocal(interface.right, right_cell, hx, hy, false, tangential_point);
		}

		ComputeInterfaceCommonFluxFromRecon2D(interface, dt, x_normal, common_flux);
	}

	void ComputeInterfaceCommonFluxFRStencil2D(
		const GKSFRMesh2D& mesh,
		int normal_index,
		int tangential_index,
		int tangential_point_id,
		double dt,
		GKSFRBoundary2D boundary,
		bool x_normal,
		double common_flux[4])
	{
		Recon2d interface;
		if (x_normal)
		{
			BuildXFaceTraceFR(
				interface,
				mesh,
				normal_index,
				tangential_index,
				boundary,
				GKSFR_GL_Point(tangential_point_id));
		}
		else
		{
			BuildYFaceTraceFRLocal(
				interface,
				mesh,
				tangential_index,
				normal_index,
				boundary,
				GKSFR_GL_Point(tangential_point_id));
		}
		BuildCenterTraceFromFaces2D(interface);
		ComputeInterfaceCommonFluxPreparedCenter2D(interface, dt, x_normal, common_flux);
	}

	void EvaluateInternalFluxPoly2D(
		double* flux,
		const double coef[6][4],
		double xbar,
		double ybar)
	{
		for (int m = 0; m < 4; ++m)
		{
			flux[m] = coef[0][m]
				+ coef[1][m] * xbar
				+ coef[2][m] * ybar
				+ 0.5 * coef[3][m] * xbar * xbar
				+ coef[4][m] * xbar * ybar
				+ 0.5 * coef[5][m] * ybar * ybar;
		}
	}

	int XFaceIndex(const GKSFRMesh2D& mesh, int face_i, int cell_j)
	{
		return face_i * mesh.cells_y + cell_j;
	}

	int YFaceIndex(const GKSFRMesh2D& mesh, int cell_i, int face_j)
	{
		return cell_i * (mesh.cells_y + 1) + face_j;
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

	if (boundary == gksfr_reflective)
	{
		ComputeReflectiveBoundaryCommonFlux(
			mesh.cell[0], mesh.dx, dt, true, face_fluxes[0].F);
	}
	else if (boundary == gksfr_transmissive_strict)
	{
		ComputeStrictTransmissiveBoundaryCommonFlux(
			mesh.cell[0], mesh.dx, dt, true, face_fluxes[0].F);
	}
	else
	{
		ComputeBoundaryCommonFlux(mesh.cell[0], mesh.dx, dt, true, face_fluxes[0].F);
	}
	for (int e = 0; e < mesh.cells - 1; ++e)
	{
		GKSFR_ComputeInterfaceCommonFlux(
			mesh.cell[e], mesh.cell[e + 1], mesh.dx, dt, face_fluxes[e + 1].F);
	}
	if (boundary == gksfr_reflective)
	{
		ComputeReflectiveBoundaryCommonFlux(
			mesh.cell[mesh.cells - 1], mesh.dx, dt, false, face_fluxes[mesh.cells].F);
	}
	else if (boundary == gksfr_transmissive_strict)
	{
		ComputeStrictTransmissiveBoundaryCommonFlux(
			mesh.cell[mesh.cells - 1], mesh.dx, dt, false, face_fluxes[mesh.cells].F);
	}
	else
	{
		ComputeBoundaryCommonFlux(
			mesh.cell[mesh.cells - 1], mesh.dx, dt, false, face_fluxes[mesh.cells].F);
	}
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

GKSFRCell2D::GKSFRCell2D() : tau(0.0)
{
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int m = 0; m < 4; ++m)
			{
				Q[i][j][m] = 0.0;
				F_gl[i][j][m] = 0.0;
				G_gl[i][j][m] = 0.0;
				residual[i][j][m] = 0.0;
			}
		}
		for (int side = 0; side < 2; ++side)
		{
			for (int m = 0; m < 4; ++m)
			{
				F_face_local[side][i][m] = 0.0;
				G_face_local[side][i][m] = 0.0;
			}
		}
	}
	for (int m = 0; m < 4; ++m)
	{
		Qc[m] = 0.0;
		Qx[m] = 0.0;
		Qy[m] = 0.0;
		Qxx[m] = 0.0;
		Qxy[m] = 0.0;
		Qyy[m] = 0.0;
		prim[m] = 0.0;
		ax[m] = 0.0;
		ay[m] = 0.0;
		A[m] = 0.0;
		bxx[m] = 0.0;
		bxy[m] = 0.0;
		byy[m] = 0.0;
		cx[m] = 0.0;
		cy[m] = 0.0;
		d[m] = 0.0;
	}
	for (int k = 0; k < 6; ++k)
	{
		for (int m = 0; m < 4; ++m)
		{
			Fcoef[k][m] = 0.0;
			Gcoef[k][m] = 0.0;
		}
	}
}

GKSFRMesh2D::GKSFRMesh2D()
	: cells_x(0), cells_y(0), x_left(0.0), x_right(0.0),
	y_bottom(0.0), y_top(0.0), dx(0.0), dy(0.0)
{
}

void GKSFR_ResizeUniformMesh2D(
	GKSFRMesh2D& mesh,
	int cells_x,
	int cells_y,
	double x_left,
	double x_right,
	double y_bottom,
	double y_top)
{
	mesh.cells_x = cells_x;
	mesh.cells_y = cells_y;
	mesh.x_left = x_left;
	mesh.x_right = x_right;
	mesh.y_bottom = y_bottom;
	mesh.y_top = y_top;
	mesh.dx = (x_right - x_left) / cells_x;
	mesh.dy = (y_top - y_bottom) / cells_y;
	mesh.cell.assign(cells_x * cells_y, GKSFRCell2D());
}

int GKSFR_CellIndex2D(const GKSFRMesh2D& mesh, int i, int j)
{
	return i * mesh.cells_y + j;
}

double GKSFR_CellCenterX2D(const GKSFRMesh2D& mesh, int i)
{
	return mesh.x_left + (i + 0.5) * mesh.dx;
}

double GKSFR_CellCenterY2D(const GKSFRMesh2D& mesh, int j)
{
	return mesh.y_bottom + (j + 0.5) * mesh.dy;
}

double GKSFR_SolutionPointX2D(const GKSFRMesh2D& mesh, int i, int p)
{
	return GKSFR_CellCenterX2D(mesh, i) + 0.5 * mesh.dx * GKSFR_GL_Point(p);
}

double GKSFR_SolutionPointY2D(const GKSFRMesh2D& mesh, int j, int q)
{
	return GKSFR_CellCenterY2D(mesh, j) + 0.5 * mesh.dy * GKSFR_GL_Point(q);
}

void GKSFR_SetBoundaryTime2D(double t)
{
	g_boundary_time_2d = t;
}

double GKSFR_GetBoundaryTime2D()
{
	return g_boundary_time_2d;
}

void GKSFR_DoubleMachPrimitive2D(double prim[4], double x, double y, double t)
{
	const double shock_x = 1.0 / 6.0 + (y + 20.0 * t) / std::sqrt(3.0);
	if (x < shock_x)
	{
		prim[0] = 8.0;
		prim[1] = 8.25 * std::cos(kPi / 6.0);
		prim[2] = -8.25 * std::sin(kPi / 6.0);
		prim[3] = 116.5;
	}
	else
	{
		prim[0] = 1.4;
		prim[1] = 0.0;
		prim[2] = 0.0;
		prim[3] = 1.0;
	}
}

void GKSFR_BoundaryGhostState2D(
	double ghost_Q[4],
	const double inner_Q[4],
	GKSFRBoundary2D boundary,
	GKSFRBoundarySide2D side,
	double x,
	double y,
	double t)
{
	if (boundary == gksfr2d_double_mach)
	{
		const bool exact_inflow =
			side == gksfr2d_left_side ||
			side == gksfr2d_top_side;
		if (exact_inflow)
		{
			double prim[4];
			GKSFR_DoubleMachPrimitive2D(prim, x, y, t);
			Primvar_to_convar_2D(ghost_Q, prim);
			return;
		}
		if (side == gksfr2d_bottom_side && x >= 1.0 / 6.0)
		{
			double prim[4];
			double inner_copy[4];
			Copy4(inner_copy, inner_Q);
			Convar_to_primvar_2D(prim, inner_copy);
			prim[2] = -prim[2];
			Primvar_to_convar_2D(ghost_Q, prim);
			return;
		}
		Copy4(ghost_Q, inner_Q);
		return;
	}

	if (boundary == gksfr2d_reflective)
	{
		double prim[4];
		double inner_copy[4];
		Copy4(inner_copy, inner_Q);
		Convar_to_primvar_2D(prim, inner_copy);
		if (side == gksfr2d_left_side || side == gksfr2d_right_side)
		{
			prim[1] = -prim[1];
		}
		else
		{
			prim[2] = -prim[2];
		}
		Primvar_to_convar_2D(ghost_Q, prim);
		return;
	}

	Copy4(ghost_Q, inner_Q);
}

void GKSFR_ComputeCellCenterData2D(GKSFRCell2D& cell, double hx, double hy)
{
	const double a_gl = kGLNode;
	const double a2 = a_gl * a_gl;
	for (int m = 0; m < 4; ++m)
	{
		cell.Qc[m] = cell.Q[1][1][m];
		cell.Qx[m] = (cell.Q[2][1][m] - cell.Q[0][1][m]) / (a_gl * hx);
		cell.Qy[m] = (cell.Q[1][2][m] - cell.Q[1][0][m]) / (a_gl * hy);
		cell.Qxx[m] = 4.0 * (cell.Q[0][1][m] - 2.0 * cell.Q[1][1][m] + cell.Q[2][1][m]) / (a2 * hx * hx);
		cell.Qyy[m] = 4.0 * (cell.Q[1][0][m] - 2.0 * cell.Q[1][1][m] + cell.Q[1][2][m]) / (a2 * hy * hy);
		cell.Qxy[m] = (cell.Q[2][2][m] - cell.Q[0][2][m] - cell.Q[2][0][m] + cell.Q[0][0][m]) / (a2 * hx * hy);
	}
}

void GKSFR_ComputeInternalGKSCoeffs2D(GKSFRCell2D& cell)
{
	double center_convar[4];
	Copy4(center_convar, cell.Qc);
	Convar_to_ULambda_2d(cell.prim, center_convar);
	cell.tau = Get_Tau_NS(cell.prim[0], cell.prim[3]);

	SolveCoeffFromMoment2D(cell.ax, cell.Qx, cell.prim);
	SolveCoeffFromMoment2D(cell.ay, cell.Qy, cell.prim);
	SolveTimeCoeff2D(cell.A, cell.ax, cell.ay, cell.prim);

	SolveCoeffFromMoment2D(cell.bxx, cell.Qxx, cell.prim);
	SolveCoeffFromMoment2D(cell.bxy, cell.Qxy, cell.prim);
	SolveCoeffFromMoment2D(cell.byy, cell.Qyy, cell.prim);
	SolveTimeCoeff2D(cell.cx, cell.bxx, cell.bxy, cell.prim);
	SolveTimeCoeff2D(cell.cy, cell.bxy, cell.byy, cell.prim);
	SolveTimeCoeff2D(cell.d, cell.cx, cell.cy, cell.prim);
}

void GKSFR_ComputeInternalTimeAveragedFluxAtGLPoints2D(GKSFRCell2D& cell, double hx, double hy, double dt)
{
	for (int k = 0; k < 6; ++k)
	{
		Zero4(cell.Fcoef[k]);
		Zero4(cell.Gcoef[k]);
	}

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			const double xbar = 0.5 * hx * GKSFR_GL_Point(i);
			const double ybar = 0.5 * hy * GKSFR_GL_Point(j);
			Point2d point_state;
			EvaluateCellStateAtXY(point_state, cell, hx, hy, xbar, ybar);

			double convar[4];
			double prim[4];
			Copy4(convar, point_state.convar);
			Convar_to_ULambda_2d(prim, convar);
			const double tau = Get_Tau_NS(prim[0], prim[3]);

			double ax[4], ay[4], Acoef[4];
			double bxx[4], bxy[4], byy[4];
			double cx[4], cy[4], dcoef[4];
			SolveCoeffFromMoment2D(ax, point_state.der1x, prim);
			SolveCoeffFromMoment2D(ay, point_state.der1y, prim);
			SolveTimeCoeff2D(Acoef, ax, ay, prim);
			SolveCoeffFromMoment2D(bxx, point_state.der2xx, prim);
			SolveCoeffFromMoment2D(bxy, point_state.der2xy, prim);
			SolveCoeffFromMoment2D(byy, point_state.der2yy, prim);
			SolveTimeCoeff2D(cx, bxx, bxy, prim);
			SolveTimeCoeff2D(cy, bxy, byy, prim);
			SolveTimeCoeff2D(dcoef, cx, cy, prim);

			double unit[4]{ 1.0, 0.0, 0.0, 0.0 };
			Zero4(cell.F_gl[i][j]);
			AddMoment2D(cell.F_gl[i][j], 1, 0, unit, prim, 1.0);
			AddMoment2D(cell.F_gl[i][j], 1, 0, Acoef, prim, 0.5 * dt);
			AddMoment2D(cell.F_gl[i][j], 1, 0, dcoef, prim, dt * dt / 6.0);
			AddMoment2D(cell.F_gl[i][j], 2, 0, ax, prim, -tau);
			AddMoment2D(cell.F_gl[i][j], 1, 1, ay, prim, -tau);
			AddMoment2D(cell.F_gl[i][j], 1, 0, Acoef, prim, -tau);
			AddMoment2D(cell.F_gl[i][j], 2, 0, cx, prim, -0.5 * tau * dt);
			AddMoment2D(cell.F_gl[i][j], 1, 1, cy, prim, -0.5 * tau * dt);
			AddMoment2D(cell.F_gl[i][j], 1, 0, dcoef, prim, -0.5 * tau * dt);

			Zero4(cell.G_gl[i][j]);
			AddMoment2D(cell.G_gl[i][j], 0, 1, unit, prim, 1.0);
			AddMoment2D(cell.G_gl[i][j], 0, 1, Acoef, prim, 0.5 * dt);
			AddMoment2D(cell.G_gl[i][j], 0, 1, dcoef, prim, dt * dt / 6.0);
			AddMoment2D(cell.G_gl[i][j], 1, 1, ax, prim, -tau);
			AddMoment2D(cell.G_gl[i][j], 0, 2, ay, prim, -tau);
			AddMoment2D(cell.G_gl[i][j], 0, 1, Acoef, prim, -tau);
			AddMoment2D(cell.G_gl[i][j], 1, 1, cx, prim, -0.5 * tau * dt);
			AddMoment2D(cell.G_gl[i][j], 0, 2, cy, prim, -0.5 * tau * dt);
			AddMoment2D(cell.G_gl[i][j], 0, 1, dcoef, prim, -0.5 * tau * dt);
		}
	}
}

void GKSFR_ExtrapolateLocalFluxToFaces2D(GKSFRCell2D& cell)
{
	for (int j = 0; j < 3; ++j)
	{
		Zero4(cell.F_face_local[0][j]);
		Zero4(cell.F_face_local[1][j]);
		for (int p = 0; p < 3; ++p)
		{
			const double left_w = LagrangeAt(p, -1.0);
			const double right_w = LagrangeAt(p, 1.0);
			AccumulateMoment4(cell.F_face_local[0][j], cell.F_gl[p][j], left_w);
			AccumulateMoment4(cell.F_face_local[1][j], cell.F_gl[p][j], right_w);
		}
	}
	for (int i = 0; i < 3; ++i)
	{
		Zero4(cell.G_face_local[0][i]);
		Zero4(cell.G_face_local[1][i]);
		for (int q = 0; q < 3; ++q)
		{
			const double bottom_w = LagrangeAt(q, -1.0);
			const double top_w = LagrangeAt(q, 1.0);
			AccumulateMoment4(cell.G_face_local[0][i], cell.G_gl[i][q], bottom_w);
			AccumulateMoment4(cell.G_face_local[1][i], cell.G_gl[i][q], top_w);
		}
	}
}

void GKSFR_PrepareCell2D(GKSFRCell2D& cell, double hx, double hy, double dt)
{
	GKSFR_ComputeCellCenterData2D(cell, hx, hy);
	GKSFR_ComputeInternalGKSCoeffs2D(cell);
	GKSFR_ComputeInternalTimeAveragedFluxAtGLPoints2D(cell, hx, hy, dt);
	GKSFR_ExtrapolateLocalFluxToFaces2D(cell);
}

void GKSFR_PrepareCells2D(GKSFRMesh2D& mesh, double dt)
{
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			GKSFR_PrepareCell2D(mesh.cell[GKSFR_CellIndex2D(mesh, i, j)], mesh.dx, mesh.dy, dt);
		}
	}
}

void GKSFR_ComputeCommonInterfaceFluxes2D(
	const GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	std::vector<GKSFRFaceFlux2D>& y_face_fluxes)
{
	x_face_fluxes.assign((mesh.cells_x + 1) * mesh.cells_y, GKSFRFaceFlux2D());
	y_face_fluxes.assign(mesh.cells_x * (mesh.cells_y + 1), GKSFRFaceFlux2D());
	if (mesh.cells_x <= 0 || mesh.cells_y <= 0)
	{
		return;
	}

	for (int face_i = 0; face_i <= mesh.cells_x; ++face_i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int q = 0; q < 3; ++q)
			{
				ComputeInterfaceCommonFluxFRStencil2D(
					mesh,
					face_i,
					j,
					q,
					dt,
					boundary,
					true,
					x_face_fluxes[XFaceIndex(mesh, face_i, j)].F[q]);
			}
		}
	}

	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int face_j = 0; face_j <= mesh.cells_y; ++face_j)
		{
			for (int p = 0; p < 3; ++p)
			{
				ComputeInterfaceCommonFluxFRStencil2D(
					mesh,
					face_j,
					i,
					p,
					dt,
					boundary,
					false,
					y_face_fluxes[YFaceIndex(mesh, i, face_j)].F[p]);
			}
		}
	}
}

void GKSFR_ComputeResiduals2D(
	GKSFRMesh2D& mesh,
	const std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
	const std::vector<GKSFRFaceFlux2D>& y_face_fluxes)
{
	for (int e_i = 0; e_i < mesh.cells_x; ++e_i)
	{
		for (int e_j = 0; e_j < mesh.cells_y; ++e_j)
		{
			GKSFRCell2D& cell = mesh.cell[GKSFR_CellIndex2D(mesh, e_i, e_j)];
			for (int i = 0; i < 3; ++i)
			{
				const double r = GKSFR_GL_Point(i);
				for (int j = 0; j < 3; ++j)
				{
					const double s = GKSFR_GL_Point(j);
					for (int m = 0; m < 4; ++m)
					{
						double dF_dr = 0.0;
						double dG_ds = 0.0;
						for (int p = 0; p < 3; ++p)
						{
							for (int q = 0; q < 3; ++q)
							{
								dF_dr += cell.F_gl[p][q][m] * LagrangeDerivAt(p, r) * LagrangeAt(q, s);
								dG_ds += cell.G_gl[p][q][m] * LagrangeAt(p, r) * LagrangeDerivAt(q, s);
							}
						}
						const double* x_left_hat = x_face_fluxes[XFaceIndex(mesh, e_i, e_j)].F[j];
						const double* x_right_hat = x_face_fluxes[XFaceIndex(mesh, e_i + 1, e_j)].F[j];
						const double* y_bottom_hat = y_face_fluxes[YFaceIndex(mesh, e_i, e_j)].F[i];
						const double* y_top_hat = y_face_fluxes[YFaceIndex(mesh, e_i, e_j + 1)].F[i];
						dF_dr += (x_left_hat[m] - cell.F_face_local[0][j][m]) * RadauLeftDeriv(r);
						dF_dr += (x_right_hat[m] - cell.F_face_local[1][j][m]) * RadauRightDeriv(r);
						dG_ds += (y_bottom_hat[m] - cell.G_face_local[0][i][m]) * RadauLeftDeriv(s);
						dG_ds += (y_top_hat[m] - cell.G_face_local[1][i][m]) * RadauRightDeriv(s);

						cell.residual[i][j][m] = -2.0 / mesh.dx * dF_dr - 2.0 / mesh.dy * dG_ds;
					}
				}
			}
		}
	}
}

void GKSFR_UpdateSolutionPoints2D(GKSFRMesh2D& mesh, double dt)
{
	for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					mesh.cell[e].Q[i][j][m] += dt * mesh.cell[e].residual[i][j][m];
				}
			}
		}
	}
}

void GKSFR_AdvanceOneStep2D(
	GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary)
{
	std::vector<GKSFRFaceFlux2D> x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> y_face_fluxes;
	GKSFR_PrepareCells2D(mesh, dt);
	GKSFR_ComputeCommonInterfaceFluxes2D(mesh, dt, boundary, x_face_fluxes, y_face_fluxes);
	GKSFR_ComputeResiduals2D(mesh, x_face_fluxes, y_face_fluxes);
	GKSFR_UpdateSolutionPoints2D(mesh, dt);
}
