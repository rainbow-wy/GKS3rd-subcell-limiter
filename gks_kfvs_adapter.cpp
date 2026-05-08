#include "gks_kfvs_adapter.h"

#include "gks_basic.h"

namespace
{
	void BuildCenterState(
		const double left_Q[3],
		const double right_Q[3],
		double center_Q[3])
	{
		double left_prim[3];
		double right_prim[3];
		double center_prim[3];
		Convar_to_primvar_1D(left_prim, const_cast<double*>(left_Q));
		Convar_to_primvar_1D(right_prim, const_cast<double*>(right_Q));

		// The existing KFVS branch in GKS(...) only uses left/right moments for
		// the returned flux, but it still constructs a center state internally.
		// Use a simple physical primitive average here to keep the wrapper
		// independent of the original large-cell data layout.
		center_prim[0] = 0.5 * (left_prim[0] + right_prim[0]);
		center_prim[1] = 0.5 * (left_prim[1] + right_prim[1]);
		center_prim[2] = 0.5 * (left_prim[2] + right_prim[2]);
		Primvar_to_convar_1D(center_Q, center_prim);
	}

	void BuildCenterState2D(
		const double left_Q[4],
		const double right_Q[4],
		double center_Q[4])
	{
		double left_prim[4];
		double right_prim[4];
		double center_prim[4];
		double left_local[4];
		double right_local[4];
		for (int m = 0; m < 4; ++m)
		{
			left_local[m] = left_Q[m];
			right_local[m] = right_Q[m];
		}
		Convar_to_primvar_2D(left_prim, left_local);
		Convar_to_primvar_2D(right_prim, right_local);
		for (int m = 0; m < 4; ++m)
		{
			center_prim[m] = 0.5 * (left_prim[m] + right_prim[m]);
		}
		Primvar_to_convar_2D(center_Q, center_prim);
	}

	void ZeroPoint2D(Point2d& point)
	{
		for (int m = 0; m < 4; ++m)
		{
			point.convar[m] = 0.0;
			point.der1x[m] = 0.0;
			point.der1y[m] = 0.0;
			point.der2xx[m] = 0.0;
			point.der2xy[m] = 0.0;
			point.der2yy[m] = 0.0;
		}
		point.x = 0.0;
		point.y = 0.0;
		point.normal[0] = 1.0;
		point.normal[1] = 0.0;
		point.is_reduce_order = false;
	}

	void RotateGlobalToYLocal(double dst[4], const double src[4])
	{
		dst[0] = src[0];
		dst[1] = src[2];
		dst[2] = -src[1];
		dst[3] = src[3];
	}

	void RotateYLocalFluxToGlobal(double dst[4], const double src[4])
	{
		dst[0] = src[0];
		dst[1] = -src[2];
		dst[2] = src[1];
		dst[3] = src[3];
	}
}

void KFVS1_TimeIntegratedFlux1D(
	const double left_Q[3],
	const double right_Q[3],
	double dt,
	double flux_int[3])
{
	Interface1d interface;
	Flux1d flux;
	double center_Q[3];
	BuildCenterState(left_Q, right_Q, center_Q);

	for (int m = 0; m < 3; ++m)
	{
		interface.left.convar[m] = left_Q[m];
		interface.right.convar[m] = right_Q[m];
		interface.center.convar[m] = center_Q[m];
		interface.left.der1[m] = 0.0;
		interface.right.der1[m] = 0.0;
		interface.center.der1[m] = 0.0;
		interface.left.der2[m] = 0.0;
		interface.right.der2[m] = 0.0;
		interface.center.der2[m] = 0.0;
		flux.F[m] = 0.0;
		flux.f[m] = 0.0;
		flux.derf[m] = 0.0;
		flux.der2f[m] = 0.0;
	}

	const GKS1d_type old_solver = gks1dsolver;
	gks1dsolver = kfvs1st;
	GKS(flux, interface, dt);
	gks1dsolver = old_solver;

	for (int m = 0; m < 3; ++m)
	{
		flux_int[m] = flux.f[m];
	}
}

void KFVS1_TimeAveragedFlux1D(
	const double left_Q[3],
	const double right_Q[3],
	double dt,
	double flux_avg[3])
{
	double flux_int[3];
	KFVS1_TimeIntegratedFlux1D(left_Q, right_Q, dt, flux_int);
	for (int m = 0; m < 3; ++m)
	{
		flux_avg[m] = flux_int[m] / dt;
	}
}

void KFVS1_TimeIntegratedFlux2D_X(
	const double left_Q[4],
	const double right_Q[4],
	double dt,
	double flux_int[4])
{
	Recon2d interface;
	Flux2d flux;
	double center_Q[4];
	BuildCenterState2D(left_Q, right_Q, center_Q);
	ZeroPoint2D(interface.left);
	ZeroPoint2D(interface.right);
	ZeroPoint2D(interface.center);
	for (int m = 0; m < 4; ++m)
	{
		interface.left.convar[m] = left_Q[m];
		interface.right.convar[m] = right_Q[m];
		interface.center.convar[m] = center_Q[m];
		flux.f[m] = 0.0;
		flux.derf[m] = 0.0;
		flux.der2f[m] = 0.0;
	}
	interface.normal[0] = 1.0;
	interface.normal[1] = 0.0;
	flux.normal[0] = 1.0;
	flux.normal[1] = 0.0;

	const GKS2d_type old_solver = gks2dsolver;
	gks2dsolver = kfvs1st_2d;
	GKS2D(flux, interface, dt);
	gks2dsolver = old_solver;
	for (int m = 0; m < 4; ++m)
	{
		flux_int[m] = flux.f[m];
	}
}

void KFVS1_TimeAveragedFlux2D_X(
	const double left_Q[4],
	const double right_Q[4],
	double dt,
	double flux_avg[4])
{
	double flux_int[4];
	KFVS1_TimeIntegratedFlux2D_X(left_Q, right_Q, dt, flux_int);
	for (int m = 0; m < 4; ++m)
	{
		flux_avg[m] = flux_int[m] / dt;
	}
}

void KFVS1_TimeIntegratedFlux2D_Y(
	const double bottom_Q[4],
	const double top_Q[4],
	double dt,
	double flux_int[4])
{
	double bottom_local[4];
	double top_local[4];
	double local_flux[4];
	RotateGlobalToYLocal(bottom_local, bottom_Q);
	RotateGlobalToYLocal(top_local, top_Q);
	KFVS1_TimeIntegratedFlux2D_X(bottom_local, top_local, dt, local_flux);
	RotateYLocalFluxToGlobal(flux_int, local_flux);
}

void KFVS1_TimeAveragedFlux2D_Y(
	const double bottom_Q[4],
	const double top_Q[4],
	double dt,
	double flux_avg[4])
{
	double flux_int[4];
	KFVS1_TimeIntegratedFlux2D_Y(bottom_Q, top_Q, dt, flux_int);
	for (int m = 0; m < 4; ++m)
	{
		flux_avg[m] = flux_int[m] / dt;
	}
}
