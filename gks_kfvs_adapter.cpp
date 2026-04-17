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
