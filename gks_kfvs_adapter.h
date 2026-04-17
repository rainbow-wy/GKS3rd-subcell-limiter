#pragma once

#include "fvm_gks1d.h"

// Wrap the existing 1D KFVS implementation in fvm_gks1d.cpp so it can be
// reused on arbitrary subcell faces. The returned flux is the time-averaged
// flux \bar F, which matches the FR/subcell framework update formulas.
void KFVS1_TimeAveragedFlux1D(
	const double left_Q[3],
	const double right_Q[3],
	double dt,
	double flux_avg[3]);

void KFVS1_TimeIntegratedFlux1D(
	const double left_Q[3],
	const double right_Q[3],
	double dt,
	double flux_int[3]);
