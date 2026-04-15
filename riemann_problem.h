#pragma once
#include"fvm_gks2d.h"
#include"output.h"
void riemann_problem_1d();
void ICfor1dRM(Fluid1d* fluids, Fluid1d zone1, Fluid1d zone2, Block1d block);
void riemann_problem_2d();
void IC_for_riemann_2d(Fluid2d* fluid, double* zone1, double* zone2, double* zone3, double* zone4, Block2d block);

