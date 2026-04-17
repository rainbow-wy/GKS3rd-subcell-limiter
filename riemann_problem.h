#pragma once
#include"fvm_gks2d.h"
#include"output.h"
#include "gks_fr_high_order.h"

struct RiemannProblem1D
{
	double x_discontinuity;
	double left_prim[3];
	double right_prim[3];
};

RiemannProblem1D RiemannProblem1D_Sod();
RiemannProblem1D RiemannProblem1D_DoubleRarefaction();

void riemann_problem_1d();
void ICfor1dRM(Fluid1d* fluids, Fluid1d zone1, Fluid1d zone2, Block1d block);
void ICfor1dRM(Fluid1d* fluids, const RiemannProblem1D& problem, Block1d block);
void ICfor1dRM(GKSFRMesh1D& mesh, const RiemannProblem1D& problem);
void riemann_problem_2d();
void IC_for_riemann_2d(Fluid2d* fluid, double* zone1, double* zone2, double* zone3, double* zone4, Block2d block);
