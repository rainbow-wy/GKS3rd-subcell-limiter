#pragma once
#include"fvm_gks2d.h"
#include"output.h"
void accuracy_sinwave_1d();
void accuracy_sinwave_1d(double& CFL, double& dt_ratio, int& mesh_number, double* error);
void ICfor_sinwave(Fluid1d* fluids, Block1d block);
void error_for_sinwave(Fluid1d* fluids, Block1d block, double tstop, double* error);
void accuracy_sinwave_2d();
void sinwave_2d(double &CFL, double &dt_ratio, int& mesh_number, double *error);
void ICfor_sinwave_2d(Fluid2d *fluids, Block2d block);
void error_for_sinwave_2d(Fluid2d *fluids, Block2d block, double tstop,double *error);
