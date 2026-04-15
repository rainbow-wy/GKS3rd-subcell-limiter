#pragma once
#include"fvm_gks2d.h"
#include"output.h"
#include"mesh_read.h"
void boundary_layer();
void boundaryforBoundary_layer(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void output_boundary_layer(Fluid2d* fluids, Block2d block, double* bcvalue, double coordx, char* label);
void output_wall_info_boundary_layer(Interface2d* xInterfaces, Interface2d* yInterfaces,
	Flux2d_gauss** xFluxes, Flux2d_gauss** yFluxes, Block2d block, double L_ref);
void output_skin_friction_boundary_layer
(Fluid2d* fluids, Interface2d* yInterfaces, Flux2d_gauss** yFluxes, Block2d block, double u_ref, double re_ref, double L_ref);