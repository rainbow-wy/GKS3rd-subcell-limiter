#pragma once
#include"fvm_gks2d.h"
#include"output.h"
#include"mesh_read.h"
void cylinder();
void boundary_for_cylinder(Fluid2d* fluids, Interface2d* faces, Block2d block, Fluid2d bcvalue);