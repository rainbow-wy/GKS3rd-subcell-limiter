#pragma once
#include"fvm_gks2d.h"
#include"output.h"
using namespace std;
string add_mesh_directory_modify_for_linux();
void Read_structured_mesh(string& filename,
	Fluid2d** fluids, Interface2d** xinterfaces, Interface2d** yinterfaces,
	Flux2d_gauss*** xfluxes, Flux2d_gauss*** yfluxes, Block2d& block);
void read_ascii_to_string(string& value, ifstream& in);
void create_mirror_ghost_cell(Block2d& block, double** node2d);
