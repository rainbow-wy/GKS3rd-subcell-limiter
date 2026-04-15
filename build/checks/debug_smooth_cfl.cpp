#include "accuracy_test.h"
#include <iostream>
#include <cmath>
int main(){
    K=4; Gamma=1.4; R_gas=1.0; tau_type=Euler; c1_euler=0.0; c2_euler=0.0;
    flux_function=GKS; gks1dsolver=gks3rd; cellreconstruction=WENO5_AO; reconstruction_variable=characteristic; wenotype=wenoz; g0reconstruction=Center_GKS3rd; timecoe_list=S1O3;
    Block1d block; block.nodex=40; block.ghost=3; block.nx=block.nodex+2*block.ghost; block.nodex_begin=block.ghost; block.nodex_end=block.nodex+block.ghost-1; block.left=0.0; block.right=2.0; block.dx=(block.right-block.left)/block.nodex; block.CFL=0.1; Initial_stages(block);
    Fluid1d* fluids = new Fluid1d[block.nx]; Interface1d* interfaces = new Interface1d[block.nx+1]; Flux1d** fluxes=Setflux_array(block); SetUniformMesh(block, fluids, interfaces, fluxes); ICfor_sinwave(fluids, block); Fluid1d* bcvalue = new Fluid1d[2];
    block.t=0; block.step=0; bool has_nan=false; while(block.t < 2.0-1e-14){ CopyFluid_new_to_old(fluids, block); block.dt=Get_CFL(block, fluids, 2.0); periodic_boundary_left(fluids, block, bcvalue[0]); periodic_boundary_right(fluids, block, bcvalue[1]); Reconstruction_within_cell(interfaces, fluids, block); Reconstruction_forg0(interfaces, fluids, block); Calculate_flux(fluxes, interfaces, block, 0); Update(fluids, fluxes, block, 0); block.t += block.dt; block.step++; for(int i=block.ghost;i<block.nx-block.ghost;i++){ if(!(fluids[i].convar[0]==fluids[i].convar[0])) { has_nan=true; break; } } if(has_nan) break; }
    std::cout << "step=" << block.step << " t=" << block.t << " has_nan=" << has_nan << std::endl; return 0; }
