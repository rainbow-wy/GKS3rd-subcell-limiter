#pragma once
#include"function.h"
#include<omp.h>
#include"gks_basic.h"
//extern type variables must be initialized (could not give value) in the corresponding .cpp file
enum GKS1d_type{nothing, kfvs1st, kfvs2nd, gks1st, gks2nd, gks3rd};
extern GKS1d_type gks1dsolver;
enum Reconstruction_variable{conservative,characteristic};
extern Reconstruction_variable reconstruction_variable;
enum WENOtype { linear, wenojs, wenoz };
extern WENOtype wenotype;
extern bool is_reduce_order_warning;
extern bool is_show_1d_timestep;

using namespace std;

// to remember the geometry information
class Block1d
{
public:
	//bool uniform;
	int ghost; //虚拟单元
	int nodex; // mesh number
	int nx;
	int nodex_begin;
	int nodex_end;
	double dx;
	double left;
	double right;
	int stages;
	double timecoefficient[5][5][3];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	int step; //start from which step
};

// to remeber the cell avg values
class Fluid1d
{
public:
	double primvar[3];
	double convar[3];
	double convar_old[3];
	double cx; //center coordinate in x direction
	double dx; //the mesh size dx
};

// to remember the fluid information in a fixed point,
// such as reconstructioned value, or middle point value
class Point1d
{
public:
	double convar[3];
	double convar_old[3]; //this one is for compact point-wise reconstruction
	double der1[3];
	double der2[3]; // second derivative in x, used by 1-stage 3rd-order GKS
	double x; // coordinate of a point
};

// remember the flux in a fixed interface,
// for RK method, we only need f
// for 2nd der method, we need derf
class Flux1d
{
public:
	double F[3]; //total flux in dt time
	double f[3]; // the f0 in t=0
	double derf[4]; // the f_t in t=0
	double der2f[3]; // the f_tt in t=0, reserved for 3rd-order GKS diagnostics
};

// every interfaces have left center and right value.
// and center value for GKS
// and a group of flux for RK method or multi-stage GKS method
class Interface1d
{
public:
	Point1d left;
	Point1d center;
	Point1d right;
	Flux1d* flux;
	double x; // coordinate of the interface, equal to point1d.x
};

// 单元重构多项式，正交 Legendre 基
// P(ξ) = q + A·ξ + B·(ξ²-1/12) + C·(ξ³-3ξ/20) + D·(ξ⁴-3ξ²/14+3/560)
// ξ = (x - x_center) / dx ∈ [-0.5, 0.5]
struct Cellpoly1d
{
	double q[3];
	double A[3];
	double B[3];
	double C[3];
	double D[3];

	void reset()
	{
		for(int m = 0; m < 3; m++)
		{
			q[m]=0.0;
			A[m]=0.0;
			B[m]=0.0;
			C[m]=0.0;
			D[m]=0.0;
		}
	}

	// 在局部坐标 xi 处求值, xi ∈ [-0.5, 0.5]
	double value(int m, double xi) const
	{
		const double xi2 = xi * xi;
		const double xi3 = xi2 * xi;
		const double xi4 = xi2 * xi2;
		return q[m]
		     + A[m] * xi
		     + B[m] * (xi2 - 1.0 / 12.0)
		     + C[m] * (xi3 - 3.0 * xi / 20.0)
		     + D[m] * (xi4 - 3.0 * xi2 / 14.0 + 3.0 / 560.0);
	}

	// 对物理坐标 x 的一阶导数
	double dxValue(int m, double xi, double dx) const
	{
		const double xi2 = xi * xi;
		return ( A[m]
			+ 2.0 * B[m] * xi
			+ C[m] * (3.0 * xi2 - 3.0 / 20.0)
			+ D[m] * (4.0 * xi2 * xi - 3.0 * xi / 7.0) ) / dx;
	}

	// 对物理坐标 x 的二阶导数
	double dxxValue(int m, double xi, double dx) const
	{
		const double xi2 = xi * xi;
		return ( 2.0 * B[m]
			+ 6.0 * C[m] * xi
			+ D[m] * (12.0 * xi2 - 3.0 / 7.0) ) / (dx * dx);
	}
};

void Convar_to_primvar(Fluid1d *fluids, Block1d block);

typedef void(*BoundaryCondition) (Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void free_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void reflection_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void reflection_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void periodic_boundary_left(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);
void periodic_boundary_right(Fluid1d *fluids, Block1d block, Fluid1d bcvalue);


void Reconstruction_within_cell(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void(*Reconstruction_within_Cell)(Point1d &left, Point1d &right, Fluid1d *fluids);
extern Reconstruction_within_Cell cellreconstruction;

void Check_Order_Reduce(Point1d &left, Point1d &right, Fluid1d &fluid);
void Check_Order_Reduce_by_Lambda_1D(bool &order_reduce, double *convar);

void Vanleer(Point1d &left, Point1d &right, Fluid1d *fluids);
void weno_3rd_left(double& var, double& der1, double& der2,
	double wn1, double w0, double wp1, double h);
void weno_3rd_right(double& var, double& der1, double& der2,
	double wn1, double w0, double wp1, double h);
void WENO3_left(double& left, double w_back, double w, double w_forward, double h);
void WENO3_right(double& right, double w_back, double w, double w_forward, double h);

void WENO5_AO(Point1d &left, Point1d &right, Fluid1d *fluids);
void weno_5th_ao_left(double &var, double &der1, double &der2, double wn2, double wn1, double w0, double wp1, double wp2, double h);
void weno_5th_ao_right(double &var, double &der1, double &der2,	double wn2, double wn1, double w0, double wp1, double wp2,	double h);
void WENO5_AO_poly(Cellpoly1d& poly, Fluid1d* fluids);
void verify_WENO5_AO_poly(Fluid1d* fluids, Block1d block);


void Reconstruction_forg0(Interface1d *interfaces, Fluid1d *fluids, Block1d block);
typedef void (*Reconstruction_forG0)(Interface1d &interfaces, Fluid1d *fluids);
extern Reconstruction_forG0 g0reconstruction;
void Center_2nd_collisionless(Interface1d& interfaces, Fluid1d* fluids);
void Center_3rd(Interface1d& interfaces, Fluid1d* fluids);
void Center_GKS3rd(Interface1d& interfaces, Fluid1d* fluids);
void Center_simple_avg(Interface1d &interfaces, Fluid1d *fluids);
void Center_collision(Interface1d &interfaces, Fluid1d *fluids);


void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d &block, int stage);
typedef void(*Flux_function)(Flux1d &flux, Interface1d& interface, double dt);
extern Flux_function flux_function;
void GKS(Flux1d &flux, Interface1d& interface, double dt);


typedef void(*TimeMarchingCoefficient)(Block1d &block);
extern TimeMarchingCoefficient timecoe_list;
void S1O1(Block1d& block);
void S1O2(Block1d& block);
void S1O3(Block1d& block);
void S2O4(Block1d& block);
void RK4(Block1d &block);

void Initial_stages(Block1d &block);
Flux1d** Setflux_array(Block1d block);
void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d *interfaces, Flux1d **fluxes);


void CopyFluid_new_to_old(Fluid1d *fluids, Block1d block);
double Get_CFL(Block1d &block, Fluid1d *fluids, double tstop);
double Dtx(double dtx, double dx, double CFL, double convar[3]);


void Update(Fluid1d *fluids, Flux1d **fluxes, Block1d block, int stage);


// two dimensional hllc and its subprogram
void ESTIME(double& SL, double& SM, double& SR, double DL, double UL, double PL, double CL, double DR, double UR, double PR, double CR);
void get_flux(double p[4], double* flux);
void ustarforHLLC(double d1, double u1, double v1, double p1, double s1, double star1, double* ustar);
void HLLC(Flux1d& flux, Interface1d& interface, double dt);
void get_Euler_flux(double p[3], double* flux);
void ustarforHLLC(double d1, double u1, double p1, double s1, double star1, double* ustar);

void LF(Flux1d& flux, Interface1d& interface, double dt);
