#pragma once
#include"fvm_gks1d.h"
extern int gausspoint;
extern double *gauss_loc;
extern double *gauss_weight;
enum GKS2d_type { nothing_2d, kfvs1st_2d, kfvs2nd_2d, gks1st_2d, gks2nd_2d, gks3rd_2d};
extern GKS2d_type gks2dsolver;



class Block2d
{
	//geometry information of 2D block
public:
	int index;
	bool uniform;
	int ghost;
	int nodex;
	int nodey;
	int nx; //total = real + ghost
	int ny; //total = real + ghost
	int xcell_begin;
	int xcell_end;
	int ycell_begin;
	int ycell_end;
	int xinterface_begin_n;
	int xinterface_begin_t;
	int xinterface_end_n;
	int xinterface_end_t;
	int yinterface_begin_n;
	int yinterface_begin_t;
	int yinterface_end_n;
	int yinterface_end_t;
	double dx;
	double dy;
	double overdx;
	double overdy;
	double left;
	double right;
	double down;
	double up;
	int stages;
	double timecoefficient[5][5][3];
	double t; //current simulation time
	double CFL; //cfl number, actually just a amplitude factor
	double dt;//current_global time size
	int step; //start from which step
	double T; //stop time
	int outputperstep;
};

// to remember the fluid information in a fixed point,
// such as reconstructioned value, or middle point value
class Point2d
{
	// Note : store values relating to one point
public:
	double convar[4];
	double der1x[4];
	double der1y[4];
	double der2xx[4];
	double der2xy[4];
	double der2yy[4];
	double x;
	double y;
	double normal[2];
	bool is_reduce_order;
};

// to remeber the cell avg values
//typedef class Fluid2d
class Fluid2d
{
	// Note : store the cell-averaged values
	// do not have the global geometry/block information
public:
	double primvar[4];
	double convar[4];
	double convar_old[4];
	double res[4];
	double xindex;
	double yindex;
	double coordx;
	double coordy;
	double dx;
	double dy;
	double area;
	bool is_hweno;
	double node[8];
	bool boundary;
};


// remember the flux in a fixed interface,
// for RK method, we only need f
// for 2nd der method, we need derf
class Flux2d
{
public:
	double x[4];
	double f[4];
	double derf[4];
	double der2f[4];
	double length;
	double normal[2];
};

// every interfaces have left center and right value.
// and center value for GKS
// and a group of flux for RK method or multi-stage GKS method
class Recon2d
{
	// Note : reconstruction can obtain left, cenrer, and right value
public:
	Point2d left;
	Point2d center;
	Point2d right;
	double x;
	double y;
	double normal[2];
};
class Interface2d
{
	// Note : for one interface, reconstruction can be for line-averagend or gauss point
	// These reconstructed values can be used for the calculation of the flux at gauss points (stored in Flux2d_gauss) of the interface
public:

	Recon2d line;
	Recon2d* gauss;
	double x;
	double y;
	double normal[2];
	double length;
};

class Flux2d_gauss
{
	// Note : Flux2d store the variables relating to flux
	// Flux2d_gauss, includes several gauss points
	// Final, one interface should have one Flux2d_gauss
public:
	Flux2d* gauss;
};

void SetGuassPoint();
double Get_CFL(Block2d &block, Fluid2d *fluids, double tstop);



void A_point(double *a, double* der, double* prim);
void G_address(int no_u, int no_v, int no_xi, double *psi, double a[4], MMDF& m);
void GL_address(int no_u, int no_v, int no_xi, double *psi, double a[4], MMDF& m);
void GR_address(int no_u, int no_v, int no_xi, double *psi, double a[4], MMDF& m);



typedef void(*BoundaryCondition2dnew) (Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void free_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_left(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_right(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_down(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void periodic_boundary_up(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_left(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary(double* w_ghost, double* w_inner, double* normal, double* local_wall_velocity);
void noslip_adiabatic_boundary_right(Fluid2d* fluids, Interface2d* face, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_down(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void noslip_adiabatic_boundary_up(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_up(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_down(Fluid2d *fluids, Block2d block, Fluid2d bcvalue);
void reflection_boundary_right(Fluid2d* fluids, Interface2d* face, Block2d block, Fluid2d bcvalue);
void reflection_boundary_left(Fluid2d *fluids, Interface2d* face, Block2d block, Fluid2d bcvalue);
void reflection_boundary(double* ghost, double* inner, double* normal);
void inflow_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void inflow_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void inflow_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue);
void ICfor_uniform_2d(Fluid2d* fluid, double* prim, Block2d block);

bool negative_density_or_pressure(double* primvar);
void Convar_to_Primvar(Fluid2d *fluids, Block2d block);
void Convar_to_primvar(Fluid2d *fluid, Block2d block);
void YchangetoX(double *fluidtmp, double *fluid);
void XchangetoY(double *fluidtmp, double *fluid);
void CopyFluid_new_to_old(Fluid2d *fluids, Block2d block);



void Reconstruction_within_cell(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
typedef void(*Reconstruction_within_Cell_2D_normal)(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
extern Reconstruction_within_Cell_2D_normal cellreconstruction_2D_normal;
void First_order_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void first_order(Point2d& left, Point2d& right, double* normal_l, double* normal_r, double* w);
void Vanleer_normal(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
void Vanleer_non_uniform(Point2d& left, Point2d& right, double* normal_l, double* normal_r, double* wn1, double* w, double* wp1, double h);
void Vanleer(Point2d &left, Point2d &right, double *wn1, double *w, double *wp1, double h);
void WENO3_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void WENO3(Point2d& left, Point2d& right, double* wn1, double* w, double* wp1, double h);
void WENO(Point2d& left, Point2d& right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h);
void WENO5_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block);
void WENO5_left(double& var, double& der1, double& der2, double wn2, double wn1, double w, double wp1, double wp2, double h);
void WENO5_right(double& var, double& der1, double& der2, double wn2, double wn1, double w, double wp1, double wp2, double h);

void WENO5_AO_normal(Interface2d &left, Interface2d &right, Interface2d &down, Interface2d &up, Fluid2d *fluids, Block2d block);
void WENO5_AO(Point2d &left, Point2d &right, double * wn2, double * wn1, double * w, double * wp1, double * wp2, double h);



typedef void(*Reconstruction_within_Cell_2D_tangent)(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
extern Reconstruction_within_Cell_2D_tangent cellreconstruction_2D_tangent;
void First_order_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);
void first_order_tangent(Point2d& gauss, Point2d& w0);
void Vanleer_tangent(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
void Vanleer(Point2d& gauss, Point2d &wn1, Point2d &w0, Point2d &wp1, double h);
void Vanleer_tangential(double *coe, double wn1, double w0, double wp1, double h);

void Polynomial_3rd(double* coefficent, double pn1, double w0, double pp1, double h);

void Polynomial_5th(double* coefficent, double wn2, double wn1, double w0, double wp1, double wp2, double h);
double Value_Polynomial(int order, double x0, double* coefficient);
double Der1_Polynomial(int order, double x0, double* coefficient);

void Polynomial_3rd_avg(double* coefficent, double wn1, double w0, double wp1, double h);
double Value_Polynomial_3rd(double x, double coefficient[3]);
double Der1_Polynomial_3rd(double x, double coefficient[3]);


void WENO5_left(double& var, double wn2, double wn1, double w, double wp1, double wp2, double h);
void WENO5_right(double& var, double wn2, double wn1, double w, double wp1, double wp2, double h);
void WENO5_tangential(double* left, double* right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h);
void WENO5_tangential_for_slope(double* left, double* right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h);

void WENO5_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double h);
void WENO5_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);

void WENO5_AO_tangent(Interface2d *left, Interface2d *right, Interface2d *down, Interface2d *up, Fluid2d *fluids, Block2d block);
void WENO5_AO_tangential(Recon2d *re, Recon2d &wn2, Recon2d &wn1, Recon2d &w0, Recon2d &wp1, Recon2d &wp2, double h);
void weno_5th_ao_2gauss(double &g1, double &g1x, double &g1xx, double &g2, double &g2x, double &g2xx, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order);
void WENO3_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block);
void WENO3_tangential(Recon2d* re, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, double h);
void WENO3_tangential(double* left, double* right, double* wn1, double* w, double* wp1, double h);
void WENO3_tangential_for_slope(double* left, double* right, double* wn1, double* w, double* wp1, double h);

void Reconstruction_forg0(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
typedef void(*Reconstruction_forG0_2D_normal)(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
extern Reconstruction_forG0_2D_normal g0reconstruction_2D_normal;
void Center_3rd_normal(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
void Center_3rd_Splitting(Interface2d& interface, double* w, double* wp1, double h);

void Center_collision_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);
void Center_do_nothing_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);

void Center_4th_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);
void Center_4th_normal(Interface2d& interface, double* wn1, double* w, double* wp1, double* wp2, double h);

void Center_5th_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);
void Center_5th_normal(Interface2d& interface, double* wn1, double* w, double* wp1, double* wp2, double h);



typedef void(*Reconstruction_forG0_2D_tangent)(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
extern Reconstruction_forG0_2D_tangent g0reconstruction_2D_tangent;
void Center_all_collision_multi(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
void Center_all_collision_2d_multi(Recon2d &gauss);
void Center_first_order_tangent(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block);
void Center_3rd_tangent(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
void Center_3rd_Multi(Point2d& gauss, Point2d& wn1, Point2d& w0, Point2d& wp1, double h);

void Center_5th_tangent(Interface2d *xinterfaces, Interface2d *yinterfaces, Fluid2d *fluids, Block2d block);
void Center_5th_Multi(Recon2d *re, Recon2d &wn2, Recon2d &wn1, Recon2d &w0, Recon2d &wp1, Recon2d &wp2, double h);




void Calculate_flux(Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Interface2d* xinterfaces, Interface2d* yinterfaces, Block2d block, int stage);
typedef void(*Flux_function_2d)(Flux2d &flux, Recon2d & re, double dt);
extern Flux_function_2d flux_function_2d;
void GKS2D_smooth(Flux2d &flux, Recon2d& interface, double dt);
void GKS2D(Flux2d &flux, Recon2d& interface, double dt);
void LF2D(Flux2d &flux, Recon2d& interface, double dt);
void HLLC2D(Flux2d& flux, Recon2d& interface, double dt);
void NS_by_central_difference_prim_2D(Flux2d& flux, Recon2d& interface, double dt);
void NS_by_central_difference_convar_2D(Flux2d& flux, Recon2d& interface, double dt);



typedef void(*TimeMarchingCoefficient_2d)(Block2d &block);
extern TimeMarchingCoefficient_2d timecoe_list_2d;


void S1O1_2D(Block2d &block);
void S1O2_2D(Block2d& block);
void S1O3_2D(Block2d& block);
void RK2_2D (Block2d& block);
void S2O4_2D(Block2d &block);
void RK4_2D(Block2d& block);
void Initial_stages(Block2d &block);



Fluid2d *Setfluid(Block2d &block);

Interface2d *Setinterface_array(Block2d block);

Flux2d_gauss** Setflux_gauss_array(Block2d block);



void SetUniformMesh(Block2d& block, Fluid2d* fluids, Interface2d *xinterfaces, Interface2d *yinterfaces, Flux2d_gauss **xfluxes, Flux2d_gauss **yfluxes);
void Copy_geo_from_interface_to_line(Interface2d& interface);
void Copy_geo_from_interface_to_flux(Interface2d& interface, Flux2d_gauss* flux, int stages);
void Set_Gauss_Coordinate(Interface2d &xinterface, Interface2d &yinterface );
void Set_Gauss_Intergation_Location_x(Point2d &xgauss, int index, double h);
void Set_Gauss_Intergation_Location_y(Point2d &ygauss, int index, double h);





void Update(Fluid2d *fluids, Flux2d_gauss **xfluxes, Flux2d_gauss **yfluxes, Block2d block, int stage);
void Update_with_gauss(Fluid2d *fluids, Flux2d_gauss **xfluxes, Flux2d_gauss **yfluxes, Block2d block, int stage);

void Check_Order_Reduce_by_Lambda_2D(bool &order_reduce, double *convar);
void Recompute_KFVS_1st(Fluid2d &fluid,	double* center, double* left, double* right, double* down, double* up, Block2d& block);
void KFVS_1st(double *flux, double* left, double* right, double dt);

void SetNonUniformMesh(Block2d& block, Fluid2d* fluids,
	Interface2d* xinterfaces, Interface2d* yinterfaces,
	Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, double** node2d);
void Set_cell_geo_from_quad_node
(Fluid2d& fluid, double* n1, double* n2, double* n3, double* n4);
void Set_interface_geo_from_two_node
(Interface2d& interface, double* n1, double* n2, int direction);
void Set_Gauss_Coordinate_general_mesh_x
(Interface2d& xinterface, double* node0, double* nodex1);
void Set_Gauss_Coordinate_general_mesh_y
(Interface2d& yinterface, double* node0, double* nodey1);

void Residual2d(Fluid2d* fluids, Block2d block, int outputstep);
