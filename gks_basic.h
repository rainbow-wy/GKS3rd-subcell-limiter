#pragma once
#include"function.h"


extern double c1_euler;
extern double c2_euler;
extern bool is_Prandtl_fix;
extern double Pr;



//basic gks function
// to store the moment
class MMDF1d
{
private:
	double u;
	double lambda;

public:
	double uwhole[10];
	double uplus[10];
	double uminus[10];
	double upxi[10][4];
	double unxi[10][4];
	double uxi[10][4];
	double xi2;
	double xi4;
	double xi6;
	MMDF1d();
	MMDF1d(double u_in, double lambda_in);
	void calcualte_MMDF1d();
};

// to calculate the microsolpe moment
void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);
void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m);

//moments of the maxwellian distribution function
class MMDF
{
private:
	double u;
	double v;
	double lambda;

public:
	double uwhole[7];
	double uplus[7];
	double uminus[7];
	double vwhole[7];
	double upvxi[7][7][3];
	double unvxi[7][7][3];
	double uvxi[7][7][3];
	double xi2;
	double xi4;
	MMDF();
	MMDF(double u_in, double v_in, double lambda_in);
	void calcualte_MMDF();
};

class MMDF1st
{
private:
	double u;
	double v;
	double lambda;

public:
	double uwhole[4];
	double uplus[4];
	double uminus[4];
	double vwhole[3];
	double xi2;
	MMDF1st();
	MMDF1st(double u_in, double v_in, double lambda_in);
	void calcualte_MMDF1st();
};


void Collision(double *w0, double left, double right, MMDF &m2, MMDF &m3);
void Collision(double *w0, double left, double right, MMDF1st &m2, MMDF1st &m3);

void A(double *a, double der[4], double prim[4]);

double Get_Tau_NS(double density0, double lambda0); // solve the smooth tau
double TauNS_Sutherland(double density0, double lambda0); // solve the smooth tau by using sutherland
double TauNS_power_law(double density0, double lambda0); //solver the smooth tau by using power law
double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt); // solve not smooth numerical tau


double Lambda(double density, double u, double densityE);
double Lambda(double density, double u, double v, double densityE);
void Convar_to_ULambda_1d(double* primvar, double convar[3]);

//计算α，α=erfc（±√λ*U）
double Alpha(double lambda, double u);
//计算β，β=e^(-λ*U^2)/(√πλ)
double Beta(double lambda, double u);
//solution of matrix equation b=Ma
void Microslope(double *a, double der[3], double prim[3]); //in one dimensional
void SolveMicroCoeff1D(double* coeff, double moment[3], double prim[3]);
void Moment_uCoeff1D(double* moment, double coeff[3], double prim[3]);
void SolveGKS_a1D(double* a, double Wx[3], double prim[3]);
void SolveGKS_A1D(double* A, double a[3], double prim[3]);
void SolveGKS_axx1D(double* axx, double Wxx[3], double prim[3]);
void SolveGKS_axt1D(double* axt, double axx[3], double prim[3]);
void SolveGKS_att1D(double* att, double axt[3], double prim[3]);

void Convar_to_ULambda_2d(double* primvar, double convar[4]);
