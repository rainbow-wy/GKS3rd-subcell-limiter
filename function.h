#pragma once
#include<iostream>
#include<cmath>
#include<assert.h>
using namespace std;

extern int K;
extern double Gamma;
extern double Mu;
extern double Nu;
extern double T_inf; //the real temperature
extern double R_gas; //the gas constant for real ideal air
enum TAU_TYPE { Euler, NS, Sutherland, Power_law };
extern TAU_TYPE tau_type;


void Copy_Array(double* target, double* origin, int dim);
void Array_zero(double* target, int dim);
void Array_subtract(double* target, double* subtractor, int dim);
void Array_scale(double* target, int dim, double value);

void Convar_to_primvar_1D(double* primvar, double convar[3]);
double Pressure(double density, double densityu, double densityE);
double Temperature(double density, double pressure);
double entropy(double density, double pressure);
double Soundspeed(double density, double pressure);
void Primvar_to_convar_1D(double* convar, double primvar[3]);
double DensityU(double density, double u);
double DensityE(double density, double u, double pressure);
void Convar_to_char1D(double* character, double primvar[3], double convar[3]);
void Char_to_convar1D(double* convar, double primvar[3], double charvar[3]);

double U(double density, double q_densityu);
double V(double density, double q_densityv);
double Pressure(double density, double q_densityu, double q_densityv, double q_densityE);

double Q_densityu(double density, double u);
double Q_densityv(double density, double v);
double Q_densityE(double density, double u, double v, double pressure);


void Primvar_to_convar_2D(double* convar, double primvar[4]);//2d version
void Convar_to_primvar_2D(double* primvar, double convar[4]);//2d version

void Convar_to_char(double* character, double* primvar, double convar[4]);
void Char_to_convar(double* convar, double* primvar, double character[4]);

void Global_to_Local(double* change, double* normal);
void Global_to_Local(double* change, double* origin, double* normal);
void Local_to_Global(double* change, double* normal);
void Local_to_Global(double* change, double* origin, double* normal);

double Dtx(double dtx, double dx, double CFL, double density, double u, double v, double pressure);
