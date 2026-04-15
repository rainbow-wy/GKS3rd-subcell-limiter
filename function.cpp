#include"function.h"

int K = 0; //initialization  //双原子，r=1.4，K=（4-2r）/（r-1）=3
double Gamma = 2.0 / (K + 2) + 1.0; //initialization

double T_inf = 285; //the real temperature, 20 Celcius degree
double R_gas = 8.31441 / (28.959e-3); //initialization //the gas constant for real ideal air
double Mu = -1.0; // mu = rho*u*L/Re
double Nu = -1.0; // Nu = u*L/Re
TAU_TYPE tau_type=Euler; //无粘还是粘性问题
//tau_type, euler==0,ns==1,sutherland==2, Power_law==3

void Copy_Array(double* target, double* origin, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		target[i] = origin[i];
	}
}
void Array_zero(double* target, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		target[i] = 0.0;
	}
}
void Array_subtract(double* target, double* subtractor, int dim)
{
	for (int i = 0; i < dim; i++)
	{
		target[i] -= subtractor[i];
	}
}

void Array_scale(double* target, int dim, double value)
{
	for (int i = 0; i < dim; i++)
	{
		target[i] *= value;
	}
}

//从守恒变量计算原始变量
void Convar_to_primvar_1D(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Pressure(convar[0], convar[1], convar[2]);

}
double Pressure(double density, double densityu, double densityE)
{
	return (Gamma - 1) * (densityE - 0.5 * densityu * densityu / density);
}
double Temperature(double density, double pressure)
{
	return pressure / R_gas / density;
}
double entropy(double density, double pressure)
{
	return log(pressure / pow(density, Gamma));
}
double Soundspeed(double density, double pressure)
{
	return sqrt(Gamma * pressure / density);
}
void Primvar_to_convar_1D(double* convar, double primvar[3])
{
	convar[0] = primvar[0];
	convar[1] = DensityU(primvar[0], primvar[1]);
	convar[2] = DensityE(primvar[0], primvar[1], primvar[2]);
}
double DensityU(double density, double u)
{
	return density * u;
}
double DensityE(double density, double u, double pressure)
{
	return density * (pressure / density / (Gamma - 1) + 0.5 * (u * u));
}
void Convar_to_char1D(double* character, double primvar[3], double convar[3])
{
	double c = sqrt(Gamma * primvar[2] / primvar[0]);
	double	alfa = (Gamma - 1.0) / (2.0 * c * c);
	double u = primvar[1];
	double s[3][3];
	s[0][0] = alfa * (0.5 * u * u + u * c / (Gamma - 1.0));
	s[0][1] = alfa * (-u - c / (Gamma - 1.0));
	s[0][2] = alfa;
	s[1][0] = alfa * (-u * u + 2.0 * c * c / (Gamma - 1.0));
	s[1][1] = alfa * 2.0 * u;
	s[1][2] = -2.0 * alfa;
	s[2][0] = alfa * (0.5 * u * u - u * c / (Gamma - 1.0));
	s[2][1] = alfa * (-u + c / (Gamma - 1.0));
	s[2][2] = alfa;

	for (int i = 0; i < 3; i++)
	{
		character[i] = 0;
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			character[i] = character[i] + s[i][j] * convar[j];
		}
	}
}
void Char_to_convar1D(double* convar, double primvar[3], double charvar[3])
{
	double c = sqrt(Gamma * primvar[2] / primvar[0]);
	double u = primvar[1];
	double	h = 0.5 * u * u + c * c / (Gamma - 1.0);
	double s[3][3];
	s[0][0] = 1.0;
	s[0][1] = 1.0;
	s[0][2] = 1.0;
	s[1][0] = u - c;
	s[1][1] = u;
	s[1][2] = u + c;
	s[2][0] = h - u * c;
	s[2][1] = u * u / 2.0;
	s[2][2] = h + u * c;

	for (int i = 0; i < 3; i++)
	{
		convar[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			convar[i] = convar[i] + s[i][j] * charvar[j];
		}
	}
}

double U(double density, double q_densityu)
{
	return q_densityu / density;
}
double V(double density, double q_densityv)
{
	return q_densityv / density;
}
double Pressure(double density, double q_densityu, double q_densityv, double q_densityE)
{
	return (Gamma - 1.0)*(q_densityE - 0.5*q_densityu*q_densityu / density - 0.5*q_densityv*q_densityv / density);
}




//从原始变量计算守恒变量
double Q_densityu(double density, double u)
{
	double q_densityu = density*u;
	return q_densityu;
}
double Q_densityv(double density, double v)
{
	double q_densityv = density*v;
	return q_densityv;
}
double Q_densityE(double density, double u, double v, double pressure)
{
	double q_densityE = density*(pressure / density / (Gamma - 1) + 0.5*(u*u + v*v));
	return q_densityE;
}



void Primvar_to_convar_2D(double *convar, double primvar[4])//2d version
{
    convar[0] = primvar[0];
    convar[1] = Q_densityu(primvar[0],primvar[1]);
    convar[2] = Q_densityv(primvar[0],primvar[2]);
    convar[3] = Q_densityE(primvar[0], primvar[1], primvar[2], primvar[3]);
}
void Convar_to_primvar_2D(double *primvar, double convar[4])//2d version //density, U, V, pressure
{
    primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = V(convar[0], convar[2]);
	primvar[3] = Pressure(convar[0], convar[1], convar[2], convar[3]);
}

void Convar_to_char(double *character, double *primvar, double convar[4])
{
	double c = sqrt(Gamma *primvar[3] / primvar[0]);
	double overc = 1.0 / c;
	double	alfa = (Gamma - 1.0)*0.5*overc*overc;
	double u = primvar[1];
	double v = primvar[2];
	double s[4][4];
	s[0][0] = 1.0 - alfa*(u*u + v*v);
	s[0][1] = 2.0*alfa*u;
	s[0][2] = 2.0*alfa*v;
	s[0][3] = -2.0*alfa;
	s[1][0] = -v;
	s[1][1] = 0.;
	s[1][2] = 1.;
	s[1][3] = 0.;
	s[2][0] = 0.5*(-u *overc + alfa*(u*u + v*v));
	s[2][1] = 0.5 *overc - alfa*u;
	s[2][2] = -alfa*v;
	s[2][3] = alfa;
	s[3][0] = 0.5*(u *overc + alfa*(u*u + v*v));
	s[3][1] = -0.5 *overc - alfa*u;
	s[3][2] = -alfa*v;
	s[3][3] = alfa;
	for (int i = 0; i < 4; i++)
	{
		character[i] = 0;
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			character[i] = character[i] + s[i][j] * convar[j];
		}
	}
}
void Char_to_convar(double *convar, double *primvar, double character[4])
{
	double c = sqrt(Gamma *primvar[3] / primvar[0]);
	double	alfa = (Gamma - 1.0) / (2.0*c*c);
	double u = primvar[1];
	double v = primvar[2];
	double s[4][4];
	s[0][0] = 1.;
	s[0][1] = 0.0;
	s[0][2] = 1.0;
	s[0][3] = 1.0;
	s[1][0] = u;
	s[1][1] = 0.;
	s[1][2] = u + c;
	s[1][3] = u - c;
	s[2][0] = v;
	s[2][1] = 1.;
	s[2][2] = v;
	s[2][3] = v;
	s[3][0] = 0.5*(u*u + v*v);
	s[3][1] = v;
	s[3][2] = 0.5*(u*u + v*v) + c*u + c*c / (Gamma - 1.0);
	s[3][3] = 0.5*(u*u + v*v) - c*u + c*c / (Gamma - 1.0);

	for (int i = 0; i < 4; i++)
	{
		convar[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			convar[i] = convar[i] + s[i][j] * character[j];
		}
	}
}

void Global_to_Local(double* change, double *normal)
{
	double temp[2];
	temp[0] = change[1];
	temp[1] = change[2];
	change[1] = temp[0] * normal[0] + temp[1] * normal[1];
	change[2] = -temp[0] * normal[1] + temp[1] * normal[0];
}

void Global_to_Local(double* change, double* origin, double* normal)
{
	change[0] = origin[0];
	change[1] = origin[1] * normal[0] + origin[2] * normal[1];
	change[2] = -origin[1] * normal[1] + origin[2] * normal[0];
	change[3] = origin[3];
}

void Local_to_Global(double *change,double *normal)
{
    double temp[2];
    temp[0] = change[1];
    temp[1] = change[2];
    change[1] = temp[0] * normal[0] - temp[1] * normal[1];
    change[2] = temp[0] * normal[1] + temp[1] * normal[0];
}

void Local_to_Global(double* change, double* origin, double *normal)
{
	change[0] = origin[0];
	change[1] = origin[1] * normal[0] - origin[2] * normal[1];
	change[2] = origin[1] * normal[1] + origin[2] * normal[0];
	change[3] = origin[3];
}

double Dtx(double dtx, double dx, double CFL, double density, double u,double v, double pressure)
{
	double tmp;
	tmp = sqrt(u*u+v*v) + sqrt(Gamma *pressure / density);
	if (tmp>CFL*dx / dtx)
	{
		dtx = CFL*dx / tmp;
	}
	if (dtx > 0.25*CFL*dx*dx / Mu&&tau_type == NS&&Mu>0)
	{
		dtx = 0.25*CFL*dx*dx / Mu;
	}
	return dtx;
}

