//this libarary includes some basic functions for gks
#include"gks_basic.h"
#define S 110.4

//artifical viscosity for euler problem
double c1_euler = -1.0; //initialization
double c2_euler = -1.0; //initialization
bool is_Prandtl_fix = false; //initialization
double Pr = 1.0; //initialization



// to prepare the basic element for moment calculation
MMDF1d::MMDF1d() { u = 1.0; lambda = 1.0; };
MMDF1d::MMDF1d(double u_in, double lambda_in)
{
	u = u_in;
	lambda = lambda_in;
	calcualte_MMDF1d();
}
void MMDF1d::calcualte_MMDF1d()
{
	uwhole[0] = 1;
	uwhole[1] = u;
	uplus[0] = 0.5 * Alpha(lambda, -u);
	uminus[0] = 0.5 * Alpha(lambda, u);
	uplus[1] = u * uplus[0] + 0.5 * Beta(lambda, u);
	uminus[1] = u * uminus[0] - 0.5 * Beta(lambda, u);
	for (int i = 2; i <= 9; i++)
	{
		uwhole[i] = u * uwhole[i - 1] + 0.5 * (i - 1) / lambda * uwhole[i - 2];
	}
	for (int i = 2; i <= 9; i++)
	{
		uplus[i] = u * uplus[i - 1] + 0.5 * (i - 1) / lambda * uplus[i - 2];
		uminus[i] = u * uminus[i - 1] + 0.5 * (i - 1) / lambda * uminus[i - 2];
	}
	xi2 = 0.5 * K / lambda;
	xi4 = 0.25 * (K * K + 2 * K) / (lambda * lambda);
	//xi6?? how to calculate
	xi6 = 0.5 * (K + 4) / lambda * xi4;

	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					upxi[i][k] = uplus[i];
					unxi[i][k] = uminus[i];
				}
				if (k == 1)
				{
					upxi[i][k] = uplus[i] * xi2;
					unxi[i][k] = uminus[i] * xi2;
				}
				if (k == 2)
				{
					upxi[i][k] = uplus[i] * xi4;
					unxi[i][k] = uminus[i] * xi4;
				}
				if (k == 3)
				{
					upxi[i][k] = uplus[i] * xi6;
					unxi[i][k] = uminus[i] * xi6;
				}
			}
		}
	}
	for (int i = 0; i < 10; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			if ((i + 2 * k) <= 9)
			{
				if (k == 0)
				{
					uxi[i][k] = uwhole[i];
				}
				if (k == 1)
				{
					uxi[i][k] = uwhole[i] * xi2;
				}
				if (k == 2)
				{
					uxi[i][k] = uwhole[i] * xi4;
				}
				if (k == 3)
				{
					uxi[i][k] = uwhole[i] * xi6;
				}
			}
		}
	}
}

//a general G function
void G(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.uxi[no_u][no_xi] + a[1] * m.uxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.uxi[no_u + 1][no_xi] + a[1] * m.uxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.uxi[no_u + 2][no_xi] + m.uxi[no_u][no_xi + 1]) +
		a[1] * (m.uxi[no_u + 3][no_xi] + m.uxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.uxi[no_u + 4][no_xi] + m.uxi[no_u][no_xi + 2] + 2 * m.uxi[no_u + 2][no_xi + 1]));

}
void GL(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.upxi[no_u][no_xi] + a[1] * m.upxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.upxi[no_u + 1][no_xi] + a[1] * m.upxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.upxi[no_u + 2][no_xi] + m.upxi[no_u][no_xi + 1]) +
		a[1] * (m.upxi[no_u + 3][no_xi] + m.upxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.upxi[no_u + 4][no_xi] + m.upxi[no_u][no_xi + 2] + 2 * m.upxi[no_u + 2][no_xi + 1]));

}
void GR(int no_u, int no_xi, double* psi, double a[3], MMDF1d m)
{
	psi[0] = a[0] * m.unxi[no_u][no_xi] + a[1] * m.unxi[no_u + 1][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]);
	psi[1] = a[0] * m.unxi[no_u + 1][no_xi] + a[1] * m.unxi[no_u + 2][no_xi] + a[2] * 0.5 * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]);
	psi[2] = 0.5 * (a[0] * (m.unxi[no_u + 2][no_xi] + m.unxi[no_u][no_xi + 1]) +
		a[1] * (m.unxi[no_u + 3][no_xi] + m.unxi[no_u + 1][no_xi + 1]) +
		a[2] * 0.5 * (m.unxi[no_u + 4][no_xi] + m.unxi[no_u][no_xi + 2] + 2 * m.unxi[no_u + 2][no_xi + 1]));

}

//moments of the maxwellian distribution function
MMDF::MMDF() { u = 1.0; v = 1.0; lambda = 1.0; };
MMDF::MMDF(double u_in, double v_in, double lambda_in)
{
	u = u_in;
	v = v_in;
	lambda = lambda_in;
	calcualte_MMDF();
}
void MMDF::calcualte_MMDF()
{
	//the moment of u in whole domain
	uwhole[0] = 1;
	uwhole[1] = u;
	//the moment of u in half domain (>0 && <0)
	uplus[0] = 0.5 * Alpha(lambda, -u);
	uminus[0] = 1.0 - uplus[0];
	uplus[1] = u * uplus[0] + 0.5 * Beta(lambda, u);
	uminus[1] = u - uplus[1];
	double overlambda = 1.0 / lambda;
	for (int i = 2; i <= 6; i++)
	{
		uwhole[i] = u * uwhole[i - 1] + 0.5 * (i - 1) * overlambda * uwhole[i - 2];
	}
	for (int i = 2; i <= 6; i++)
	{
		uplus[i] = u * uplus[i - 1] + 0.5 * (i - 1) * overlambda * uplus[i - 2];
		uminus[i] = uwhole[i] - uplus[i];
	}
	//the moment of v in whole domain
	//the moment of v in half domain is useless
	vwhole[0] = 1;
	vwhole[1] = v;
	for (int i = 2; i <= 6; i++)
	{
		vwhole[i] = v * vwhole[i - 1] + 0.5 * (i - 1) * overlambda * vwhole[i - 2];
	}
	//the moment of xi (kesi)
	xi2 = 0.5 * K * overlambda;
	xi4 = 0.25 * (K * K + 2 * K) * overlambda * overlambda;
	//xi6?? how to calculate
	//xi6 = 0.5*(K + 4) / lambda*xi4; // no use

	//get other variables by the above mement variables
	//i, j, k, means i power of u, j power of v, and 2k power of xi
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if ((i + j + 2 * k) <= 6)
				{
					if (k == 0)
					{
						// u, half domain > 0; v, whole domain; xi, equals 0
						upvxi[i][j][k] = uplus[i] * vwhole[j];
						// u, half domain < 0; v, whole domain; xi, equals 0
						unvxi[i][j][k] = uminus[i] * vwhole[j];
						// u, whole domain ; v, whole domain; xi, equals 0
						uvxi[i][j][k] = upvxi[i][j][k] + unvxi[i][j][k];
						// = (uplus[i] + uminus[i]) * vwhole[j] = uwhole[i] * vwhole[j]
					}
					if (k == 1)
					{
						// u, half domain > 0; v, whole domain; xi, equals 2
						upvxi[i][j][k] = uplus[i] * vwhole[j] * xi2;
						// u, half domain < 0; v, whole domain; xi, equals 2
						unvxi[i][j][k] = uminus[i] * vwhole[j] * xi2;
						// u, whole domain ; v, whole domain; xi, equals 2
						uvxi[i][j][k] = upvxi[i][j][k] + unvxi[i][j][k];
					}
					if (k == 2)
					{
						// u, half domain > 0; v, whole domain; xi, equals 4
						upvxi[i][j][k] = uplus[i] * vwhole[j] * xi4;
						// u, half domain < 0; v, whole domain; xi, equals 4
						unvxi[i][j][k] = uminus[i] * vwhole[j] * xi4;
						// u, whole domain ; v, whole domain; xi, equals 4
						uvxi[i][j][k] = upvxi[i][j][k] + unvxi[i][j][k];
					}

				}
			}
		}
	}
}

MMDF1st::MMDF1st() { u = 1.0; v = 1.0; lambda = 1.0; };
MMDF1st::MMDF1st(double u_in, double v_in, double lambda_in)
{
	u = u_in;
	v = v_in;
	lambda = lambda_in;
	calcualte_MMDF1st();
}
void MMDF1st::calcualte_MMDF1st()
{
	uwhole[0] = 1.0;
	uwhole[1] = u;

	//return erfc(sqrt(lambda)*u);

	//return exp(-lambda*u*u) / sqrt(pi*lambda);
	uminus[0] = 0.5 * Alpha(lambda, u);
	uplus[0] = 1.0 - uminus[0];
	double beta = Beta(lambda, u);
	//uplus[1] = u*uplus[0] + 0.5*Beta(lambda, u);
	//uminus[1] = u*uminus[0] - 0.5*Beta(lambda, u);
	uplus[1] = u * uplus[0] + 0.5 * beta;
	uminus[1] = u * uminus[0] - 0.5 * beta;
	double overlambda = 1.0 / lambda;
	for (int i = 2; i <= 3; i++)
	{
		uwhole[i] = u * uwhole[i - 1] + 0.5 * (i - 1) * overlambda * uwhole[i - 2];
	}
	for (int i = 2; i <= 3; i++)
	{
		uplus[i] = u * uplus[i - 1] + 0.5 * (i - 1) * overlambda * uplus[i - 2];
		uminus[i] = uwhole[i] - uplus[i];
	}
	vwhole[0] = 1;
	vwhole[1] = v;
	for (int i = 2; i <= 2; i++)
	{
		vwhole[i] = v * vwhole[i - 1] + 0.5 * (i - 1) * overlambda * vwhole[i - 2];
	}
	xi2 = 0.5 * K * overlambda;
}



void Collision(double* w0, double left, double right, MMDF& m2, MMDF& m3)
{
	// get the equilibrium variables by collision
	w0[0] = left * m2.uplus[0] + right * m3.uminus[0];
	w0[1] = left * m2.uplus[1] + right * m3.uminus[1];
	w0[2] = left * m2.vwhole[1] * m2.uplus[0] + right * m3.vwhole[1] * m3.uminus[0];
	w0[3] = 0.5 * left * (m2.uplus[2] + m2.uplus[0] * m2.vwhole[2] + m2.uplus[0] * m2.xi2) +
		0.5 * right * (m3.uminus[2] + m3.uminus[0] * m3.vwhole[2] + m3.uminus[0] * m3.xi2);
}
void Collision(double* w0, double left, double right, MMDF1st& m2, MMDF1st& m3)
{
	w0[0] = left * m2.uplus[0] + right * m3.uminus[0];
	w0[1] = left * m2.uplus[1] + right * m3.uminus[1];
	w0[2] = left *  m2.vwhole[1] * m2.uplus[0] + right * m3.vwhole[1] * m3.uminus[0];
	w0[3] = 0.5 * left * (m2.uplus[2] + m2.uplus[0] * m2.vwhole[2] + m2.uplus[0] * m2.xi2) +
		0.5 * right * (m3.uminus[2] + m3.uminus[0] * m3.vwhole[2] + m3.uminus[0] * m3.xi2);
}



void A(double* a, double der[4], double prim[4])
{
	//Note: solve the a coefficient (slope information), equivalent to solving the matrix to get a
	//input: der[4], the gradient of conservative variables, density, momentumX, momentumY, energy, partialW/patialX
	//input: prim[4], the primary variables (density, u, v, lambda) at one point or one interface
	//output: a, solve the results

	double R4, R3, R2;
	double overden = 1.0 / prim[0];
	//double overlambda = 1.0 / prim[3];
	R4 = der[3] * overden - 0.5 * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]) * der[0] * overden;
	R3 = (der[2] - prim[2] * der[0]) * overden;
	R2 = (der[1] - prim[1] * der[0]) * overden;
	a[3] = (4.0 / (K + 2)) * prim[3] * prim[3] * (2 * R4 - 2 * prim[1] * R2 - 2 * prim[2] * R3);
	a[2] = 2 * prim[3] * R3 - prim[2] * a[3];
	a[1] = 2 * prim[3] * R2 - prim[1] * a[3];
	a[0] = der[0] * overden - prim[1] * a[1] - prim[2] * a[2] - 0.5 * a[3] * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]);

	//R4 = der[3] / prim[0] - 0.5*(prim[1] * prim[1] + prim[2] * prim[2] + 0.5*(K + 2) / prim[3])*der[0] / prim[0];
	//R3 = (der[2] - prim[2]*der[0]) / prim[0];
	//R2 = (der[1] - prim[1]*der[0]) / prim[0];
	//a[3] = 4.0 * prim[3]*prim[3] / (K + 2)*(2.0 * R4 - 2.0 * prim[1]*R2 - 2.0 * prim[2]*R3);
	//a[2] = 2.0 * prim[3]*R3 - prim[2]*a[3];
	//a[1] = 2.0 * prim[3]*R2 - prim[1]*a[3];
	//a[0] = der[0] / prim[0] - prim[1]*a[1] - prim[2]*a[2] - 0.5*a[3] * (prim[1]*prim[1] + prim[2]*prim[2] + 0.5*(K + 2) / prim[3]);

}



double Get_Tau_NS(double density0, double lambda0)
{
	if (tau_type == Euler)
	{
		return 0.0;
	}
	else if (tau_type == NS)
	{
		if (Mu > 0.0)
		{
			//cout << "here" << endl;
			return 2.0 * Mu * lambda0 / density0;
		}
		else if (Nu > 0.0)
		{
			return 2 * Nu * lambda0;
		}
		else
		{
			return 0.0;
		}
	}
	else if (tau_type == Sutherland)
	{
		return TauNS_Sutherland(density0, lambda0);
	}
	else
	{
		return TauNS_power_law(density0, lambda0);
	}

}
double Get_Tau(double density_left, double density_right, double density0, double lambda_left, double lambda_right, double lambda0, double dt)
{
	if (tau_type == Euler)
	{
		if (c1_euler <= 0 && c2_euler <= 0)
		{
			return 0.0;
		}
		else
		{
			double C = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right);
			//C+= c2_euler * abs(density_left - density_right) / abs(density_left + density_right);
			//if (C < 10)
			//{
			return c1_euler * dt + dt * C;
			//}
			//else
			//{
			//	return c1_euler*dt + 10*dt;
			//	//return c1_euler*dt + dt*C;
			//}
		}
	}
	else if (tau_type == NS)
	{
		double tau_n = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right) * dt;
		if (tau_n != tau_n)
		{
			tau_n = 0.0;
		}
		if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
		{
			return 0.0;
		}
		else
		{
			if (Mu > 0.0)
			{
				return tau_n + 2.0 * Mu * lambda0 / density0;
			}
			else if (Nu > 0.0)
			{
				return tau_n + 2 * Nu * lambda0;
			}
			else
			{
				return 0.0;
			}
		}
	}
	else if (tau_type == Sutherland)
	{
		double tau_n = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right) * dt;
		if (tau_n != tau_n)
		{
			tau_n = 0.0;
		}
		if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
		{

			return 0.0;
		}
		else
		{
			if (Mu > 0.0)
			{
				return tau_n + TauNS_Sutherland(density0, lambda0);

			}
			else if (Nu > 0.0)
			{
				return tau_n + 2 * Nu * lambda0;
			}
			else
			{
				return 0.0;
			}
		}
	}
	else if (tau_type == Power_law)
	{

		double tau_n = c2_euler * abs(density_left / lambda_left - density_right / lambda_right) / abs(density_left / lambda_left + density_right / lambda_right) * dt;
		if (tau_n != tau_n)
		{
			tau_n = 0.0;
		}
		if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
		{
			return 0.0;
		}
		else
		{
			if (Mu > 0.0)
			{
				return tau_n + TauNS_power_law(density0, lambda0);
			}
			else if (Nu > 0.0)
			{
				return tau_n + 2 * Nu * lambda0;
			}
			else
			{
				return 0.0;
			}
		}

	}
	else
	{
		return 0.0;
	}
}
double TauNS_Sutherland(double density0, double lambda0)
{
	double T0 = 1.0 / 2.0 / (lambda0 * R_gas);
	double mu = Mu * pow(T0 / T_inf, 1.5) * (T_inf + S) / (T0 + S);
	return 2.0 * mu * lambda0 / density0;
}
double TauNS_power_law(double density0, double lambda0)
{
	double T0 = 1.0 / 2.0 / (lambda0 * R_gas);
	double mu = Mu * pow(T0 / T_inf, 0.76);
	return 2.0 * mu * lambda0 / density0;
}



//计算λ
double Lambda(double density, double u, double densityE)
{
	return (K + 1.0) * 0.25 * (density / (densityE - 0.5 * density * (u * u)));
}
double Lambda(double density, double u, double v, double densityE)
{
	return (K + 2.0) * 0.25 * (density / (densityE - 0.5 * density * (u * u + v * v)));
}

void Convar_to_ULambda_1d(double* primvar, double convar[3])
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = Lambda(convar[0], primvar[1], convar[2]);
}

//计算α，α=erfc（±√λ*U）
double Alpha(double lambda, double u)
{
	return erfc(sqrt(lambda) * u);
	//return erfcmy(sqrt(lambda)*u);
}
//计算β，β=e^(-λ*U^2)/(√πλ)
double Beta(double lambda, double u)
{
	double pi = 3.14159265358979323846;
	return exp(-lambda * u * u) / sqrt(pi * lambda);
}


//solution of matrix equation b=Ma 1D
void Microslope(double* a, double der[3], double prim[3])
{
	SolveMicroCoeff1D(a, der, prim);
}

void SolveMicroCoeff1D(double* coeff, double moment[3], double prim[3])
{
	// Solve c in <c> = moment with
	// c(u,xi) = c0 + c1*u + c2*0.5*(u^2 + xi^2).
	// This is the 1D Euler coefficient form used by paper Eq. (24)-(26).
	double R4, R2;
	R4 = moment[2] / prim[0] - 0.5 * (prim[1] * prim[1] + 0.5 * (K + 1) / prim[2]) * moment[0] / prim[0];
	R2 = (moment[1] - prim[1] * moment[0]) / prim[0];
	coeff[2] = 4 * prim[2] * prim[2] / (K + 1) * (2 * R4 - 2 * prim[1] * R2);
	coeff[1] = 2 * prim[2] * R2 - prim[1] * coeff[2];
	coeff[0] = moment[0] / prim[0] - prim[1] * coeff[1] - 0.5 * coeff[2] * (prim[1] * prim[1] + 0.5 * (K + 1) / prim[2]);
}

void Moment_uCoeff1D(double* moment, double coeff[3], double prim[3])
{
	// Return rho * <u * coeff * psi> for the local Maxwellian.
	MMDF1d m(prim[1], prim[2]);
	G(1, 0, moment, coeff, m);
	for (int i = 0; i < 3; ++i)
	{
		moment[i] *= prim[0];
	}
}

void SolveGKS_a1D(double* a, double Wx[3], double prim[3])
{
	// Paper Eq. (25) in 1D: <a> = dW/dx.
	SolveMicroCoeff1D(a, Wx, prim);
}

void SolveGKS_A1D(double* A, double a[3], double prim[3])
{
	// Paper Eq. (25) in 1D: <A + a*u> = 0.
	double rhs[3];
	Moment_uCoeff1D(rhs, a, prim);
	for (int i = 0; i < 3; ++i)
	{
		rhs[i] = -rhs[i];
	}
	SolveMicroCoeff1D(A, rhs, prim);
}

void SolveGKS_axx1D(double* axx, double Wxx[3], double prim[3])
{
	// Paper Eq. (26) in 1D: <a^2 + d> = d^2W/dx^2.
	// In code we solve the combined coefficient axx := a^2 + d directly.
	SolveMicroCoeff1D(axx, Wxx, prim);
}

void SolveGKS_axt1D(double* axt, double axx[3], double prim[3])
{
	// Paper Eq. (26) in 1D: <(A*a + b) + (a^2 + d)u> = 0.
	// In code we solve the combined coefficient axt := A*a + b directly.
	double rhs[3];
	Moment_uCoeff1D(rhs, axx, prim);
	for (int i = 0; i < 3; ++i)
	{
		rhs[i] = -rhs[i];
	}
	SolveMicroCoeff1D(axt, rhs, prim);
}

void SolveGKS_att1D(double* att, double axt[3], double prim[3])
{
	// Paper Eq. (26) in 1D: <(A^2 + B) + (A*a + b)u> = 0.
	// In code we solve the combined coefficient att := A^2 + B directly.
	double rhs[3];
	Moment_uCoeff1D(rhs, axt, prim);
	for (int i = 0; i < 3; ++i)
	{
		rhs[i] = -rhs[i];
	}
	SolveMicroCoeff1D(att, rhs, prim);
}

void Convar_to_ULambda_2d(double* primvar, double convar[4]) //density, U, V, lambda
{
	primvar[0] = convar[0];
	primvar[1] = U(convar[0], convar[1]);
	primvar[2] = V(convar[0], convar[2]);
	//here primvar[3] refers lambda
	primvar[3] = Lambda(convar[0], primvar[1], primvar[2], convar[3]);
}
