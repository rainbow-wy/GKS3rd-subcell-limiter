#include "fvm_gks2d.h"
#include "gks_smooth_indicator.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace
{
	void Configure2DBase(double mu, double c1, double c2, TAU_TYPE tau_mode)
	{
		K = 3;
		Gamma = 1.4;
		R_gas = 1.0;
		Pr = 0.73;
		Mu = mu;
		Nu = -1.0;
		tau_type = tau_mode;
		c1_euler = c1;
		c2_euler = c2;
		gks2dsolver = gks3rd_2d;
	}

	void Configure1DBase(double c1, double c2)
	{
		K = 4;
		Gamma = 1.4;
		R_gas = 1.0;
		tau_type = Euler;
		Mu = 0.0;
		Nu = -1.0;
		c1_euler = c1;
		c2_euler = c2;
		gks1dsolver = gks3rd;
	}

	void PrimToConservative2D(double* convar, double rho, double u, double v, double p)
	{
		double prim[4] = { rho, u, v, p };
		Primvar_to_convar_2D(convar, prim);
	}

	void PrimToConservative1D(double* convar, double rho, double u, double p)
	{
		double prim[3] = { rho, u, p };
		Primvar_to_convar_1D(convar, prim);
	}

	void ClearPoint2D(Point2d& point)
	{
		for (int m = 0; m < 4; ++m)
		{
			point.convar[m] = 0.0;
			point.der1x[m] = 0.0;
			point.der1y[m] = 0.0;
			point.der2xx[m] = 0.0;
			point.der2xy[m] = 0.0;
			point.der2yy[m] = 0.0;
		}
		point.x = 0.0;
		point.y = 0.0;
		point.normal[0] = 1.0;
		point.normal[1] = 0.0;
		point.is_reduce_order = false;
	}

	void SetPoint2D(
		Point2d& point,
		double rho,
		double u,
		double v,
		double p,
		const double* der1x,
		const double* der1y,
		const double* der2xx,
		const double* der2xy,
		const double* der2yy)
	{
		ClearPoint2D(point);
		PrimToConservative2D(point.convar, rho, u, v, p);
		for (int m = 0; m < 4; ++m)
		{
			point.der1x[m] = der1x[m];
			point.der1y[m] = der1y[m];
			point.der2xx[m] = der2xx[m];
			point.der2xy[m] = der2xy[m];
			point.der2yy[m] = der2yy[m];
		}
	}

	void SetPoint1D(
		Point1d& point,
		double rho,
		double u,
		double p,
		const double* der1,
		const double* der2)
	{
		for (int m = 0; m < 3; ++m)
		{
			point.convar[m] = 0.0;
			point.der1[m] = der1[m];
			point.der2[m] = der2[m];
		}
		PrimToConservative1D(point.convar, rho, u, p);
		point.x = 0.0;
	}

	double MaxAbsDiff(const double* a, const double* b, int n)
	{
		double err = 0.0;
		for (int i = 0; i < n; ++i)
		{
			err = std::max(err, std::fabs(a[i] - b[i]));
		}
		return err;
	}

	bool UniformFlowViscousExactTest()
	{
		Configure2DBase(0.02, 0.0, 0.0, NS);
		Recon2d interface;
		double zero[4] = { 0.0, 0.0, 0.0, 0.0 };
		SetPoint2D(interface.left, 1.0, 0.8, 0.3, 1.0, zero, zero, zero, zero, zero);
		SetPoint2D(interface.right, 1.0, 0.8, 0.3, 1.0, zero, zero, zero, zero, zero);
		SetPoint2D(interface.center, 1.0, 0.8, 0.3, 1.0, zero, zero, zero, zero, zero);

		const double dt = 0.02;
		Flux2d flux;
		GKS2D(flux, interface, dt);

		double expect[4];
		double convar[4];
		PrimToConservative2D(convar, 1.0, 0.8, 0.3, 1.0);
		const double rho_u = convar[1];
		const double rho_v = convar[2];
		const double rho_e = convar[3];
		expect[0] = rho_u * dt;
		expect[1] = (rho_u * 0.8 + 1.0) * dt;
		expect[2] = rho_v * 0.8 * dt;
		expect[3] = (rho_e + 1.0) * 0.8 * dt;

		const double err = MaxAbsDiff(flux.f, expect, 4);
		std::cout << "uniform_viscous_exact max_err=" << err << std::endl;
		return err < 1.0e-10;
	}

	bool OneDimensionalReductionTest()
	{
		const double dt = 0.01;
		const double c1 = 0.05;
		const double c2 = 0.10;

		double der1_left_1d[3] = { 0.02, 0.01, 0.03 };
		double der1_right_1d[3] = { 0.01, -0.01, 0.02 };
		double der1_bar_1d[3] = { 0.015, 0.0, 0.025 };
		double der2_left_1d[3] = { 0.003, -0.002, 0.004 };
		double der2_right_1d[3] = { -0.002, 0.001, 0.003 };
		double der2_bar_1d[3] = { 0.001, -0.001, 0.002 };

		Configure1DBase(c1, c2);
		Interface1d interface1d;
		SetPoint1D(interface1d.left, 1.0, 0.35, 1.0, der1_left_1d, der2_left_1d);
		SetPoint1D(interface1d.right, 0.9, 0.20, 0.95, der1_right_1d, der2_right_1d);
		SetPoint1D(interface1d.center, 0.95, 0.28, 0.98, der1_bar_1d, der2_bar_1d);
		Flux1d flux1d;
		GKS(flux1d, interface1d, dt);

		double zero4[4] = { 0.0, 0.0, 0.0, 0.0 };
		double der1_left_2d[4] = { der1_left_1d[0], der1_left_1d[1], 0.0, der1_left_1d[2] };
		double der1_right_2d[4] = { der1_right_1d[0], der1_right_1d[1], 0.0, der1_right_1d[2] };
		double der1_bar_2d[4] = { der1_bar_1d[0], der1_bar_1d[1], 0.0, der1_bar_1d[2] };
		double der2_left_2d[4] = { der2_left_1d[0], der2_left_1d[1], 0.0, der2_left_1d[2] };
		double der2_right_2d[4] = { der2_right_1d[0], der2_right_1d[1], 0.0, der2_right_1d[2] };
		double der2_bar_2d[4] = { der2_bar_1d[0], der2_bar_1d[1], 0.0, der2_bar_1d[2] };

		Configure2DBase(0.0, c1, c2, Euler);
		Recon2d interface2d;
		SetPoint2D(interface2d.left, 1.0, 0.35, 0.0, 1.0, der1_left_2d, zero4, der2_left_2d, zero4, zero4);
		SetPoint2D(interface2d.right, 0.9, 0.20, 0.0, 0.95, der1_right_2d, zero4, der2_right_2d, zero4, zero4);
		SetPoint2D(interface2d.center, 0.95, 0.28, 0.0, 0.98, der1_bar_2d, zero4, der2_bar_2d, zero4, zero4);
		Flux2d flux2d;
		GKS2D(flux2d, interface2d, dt);

		double flux2d_1dmap[3] = { flux2d.f[0], flux2d.f[1], flux2d.f[3] };
		const double err = MaxAbsDiff(flux1d.f, flux2d_1dmap, 3);
		std::cout << "one_dimensional_reduction max_err=" << err
			<< " tangential_flux=" << std::fabs(flux2d.f[2]) << std::endl;
		return err < 2.0e-6 && std::fabs(flux2d.f[2]) < 1.0e-10;
	}

	bool ShearViscousDirectionTest()
	{
		const double dt = 0.01;
		double zero[4] = { 0.0, 0.0, 0.0, 0.0 };
		double der1x[4] = { 0.0, 0.0, 0.40, 0.0 };

		Recon2d interface;
		SetPoint2D(interface.left, 1.0, 0.0, 0.0, 1.0, der1x, zero, zero, zero, zero);
		SetPoint2D(interface.right, 1.0, 0.0, 0.0, 1.0, der1x, zero, zero, zero, zero);
		SetPoint2D(interface.center, 1.0, 0.0, 0.0, 1.0, der1x, zero, zero, zero, zero);

		Configure2DBase(0.0, 0.0, 0.0, Euler);
		Flux2d flux_inv;
		GKS2D(flux_inv, interface, dt);

		Configure2DBase(0.05, 0.0, 0.0, NS);
		Flux2d flux_vis;
		GKS2D(flux_vis, interface, dt);

		std::cout << "shear_toggle inv_tangential=" << flux_inv.f[2]
			<< " visc_tangential=" << flux_vis.f[2]
			<< " mass_delta=" << std::fabs(flux_vis.f[0] - flux_inv.f[0]) << std::endl;
		return std::fabs(flux_inv.f[2]) < 1.0e-12 &&
			flux_vis.f[2] < 0.0 &&
			std::fabs(flux_vis.f[0] - flux_inv.f[0]) < 1.0e-12;
	}

	double OrthonormalPhi2(double x)
	{
		return std::sqrt(2.5) * 0.5 * (3.0 * x * x - 1.0);
	}

	void BuildSmoothIndicatorState(int profile, double amplitude, double point_Q[3][3][4])
	{
		const double node[3] = {
			-std::sqrt(3.0 / 5.0),
			0.0,
			std::sqrt(3.0 / 5.0)
		};
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				double perturbation = 0.0;
				if (profile == 1)
				{
					perturbation = node[i];
				}
				else if (profile == 2)
				{
					perturbation = node[j];
				}
				else if (profile == 3)
				{
					perturbation = OrthonormalPhi2(node[i]);
				}
				else if (profile == 4)
				{
					perturbation = OrthonormalPhi2(node[j]);
				}
				else if (profile == 5)
				{
					perturbation = OrthonormalPhi2(node[i]) * OrthonormalPhi2(node[j]);
				}
				PrimToConservative2D(point_Q[i][j], 1.0, 0.0, 0.0, 1.0 + amplitude * perturbation);
			}
		}
	}

	bool SmoothIndicatorTensorModalSelfTest()
	{
		Gamma = 1.4;
		const double amplitude = 0.1;
		const double tolerance = 2.0e-13;
		double max_error = 0.0;
		bool ok = true;
		for (int profile = 0; profile < 6; ++profile)
		{
			double point_Q[3][3][4];
			BuildSmoothIndicatorState(profile, amplitude, point_Q);
			GKSSmoothIndicatorParam2D param;
			GKSSmoothIndicatorCellDiag2D diag;
			GKSSmoothIndicatorCellDiagnostics2D(point_Q, param, diag);

			double expected[3][3]{};
			expected[0][0] = 2.0;
			if (profile == 1)
			{
				expected[1][0] = 2.0 * amplitude / std::sqrt(3.0);
			}
			else if (profile == 2)
			{
				expected[0][1] = 2.0 * amplitude / std::sqrt(3.0);
			}
			else if (profile == 3)
			{
				expected[2][0] = std::sqrt(2.0) * amplitude;
			}
			else if (profile == 4)
			{
				expected[0][2] = std::sqrt(2.0) * amplitude;
			}
			else if (profile == 5)
			{
				expected[2][2] = amplitude;
			}

			for (int m = 0; m < 3; ++m)
			{
				for (int n = 0; n < 3; ++n)
				{
					max_error = std::max(max_error, std::fabs(diag.qhat[m][n] - expected[m][n]));
				}
			}

			double expected_E1 = 0.0;
			double expected_E2 = 0.0;
			if (profile == 1 || profile == 2)
			{
				expected_E1 = (amplitude * amplitude / 3.0)
					/ (1.0 + amplitude * amplitude / 3.0);
			}
			else if (profile == 3 || profile == 4)
			{
				expected_E2 = amplitude * amplitude / (2.0 + amplitude * amplitude);
			}
			else if (profile == 5)
			{
				expected_E2 = amplitude * amplitude / (4.0 + amplitude * amplitude);
			}
			max_error = std::max(max_error, std::fabs(diag.E1 - expected_E1));
			max_error = std::max(max_error, std::fabs(diag.E2 - expected_E2));
			max_error = std::max(max_error, std::fabs(diag.E - std::max(expected_E1, expected_E2)));
			max_error = std::max(max_error, std::fabs(diag.S2 - (diag.S1
				+ diag.qhat[2][0] * diag.qhat[2][0]
				+ diag.qhat[2][1] * diag.qhat[2][1]
				+ diag.qhat[2][2] * diag.qhat[2][2]
				+ diag.qhat[0][2] * diag.qhat[0][2]
				+ diag.qhat[1][2] * diag.qhat[1][2])));
			if (profile == 0)
			{
				ok = diag.alpha_raw == 0.0 && ok;
			}
		}
		ok = max_error < tolerance && ok;
		std::cout << "smooth_indicator_tensor_modal max_err=" << max_error << std::endl;
		return ok;
	}
}

int fluxtest_main()
{
	std::cout << std::setprecision(16);
	bool ok = true;
	ok = UniformFlowViscousExactTest() && ok;
	ok = OneDimensionalReductionTest() && ok;
	ok = ShearViscousDirectionTest() && ok;
	ok = SmoothIndicatorTensorModalSelfTest() && ok;
	std::cout << (ok ? "ALL_GKS3RD_2D_FLUX_TESTS_PASS" : "GKS3RD_2D_FLUX_TESTS_FAIL") << std::endl;
	return ok ? 0 : 1;
}
