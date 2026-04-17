#include "gks3rd_demo.h"

#include "accuracy_test.h"
#include "riemann_problem.h"

#include <fstream>
#include <iomanip>
#include <string>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace
{
	void Ensure_GKS3rd_Result_Directory()
	{
		// Keep the formal numerical outputs under build/result so they are
		// separated from temporary debug artifacts in build/checks.
#if defined(_WIN32)
		_mkdir("build");
		_mkdir("build\\result");
#else
		mkdir("build", 0777);
		mkdir("build/result", 0777);
#endif
	}

	void Configure_GKS3rd_1D(double c1, double c2)
	{
		K = 4;
		Gamma = 1.4;
		R_gas = 1.0;

		tau_type = Euler;
		c1_euler = c1;
		c2_euler = c2;

		flux_function = GKS;
		gks1dsolver = gks3rd;
		cellreconstruction = WENO5_AO;
		reconstruction_variable = characteristic;
		wenotype = wenoz;
		g0reconstruction = Center_GKS3rd;
		timecoe_list = S1O3;
		is_reduce_order_warning = false;
		is_show_1d_timestep = true;
	}

	void Prepare_Block_1D(Block1d& block, int mesh_number, double left, double right, double CFL)
	{
		block.nodex = mesh_number;
		block.ghost = 3;
		block.nx = block.nodex + 2 * block.ghost;
		block.nodex_begin = block.ghost;
		block.nodex_end = block.nodex + block.ghost - 1;
		block.left = left;
		block.right = right;
		block.dx = (block.right - block.left) / block.nodex;
		block.CFL = CFL;
		Initial_stages(block);
	}

	void Write_Density_Output_Tecplot(Fluid1d* fluids, Block1d block, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,density\n";
		out << "zone i = " << block.nodex << ", F=POINT\n";
		for (int i = block.ghost; i < block.nx - block.ghost; ++i)
		{
			out << fluids[i].cx << " " << fluids[i].convar[0] << "\n";
		}
	}

	std::string StepTaggedPath(const std::string& prefix, int final_step, const char* ext)
	{
		return "build/result/" + prefix + "_" + std::to_string(final_step) + ext;
	}

	void Advance_1D_Case(
		Fluid1d* fluids,
		Interface1d* interfaces,
		Flux1d** fluxes,
		Block1d& block,
		BoundaryCondition leftboundary,
		BoundaryCondition rightboundary,
		Fluid1d* bcvalue,
		double tstop)
	{
		block.t = 0.0;
		block.step = 0;
		while (block.t < tstop - 1e-14)
		{
			CopyFluid_new_to_old(fluids, block);
			block.dt = Get_CFL(block, fluids, tstop);
			for (int stage = 0; stage < block.stages; ++stage)
			{
				leftboundary(fluids, block, bcvalue[0]);
				rightboundary(fluids, block, bcvalue[1]);
				Reconstruction_within_cell(interfaces, fluids, block);
				Reconstruction_forg0(interfaces, fluids, block);
				Calculate_flux(fluxes, interfaces, block, stage);
				Update(fluids, fluxes, block, stage);
			}
			block.t += block.dt;
			block.step++;
		}
	}

	void Compute_Sinwave_Error(Fluid1d* fluids, Block1d block, double* error)
	{
		double error1 = 0.0;
		double error2 = 0.0;
		double errorinf = 0.0;
		for (int i = block.ghost; i < block.nx - block.ghost; ++i)
		{
			const double pi = 3.14159265358979323846;
			const double exact = 1.0 - 0.2 / pi / block.dx *
				(cos(pi * (i + 1 - block.ghost) * block.dx) - cos(pi * (i - block.ghost) * block.dx));
			const double diff = fabs(fluids[i].convar[0] - exact);
			error1 += diff;
			error2 += diff * diff;
			errorinf = (diff > errorinf) ? diff : errorinf;
		}
		error[0] = error1 / block.nodex;
		error[1] = sqrt(error2 / block.nodex);
		error[2] = errorinf;
	}

	void Write_Sod_Output(Fluid1d* fluids, Block1d block, const char* path)
	{
		Convar_to_primvar(fluids, block);
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,density,u,pressure,temperature,entropy,Ma\n";
		out << "zone i = " << block.nodex << ", F=POINT\n";
		for (int i = block.ghost; i < block.nx - block.ghost; ++i)
		{
			out << fluids[i].cx << " "
				<< fluids[i].primvar[0] << " "
				<< fluids[i].primvar[1] << " "
				<< fluids[i].primvar[2] << " "
				<< Temperature(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
				<< entropy(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
				<< sqrt(fluids[i].primvar[1] * fluids[i].primvar[1]) /
					Soundspeed(fluids[i].primvar[0], fluids[i].primvar[2]) << "\n";
		}
	}
}

void accuracy_sinwave_1d_gks3rd()
{
	Ensure_GKS3rd_Result_Directory();
	Configure_GKS3rd_1D(0.0, 0.0);
	const double accuracy_cfl = 0.05;

	const int mesh_set = 5;
	const int mesh_number[mesh_set] = { 20, 40, 80, 160, 320};
	double error[mesh_set][3]{};

	for (int imesh = 0; imesh < mesh_set; ++imesh)
	{
		Block1d block;
		Prepare_Block_1D(block, mesh_number[imesh], 0.0, 2.0, accuracy_cfl);

		Fluid1d* fluids = new Fluid1d[block.nx];
		Interface1d* interfaces = new Interface1d[block.nx + 1];
		Flux1d** fluxes = Setflux_array(block);
		SetUniformMesh(block, fluids, interfaces, fluxes);
		ICfor_sinwave(fluids, block);

		Fluid1d* bcvalue = new Fluid1d[2];
		Advance_1D_Case(
			fluids, interfaces, fluxes, block,
			periodic_boundary_left, periodic_boundary_right, bcvalue, 2.0);
		Compute_Sinwave_Error(fluids, block, error[imesh]);
		std::string path = "build/result/gks3rd_sinwave_mesh" + std::to_string(mesh_number[imesh]) + ".plt";
		Write_Density_Output_Tecplot(fluids, block, path.c_str());

		delete[] bcvalue;
		for (int i = 0; i <= block.nx; ++i)
		{
			delete[] fluxes[i];
		}
		delete[] fluxes;
		delete[] interfaces;
		delete[] fluids;
	}

	std::cout << "GKS3rd 1D sinwave errors" << std::endl;
	for (int i = 0; i < mesh_set; ++i)
	{
		double order_l1 = 0.0;
		double order_l2 = 0.0;
		double order_linf = 0.0;
		if (i > 0)
		{
			order_l1 = log(error[i - 1][0] / error[i][0]) / log(2.0);
			order_l2 = log(error[i - 1][1] / error[i][1]) / log(2.0);
			order_linf = log(error[i - 1][2] / error[i][2]) / log(2.0);
		}
		std::cout << "mesh=" << mesh_number[i]
			<< " L1=" << error[i][0]
			<< " L2=" << error[i][1]
			<< " Linf=" << error[i][2];
		if (i > 0)
		{
			std::cout << " order(L1)=" << order_l1
				<< " order(L2)=" << order_l2
				<< " order(Linf)=" << order_linf;
		}
		std::cout << std::endl;
	}
}

void riemann_problem_1d_gks3rd()
{
	Ensure_GKS3rd_Result_Directory();
	// Paper Sec. 4 suggests C1 << 1 and C2 = O(1) for Euler computations.
	Configure_GKS3rd_1D(0.1, 1.0);

	Block1d block;
	Prepare_Block_1D(block, 400, 0.0, 1.0, 0.5);

	Fluid1d* fluids = new Fluid1d[block.nx];
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	Flux1d** fluxes = Setflux_array(block);
		SetUniformMesh(block, fluids, interfaces, fluxes);

	ICfor1dRM(fluids, RiemannProblem1D_Sod(), block);

	Fluid1d* bcvalue = new Fluid1d[2];
	Advance_1D_Case(
		fluids, interfaces, fluxes, block,
		free_boundary_left, free_boundary_right, bcvalue, 0.2);
	Write_Sod_Output(fluids, block, StepTaggedPath("gks3rd", block.step, ".plt").c_str());

	delete[] bcvalue;
	for (int i = 0; i <= block.nx; ++i)
	{
		delete[] fluxes[i];
	}
	delete[] fluxes;
	delete[] interfaces;
	delete[] fluids;
}
