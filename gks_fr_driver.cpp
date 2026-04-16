#include "gks_fr_driver.h"

#include "gks_fr_high_order.h"
#include "function.h"
#include "gks_basic.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace
{
	const double kPi = 3.14159265358979323846;

	void Ensure_GKSFR_Result_Directory()
	{
#if defined(_WIN32)
		_mkdir("build");
		_mkdir("build\\result");
#else
		mkdir("build", 0777);
		mkdir("build/result", 0777);
#endif
	}

	void Configure_GKSFR_1D(double c1, double c2)
	{
		K = 4;
		Gamma = 1.4;
		R_gas = 1.0;

		tau_type = Euler;
		c1_euler = c1;
		c2_euler = c2;

		flux_function = GKS;
		gks1dsolver = gks3rd;
	}

	double CellCenterX(const GKSFRMesh1D& mesh, int e)
	{
		return mesh.x_left + (e + 0.5) * mesh.dx;
	}

	double SolutionPointX(const GKSFRMesh1D& mesh, int e, int i)
	{
		return CellCenterX(mesh, e) + 0.5 * mesh.dx * GKSFR_GL_Point(i);
	}

	void SetPrimitiveAtPoint(double* Q, double rho, double u, double p)
	{
		double prim[3];
		double convar[3];
		prim[0] = rho;
		prim[1] = u;
		prim[2] = p;
		Primvar_to_convar_1D(convar, prim);
		for (int m = 0; m < 3; ++m)
		{
			Q[m] = convar[m];
		}
	}

	void InitializeSinwave(GKSFRMesh1D& mesh)
	{
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				const double x = SolutionPointX(mesh, e, i);
				const double rho = 1.0 + 0.2 * std::sin(kPi * x);
				SetPrimitiveAtPoint(mesh.cell[e].Q[i], rho, 1.0, 1.0);
			}
		}
	}

	void InitializeSod(GKSFRMesh1D& mesh)
	{
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				const double x = SolutionPointX(mesh, e, i);
				if (x < 0.5)
				{
					SetPrimitiveAtPoint(mesh.cell[e].Q[i], 1.0, 0.0, 1.0);
				}
				else
				{
					SetPrimitiveAtPoint(mesh.cell[e].Q[i], 0.125, 0.0, 0.1);
				}
			}
		}
	}

	double GetFRTimeStep(const GKSFRMesh1D& mesh, double CFL, double t, double tstop)
	{
		double dt = mesh.dx;
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				double convar[3];
				for (int m = 0; m < 3; ++m)
				{
					convar[m] = mesh.cell[e].Q[i][m];
				}
				dt = Dtx(dt, mesh.dx, CFL, convar);
			}
		}
		if (t + dt > tstop)
		{
			dt = tstop - t;
		}
		return dt;
	}

	bool CheckPhysicalState(const GKSFRMesh1D& mesh, int& bad_e, int& bad_i)
	{
		double prim[3];
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				double convar[3];
				for (int m = 0; m < 3; ++m)
				{
					convar[m] = mesh.cell[e].Q[i][m];
				}
				Convar_to_primvar_1D(prim, convar);
				if (!(convar[0] == convar[0]) || !(prim[2] == prim[2]) || convar[0] <= 0.0 || prim[2] <= 0.0)
				{
					bad_e = e;
					bad_i = i;
					return false;
				}
			}
		}
		bad_e = -1;
		bad_i = -1;
		return true;
	}

	void WriteCellCenterDensity(const GKSFRMesh1D& mesh, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "# x rho\n";
		for (int e = 0; e < mesh.cells; ++e)
		{
			GKSFRCell1D cell = mesh.cell[e];
			GKSFR_ComputeCellCenterData(cell, mesh.dx);
			out << CellCenterX(mesh, e) << " " << cell.Qc[0] << "\n";
		}
	}

	void WriteSodTecplot(const GKSFRMesh1D& mesh, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,density,u,pressure,temperature,entropy,Ma\n";
		out << "zone i = " << mesh.cells << ", F=POINT\n";

		double prim[3];
		for (int e = 0; e < mesh.cells; ++e)
		{
			GKSFRCell1D cell = mesh.cell[e];
			GKSFR_ComputeCellCenterData(cell, mesh.dx);
			double convar[3];
			for (int m = 0; m < 3; ++m)
			{
				convar[m] = cell.Qc[m];
			}
			Convar_to_primvar_1D(prim, convar);
			out << CellCenterX(mesh, e) << " "
				<< prim[0] << " "
				<< prim[1] << " "
				<< prim[2] << " "
				<< Temperature(prim[0], prim[2]) << " "
				<< entropy(prim[0], prim[2]) << " "
				<< std::fabs(prim[1]) / Soundspeed(prim[0], prim[2]) << "\n";
		}
	}

	void ComputeSinwaveError(const GKSFRMesh1D& mesh, double t, double* error)
	{
		double e1 = 0.0;
		double e2 = 0.0;
		double einf = 0.0;
		int count = 0;
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				const double x = SolutionPointX(mesh, e, i);
				const double exact = 1.0 + 0.2 * std::sin(kPi * (x - t));
				const double diff = std::fabs(mesh.cell[e].Q[i][0] - exact);
				e1 += diff;
				e2 += diff * diff;
				if (diff > einf)
				{
					einf = diff;
				}
				++count;
			}
		}
		error[0] = e1 / count;
		error[1] = std::sqrt(e2 / count);
		error[2] = einf;
	}

	bool AdvanceCase(
		GKSFRMesh1D& mesh,
		double CFL,
		double tstop,
		GKSFRBoundary1D boundary,
		bool show_step)
	{
		int step = 0;
		double t = 0.0;
		while (t < tstop - 1e-14)
		{
			const double dt = GetFRTimeStep(mesh, CFL, t, tstop);
			if (show_step)
			{
				std::cout << "step= " << step
					<< " time size is " << dt
					<< " time= " << t << std::endl;
			}
			GKSFR_AdvanceOneStep(mesh, dt, boundary);
			t += dt;
			++step;

			int bad_e = -1;
			int bad_i = -1;
			if (!CheckPhysicalState(mesh, bad_e, bad_i))
			{
				std::cout << "GKS-FR invalid state at cell=" << bad_e
					<< " point=" << bad_i
					<< " after step=" << step
					<< " time=" << t << std::endl;
				return false;
			}
		}
		return true;
	}
}

void accuracy_sinwave_1d_gksfr()
{
	Ensure_GKSFR_Result_Directory();
	Configure_GKSFR_1D(0.0, 0.0);

	const double CFL = 0.02;
	const double tstop = 2.0;
	const int mesh_set = 5;
	const int mesh_number[mesh_set] = { 20, 40, 80, 160, 320 };
	double error[mesh_set][3]{};

	for (int imesh = 0; imesh < mesh_set; ++imesh)
	{
		GKSFRMesh1D mesh;
		GKSFR_ResizeUniformMesh(mesh, mesh_number[imesh], 0.0, 2.0);
		InitializeSinwave(mesh);
		const bool ok = AdvanceCase(mesh, CFL, tstop, gksfr_periodic, false);
		if (!ok)
		{
			std::cout << "GKS-FR sinwave failed on mesh " << mesh_number[imesh] << std::endl;
			return;
		}
		ComputeSinwaveError(mesh, tstop, error[imesh]);

		std::string path = "build/result/gksfr_sinwave_mesh" + std::to_string(mesh_number[imesh]) + ".dat";
		WriteCellCenterDensity(mesh, path.c_str());
	}

	std::cout << "GKS-FR 1D sinwave errors" << std::endl;
	for (int i = 0; i < mesh_set; ++i)
	{
		std::cout << "mesh=" << mesh_number[i]
			<< " L1=" << error[i][0]
			<< " L2=" << error[i][1]
			<< " Linf=" << error[i][2];
		if (i > 0)
		{
			std::cout << " order(L1)=" << std::log(error[i - 1][0] / error[i][0]) / std::log(2.0)
				<< " order(L2)=" << std::log(error[i - 1][1] / error[i][1]) / std::log(2.0)
				<< " order(Linf)=" << std::log(error[i - 1][2] / error[i][2]) / std::log(2.0);
		}
		std::cout << std::endl;
	}
}

void riemann_problem_1d_gksfr()
{
	Ensure_GKSFR_Result_Directory();
	Configure_GKSFR_1D(0.1, 1.0);

	GKSFRMesh1D mesh;
	GKSFR_ResizeUniformMesh(mesh, 200, 0.0, 1.0);
	InitializeSod(mesh);

	const bool ok = AdvanceCase(mesh, 0.02, 0.2, gksfr_free, true);
	if (!ok)
	{
		std::cout << "GKS-FR Sod run stopped because of an invalid state." << std::endl;
		return;
	}

	WriteSodTecplot(mesh, "build/result/gksfr_sod_t020.plt");
}
