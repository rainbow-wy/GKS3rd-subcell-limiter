#include "gks_fr_driver.h"

#include "gks_fr_high_order.h"
#include "function.h"
#include "gks_basic.h"
#include "riemann_problem.h"

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

	void Configure_GKSFR_2D(double c1, double c2)
	{
		K = 3;
		Gamma = 1.4;
		R_gas = 1.0;
		Pr = 0.73;

		tau_type = Euler;
		c1_euler = c1;
		c2_euler = c2;
		Mu = 0.0;
		Nu = -1.0;

		gks2dsolver = gks3rd_2d;
		reconstruction_variable = conservative;
		wenotype = linear;
		is_reduce_order_warning = false;
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

	void SetPrimitiveAtPoint2D(double* Q, double rho, double u, double v, double p)
	{
		double prim[4];
		double convar[4];
		prim[0] = rho;
		prim[1] = u;
		prim[2] = v;
		prim[3] = p;
		Primvar_to_convar_2D(convar, prim);
		for (int m = 0; m < 4; ++m)
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

	void InitializeSinwave2D(GKSFRMesh2D& mesh)
	{
		for (int i = 0; i < mesh.cells_x; ++i)
		{
			for (int j = 0; j < mesh.cells_y; ++j)
			{
				GKSFRCell2D& cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
				for (int p = 0; p < 3; ++p)
				{
					for (int q = 0; q < 3; ++q)
					{
						const double x = GKSFR_SolutionPointX2D(mesh, i, p);
						const double y = GKSFR_SolutionPointY2D(mesh, j, q);
						const double rho = 1.0 + 0.2 * std::sin(kPi * (x + y));
						SetPrimitiveAtPoint2D(cell.Q[p][q], rho, 1.0, 1.0, 1.0);
					}
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

	double GetFRTimeStep2D(const GKSFRMesh2D& mesh, double CFL, double t, double tstop)
	{
		double dt = (mesh.dx < mesh.dy) ? mesh.dx : mesh.dy;
		const double h = dt;
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					double prim[4];
					double convar[4];
					for (int m = 0; m < 4; ++m)
					{
						convar[m] = mesh.cell[e].Q[i][j][m];
					}
					Convar_to_primvar_2D(prim, convar);
					dt = Dtx(dt, h, CFL, prim[0], prim[1], prim[2], prim[3]);
				}
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

	bool CheckPhysicalState2D(const GKSFRMesh2D& mesh, int& bad_e, int& bad_i, int& bad_j)
	{
		double prim[4];
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					double convar[4];
					for (int m = 0; m < 4; ++m)
					{
						convar[m] = mesh.cell[e].Q[i][j][m];
					}
					Convar_to_primvar_2D(prim, convar);
					if (!(convar[0] == convar[0]) || !(prim[3] == prim[3]) || convar[0] <= 0.0 || prim[3] <= 0.0)
					{
						bad_e = e;
						bad_i = i;
						bad_j = j;
						return false;
					}
				}
			}
		}
		bad_e = -1;
		bad_i = -1;
		bad_j = -1;
		return true;
	}

	void WriteCellCenterDensityTecplot(const GKSFRMesh1D& mesh, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,density\n";
		out << "zone i = " << mesh.cells << ", F=POINT\n";
		for (int e = 0; e < mesh.cells; ++e)
		{
			GKSFRCell1D cell = mesh.cell[e];
			GKSFR_ComputeCellCenterData(cell, mesh.dx);
			out << CellCenterX(mesh, e) << " " << cell.Qc[0] << "\n";
		}
	}

	void WriteCellCenterDensityTecplot2D(const GKSFRMesh2D& mesh, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,y,density\n";
		out << "zone i = " << mesh.cells_x << ", j = " << mesh.cells_y << ", F=POINT\n";
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int i = 0; i < mesh.cells_x; ++i)
			{
				GKSFRCell2D cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
				GKSFR_ComputeCellCenterData2D(cell, mesh.dx, mesh.dy);
				out << GKSFR_CellCenterX2D(mesh, i) << " "
					<< GKSFR_CellCenterY2D(mesh, j) << " "
					<< cell.Qc[0] << "\n";
			}
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

	void ComputeSinwaveError2D(const GKSFRMesh2D& mesh, double t, double* error)
	{
		double e1 = 0.0;
		double e2 = 0.0;
		double einf = 0.0;
		int count = 0;
		for (int ei = 0; ei < mesh.cells_x; ++ei)
		{
			for (int ej = 0; ej < mesh.cells_y; ++ej)
			{
				const GKSFRCell2D& cell = mesh.cell[GKSFR_CellIndex2D(mesh, ei, ej)];
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						const double x = GKSFR_SolutionPointX2D(mesh, ei, i);
						const double y = GKSFR_SolutionPointY2D(mesh, ej, j);
						const double exact = 1.0 + 0.2 * std::sin(kPi * (x + y - 2.0 * t));
						const double diff = std::fabs(cell.Q[i][j][0] - exact);
						e1 += diff;
						e2 += diff * diff;
						if (diff > einf)
						{
							einf = diff;
						}
						++count;
					}
				}
			}
		}
		error[0] = e1 / count;
		error[1] = std::sqrt(e2 / count);
		error[2] = einf;
	}

	std::string StepTaggedPath(const std::string& prefix, int final_step, const char* ext)
	{
		return "build/result/" + prefix + "_" + std::to_string(final_step) + ext;
	}

	bool AdvanceCase(
		GKSFRMesh1D& mesh,
		double CFL,
		double tstop,
		GKSFRBoundary1D boundary,
		bool show_step,
		int* final_step = nullptr)
	{
		int step = 0;
		double t = 0.0;
		while (t < tstop - 1e-14)
		{
			const double dt = GetFRTimeStep(mesh, CFL, t, tstop);
			if (show_step)
			{
				std::cout << "step= " << step
					<< " dt= " << dt
					<< " t= " << t << std::endl;
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
				if (final_step != nullptr)
				{
					*final_step = step;
				}
				return false;
			}
		}
		if (final_step != nullptr)
		{
			*final_step = step;
		}
		if (show_step)
		{
			std::cout << "Numerical simulation completed successfully." << std::endl;
		}
		return true;
	}

	bool AdvanceCase2D(
		GKSFRMesh2D& mesh,
		double CFL,
		double tstop,
		GKSFRBoundary2D boundary,
		bool show_step)
	{
		int step = 0;
		double t = 0.0;
		while (t < tstop - 1e-14)
		{
			const double dt = GetFRTimeStep2D(mesh, CFL, t, tstop);
			if (show_step)
			{
				std::cout << "step= " << step
					<< " time size is " << dt
					<< " time= " << t << std::endl;
			}
			GKSFR_AdvanceOneStep2D(mesh, dt, boundary);
			t += dt;
			++step;

			int bad_e = -1;
			int bad_i = -1;
			int bad_j = -1;
			if (!CheckPhysicalState2D(mesh, bad_e, bad_i, bad_j))
			{
				const int bad_cell_i = bad_e / mesh.cells_y;
				const int bad_cell_j = bad_e - bad_cell_i * mesh.cells_y;
				double bad_convar[4];
				double bad_prim[4];
				for (int m = 0; m < 4; ++m)
				{
					bad_convar[m] = mesh.cell[bad_e].Q[bad_i][bad_j][m];
				}
				Convar_to_primvar_2D(bad_prim, bad_convar);
				std::cout << "GKS-FR 2D invalid state at cell=" << bad_e
					<< " cell_ij=(" << bad_cell_i << "," << bad_cell_j << ")"
					<< " point=(" << bad_i << "," << bad_j << ")"
					<< " x=" << GKSFR_SolutionPointX2D(mesh, bad_cell_i, bad_i)
					<< " y=" << GKSFR_SolutionPointY2D(mesh, bad_cell_j, bad_j)
					<< " after step=" << step
					<< " time=" << t << std::endl;
				std::cout << "  primitive rho=" << bad_prim[0]
					<< " u=" << bad_prim[1]
					<< " v=" << bad_prim[2]
					<< " p=" << bad_prim[3] << std::endl;
				std::cout << "  conservative rho=" << bad_convar[0]
					<< " rhou=" << bad_convar[1]
					<< " rhov=" << bad_convar[2]
					<< " rhoE=" << bad_convar[3] << std::endl;
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

		std::string path = "build/result/gksfr_sinwave_mesh" + std::to_string(mesh_number[imesh]) + ".plt";
		WriteCellCenterDensityTecplot(mesh, path.c_str());
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

void accuracy_sinwave_2d_gksfr()
{
	Ensure_GKSFR_Result_Directory();
	Configure_GKSFR_2D(0.0, 0.0);

	const double CFL = 0.05;
	const double tstop = 0.1;
	const int mesh_set = 5;
	const int mesh_number[mesh_set] = { 20, 40, 80, 160, 320};
	double error[mesh_set][5]{};

	for (int imesh = 0; imesh < mesh_set; ++imesh)
	{
		GKSFRMesh2D mesh;
		GKSFR_ResizeUniformMesh2D(mesh, mesh_number[imesh], mesh_number[imesh], 0.0, 2.0, 0.0, 2.0);
		InitializeSinwave2D(mesh);
		std::cout << "GKS-FR 2D sinwave mesh="
			<< mesh_number[imesh] << "x" << mesh_number[imesh] << std::endl;
		const bool ok = AdvanceCase2D(mesh, CFL, tstop, gksfr2d_periodic, true);
		if (!ok)
		{
			std::cout << "GKS-FR 2D sinwave failed on mesh "
				<< mesh_number[imesh] << "x" << mesh_number[imesh] << std::endl;
			return;
		}
		ComputeSinwaveError2D(mesh, tstop, error[imesh]);

		std::string path = "build/result/gksfr_2d_sinwave_mesh" + std::to_string(mesh_number[imesh]) + ".plt";
		WriteCellCenterDensityTecplot2D(mesh, path.c_str());
	}

	std::cout << "GKS-FR 2D sinwave errors" << std::endl;
	std::cout << "Nx, Ny, L1 error, L2 error, Linf error, observed order" << std::endl;
	for (int i = 0; i < mesh_set; ++i)
	{
		double order_l1 = 0.0;
		double order_l2 = 0.0;
		double order_linf = 0.0;
		if (i > 0)
		{
			order_l1 = std::log(error[i - 1][0] / error[i][0]) / std::log(2.0);
			order_l2 = std::log(error[i - 1][1] / error[i][1]) / std::log(2.0);
			order_linf = std::log(error[i - 1][2] / error[i][2]) / std::log(2.0);
		}
		std::cout << mesh_number[i] << ", "
			<< mesh_number[i] << ", "
			<< error[i][0] << ", "
			<< error[i][1] << ", "
			<< error[i][2] << ", ";
		if (i == 0)
		{
			std::cout << "-";
		}
		else
		{
			std::cout << "L1=" << order_l1
				<< " L2=" << order_l2
				<< " Linf=" << order_linf;
		}
		std::cout << std::endl;
	}
}

void riemann_problem_1d_gksfr()
{
	Ensure_GKSFR_Result_Directory();
	Configure_GKSFR_1D(0.1, 1.0);

	GKSFRMesh1D mesh;
	GKSFR_ResizeUniformMesh(mesh, 400, -5.0, 5.0);
	ICfor1dRM(mesh, RiemannProblem1D_Sod());
	//ICfor1dRM(mesh, RiemannProblem1D_ShuOsher());

	int final_step = 0;
	const bool ok = AdvanceCase(mesh, 0.02, 1.8, gksfr_free, true, &final_step);
	if (!ok)
	{
		std::cout << "GKS-FR Sod run stopped because of an invalid state." << std::endl;
		return;
	}

	WriteSodTecplot(mesh, StepTaggedPath("gksfr", final_step, ".plt").c_str());
}
