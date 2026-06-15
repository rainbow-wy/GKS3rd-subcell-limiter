#include "gks_subcell_blend_driver.h"

#include "function.h"
#include "gks_basic.h"
#include "gks_fr_adapter.h"
#include "gks_kfvs_adapter.h"
#include "riemann_problem.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#endif

GKSSubcellFrameworkConfig1D::GKSSubcellFrameworkConfig1D()
	: blend_mode(gks_subcell_hybrid),
	  low_mode(KFVS1),
	  use_flux_limiter(true),
	  use_scaling_limiter(true)
{
}

GKSSubcellFrameworkDiag1D::GKSSubcellFrameworkDiag1D()
	: min_rho(1.0e300),
	  min_p(1.0e300),
	  min_rho_cell(-1),
	  min_rho_point(-1),
	  min_p_cell(-1),
	  min_p_point(-1)
{
}

GKSSubcellFrameworkConfig2D::GKSSubcellFrameworkConfig2D()
	: blend_mode(gks_subcell2d_hybrid),
	  low_mode(KFVS1),
	  use_flux_limiter(true),
	  use_scaling_limiter(true)
{
}

GKSSubcellMask2D::GKSSubcellMask2D()
	: enabled(false),
	  cells_x(0),
	  cells_y(0)
{
}

GKSSubcellFrameworkDiag2D::GKSSubcellFrameworkDiag2D()
	: min_rho(1.0e300),
	  min_p(1.0e300),
	  max_alpha(0.0),
	  troubled_cells(0)
{
}

namespace
{
	const double kPi = 3.14159265358979323846;

	void Ensure_Result_Directory()
	{
#if defined(_WIN32)
		_mkdir("build");
		_mkdir("build\\result");
#else
		mkdir("build", 0777);
		mkdir("build/result", 0777);
#endif
	}

	void Ensure_Debug_Directory(const std::string& dir)
	{
		if (dir.empty())
		{
			return;
		}
#if defined(_WIN32)
		_mkdir(dir.c_str());
#else
		mkdir(dir.c_str(), 0777);
#endif
	}

	std::string GetEnvString(const char* name, const std::string& default_value = "")
	{
		const char* value = std::getenv(name);
		return value ? std::string(value) : default_value;
	}

	int GetEnvInt(const char* name, int default_value)
	{
		const char* value = std::getenv(name);
		return value ? std::atoi(value) : default_value;
	}

	double GetEnvDouble(const char* name, double default_value)
	{
		const char* value = std::getenv(name);
		return value ? std::atof(value) : default_value;
	}

	std::vector<int> ParseStepList(const std::string& text)
	{
		std::vector<int> steps;
		std::stringstream ss(text);
		std::string item;
		while (std::getline(ss, item, ','))
		{
			if (!item.empty())
			{
				steps.push_back(std::atoi(item.c_str()));
			}
		}
		return steps;
	}

	bool StepInList(const std::vector<int>& steps, int step)
	{
		return std::find(steps.begin(), steps.end(), step) != steps.end();
	}

	std::string SanitizedTimeTag(double t)
	{
		std::ostringstream ss;
		ss << std::setprecision(10) << t;
		std::string s = ss.str();
		for (std::size_t i = 0; i < s.size(); ++i)
		{
			if (s[i] == '.')
			{
				s[i] = 'p';
			}
			else if (s[i] == '-')
			{
				s[i] = 'm';
			}
			else if (s[i] == '+')
			{
				s[i] = 'p';
			}
		}
		return s;
	}

	std::string JetCheckpointPath(const std::string& dir, int step, double t)
	{
		std::ostringstream path;
		path << dir << "/debug_jet_chk_step" << step << "_t" << SanitizedTimeTag(t) << ".bin";
		return path.str();
	}

	bool SaveJetCheckpoint2D(const GKSFRMesh2D& mesh, int step, double t, double last_dt, const std::string& path)
	{
		std::ofstream out(path.c_str(), std::ios::binary);
		if (!out)
		{
			std::cout << "cannot write jet checkpoint: " << path << std::endl;
			return false;
		}
		const char magic[16] = { 'G','K','S','J','E','T','C','H','K','2','D','0','0','0','1',0 };
		const int version = 1;
		out.write(magic, sizeof(magic));
		out.write(reinterpret_cast<const char*>(&version), sizeof(version));
		out.write(reinterpret_cast<const char*>(&mesh.cells_x), sizeof(mesh.cells_x));
		out.write(reinterpret_cast<const char*>(&mesh.cells_y), sizeof(mesh.cells_y));
		out.write(reinterpret_cast<const char*>(&mesh.x_left), sizeof(mesh.x_left));
		out.write(reinterpret_cast<const char*>(&mesh.x_right), sizeof(mesh.x_right));
		out.write(reinterpret_cast<const char*>(&mesh.y_bottom), sizeof(mesh.y_bottom));
		out.write(reinterpret_cast<const char*>(&mesh.y_top), sizeof(mesh.y_top));
		out.write(reinterpret_cast<const char*>(&mesh.dx), sizeof(mesh.dx));
		out.write(reinterpret_cast<const char*>(&mesh.dy), sizeof(mesh.dy));
		out.write(reinterpret_cast<const char*>(&step), sizeof(step));
		out.write(reinterpret_cast<const char*>(&t), sizeof(t));
		out.write(reinterpret_cast<const char*>(&last_dt), sizeof(last_dt));
		const int vars = 4;
		const int points = 9;
		out.write(reinterpret_cast<const char*>(&vars), sizeof(vars));
		out.write(reinterpret_cast<const char*>(&points), sizeof(points));
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			out.write(reinterpret_cast<const char*>(mesh.cell[e].Q), sizeof(mesh.cell[e].Q));
		}
		if (!out)
		{
			std::cout << "failed while writing jet checkpoint: " << path << std::endl;
			return false;
		}
		std::cout << "wrote jet checkpoint: " << path << std::endl;
		return true;
	}

	bool LoadJetCheckpoint2D(GKSFRMesh2D& mesh, int& step, double& t, double& last_dt, const std::string& path)
	{
		std::ifstream in(path.c_str(), std::ios::binary);
		if (!in)
		{
			std::cout << "cannot read jet checkpoint: " << path << std::endl;
			return false;
		}
		char magic[16];
		int version = 0, nx = 0, ny = 0, vars = 0, points = 0;
		double x_left = 0.0, x_right = 0.0, y_bottom = 0.0, y_top = 0.0, dx = 0.0, dy = 0.0;
		in.read(magic, sizeof(magic));
		in.read(reinterpret_cast<char*>(&version), sizeof(version));
		in.read(reinterpret_cast<char*>(&nx), sizeof(nx));
		in.read(reinterpret_cast<char*>(&ny), sizeof(ny));
		in.read(reinterpret_cast<char*>(&x_left), sizeof(x_left));
		in.read(reinterpret_cast<char*>(&x_right), sizeof(x_right));
		in.read(reinterpret_cast<char*>(&y_bottom), sizeof(y_bottom));
		in.read(reinterpret_cast<char*>(&y_top), sizeof(y_top));
		in.read(reinterpret_cast<char*>(&dx), sizeof(dx));
		in.read(reinterpret_cast<char*>(&dy), sizeof(dy));
		in.read(reinterpret_cast<char*>(&step), sizeof(step));
		in.read(reinterpret_cast<char*>(&t), sizeof(t));
		in.read(reinterpret_cast<char*>(&last_dt), sizeof(last_dt));
		in.read(reinterpret_cast<char*>(&vars), sizeof(vars));
		in.read(reinterpret_cast<char*>(&points), sizeof(points));
		const char expected[16] = { 'G','K','S','J','E','T','C','H','K','2','D','0','0','0','1',0 };
		if (!in || std::string(magic, magic + sizeof(magic)) != std::string(expected, expected + sizeof(expected)) ||
			version != 1 || nx != mesh.cells_x || ny != mesh.cells_y || vars != 4 || points != 9 ||
			std::fabs(x_left - mesh.x_left) > 1.0e-14 || std::fabs(x_right - mesh.x_right) > 1.0e-14 ||
			std::fabs(y_bottom - mesh.y_bottom) > 1.0e-14 || std::fabs(y_top - mesh.y_top) > 1.0e-14 ||
			std::fabs(dx - mesh.dx) > 1.0e-14 || std::fabs(dy - mesh.dy) > 1.0e-14)
		{
			std::cout << "jet checkpoint header mismatch: " << path << std::endl;
			return false;
		}
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			in.read(reinterpret_cast<char*>(mesh.cell[e].Q), sizeof(mesh.cell[e].Q));
		}
		if (!in)
		{
			std::cout << "failed while reading jet checkpoint: " << path << std::endl;
			return false;
		}
		std::cout << "loaded jet checkpoint: " << path
			<< " step=" << step << " t=" << t << " last_dt=" << last_dt << std::endl;
		return true;
	}

	void Configure_GKS_Subcell_1D(double c1, double c2)
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

	void Configure_GKS_Subcell_2D(double c1, double c2)
	{
		K = 3;
		Gamma = 1.4;
		R_gas = 1.0;
		Pr = 0.73;
		tau_type = Euler;
		Mu = 0.0;
		Nu = -1.0;
		c1_euler = c1;
		c2_euler = c2;
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
		double prim[3]{ rho, u, p };
		Primvar_to_convar_1D(Q, prim);
	}

	void SetPrimitiveAtPoint2D(double* Q, double rho, double u, double v, double p)
	{
		double prim[4]{ rho, u, v, p };
		Primvar_to_convar_2D(Q, prim);
	}

	void InitializeSinwave(GKSFRMesh1D& mesh)
	{
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				const double x = SolutionPointX(mesh, e, i);
				SetPrimitiveAtPoint(mesh.cell[e].Q[i], 1.0 + 0.2 * std::sin(kPi * x), 1.0, 1.0);
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
						SetPrimitiveAtPoint2D(cell.Q[p][q], 1.0 + 0.2 * std::sin(kPi * (x + y)), 1.0, 1.0, 1.0);
					}
				}
			}
		}
	}

	double GetTimeStep(const GKSFRMesh1D& mesh, double CFL, double t, double tstop)
	{
		double dt = mesh.dx;
		for (int e = 0; e < mesh.cells; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				dt = Dtx(dt, mesh.dx, CFL, const_cast<double*>(mesh.cell[e].Q[i]));
			}
		}
		if (t + dt > tstop)
		{
			dt = tstop - t;
		}
		return dt;
	}

	double GetTimeStep2D(const GKSFRMesh2D& mesh, double CFL, double t, double tstop)
	{
		double dt = std::min(mesh.dx, mesh.dy);
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

	double GetTimeStepAstrophysicalJet2D(const GKSFRMesh2D& mesh, double CFL, double t, double tstop)
	{
		double dt = GetTimeStep2D(mesh, CFL, t, tstop);
		const double subcell_weight_min = 5.0 / 18.0;
		const double kx = 0.5;
		const double ky = 0.5;
		const double effective_h = std::min(kx * subcell_weight_min * mesh.dx,ky * subcell_weight_min * mesh.dy);
		double jet_prim[4];
		GKSFR_AstrophysicalJetPrimitive2D(jet_prim, mesh.x_left, 0.0, t);
		const double max_jet_speed =
			std::sqrt(jet_prim[1] * jet_prim[1] + jet_prim[2] * jet_prim[2]) +
			Soundspeed(jet_prim[0], jet_prim[3]);
		if (max_jet_speed > 0.0)
		{
			const double jet_safety_factor = 0.5;
			dt = std::min(dt, jet_safety_factor * CFL * effective_h / max_jet_speed);
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
				Convar_to_primvar_1D(prim, const_cast<double*>(mesh.cell[e].Q[i]));
				if (!(mesh.cell[e].Q[i][0] == mesh.cell[e].Q[i][0]) ||
					!(prim[2] == prim[2]) ||
					mesh.cell[e].Q[i][0] <= 0.0 ||
					prim[2] <= 0.0)
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

	double Pressure2DLocal(const double U[4])
	{
		if (U[0] <= 0.0)
		{
			return -1.0;
		}
		return (Gamma - 1.0) *
			(U[3] - 0.5 * (U[1] * U[1] + U[2] * U[2]) / U[0]);
	}

	int g_jet_trace_step = -1;
	double g_jet_trace_time = 0.0;
	double g_jet_trace_dt = 0.0;

	void SetJetTraceContext(int step, double t, double dt)
	{
		g_jet_trace_step = step;
		g_jet_trace_time = t;
		g_jet_trace_dt = dt;
	}

	bool JetTraceEnabled()
	{
		return GetEnvInt("JET_TRACE_CELL", -1) >= 0;
	}

	int JetTraceCell()
	{
		return GetEnvInt("JET_TRACE_CELL", -1);
	}

	int JetTraceRadius()
	{
		return std::max(0, GetEnvInt("JET_TRACE_RADIUS", 2));
	}

	bool IsFinite4(const double U[4])
	{
		for (int m = 0; m < 4; ++m)
		{
			if (!std::isfinite(U[m]))
			{
				return false;
			}
		}
		return true;
	}

	bool StateBadForTrace(const double U[4])
	{
		const double p = Pressure2DLocal(U);
		return !IsFinite4(U) || !std::isfinite(p) || U[0] <= 0.0 || p <= 0.0;
	}

	void ComputeAverageFromPoints2D(const GKSFRCell2D& cell, double avg[4])
	{
		const double w[3] = { 5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0 };
		for (int m = 0; m < 4; ++m)
		{
			avg[m] = 0.0;
		}
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					avg[m] += w[i] * w[j] * cell.Q[i][j][m];
				}
			}
		}
	}

	bool JetCoreCell2D(const GKSFRMesh2D& mesh, int i, int j)
	{
		const double x = GKSFR_CellCenterX2D(mesh, i);
		const double y = GKSFR_CellCenterY2D(mesh, j);
		return x - mesh.x_left <= 0.15 && std::fabs(y) <= 0.12;
	}

	void PrintTraceState2D(const char* phase, const GKSFRMesh2D& mesh, int e, int p_id, int q_id, const double U[4])
	{
		const int ci = e / mesh.cells_y;
		const int cj = e % mesh.cells_y;
		const double x = (p_id >= 0) ? GKSFR_SolutionPointX2D(mesh, ci, p_id) : GKSFR_CellCenterX2D(mesh, ci);
		const double y = (q_id >= 0) ? GKSFR_SolutionPointY2D(mesh, cj, q_id) : GKSFR_CellCenterY2D(mesh, cj);
		const double rho = U[0];
		const double u = (std::isfinite(rho) && std::fabs(rho) > 0.0) ? U[1] / rho : 0.0;
		const double v = (std::isfinite(rho) && std::fabs(rho) > 0.0) ? U[2] / rho : 0.0;
		const double p = Pressure2DLocal(U);
		std::cout << "JET_TRACE phase=" << phase
			<< " step=" << g_jet_trace_step << " t=" << g_jet_trace_time << " dt=" << g_jet_trace_dt
			<< " cell=" << e << " i=" << ci << " j=" << cj
			<< " point=(" << p_id << "," << q_id << ")"
			<< " x=" << x << " y=" << y
			<< " core=" << (JetCoreCell2D(mesh, ci, cj) ? 1 : 0)
			<< " rho=" << U[0] << " rhoU=" << U[1] << " rhoV=" << U[2] << " rhoE=" << U[3]
			<< " u=" << u << " v=" << v << " p=" << p << std::endl;
	}

	void PrintNeighborAverages2D(const GKSFRMesh2D& mesh, int center_i, int center_j)
	{
		const int di[5] = { 0, -1, 1, 0, 0 };
		const int dj[5] = { 0, 0, 0, -1, 1 };
		for (int k = 0; k < 5; ++k)
		{
			const int i = center_i + di[k];
			const int j = center_j + dj[k];
			if (i < 0 || i >= mesh.cells_x || j < 0 || j >= mesh.cells_y)
			{
				continue;
			}
			double avg[4];
			const int e = GKSFR_CellIndex2D(mesh, i, j);
			ComputeAverageFromPoints2D(mesh.cell[e], avg);
			std::cout << "JET_TRACE neighbor cell=" << e << " i=" << i << " j=" << j
				<< " avg_rho=" << avg[0] << " avg_p=" << Pressure2DLocal(avg) << " avg_E=" << avg[3] << std::endl;
		}
	}

	bool JetTraceCellInNeighborhood(const GKSFRMesh2D& mesh, int e)
	{
		const int target = JetTraceCell();
		if (target < 0 || target >= mesh.cells_x * mesh.cells_y)
		{
			return false;
		}
		const int ti = target / mesh.cells_y;
		const int tj = target % mesh.cells_y;
		const int ci = e / mesh.cells_y;
		const int cj = e % mesh.cells_y;
		const int r = JetTraceRadius();
		return std::abs(ci - ti) <= r && std::abs(cj - tj) <= r;
	}

	void JetTraceMeshBadStates2D(const char* phase, const GKSFRMesh2D& mesh, const std::vector<double>* alpha = nullptr)
	{
		if (!JetTraceEnabled())
		{
			return;
		}
		bool printed_neighbor = false;
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			if (!JetTraceCellInNeighborhood(mesh, e))
			{
				continue;
			}
			double min_rho = 1.0e300, min_p = 1.0e300;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					const double* U = mesh.cell[e].Q[i][j];
					min_rho = std::min(min_rho, U[0]);
					min_p = std::min(min_p, Pressure2DLocal(U));
					if (StateBadForTrace(U))
					{
						PrintTraceState2D(phase, mesh, e, i, j, U);
						std::cout << "JET_TRACE cell_min rho=" << min_rho << " p=" << min_p;
						if (alpha != nullptr && e < static_cast<int>(alpha->size()))
						{
							std::cout << " alpha=" << (*alpha)[e];
						}
						std::cout << std::endl;
						if (!printed_neighbor)
						{
							PrintNeighborAverages2D(mesh, e / mesh.cells_y, e % mesh.cells_y);
							printed_neighbor = true;
						}
					}
				}
			}
		}
	}

	void JetTraceBranchBadStates2D(const char* phase, const GKSFRMesh2D& mesh, const GKSSubcellBranch2D& branch, bool check_new)
	{
		if (!JetTraceEnabled())
		{
			return;
		}
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			if (!JetTraceCellInNeighborhood(mesh, e))
			{
				continue;
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					const double* U = check_new ? branch.cell[e].low_dof_new[i][j] : branch.cell[e].low_dof[i][j];
					if (StateBadForTrace(U))
					{
						PrintTraceState2D(phase, mesh, e, i, j, U);
					}
				}
			}
			if (check_new && StateBadForTrace(branch.cell[e].cell_avg_new))
			{
				PrintTraceState2D("safe_average", mesh, e, -1, -1, branch.cell[e].cell_avg_new);
			}
		}
	}

	void JetTraceFaceFluxes2D(
		const char* phase,
		const GKSFRMesh2D& mesh,
		const std::vector<GKSFRFaceFlux2D>& x_fluxes,
		const std::vector<GKSFRFaceFlux2D>& y_fluxes,
		const GKSFluxLimiterDiag2D* flux_diag = nullptr)
	{
		if (!JetTraceEnabled())
		{
			return;
		}
		const int target = JetTraceCell();
		if (target < 0 || target >= mesh.cells_x * mesh.cells_y)
		{
			return;
		}
		const int ti = target / mesh.cells_y;
		const int tj = target % mesh.cells_y;
		if (static_cast<int>(x_fluxes.size()) < (mesh.cells_x + 1) * mesh.cells_y ||
			static_cast<int>(y_fluxes.size()) < mesh.cells_x * (mesh.cells_y + 1))
		{
			return;
		}
		for (int face_i = std::max(0, ti - JetTraceRadius()); face_i <= std::min(mesh.cells_x, ti + JetTraceRadius() + 1); ++face_i)
		{
			for (int j = std::max(0, tj - JetTraceRadius()); j <= std::min(mesh.cells_y - 1, tj + JetTraceRadius()); ++j)
			{
				const int face = face_i * mesh.cells_y + j;
				for (int q = 0; q < 3; ++q)
				{
					if (!IsFinite4(x_fluxes[face].F[q]))
					{
						const int d = face * 3 + q;
						std::cout << "JET_TRACE phase=" << phase << " nonfinite x_flux face_i=" << face_i
							<< " j=" << j << " q=" << q;
						if (flux_diag != nullptr && d < static_cast<int>(flux_diag->alpha_face_x.size()))
						{
							std::cout << " alpha_face=" << flux_diag->alpha_face_x[d] << " theta=" << flux_diag->theta_x[d];
						}
						std::cout << std::endl;
					}
				}
			}
		}
		for (int i = std::max(0, ti - JetTraceRadius()); i <= std::min(mesh.cells_x - 1, ti + JetTraceRadius()); ++i)
		{
			for (int face_j = std::max(0, tj - JetTraceRadius()); face_j <= std::min(mesh.cells_y, tj + JetTraceRadius() + 1); ++face_j)
			{
				const int face = i * (mesh.cells_y + 1) + face_j;
				for (int p = 0; p < 3; ++p)
				{
					if (!IsFinite4(y_fluxes[face].F[p]))
					{
						const int d = face * 3 + p;
						std::cout << "JET_TRACE phase=" << phase << " nonfinite y_flux i=" << i
							<< " face_j=" << face_j << " p=" << p;
						if (flux_diag != nullptr && d < static_cast<int>(flux_diag->alpha_face_y.size()))
						{
							std::cout << " alpha_face=" << flux_diag->alpha_face_y[d] << " theta=" << flux_diag->theta_y[d];
						}
						std::cout << std::endl;
					}
				}
			}
		}
	}

	bool MaskCellActive2D(const GKSSubcellMask2D* mask, const GKSFRMesh2D& mesh, int i, int j)
	{
		if (mask == nullptr || !mask->enabled)
		{
			return true;
		}
		if (i < 0 || i >= mesh.cells_x || j < 0 || j >= mesh.cells_y)
		{
			return false;
		}
		const int e = GKSFR_CellIndex2D(mesh, i, j);
		if (e < 0 || e >= static_cast<int>(mask->active.size()))
		{
			return false;
		}
		return mask->active[e] != 0;
	}

	bool MaskCellActive2D(const GKSSubcellMask2D* mask, const GKSSubcellBranch2D& branch, int i, int j)
	{
		if (mask == nullptr || !mask->enabled)
		{
			return true;
		}
		if (i < 0 || i >= branch.cells_x || j < 0 || j >= branch.cells_y)
		{
			return false;
		}
		const int e = i * branch.cells_y + j;
		if (e < 0 || e >= static_cast<int>(mask->active.size()))
		{
			return false;
		}
		return mask->active[e] != 0;
	}

	void ReflectState2D(double out_Q[4], const double in_Q[4], bool x_normal)
	{
		double prim[4];
		double local_Q[4];
		for (int m = 0; m < 4; ++m)
		{
			local_Q[m] = in_Q[m];
		}
		Convar_to_primvar_2D(prim, local_Q);
		if (x_normal)
		{
			prim[1] = -prim[1];
		}
		else
		{
			prim[2] = -prim[2];
		}
		Primvar_to_convar_2D(out_Q, prim);
	}

	void OverrideMaskedWallFaceFluxes2D(
		const GKSSubcellMask2D* mask,
		const GKSSubcellBranch2D& branch,
		double dt,
		std::vector<GKSFRFaceFlux2D>& x_face_fluxes,
		std::vector<GKSFRFaceFlux2D>& y_face_fluxes)
	{
		if (mask == nullptr || !mask->enabled)
		{
			return;
		}

		for (int face_i = 1; face_i < branch.cells_x; ++face_i)
		{
			for (int j = 0; j < branch.cells_y; ++j)
			{
				const bool left_active = MaskCellActive2D(mask, branch, face_i - 1, j);
				const bool right_active = MaskCellActive2D(mask, branch, face_i, j);
				if (left_active == right_active)
				{
					continue;
				}
				const int face_index = face_i * branch.cells_y + j;
				for (int q = 0; q < 3; ++q)
				{
					double ghost_Q[4];
					if (left_active)
					{
						const GKSSubcellCell2D& left = branch.cell[GKSSubcellCellIndex2D(branch, face_i - 1, j)];
						ReflectState2D(ghost_Q, left.low_dof[2][q], true);
						KFVS1_TimeAveragedFlux2D_X(left.low_dof[2][q], ghost_Q, dt, x_face_fluxes[face_index].F[q]);
					}
					else
					{
						const GKSSubcellCell2D& right = branch.cell[GKSSubcellCellIndex2D(branch, face_i, j)];
						ReflectState2D(ghost_Q, right.low_dof[0][q], true);
						KFVS1_TimeAveragedFlux2D_X(ghost_Q, right.low_dof[0][q], dt, x_face_fluxes[face_index].F[q]);
					}
				}
			}
		}

		for (int i = 0; i < branch.cells_x; ++i)
		{
			for (int face_j = 1; face_j < branch.cells_y; ++face_j)
			{
				const bool bottom_active = MaskCellActive2D(mask, branch, i, face_j - 1);
				const bool top_active = MaskCellActive2D(mask, branch, i, face_j);
				if (bottom_active == top_active)
				{
					continue;
				}
				const int face_index = i * (branch.cells_y + 1) + face_j;
				for (int p = 0; p < 3; ++p)
				{
					double ghost_Q[4];
					if (bottom_active)
					{
						const GKSSubcellCell2D& bottom = branch.cell[GKSSubcellCellIndex2D(branch, i, face_j - 1)];
						ReflectState2D(ghost_Q, bottom.low_dof[p][2], false);
						KFVS1_TimeAveragedFlux2D_Y(bottom.low_dof[p][2], ghost_Q, dt, y_face_fluxes[face_index].F[p]);
					}
					else
					{
						const GKSSubcellCell2D& top = branch.cell[GKSSubcellCellIndex2D(branch, i, face_j)];
						ReflectState2D(ghost_Q, top.low_dof[p][0], false);
						KFVS1_TimeAveragedFlux2D_Y(ghost_Q, top.low_dof[p][0], dt, y_face_fluxes[face_index].F[p]);
					}
				}
			}
		}
	}

	void RestoreInactiveCells2D(
		const GKSSubcellMask2D* mask,
		const GKSFRMesh2D& old_mesh,
		GKSFRMesh2D& new_mesh)
	{
		if (mask == nullptr || !mask->enabled)
		{
			return;
		}
		for (int i = 0; i < new_mesh.cells_x; ++i)
		{
			for (int j = 0; j < new_mesh.cells_y; ++j)
			{
				if (MaskCellActive2D(mask, new_mesh, i, j))
				{
					continue;
				}
				const int e = GKSFR_CellIndex2D(new_mesh, i, j);
				new_mesh.cell[e] = old_mesh.cell[e];
			}
		}
	}

	void ApplyMaskToAlpha2D(
		const GKSSubcellMask2D* mask,
		const GKSFRMesh2D& mesh,
		std::vector<double>& alpha)
	{
		if (mask == nullptr || !mask->enabled)
		{
			return;
		}
		for (int i = 0; i < mesh.cells_x; ++i)
		{
			for (int j = 0; j < mesh.cells_y; ++j)
			{
				const int e = GKSFR_CellIndex2D(mesh, i, j);
				if (!MaskCellActive2D(mask, mesh, i, j))
				{
					alpha[e] = 1.0;
					continue;
				}
				const bool near_hole =
					(i > 0 && !MaskCellActive2D(mask, mesh, i - 1, j)) ||
					(i + 1 < mesh.cells_x && !MaskCellActive2D(mask, mesh, i + 1, j)) ||
					(j > 0 && !MaskCellActive2D(mask, mesh, i, j - 1)) ||
					(j + 1 < mesh.cells_y && !MaskCellActive2D(mask, mesh, i, j + 1));
				if (near_hole)
				{
					alpha[e] = 1.0;
				}
			}
		}
	}

	bool CheckPhysicalState2D(const GKSFRMesh2D& mesh, int& bad_e, int& bad_i, int& bad_j)
	{
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					const double* q = mesh.cell[e].Q[i][j];
					if (!std::isfinite(q[0]) || !std::isfinite(q[1]) ||
						!std::isfinite(q[2]) || !std::isfinite(q[3]) ||
						q[0] <= 0.0 || Pressure2DLocal(q) <= 0.0)
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

	void PrintBadStateSummary2D(
		const GKSFRMesh2D& mesh,
		int bad_e,
		int bad_i,
		int bad_j)
	{
		if (bad_e < 0 || bad_i < 0 || bad_j < 0)
		{
			return;
		}
		const int cell_i = bad_e / mesh.cells_y;
		const int cell_j = bad_e - cell_i * mesh.cells_y;
		const double x = GKSFR_SolutionPointX2D(mesh, cell_i, bad_i);
		const double y = GKSFR_SolutionPointY2D(mesh, cell_j, bad_j);
		const double* q = mesh.cell[bad_e].Q[bad_i][bad_j];
		const double p = Pressure2DLocal(q);
		std::cout << "bad point x=" << x
			<< " y=" << y
			<< " rho=" << q[0]
			<< " rhoU=" << q[1]
			<< " rhoV=" << q[2]
			<< " rhoE=" << q[3]
			<< " p=" << p << std::endl;
	}

	void EnforceDoubleMachBoundaryBuffer2D(
		const GKSFRMesh2D& mesh,
		GKSFRBoundary2D boundary,
		std::vector<double>& alpha)
	{
		if (boundary != gksfr2d_double_mach ||
			static_cast<int>(alpha.size()) != mesh.cells_x * mesh.cells_y)
		{
			return;
		}

		// The moving top inflow shock is a strong boundary discontinuity.  Use a
		// low-order buffer next to the physical boundary so the high-order FR
		// boundary trace does not inject small spurious waves into the domain.
		const int top_buffer = std::min(2, mesh.cells_y);
		for (int i = 0; i < mesh.cells_x; ++i)
		{
			for (int jj = 0; jj < top_buffer; ++jj)
			{
				const int j = mesh.cells_y - 1 - jj;
				alpha[GKSFR_CellIndex2D(mesh, i, j)] = 1.0;
			}
		}
	}

	void EnforceAstrophysicalJetBoundaryBuffer2D(
		const GKSFRMesh2D& mesh,
		GKSFRBoundary2D boundary,
		std::vector<double>& alpha)
	{
		if (boundary != gksfr2d_astrophysical_jet ||
			static_cast<int>(alpha.size()) != mesh.cells_x * mesh.cells_y)
		{
			return;
		}

		// The jet is injected as a very cold hypersonic inflow.  Keep a low-order
		// inlet core and then relax the blending coefficient gradually downstream;
		// a hard cutoff at the buffer exit tends to generate non-finite states.
		const double jet_core_length = 0.15;
		const double jet_relax_length = 0.35;
		const double jet_buffer_half_width = 0.12;
		for (int i = 0; i < mesh.cells_x; ++i)
		{
			const double x = GKSFR_CellCenterX2D(mesh, i);
			const double xi = x - mesh.x_left;
			if (xi > jet_relax_length)
			{
				continue;
			}
			double alpha_floor = 0.0;
			if (xi <= jet_core_length)
			{
				alpha_floor = 1.0;
			}
			else if (xi <= 0.22)
			{
				alpha_floor = 0.75;
			}
			else if (xi <= 0.28)
			{
				alpha_floor = 0.50;
			}
			else
			{
				alpha_floor = 0.25;
			}
			for (int j = 0; j < mesh.cells_y; ++j)
			{
				const double y = GKSFR_CellCenterY2D(mesh, j);
				if (std::fabs(y) <= jet_buffer_half_width)
				{
					const int e = GKSFR_CellIndex2D(mesh, i, j);
					alpha[e] = std::max(alpha[e], alpha_floor);
				}
			}
		}
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
			Convar_to_primvar_1D(prim, cell.Qc);
			out << CellCenterX(mesh, e) << " "
				<< prim[0] << " "
				<< prim[1] << " "
				<< prim[2] << " "
				<< Temperature(prim[0], prim[2]) << " "
				<< entropy(prim[0], prim[2]) << " "
				<< std::fabs(prim[1]) / Soundspeed(prim[0], prim[2]) << "\n";
		}
	}

	void WriteCellCenterDensityTecplot2D(const GKSFRMesh2D& mesh, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,y,density,u,v,pressure\n";
		out << "zone i = " << mesh.cells_x << ", j = " << mesh.cells_y << ", F=POINT\n";
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int i = 0; i < mesh.cells_x; ++i)
			{
				GKSFRCell2D cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
				GKSFR_ComputeCellCenterData2D(cell, mesh.dx, mesh.dy);
				double prim[4];
				double convar[4];
				for (int m = 0; m < 4; ++m)
				{
					convar[m] = cell.Qc[m];
				}
				Convar_to_primvar_2D(prim, convar);
				out << GKSFR_CellCenterX2D(mesh, i) << " "
					<< GKSFR_CellCenterY2D(mesh, j) << " "
					<< prim[0] << " "
					<< prim[1] << " "
					<< prim[2] << " "
					<< prim[3] << "\n";
			}
		}
	}

	void WriteFRSolutionPointDensityTecplot2D(const GKSFRMesh2D& mesh, const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,y,density,u,v,pressure\n";
		out << "zone i = " << 3 * mesh.cells_x << ", j = " << 3 * mesh.cells_y << ", F=POINT\n";
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int q = 0; q < 3; ++q)
			{
				for (int i = 0; i < mesh.cells_x; ++i)
				{
					const GKSFRCell2D& cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
					for (int p = 0; p < 3; ++p)
					{
						double prim[4];
						double convar[4];
						for (int m = 0; m < 4; ++m)
						{
							convar[m] = cell.Q[p][q][m];
						}
						Convar_to_primvar_2D(prim, convar);
						out << GKSFR_SolutionPointX2D(mesh, i, p) << " "
							<< GKSFR_SolutionPointY2D(mesh, j, q) << " "
							<< prim[0] << " "
							<< prim[1] << " "
							<< prim[2] << " "
							<< prim[3] << "\n";
					}
				}
			}
		}
	}

	void WriteCellCenterDensityTecplot2DMasked(
		const GKSFRMesh2D& mesh,
		const GKSSubcellMask2D& mask,
		const char* path)
	{
		std::ofstream out(path);
		out << std::setprecision(16);
		out << "variables = x,y,density,u,v,pressure,active,Ma\n";
		out << "zone i = " << mesh.cells_x << ", j = " << mesh.cells_y << ", F=POINT\n";
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int i = 0; i < mesh.cells_x; ++i)
			{
				const bool active = MaskCellActive2D(&mask, mesh, i, j);
				if (!active)
				{
					out << GKSFR_CellCenterX2D(mesh, i) << " "
						<< GKSFR_CellCenterY2D(mesh, j) << " "
						<< 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0 << " " << 0.0 << "\n";
					continue;
				}
				GKSFRCell2D cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
				GKSFR_ComputeCellCenterData2D(cell, mesh.dx, mesh.dy);
				double prim[4];
				double convar[4];
				for (int m = 0; m < 4; ++m)
				{
					convar[m] = cell.Qc[m];
				}
				Convar_to_primvar_2D(prim, convar);
				const double velocity_magnitude = std::sqrt(prim[1] * prim[1] + prim[2] * prim[2]);
				const double mach_number = velocity_magnitude / Soundspeed(prim[0], prim[3]);
				out << GKSFR_CellCenterX2D(mesh, i) << " "
					<< GKSFR_CellCenterY2D(mesh, j) << " "
					<< prim[0] << " "
					<< prim[1] << " "
					<< prim[2] << " "
					<< prim[3] << " "
					<< 1 << " "
					<< mach_number << "\n";
			}
		}
	}

	void InitializeDetonationShockDiffraction2D(GKSFRMesh2D& mesh)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int i = 0; i < mesh.cells_x; ++i)
			{
				GKSFRCell2D& cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
				for (int q = 0; q < 3; ++q)
				{
					for (int p = 0; p < 3; ++p)
					{
						const double x = GKSFR_SolutionPointX2D(mesh, i, p);
						const double y = GKSFR_SolutionPointY2D(mesh, j, q);
						double prim[4];
						GKSFR_DetonationShockPrimitive2D(prim, x, y, 0.0);
						Primvar_to_convar_2D(cell.Q[p][q], prim);
					}
				}
			}
		}
	}

	void InitializeAstrophysicalJet2D(GKSFRMesh2D& mesh)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			for (int i = 0; i < mesh.cells_x; ++i)
			{
				GKSFRCell2D& cell = mesh.cell[GKSFR_CellIndex2D(mesh, i, j)];
				for (int q = 0; q < 3; ++q)
				{
					for (int p = 0; p < 3; ++p)
					{
						double prim[4] = { 0.5, 0.0, 0.0, 0.4127 };
						Primvar_to_convar_2D(cell.Q[p][q], prim);
					}
				}
			}
		}
	}

	bool CheckPhysicalState2DMasked(
		const GKSFRMesh2D& mesh,
		const GKSSubcellMask2D& mask,
		int& bad_e,
		int& bad_i,
		int& bad_j)
	{
		for (int i_cell = 0; i_cell < mesh.cells_x; ++i_cell)
		{
			for (int j_cell = 0; j_cell < mesh.cells_y; ++j_cell)
			{
				if (!MaskCellActive2D(&mask, mesh, i_cell, j_cell))
				{
					continue;
				}
				const int e = GKSFR_CellIndex2D(mesh, i_cell, j_cell);
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						const double* q = mesh.cell[e].Q[i][j];
						if (!std::isfinite(q[0]) || !std::isfinite(q[1]) ||
							!std::isfinite(q[2]) || !std::isfinite(q[3]) ||
							q[0] <= 0.0 || Pressure2DLocal(q) <= 0.0)
						{
							bad_e = e;
							bad_i = i;
							bad_j = j;
							return false;
						}
					}
				}
			}
		}
		bad_e = -1;
		bad_i = -1;
		bad_j = -1;
		return true;
	}
	
	void PrintBadStateSummary(
		const GKSFRMesh1D& mesh,
		int bad_e,
		int bad_i)
	{
		double prim[3];
		Convar_to_primvar_1D(prim, const_cast<double*>(mesh.cell[bad_e].Q[bad_i]));
		std::cout << "bad point x=" << SolutionPointX(mesh, bad_e, bad_i)
			<< " rho=" << prim[0]
			<< " u=" << prim[1]
			<< " p=" << prim[2] << std::endl;
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
				einf = std::max(einf, diff);
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
						einf = std::max(einf, diff);
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

	void ResetFluxDiagStateVectors(int faces, GKSFluxLimiterDiag1D& flux_diag)
	{
		flux_diag.left_safe_rho.assign(faces, 0.0);
		flux_diag.left_safe_u.assign(faces, 0.0);
		flux_diag.left_safe_p.assign(faces, 0.0);
		flux_diag.right_safe_rho.assign(faces, 0.0);
		flux_diag.right_safe_u.assign(faces, 0.0);
		flux_diag.right_safe_p.assign(faces, 0.0);
		flux_diag.left_candidate_rho.assign(faces, 0.0);
		flux_diag.left_candidate_u.assign(faces, 0.0);
		flux_diag.left_candidate_p.assign(faces, 0.0);
		flux_diag.right_candidate_rho.assign(faces, 0.0);
		flux_diag.right_candidate_u.assign(faces, 0.0);
		flux_diag.right_candidate_p.assign(faces, 0.0);
		flux_diag.left_candidate_rho_bad.assign(faces, 0);
		flux_diag.left_candidate_p_bad.assign(faces, 0);
		flux_diag.right_candidate_rho_bad.assign(faces, 0);
		flux_diag.right_candidate_p_bad.assign(faces, 0);
		flux_diag.left_safe_rho_bad.assign(faces, 0);
		flux_diag.left_safe_p_bad.assign(faces, 0);
		flux_diag.right_safe_rho_bad.assign(faces, 0);
		flux_diag.right_safe_p_bad.assign(faces, 0);
	}

	void MixedFaceFluxWithoutLimiter(
		const std::vector<double>& alpha_final,
		const std::vector<GKSFRFaceFlux1D>& high_face_fluxes,
		const std::vector<GKSFRFaceFlux1D>& low_face_fluxes,
		GKSFRBoundary1D boundary,
		std::vector<GKSFRFaceFlux1D>& final_face_fluxes,
		GKSFluxLimiterDiag1D& flux_diag)
	{
		const int faces = static_cast<int>(high_face_fluxes.size());
		const int cells = faces - 1;
		final_face_fluxes.assign(faces, GKSFRFaceFlux1D());
		flux_diag.alpha_face.assign(faces, 0.0);
		flux_diag.theta_rho.assign(faces, 1.0);
		flux_diag.theta_p.assign(faces, 1.0);
		flux_diag.theta.assign(faces, 1.0);
		ResetFluxDiagStateVectors(faces, flux_diag);

		for (int f = 0; f < faces; ++f)
		{
			double alpha_face = 0.0;
			if (f == 0)
			{
				alpha_face = (boundary == gksfr_periodic)
					? 0.5 * (alpha_final[cells - 1] + alpha_final[0])
					: alpha_final[0];
			}
			else if (f == cells)
			{
				alpha_face = (boundary == gksfr_periodic)
					? flux_diag.alpha_face[0]
					: alpha_final[cells - 1];
			}
			else
			{
				alpha_face = 0.5 * (alpha_final[f - 1] + alpha_final[f]);
			}
			flux_diag.alpha_face[f] = alpha_face;
			if (alpha_face <= 1.0e-14)
			{
				for (int m = 0; m < 3; ++m)
				{
					final_face_fluxes[f].F[m] = high_face_fluxes[f].F[m];
				}
			}
			else if (alpha_face >= 1.0 - 1.0e-14)
			{
				for (int m = 0; m < 3; ++m)
				{
					final_face_fluxes[f].F[m] = low_face_fluxes[f].F[m];
				}
			}
			else
			{
				for (int m = 0; m < 3; ++m)
				{
					final_face_fluxes[f].F[m] =
						(1.0 - alpha_face) * high_face_fluxes[f].F[m]
						+ alpha_face * low_face_fluxes[f].F[m];
				}
			}
		}
		if (boundary == gksfr_periodic)
		{
			for (int m = 0; m < 3; ++m)
			{
				final_face_fluxes[cells].F[m] = final_face_fluxes[0].F[m];
			}
		}
	}

	int XFaceIndex2DLocal(int cells_y, int face_i, int cell_j)
	{
		return face_i * cells_y + cell_j;
	}

	int YFaceIndex2DLocal(int cells_y, int cell_i, int face_j)
	{
		return cell_i * (cells_y + 1) + face_j;
	}

	void MixedFaceFluxWithoutLimiter2D(
		const GKSFRMesh2D& mesh,
		const std::vector<double>& alpha_final,
		const std::vector<GKSFRFaceFlux2D>& high_x_face_fluxes,
		const std::vector<GKSFRFaceFlux2D>& high_y_face_fluxes,
		const std::vector<GKSFRFaceFlux2D>& low_x_face_fluxes,
		const std::vector<GKSFRFaceFlux2D>& low_y_face_fluxes,
		std::vector<GKSFRFaceFlux2D>& final_x_face_fluxes,
		std::vector<GKSFRFaceFlux2D>& final_y_face_fluxes,
		GKSFluxLimiterDiag2D& flux_diag)
	{
		final_x_face_fluxes.assign((mesh.cells_x + 1) * mesh.cells_y, GKSFRFaceFlux2D());
		final_y_face_fluxes.assign(mesh.cells_x * (mesh.cells_y + 1), GKSFRFaceFlux2D());
		flux_diag = GKSFluxLimiterDiag2D();
		flux_diag.alpha_face_x.assign(final_x_face_fluxes.size() * 3, 0.0);
		flux_diag.alpha_face_y.assign(final_y_face_fluxes.size() * 3, 0.0);
		flux_diag.theta_x.assign(final_x_face_fluxes.size() * 3, 1.0);
		flux_diag.theta_y.assign(final_y_face_fluxes.size() * 3, 1.0);
		flux_diag.theta_rho_x.assign(final_x_face_fluxes.size() * 3, 1.0);
		flux_diag.theta_rho_y.assign(final_y_face_fluxes.size() * 3, 1.0);
		flux_diag.theta_p_x.assign(final_x_face_fluxes.size() * 3, 1.0);
		flux_diag.theta_p_y.assign(final_y_face_fluxes.size() * 3, 1.0);

		for (int face_i = 0; face_i <= mesh.cells_x; ++face_i)
		{
			const int left_i = (face_i + mesh.cells_x - 1) % mesh.cells_x;
			const int right_i = face_i % mesh.cells_x;
			for (int j = 0; j < mesh.cells_y; ++j)
			{
				const int left_e = GKSFR_CellIndex2D(mesh, left_i, j);
				const int right_e = GKSFR_CellIndex2D(mesh, right_i, j);
				const double alpha_face = 0.5 * (alpha_final[left_e] + alpha_final[right_e]);
				const int f = XFaceIndex2DLocal(mesh.cells_y, face_i, j);
				for (int q = 0; q < 3; ++q)
				{
					flux_diag.alpha_face_x[f * 3 + q] = alpha_face;
					for (int m = 0; m < 4; ++m)
					{
						final_x_face_fluxes[f].F[q][m] =
							(1.0 - alpha_face) * high_x_face_fluxes[f].F[q][m]
							+ alpha_face * low_x_face_fluxes[f].F[q][m];
					}
				}
			}
		}
		for (int i = 0; i < mesh.cells_x; ++i)
		{
			for (int face_j = 0; face_j <= mesh.cells_y; ++face_j)
			{
				const int bottom_j = (face_j + mesh.cells_y - 1) % mesh.cells_y;
				const int top_j = face_j % mesh.cells_y;
				const int bottom_e = GKSFR_CellIndex2D(mesh, i, bottom_j);
				const int top_e = GKSFR_CellIndex2D(mesh, i, top_j);
				const double alpha_face = 0.5 * (alpha_final[bottom_e] + alpha_final[top_e]);
				const int f = YFaceIndex2DLocal(mesh.cells_y, i, face_j);
				for (int p = 0; p < 3; ++p)
				{
					flux_diag.alpha_face_y[f * 3 + p] = alpha_face;
					for (int m = 0; m < 4; ++m)
					{
						final_y_face_fluxes[f].F[p][m] =
							(1.0 - alpha_face) * high_y_face_fluxes[f].F[p][m]
							+ alpha_face * low_y_face_fluxes[f].F[p][m];
					}
				}
			}
		}
	}

	void BuildSafeCellAverage(
		const GKSSubcellBranch1D& low_branch,
		std::vector<GKSCellAverage1D>& safe_avg)
	{
		safe_avg.resize(low_branch.cell.size());
		for (int e = 0; e < static_cast<int>(low_branch.cell.size()); ++e)
		{
			for (int m = 0; m < 3; ++m)
			{
				safe_avg[e].U[m] = low_branch.cell[e].cell_avg_new[m];
			}
		}
	}

	bool AdvanceCase(
		GKSFRMesh1D& mesh,
		double CFL,
		double tstop,
		GKSFRBoundary1D boundary,
		const GKSSubcellFrameworkConfig1D& config,
		GKSSubcellFrameworkDiag1D& final_diag,
		bool show_step,
		int* final_step = nullptr)
	{
		int step = 0;
		double t = 0.0;
		while (t < tstop - 1.0e-14)
		{
			const double dt = GetTimeStep(mesh, CFL, t, tstop);
			if (show_step)
			{
				std::cout << "step= " << step
					<< " time size is " << dt
					<< " time= " << t << std::endl;
			}
			GKSSubcellFrameworkDiag1D diag;
			GKSSubcellAdvanceOneStep1D(mesh, dt, boundary, config, diag);
			final_diag = diag;
			t += dt;
			++step;

			int bad_e = -1;
			int bad_i = -1;
			if (!CheckPhysicalState(mesh, bad_e, bad_i))
			{
				std::cout << "GKS-subcell invalid state at cell=" << bad_e
					<< " point=" << bad_i
					<< " after step=" << step
					<< " time=" << t << std::endl;
				PrintBadStateSummary(mesh, bad_e, bad_i);
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

}

void GKSSubcellAdvanceOneStep1D(
	GKSFRMesh1D& mesh,
	double dt,
	GKSFRBoundary1D boundary,
	const GKSSubcellFrameworkConfig1D& config,
	GKSSubcellFrameworkDiag1D& diag)
{
	low_order_type = config.low_mode;

	GKSSmoothIndicatorAllCells1D(mesh, config.smooth_param, diag.alpha_raw, diag.alpha_final);
	if (config.blend_mode == gks_subcell_pure_high)
	{
		for (int e = 0; e < mesh.cells; ++e)
		{
			diag.alpha_final[e] = 0.0;
		}
	}
	else if (config.blend_mode == gks_subcell_pure_low)
	{
		for (int e = 0; e < mesh.cells; ++e)
		{
			diag.alpha_final[e] = 1.0;
		}
	}

	std::vector<GKSFRFaceFlux1D> high_face_fluxes;
	GKSFRAdapterComputeHighFaceFluxes(mesh, dt, boundary, high_face_fluxes);

	GKSSubcellBranch1D low_branch;
	GKSSubcellInitializeFromCurrentDofs1D(mesh, low_branch);
	GKSSubcellComputeInternalFluxes1D(low_branch, dt);

	GKSSubcellBranch1D low_branch_kfvs_safe = low_branch;
	if (config.low_mode == MUSCL_HANCOCK)
	{
		GKSSubcellComputeInternalFluxes1D(low_branch_kfvs_safe, dt, KFVS1);
	}

	std::vector<GKSFRFaceFlux1D> low_face_fluxes;
	GKSSubcellComputeLowFaceFluxes1D(low_branch_kfvs_safe, dt, boundary, low_face_fluxes);

	std::vector<GKSFRFaceFlux1D> final_face_fluxes;
	if (config.blend_mode == gks_subcell_pure_high)
	{
		final_face_fluxes = high_face_fluxes;
		diag.flux_diag.alpha_face.assign(high_face_fluxes.size(), 0.0);
		diag.flux_diag.theta_rho.assign(high_face_fluxes.size(), 1.0);
		diag.flux_diag.theta_p.assign(high_face_fluxes.size(), 1.0);
		diag.flux_diag.theta.assign(high_face_fluxes.size(), 1.0);
		ResetFluxDiagStateVectors(static_cast<int>(high_face_fluxes.size()), diag.flux_diag);
	}
	else if (config.blend_mode == gks_subcell_pure_low)
	{
		final_face_fluxes = low_face_fluxes;
		diag.flux_diag.alpha_face.assign(low_face_fluxes.size(), 1.0);
		diag.flux_diag.theta_rho.assign(low_face_fluxes.size(), 1.0);
		diag.flux_diag.theta_p.assign(low_face_fluxes.size(), 1.0);
		diag.flux_diag.theta.assign(low_face_fluxes.size(), 1.0);
		ResetFluxDiagStateVectors(static_cast<int>(low_face_fluxes.size()), diag.flux_diag);
	}
	else if (config.use_flux_limiter)
	{
		GKSFluxLimiterApply1D(
			low_branch_kfvs_safe,
			diag.alpha_final,
			high_face_fluxes,
			low_face_fluxes,
			dt,
			boundary,
			config.flux_param,
			final_face_fluxes,
			diag.flux_diag);
	}
	else
	{
		MixedFaceFluxWithoutLimiter(
			diag.alpha_final,
			high_face_fluxes,
			low_face_fluxes,
			boundary,
			final_face_fluxes,
			diag.flux_diag);
	}

	GKSFRMesh1D high_residual;
	GKSFRAdapterAssembleBlendedHighResiduals(
		mesh,
		dt,
		final_face_fluxes,
		diag.alpha_final,
		high_residual);

	GKSSubcellBranch1D low_branch_safe = low_branch;
	if (config.low_mode == MUSCL_HANCOCK)
	{
		low_branch_safe = low_branch_kfvs_safe;
	}
	GKSSubcellComputeElementResiduals1D(low_branch_safe);
	GKSSubcellAddFaceResiduals1D(low_branch_safe, final_face_fluxes);
	GKSSubcellUpdateLowOrderDofs1D(low_branch_safe, dt);

	GKSSubcellComputeElementResiduals1D(low_branch);
	GKSSubcellAddFaceResiduals1D(low_branch, final_face_fluxes);
	GKSSubcellFallbackBadMUSCLUpdates1D(low_branch, low_branch_safe, dt);
	GKSSubcellScaleResiduals1D(low_branch, diag.alpha_final);
	diag.muscl_stats = low_branch.muscl_stats;

	for (int e = 0; e < mesh.cells; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int m = 0; m < 3; ++m)
			{
				mesh.cell[e].Q[i][m] += dt * (
					high_residual.cell[e].residual[i][m]
					+ low_branch.cell[e].residual_subcell[i][m]);
			}
		}
	}

	std::vector<GKSCellAverage1D> safe_avg;
	BuildSafeCellAverage(low_branch_safe, safe_avg);
	if (config.use_scaling_limiter)
	{
		GKSScalingLimiterApply1D(mesh, safe_avg, config.scaling_param, diag.scaling_diag);
	}
	else
	{
		diag.scaling_diag.theta_rho.assign(mesh.cells, 1.0);
		diag.scaling_diag.theta_p.assign(mesh.cells, 1.0);
		diag.scaling_diag.theta.assign(mesh.cells, 1.0);
	}

	double prim[3];
	diag.min_rho = 1.0e300;
	diag.min_p = 1.0e300;
	diag.min_rho_cell = -1;
	diag.min_rho_point = -1;
	diag.min_p_cell = -1;
	diag.min_p_point = -1;
	for (int e = 0; e < mesh.cells; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			Convar_to_primvar_1D(prim, mesh.cell[e].Q[i]);
			if (prim[0] < diag.min_rho)
			{
				diag.min_rho = prim[0];
				diag.min_rho_cell = e;
				diag.min_rho_point = i;
			}
			if (prim[2] < diag.min_p)
			{
				diag.min_p = prim[2];
				diag.min_p_cell = e;
				diag.min_p_point = i;
			}
		}
	}
}

void GKSSubcellAdvanceOneStep2D(
	GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	const GKSSubcellFrameworkConfig2D& config,
	GKSSubcellFrameworkDiag2D& diag)
{
	low_order_type = config.low_mode;
	JetTraceMeshBadStates2D("step_start_solution", mesh);

	GKSSmoothIndicatorAllCells2D(mesh, config.smooth_param, diag.alpha_raw, diag.alpha_final);
	if (config.blend_mode == gks_subcell2d_pure_high)
	{
		std::fill(diag.alpha_final.begin(), diag.alpha_final.end(), 0.0);
	}
	else if (config.blend_mode == gks_subcell2d_pure_low)
	{
		std::fill(diag.alpha_final.begin(), diag.alpha_final.end(), 1.0);
	}
	EnforceDoubleMachBoundaryBuffer2D(mesh, boundary, diag.alpha_final);
	EnforceAstrophysicalJetBoundaryBuffer2D(mesh, boundary, diag.alpha_final);
	JetTraceMeshBadStates2D("after_alpha_buffer", mesh, &diag.alpha_final);

	std::vector<GKSFRFaceFlux2D> high_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> high_y_face_fluxes;
	GKSFRAdapterComputeHighFaceFluxes2D(mesh, dt, boundary, high_x_face_fluxes, high_y_face_fluxes);

	GKSSubcellBranch2D low_branch;
	GKSSubcellInitializeFromCurrentDofs2D(mesh, low_branch);
	JetTraceBranchBadStates2D("low_flux_input", mesh, low_branch, false);
	GKSSubcellComputeInternalFluxes2D(low_branch, dt, config.low_mode);

	GKSSubcellBranch2D low_branch_kfvs_safe = low_branch;
	if (config.low_mode == MUSCL_HANCOCK_2d)
	{
		GKSSubcellComputeInternalFluxes2D(low_branch_kfvs_safe, dt, KFVS1);
	}

	std::vector<GKSFRFaceFlux2D> low_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> low_y_face_fluxes;
	GKSSubcellComputeLowFaceFluxes2D(low_branch_kfvs_safe, dt, boundary, low_x_face_fluxes, low_y_face_fluxes);
	JetTraceFaceFluxes2D("low_flux_output", mesh, low_x_face_fluxes, low_y_face_fluxes);

	std::vector<GKSFRFaceFlux2D> final_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> final_y_face_fluxes;
	if (config.blend_mode == gks_subcell2d_pure_high)
	{
		final_x_face_fluxes = high_x_face_fluxes;
		final_y_face_fluxes = high_y_face_fluxes;
		diag.flux_diag = GKSFluxLimiterDiag2D();
	}
	else if (config.blend_mode == gks_subcell2d_pure_low)
	{
		final_x_face_fluxes = low_x_face_fluxes;
		final_y_face_fluxes = low_y_face_fluxes;
		diag.flux_diag = GKSFluxLimiterDiag2D();
	}
	else if (config.use_flux_limiter)
	{
		GKSFluxLimiterApply2D(
			low_branch_kfvs_safe,
			diag.alpha_final,
			high_x_face_fluxes,
			high_y_face_fluxes,
			low_x_face_fluxes,
			low_y_face_fluxes,
			dt,
			boundary,
			config.flux_param,
			final_x_face_fluxes,
			final_y_face_fluxes,
			diag.flux_diag);
	}
	else
	{
		MixedFaceFluxWithoutLimiter2D(
			mesh,
			diag.alpha_final,
			high_x_face_fluxes,
			high_y_face_fluxes,
			low_x_face_fluxes,
			low_y_face_fluxes,
			final_x_face_fluxes,
			final_y_face_fluxes,
			diag.flux_diag);
	}
	JetTraceFaceFluxes2D("final_flux_assembly", mesh, final_x_face_fluxes, final_y_face_fluxes, &diag.flux_diag);

	GKSFRMesh2D high_new;
	GKSFRAdapterAdvanceWithFaceFluxes2D(mesh, dt, final_x_face_fluxes, final_y_face_fluxes, high_new);

	GKSSubcellBranch2D low_safe = low_branch_kfvs_safe;
	GKSSubcellAdvanceWithFaceFluxes2D(low_safe, final_x_face_fluxes, final_y_face_fluxes, dt);
	JetTraceBranchBadStates2D("low_update", mesh, low_safe, true);

	GKSSubcellBranch2D low_candidate = low_branch;
	if (config.low_mode == MUSCL_HANCOCK_2d)
	{
		GKSSubcellComputeElementResiduals2D(low_candidate);
		GKSSubcellAddFaceResiduals2D(low_candidate, final_x_face_fluxes, final_y_face_fluxes);
		GKSSubcellFallbackBadMUSCLUpdates2D(low_candidate, low_safe, dt);
		GKSSubcellUpdateLowOrderDofs2D(low_candidate, dt);
	}
	else
	{
		low_candidate = low_safe;
	}
	diag.muscl_stats = low_candidate.muscl_stats;

	GKSFRMesh2D mixed = mesh;
	diag.max_alpha = 0.0;
	diag.troubled_cells = 0;
	for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
	{
		const double alpha = diag.alpha_final[e];
		diag.max_alpha = std::max(diag.max_alpha, alpha);
		if (alpha > 0.0)
		{
			diag.troubled_cells++;
		}
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					if (alpha <= 1.0e-14)
					{
						mixed.cell[e].Q[i][j][m] = high_new.cell[e].Q[i][j][m];
					}
					else if (alpha >= 1.0 - 1.0e-14)
					{
						mixed.cell[e].Q[i][j][m] = low_candidate.cell[e].low_dof_new[i][j][m];
					}
					else if (!std::isfinite(high_new.cell[e].Q[i][j][m]) &&
						std::isfinite(low_candidate.cell[e].low_dof_new[i][j][m]))
					{
						mixed.cell[e].Q[i][j][m] = low_candidate.cell[e].low_dof_new[i][j][m];
					}
					else
					{
						mixed.cell[e].Q[i][j][m] =
							(1.0 - alpha) * high_new.cell[e].Q[i][j][m]
							+ alpha * low_candidate.cell[e].low_dof_new[i][j][m];
					}
				}
			}
		}
	}
	JetTraceMeshBadStates2D("scaling_limiter_before", mixed, &diag.alpha_final);

	if (config.use_scaling_limiter)
	{
		std::vector<GKSCellAverage2D> safe_average(mesh.cells_x * mesh.cells_y);
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			for (int m = 0; m < 4; ++m)
			{
				safe_average[e].U[m] = low_safe.cell[e].cell_avg_new[m];
			}
		}
		GKSScalingLimiterApply2D(mixed, safe_average, config.scaling_param, diag.scaling_diag);
		if (JetTraceEnabled())
		{
			const int target = JetTraceCell();
			if (target >= 0 && target < mesh.cells_x * mesh.cells_y &&
				target < static_cast<int>(diag.scaling_diag.theta.size()))
			{
				std::cout << "JET_TRACE scaling theta cell=" << target
					<< " theta_rho=" << diag.scaling_diag.theta_rho[target]
					<< " theta_p=" << diag.scaling_diag.theta_p[target]
					<< " theta=" << diag.scaling_diag.theta[target] << std::endl;
			}
		}
	}
	else
	{
		diag.scaling_diag.theta_rho.assign(mesh.cells_x * mesh.cells_y, 1.0);
		diag.scaling_diag.theta_p.assign(mesh.cells_x * mesh.cells_y, 1.0);
		diag.scaling_diag.theta.assign(mesh.cells_x * mesh.cells_y, 1.0);
		diag.scaling_diag.limited_cells = 0;
	}
	JetTraceMeshBadStates2D("scaling_limiter_after", mixed, &diag.alpha_final);

	diag.min_rho = 1.0e300;
	diag.min_p = 1.0e300;
	for (int e = 0; e < mixed.cells_x * mixed.cells_y; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				const double* q = mixed.cell[e].Q[i][j];
				diag.min_rho = std::min(diag.min_rho, q[0]);
				diag.min_p = std::min(diag.min_p, Pressure2DLocal(q));
			}
		}
	}
	mesh = mixed;
	JetTraceMeshBadStates2D("final_solution_check", mesh, &diag.alpha_final);
}

void GKSSubcellBuildRectangularCutoutMask2D(
	const GKSFRMesh2D& mesh,
	double x_min,
	double x_max,
	double y_min,
	double y_max,
	GKSSubcellMask2D& mask)
{
	mask.enabled = true;
	mask.cells_x = mesh.cells_x;
	mask.cells_y = mesh.cells_y;
	mask.active.assign(mesh.cells_x * mesh.cells_y, 1);
	for (int i = 0; i < mesh.cells_x; ++i)
	{
		for (int j = 0; j < mesh.cells_y; ++j)
		{
			const double x = GKSFR_CellCenterX2D(mesh, i);
			const double y = GKSFR_CellCenterY2D(mesh, j);
			if (x >= x_min && x <= x_max && y >= y_min && y <= y_max)
			{
				mask.active[GKSFR_CellIndex2D(mesh, i, j)] = 0;
			}
		}
	}
}

void GKSSubcellAdvanceOneStep2DMasked(
	GKSFRMesh2D& mesh,
	double dt,
	GKSFRBoundary2D boundary,
	const GKSSubcellMask2D& mask,
	const GKSSubcellFrameworkConfig2D& config,
	GKSSubcellFrameworkDiag2D& diag)
{
	if (!mask.enabled)
	{
		GKSSubcellAdvanceOneStep2D(mesh, dt, boundary, config, diag);
		return;
	}

	low_order_type = config.low_mode;
	GKSSmoothIndicatorAllCells2D(mesh, config.smooth_param, diag.alpha_raw, diag.alpha_final);
	if (config.blend_mode == gks_subcell2d_pure_high)
	{
		std::fill(diag.alpha_final.begin(), diag.alpha_final.end(), 0.0);
	}
	else if (config.blend_mode == gks_subcell2d_pure_low)
	{
		std::fill(diag.alpha_final.begin(), diag.alpha_final.end(), 1.0);
	}
	ApplyMaskToAlpha2D(&mask, mesh, diag.alpha_final);

	std::vector<GKSFRFaceFlux2D> high_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> high_y_face_fluxes;
	GKSFRAdapterComputeHighFaceFluxes2D(mesh, dt, boundary, high_x_face_fluxes, high_y_face_fluxes);

	GKSSubcellBranch2D low_branch;
	GKSSubcellInitializeFromCurrentDofs2D(mesh, low_branch);
	GKSSubcellComputeInternalFluxes2D(low_branch, dt, config.low_mode);

	GKSSubcellBranch2D low_branch_kfvs_safe = low_branch;
	if (config.low_mode == MUSCL_HANCOCK_2d)
	{
		GKSSubcellComputeInternalFluxes2D(low_branch_kfvs_safe, dt, KFVS1);
	}
	OverrideMaskedWallFaceFluxes2D(&mask, low_branch_kfvs_safe, dt, high_x_face_fluxes, high_y_face_fluxes);

	std::vector<GKSFRFaceFlux2D> low_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> low_y_face_fluxes;
	GKSSubcellComputeLowFaceFluxes2D(low_branch_kfvs_safe, dt, boundary, low_x_face_fluxes, low_y_face_fluxes);
	OverrideMaskedWallFaceFluxes2D(&mask, low_branch_kfvs_safe, dt, low_x_face_fluxes, low_y_face_fluxes);

	std::vector<GKSFRFaceFlux2D> final_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> final_y_face_fluxes;
	if (config.blend_mode == gks_subcell2d_pure_high)
	{
		final_x_face_fluxes = high_x_face_fluxes;
		final_y_face_fluxes = high_y_face_fluxes;
		diag.flux_diag = GKSFluxLimiterDiag2D();
	}
	else if (config.blend_mode == gks_subcell2d_pure_low)
	{
		final_x_face_fluxes = low_x_face_fluxes;
		final_y_face_fluxes = low_y_face_fluxes;
		diag.flux_diag = GKSFluxLimiterDiag2D();
	}
	else if (config.use_flux_limiter)
	{
		GKSFluxLimiterApply2D(
			low_branch_kfvs_safe,
			diag.alpha_final,
			high_x_face_fluxes,
			high_y_face_fluxes,
			low_x_face_fluxes,
			low_y_face_fluxes,
			dt,
			boundary,
			config.flux_param,
			final_x_face_fluxes,
			final_y_face_fluxes,
			diag.flux_diag);
	}
	else
	{
		MixedFaceFluxWithoutLimiter2D(
			mesh,
			diag.alpha_final,
			high_x_face_fluxes,
			high_y_face_fluxes,
			low_x_face_fluxes,
			low_y_face_fluxes,
			final_x_face_fluxes,
			final_y_face_fluxes,
			diag.flux_diag);
	}
	OverrideMaskedWallFaceFluxes2D(&mask, low_branch_kfvs_safe, dt, final_x_face_fluxes, final_y_face_fluxes);

	GKSFRMesh2D high_new;
	GKSFRAdapterAdvanceWithFaceFluxes2D(mesh, dt, final_x_face_fluxes, final_y_face_fluxes, high_new);
	RestoreInactiveCells2D(&mask, mesh, high_new);

	GKSSubcellBranch2D low_safe = low_branch_kfvs_safe;
	GKSSubcellAdvanceWithFaceFluxes2D(low_safe, final_x_face_fluxes, final_y_face_fluxes, dt);

	GKSSubcellBranch2D low_candidate = low_branch;
	if (config.low_mode == MUSCL_HANCOCK_2d)
	{
		GKSSubcellComputeElementResiduals2D(low_candidate);
		GKSSubcellAddFaceResiduals2D(low_candidate, final_x_face_fluxes, final_y_face_fluxes);
		GKSSubcellFallbackBadMUSCLUpdates2D(low_candidate, low_safe, dt);
		GKSSubcellUpdateLowOrderDofs2D(low_candidate, dt);
	}
	else
	{
		low_candidate = low_safe;
	}
	diag.muscl_stats = low_candidate.muscl_stats;

	GKSFRMesh2D mixed = mesh;
	diag.max_alpha = 0.0;
	diag.troubled_cells = 0;
	for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
	{
		const int i_cell = e / mesh.cells_y;
		const int j_cell = e - i_cell * mesh.cells_y;
		if (!MaskCellActive2D(&mask, mesh, i_cell, j_cell))
		{
			mixed.cell[e] = mesh.cell[e];
			continue;
		}
		const double alpha = diag.alpha_final[e];
		diag.max_alpha = std::max(diag.max_alpha, alpha);
		if (alpha > 0.0)
		{
			diag.troubled_cells++;
		}
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int m = 0; m < 4; ++m)
				{
					if (alpha <= 1.0e-14)
					{
						mixed.cell[e].Q[i][j][m] = high_new.cell[e].Q[i][j][m];
					}
					else if (alpha >= 1.0 - 1.0e-14)
					{
						mixed.cell[e].Q[i][j][m] = low_candidate.cell[e].low_dof_new[i][j][m];
					}
					else if (!std::isfinite(high_new.cell[e].Q[i][j][m]) &&
						std::isfinite(low_candidate.cell[e].low_dof_new[i][j][m]))
					{
						mixed.cell[e].Q[i][j][m] = low_candidate.cell[e].low_dof_new[i][j][m];
					}
					else
					{
						mixed.cell[e].Q[i][j][m] =
							(1.0 - alpha) * high_new.cell[e].Q[i][j][m]
							+ alpha * low_candidate.cell[e].low_dof_new[i][j][m];
					}
				}
			}
		}
	}
	if (config.use_scaling_limiter)
	{
		std::vector<GKSCellAverage2D> safe_average(mesh.cells_x * mesh.cells_y);
		for (int e = 0; e < mesh.cells_x * mesh.cells_y; ++e)
		{
			for (int m = 0; m < 4; ++m)
			{
				safe_average[e].U[m] = low_safe.cell[e].cell_avg_new[m];
			}
		}
		GKSScalingLimiterApply2D(mixed, safe_average, config.scaling_param, diag.scaling_diag);
	}
	else
	{
		diag.scaling_diag.theta_rho.assign(mesh.cells_x * mesh.cells_y, 1.0);
		diag.scaling_diag.theta_p.assign(mesh.cells_x * mesh.cells_y, 1.0);
		diag.scaling_diag.theta.assign(mesh.cells_x * mesh.cells_y, 1.0);
		diag.scaling_diag.limited_cells = 0;
	}
	RestoreInactiveCells2D(&mask, mesh, mixed);

	diag.min_rho = 1.0e300;
	diag.min_p = 1.0e300;
	for (int e = 0; e < mixed.cells_x * mixed.cells_y; ++e)
	{
		const int i_cell = e / mixed.cells_y;
		const int j_cell = e - i_cell * mixed.cells_y;
		if (!MaskCellActive2D(&mask, mixed, i_cell, j_cell))
		{
			continue;
		}
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				const double* q = mixed.cell[e].Q[i][j];
				diag.min_rho = std::min(diag.min_rho, q[0]);
				diag.min_p = std::min(diag.min_p, Pressure2DLocal(q));
			}
		}
	}
	mesh = mixed;
}

void accuracy_sinwave_1d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_1D(0.0, 0.0);

	GKSSubcellFrameworkConfig1D config;
	config.blend_mode = gks_subcell_hybrid;
	config.low_mode = KFVS1;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;

	const double CFL = 0.02;
	const double tstop = 2.0;
	const int mesh_set = 5;
	const int mesh_number[mesh_set] = { 20, 40, 80, 160, 320 };
	double error[mesh_set][3]{};
	GKSSubcellFrameworkDiag1D last_diag;

	for (int imesh = 0; imesh < mesh_set; ++imesh)
	{
		GKSFRMesh1D mesh;
		GKSFR_ResizeUniformMesh(mesh, mesh_number[imesh], 0.0, 2.0);
		InitializeSinwave(mesh);
		if (!AdvanceCase(mesh, CFL, tstop, gksfr_periodic, config, last_diag, false))
		{
			std::cout << "GKS-subcell sinwave failed on mesh " << mesh_number[imesh] << std::endl;
			return;
		}
		ComputeSinwaveError(mesh, tstop, error[imesh]);
		const std::string out_path = "build/result/gks_subcell_sinwave_mesh" + std::to_string(mesh_number[imesh]) + ".plt";
		WriteCellCenterDensityTecplot(mesh, out_path.c_str());
	}

	std::cout << "GKS-subcell 1D sinwave errors" << std::endl;
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

void riemann_problem_1d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_1D(0.1, 1.0);

	GKSSubcellFrameworkConfig1D config;
	config.blend_mode = gks_subcell_hybrid;
	config.low_mode = MUSCL_HANCOCK;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;

	GKSFRMesh1D mesh;
	GKSFR_ResizeUniformMesh(mesh, 2000, -5.0, 5.0);
	//ICfor1dRM(mesh, RiemannProblem1D_Sod());
	//ICfor1dRM(mesh,RiemannProblem1D_DoubleRarefaction());
	//ICfor1dRM(mesh,RiemannProblem1D_Leblanc());
	ICfor1dRM(mesh,RiemannProblem1D_ShuOsher());//shu-osher问题
	//ICfor1dRM(mesh,RiemannProblem1D_BlastWave());
	//ICfor1dRM(mesh, RiemannProblem1D_SedovBlastWave(-2.0, 4.0 / 201.0));

	GKSSubcellFrameworkDiag1D diag;
	int final_step = 0;
	const bool ok = AdvanceCase(mesh, 0.02, 1.8, gksfr_shu_osher, config, diag, true, &final_step);
	if (!ok)
	{
		std::cout << "GKS-subcell run stopped because of an invalid state." << std::endl;
		return;
	}

	WriteSodTecplot(mesh, StepTaggedPath("gks_subcell", final_step, ".plt").c_str());
}

void accuracy_sinwave_2d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_2D(0.0, 0.0);
	GKSSubcellFrameworkConfig2D config;
	config.blend_mode = gks_subcell2d_hybrid;
	config.low_mode = KFVS1;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;

	const double CFL = 0.05;
	const double tstop = 0.1;
	const int mesh_set = 3;
	const int mesh_number[mesh_set] = { 20, 40, 80 };
	double error[mesh_set][3]{};
	GKSSubcellFrameworkDiag2D last_diag;

	for (int imesh = 0; imesh < mesh_set; ++imesh)
	{
		GKSFRMesh2D mesh;
		GKSFR_ResizeUniformMesh2D(mesh, mesh_number[imesh], mesh_number[imesh], 0.0, 2.0, 0.0, 2.0);
		InitializeSinwave2D(mesh);
		double t = 0.0;
		int step = 0;
		while (t < tstop - 1.0e-14)
		{
			const double dt = GetTimeStep2D(mesh, CFL, t, tstop);
			GKSSubcellAdvanceOneStep2D(mesh, dt, gksfr2d_periodic, config, last_diag);
			t += dt;
			++step;
			int bad_e = -1, bad_i = -1, bad_j = -1;
			if (!CheckPhysicalState2D(mesh, bad_e, bad_i, bad_j))
			{
				std::cout << "GKS-subcell 2D sinwave failed at mesh=" << mesh_number[imesh]
					<< " step=" << step
					<< " cell=" << bad_e
					<< " point=(" << bad_i << "," << bad_j << ")" << std::endl;
				return;
			}
		}
		ComputeSinwaveError2D(mesh, tstop, error[imesh]);
		WriteCellCenterDensityTecplot2D(
			mesh,
			("build/result/gks_subcell_2d_sinwave_mesh" + std::to_string(mesh_number[imesh]) + ".plt").c_str());
		std::cout << "smooth limiter stats mesh=" << mesh_number[imesh]
			<< " max alpha=" << last_diag.max_alpha
			<< " troubled cells=" << last_diag.troubled_cells
			<< " flux-limited x/y=" << last_diag.flux_diag.limited_faces_x << "/" << last_diag.flux_diag.limited_faces_y
			<< " scaling cells=" << last_diag.scaling_diag.limited_cells << std::endl;
	}

	std::cout << "GKS-subcell 2D sinwave errors" << std::endl;
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

void riemann_problem_2d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_2D(0.1, 1.0);
	GKSSubcellFrameworkConfig2D config;
	config.blend_mode = gks_subcell2d_hybrid;
	config.low_mode = KFVS1;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;
	config.flux_param.kx = 0.5;
	config.flux_param.ky = 0.5;
	config.flux_param.rho_floor = 1.0e-10;
	config.flux_param.p_floor = 1.0e-8;
	config.scaling_param.rho_floor = 1.0e-10;
	config.scaling_param.p_floor = 1.0e-8;

	GKSFRMesh2D mesh;
	GKSFR_ResizeUniformMesh2D(mesh, 40, 40, 0.0, 1.0, 0.0, 1.0);
	ICfor2dRM(mesh, RiemannProblem2D_SubcellLimiterReference());

	const double CFL = 0.02;
	const double tstop = 0.25;
	double t = 0.0;
	int step = 0;
	GKSSubcellFrameworkDiag2D diag;
	while (t < tstop - 1.0e-14)
	{
		const double dt = GetTimeStep2D(mesh, CFL, t, tstop);
		std::cout << "step=" << step << " dt=" << dt << " t=" << t << std::endl;
		GKSFR_SetBoundaryTime2D(t + 0.5 * dt);
		GKSSubcellAdvanceOneStep2D(mesh, dt, gksfr2d_transmissive, config, diag);
		t += dt;
		++step;

		int bad_e = -1, bad_i = -1, bad_j = -1;
		if (!CheckPhysicalState2D(mesh, bad_e, bad_i, bad_j))
		{
			std::cout << "GKS-subcell 2D Riemann failed at step=" << step
				<< " cell=" << bad_e
				<< " point=(" << bad_i << "," << bad_j << ")" << std::endl;
			PrintBadStateSummary2D(mesh, bad_e, bad_i, bad_j);
			break;
		}
	}

	std::cout << "GKS-subcell 2D Riemann limiter statistics" << std::endl;
	std::cout << "min rho=" << diag.min_rho << std::endl;
	std::cout << "min p=" << diag.min_p << std::endl;
	std::cout << "max alpha=" << diag.max_alpha << std::endl;
	std::cout << "number of troubled cells alpha > 0=" << diag.troubled_cells << std::endl;
	std::cout << "number of flux-limited faces x/y="
		<< diag.flux_diag.limited_faces_x << "/" << diag.flux_diag.limited_faces_y << std::endl;
	std::cout << "min theta_F=" << diag.flux_diag.min_theta_x << std::endl;
	std::cout << "min theta_G=" << diag.flux_diag.min_theta_y << std::endl;
	std::cout << "number of scaling-limited cells=" << diag.scaling_diag.limited_cells << std::endl;
	WriteCellCenterDensityTecplot2D(mesh, StepTaggedPath("gks_subcell_2d_riemann", step, ".plt").c_str());
	WriteFRSolutionPointDensityTecplot2D(mesh, StepTaggedPath("gks_subcell_2d_riemann_frpoints", step, ".plt").c_str());
}

void double_mach_reflection_2d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_2D(0.1, 1.0);
	GKSSubcellFrameworkConfig2D config;
	config.blend_mode = gks_subcell2d_hybrid;
	config.low_mode = MUSCL_HANCOCK_2d;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;
	config.flux_param.kx = 0.5;
	config.flux_param.ky = 0.5;

	GKSFRMesh2D mesh;
	GKSFR_ResizeUniformMesh2D(mesh, 600, 150, 0.0, 4.0, 0.0, 1.0);
	ICforDoubleMachReflection2D(mesh);

	const double CFL = 0.2;  
	const double tstop = 0.2;
	double t = 0.0;
	int step = 0;
	GKSSubcellFrameworkDiag2D diag;
	while (t < tstop - 1.0e-14)
	{
		const double dt = GetTimeStep2D(mesh, CFL, t, tstop);
		std::cout << "step=" << step << " dt=" << dt << " t=" << t << std::endl;
		GKSFR_SetBoundaryTime2D(t + 0.5 * dt);
		GKSSubcellAdvanceOneStep2D(mesh, dt, gksfr2d_double_mach, config, diag);
		t += dt;
		++step;

		int bad_e = -1, bad_i = -1, bad_j = -1;
		if (!CheckPhysicalState2D(mesh, bad_e, bad_i, bad_j))
		{
			std::cout << "Double Mach reflection stopped at step=" << step
				<< " cell=" << bad_e
				<< " point=(" << bad_i << "," << bad_j << ")" << std::endl;
			PrintBadStateSummary2D(mesh, bad_e, bad_i, bad_j);
			break;
		}
	}

	std::cout << "Double Mach reflection limiter statistics" << std::endl;
	std::cout << "min rho=" << diag.min_rho << std::endl;
	std::cout << "min p=" << diag.min_p << std::endl;
	std::cout << "max alpha=" << diag.max_alpha << std::endl;
	std::cout << "number of troubled cells alpha > 0=" << diag.troubled_cells << std::endl;
	std::cout << "number of flux-limited faces x/y="
		<< diag.flux_diag.limited_faces_x << "/" << diag.flux_diag.limited_faces_y << std::endl;
	std::cout << "min theta_F=" << diag.flux_diag.min_theta_x << std::endl;
	std::cout << "min theta_G=" << diag.flux_diag.min_theta_y << std::endl;
	std::cout << "number of scaling-limited cells=" << diag.scaling_diag.limited_cells << std::endl;
	WriteCellCenterDensityTecplot2D(mesh, StepTaggedPath("double_mach_2d", step, ".plt").c_str());
	WriteFRSolutionPointDensityTecplot2D(mesh, StepTaggedPath("double_mach_2d_frpoints", step, ".plt").c_str());
}

void detonation_shock_diffraction_2d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_2D(0.1, 1.0);
	GKSSubcellFrameworkConfig2D config;
	config.blend_mode = gks_subcell2d_hybrid;
	config.low_mode = MUSCL_HANCOCK_2d;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;
	config.flux_param.kx = 0.5;
	config.flux_param.ky = 0.5;

	const bool use_rectangular_cutout_domain = true;
	const int cells_per_unit = 200;
	GKSFRMesh2D mesh;
	GKSFR_ResizeUniformMesh2D(mesh , 2 * cells_per_unit , 2 * cells_per_unit , 0.0 , 2.0 , 0.0 ,2.0);
	InitializeDetonationShockDiffraction2D(mesh);

	GKSSubcellMask2D mask;
	if (use_rectangular_cutout_domain)
	{
		GKSSubcellBuildRectangularCutoutMask2D(mesh, 0.0, 0.5, 0.0, 1.0, mask);
	}

	const double CFL = 0.1;
	const double tstop = 0.01;
	double t = 0.0;
	int step = 0;
	GKSSubcellFrameworkDiag2D diag;
	while (t < tstop - 1.0e-14)
	{
		const double dt = GetTimeStep2D(mesh, CFL, t, tstop);
		std::cout << "step=" << step << " dt=" << dt << " t=" << t << std::endl;
		GKSFR_SetBoundaryTime2D(t + 0.5 * dt);
		if (use_rectangular_cutout_domain)
		{
			GKSSubcellAdvanceOneStep2DMasked(mesh,dt,gksfr2d_detonation_diffraction,mask,config,diag);
		}
		else
		{
			GKSSubcellAdvanceOneStep2D(mesh,dt,gksfr2d_detonation_diffraction,config,diag);
		}
		t += dt;
		++step;

		int bad_e = -1, bad_i = -1, bad_j = -1;
		const bool ok = use_rectangular_cutout_domain
			? CheckPhysicalState2DMasked(mesh, mask, bad_e, bad_i, bad_j)
			: CheckPhysicalState2D(mesh, bad_e, bad_i, bad_j);
		if (!ok)
		{
			std::cout << "Detonation Shock Diffraction stopped at step=" << step
				<< " cell=" << bad_e
				<< " point=(" << bad_i << "," << bad_j << ")" << std::endl;
			PrintBadStateSummary2D(mesh, bad_e, bad_i, bad_j);
			break;
		}
	}

	std::cout << "Detonation Shock Diffraction limiter statistics" << std::endl;
	std::cout << "min rho=" << diag.min_rho << std::endl;
	std::cout << "min p=" << diag.min_p << std::endl;
	std::cout << "max alpha=" << diag.max_alpha << std::endl;
	std::cout << "number of troubled cells alpha > 0=" << diag.troubled_cells << std::endl;
	std::cout << "number of flux-limited faces x/y="
		<< diag.flux_diag.limited_faces_x << "/" << diag.flux_diag.limited_faces_y << std::endl;
	std::cout << "min theta_F=" << diag.flux_diag.min_theta_x << std::endl;
	std::cout << "min theta_G=" << diag.flux_diag.min_theta_y << std::endl;
	std::cout << "number of scaling-limited cells=" << diag.scaling_diag.limited_cells << std::endl;
	if (use_rectangular_cutout_domain)
	{
		WriteCellCenterDensityTecplot2DMasked(
			mesh,
			mask,
			StepTaggedPath("detonation_diffraction_2d_masked", step, ".plt").c_str());
	}
	else
	{
		WriteCellCenterDensityTecplot2D(
			mesh,
			StepTaggedPath("detonation_diffraction_2d_full", step, ".plt").c_str());
	}
}

void astrophysical_jet_2d_gks_subcell()
{
	Ensure_Result_Directory();
	Configure_GKS_Subcell_2D(0.1, 1.0);
	GKSSubcellFrameworkConfig2D config;
	config.blend_mode = gks_subcell2d_hybrid;
	config.low_mode = KFVS1;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;
	config.flux_param.kx = 0.5;
	config.flux_param.ky = 0.5;

	GKSFRMesh2D mesh;
	const int jet_nx = GetEnvInt("JET_NX", 400);
	const int jet_ny = GetEnvInt("JET_NY", 400);
	GKSFR_ResizeUniformMesh2D(mesh, jet_nx, jet_ny, 0.0, 1.0, -0.5, 0.5);
	InitializeAstrophysicalJet2D(mesh);

	const double CFL = GetEnvDouble("JET_CFL", 0.5);
	const double tstop = GetEnvDouble("JET_TSTOP", 0.001);
	double t = 0.0;
	int step = 0;
	double last_dt = 0.0;
	const std::string checkpoint_dir = GetEnvString("JET_CHK_DIR", "debug_jet_checkpoint");
	const std::vector<int> checkpoint_steps = ParseStepList(GetEnvString("JET_CHK_STEPS", ""));
	const int checkpoint_every = GetEnvInt("JET_CHK_EVERY", 0);
	const int max_debug_steps = GetEnvInt("JET_MAX_STEPS", -1);
	const int stop_step = max_debug_steps >= 0 ? step + max_debug_steps : -1;
	const std::string restart_file = GetEnvString("JET_RESTART_FILE", "");
	if (!restart_file.empty())
	{
		if (!LoadJetCheckpoint2D(mesh, step, t, last_dt, restart_file))
		{
			return;
		}
	}
	const int restart_stop_step = max_debug_steps >= 0 ? step + max_debug_steps : -1;
	if (!checkpoint_steps.empty() || checkpoint_every > 0)
	{
		Ensure_Debug_Directory(checkpoint_dir);
	}
	GKSSubcellFrameworkDiag2D diag;
	while (t < tstop - 1.0e-14)
	{
		if ((restart_file.empty() && stop_step >= 0 && step >= stop_step) ||
			(!restart_file.empty() && restart_stop_step >= 0 && step >= restart_stop_step))
		{
			break;
		}
		const double dt = GetTimeStepAstrophysicalJet2D(mesh, CFL, t, tstop);
		std::cout << "step=" << step << " dt=" << dt << " t=" << t << std::endl;
		GKSFR_SetBoundaryTime2D(t + 0.5 * dt);
		SetJetTraceContext(step, t, dt);
		GKSSubcellAdvanceOneStep2D(mesh, dt, gksfr2d_astrophysical_jet, config, diag);
		t += dt;
		++step;
		last_dt = dt;
		if (StepInList(checkpoint_steps, step) || (checkpoint_every > 0 && step % checkpoint_every == 0))
		{
			SaveJetCheckpoint2D(mesh, step, t, last_dt, JetCheckpointPath(checkpoint_dir, step, t));
		}

		int bad_e = -1, bad_i = -1, bad_j = -1;
		if (!CheckPhysicalState2D(mesh, bad_e, bad_i, bad_j))
		{
			std::cout << "Astrophysical jet stopped at step=" << step
				<< " cell=" << bad_e
				<< " point=(" << bad_i << "," << bad_j << ")" << std::endl;
			PrintBadStateSummary2D(mesh, bad_e, bad_i, bad_j);
			break;
		}
	}

	std::cout << "Astrophysical jet limiter statistics" << std::endl;
	std::cout << "min rho=" << diag.min_rho << std::endl;
	std::cout << "min p=" << diag.min_p << std::endl;
	std::cout << "max alpha=" << diag.max_alpha << std::endl;
	std::cout << "number of troubled cells alpha > 0=" << diag.troubled_cells << std::endl;
	std::cout << "number of flux-limited faces x/y="
		<< diag.flux_diag.limited_faces_x << "/" << diag.flux_diag.limited_faces_y << std::endl;
	std::cout << "min theta_F=" << diag.flux_diag.min_theta_x << std::endl;
	std::cout << "min theta_G=" << diag.flux_diag.min_theta_y << std::endl;
	std::cout << "number of scaling-limited cells=" << diag.scaling_diag.limited_cells << std::endl;
WriteCellCenterDensityTecplot2D(mesh, StepTaggedPath("astrophysical_jet_2d", step, ".plt").c_str());
WriteFRSolutionPointDensityTecplot2D(mesh, StepTaggedPath("astrophysical_jet_2d_frpoints", step, ".plt").c_str());
}
