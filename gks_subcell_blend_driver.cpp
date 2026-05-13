#include "gks_subcell_blend_driver.h"

#include "function.h"
#include "gks_basic.h"
#include "gks_fr_adapter.h"
#include "riemann_problem.h"

#include <algorithm>
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

	std::vector<GKSFRFaceFlux2D> low_x_face_fluxes;
	std::vector<GKSFRFaceFlux2D> low_y_face_fluxes;
	GKSSubcellComputeLowFaceFluxes2D(low_branch_kfvs_safe, dt, boundary, low_x_face_fluxes, low_y_face_fluxes);

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

	GKSFRMesh2D high_new;
	GKSFRAdapterAdvanceWithFaceFluxes2D(mesh, dt, final_x_face_fluxes, final_y_face_fluxes, high_new);

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
	GKSFR_ResizeUniformMesh(mesh, 400, 0, 1.0);
	//ICfor1dRM(mesh, RiemannProblem1D_Sod());
	//ICfor1dRM(mesh,RiemannProblem1D_DoubleRarefaction());
	//ICfor1dRM(mesh,RiemannProblem1D_Leblanc());
	//ICfor1dRM(mesh,RiemannProblem1D_ShuOsher());//shu-osher问题
	ICfor1dRM(mesh,RiemannProblem1D_BlastWave());
	//ICfor1dRM(mesh, RiemannProblem1D_SedovBlastWave(-2.0, 4.0 / 201.0));

	GKSSubcellFrameworkDiag1D diag;
	int final_step = 0;
	const bool ok = AdvanceCase(mesh, 0.02, 0.038, gksfr_reflective, config, diag, true, &final_step);
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
	config.low_mode = KFVS1;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;
	config.flux_param.kx = 0.5;
	config.flux_param.ky = 0.5;

	GKSFRMesh2D mesh;
	GKSFR_ResizeUniformMesh2D(mesh, 240, 60, 0.0, 4.0, 0.0, 1.0);
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
