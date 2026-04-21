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

	std::vector<GKSFRFaceFlux1D> low_face_fluxes;
	GKSSubcellComputeLowFaceFluxes1D(low_branch, dt, boundary, low_face_fluxes);

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
			low_branch,
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
	GKSSubcellComputeElementResiduals1D(low_branch_safe);
	GKSSubcellAddFaceResiduals1D(low_branch_safe, final_face_fluxes);
	GKSSubcellUpdateLowOrderDofs1D(low_branch_safe, dt);

	GKSSubcellComputeElementResiduals1D(low_branch);
	GKSSubcellAddFaceResiduals1D(low_branch, final_face_fluxes);
	GKSSubcellScaleResiduals1D(low_branch, diag.alpha_final);

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
	config.low_mode = KFVS1;
	config.use_flux_limiter = true;
	config.use_scaling_limiter = true;

	GKSFRMesh1D mesh;
	GKSFR_ResizeUniformMesh(mesh, 400, -5.0, 5.0);
	//ICfor1dRM(mesh, RiemannProblem1D_Sod());
	//ICfor1dRM(mesh,RiemannProblem1D_DoubleRarefaction());
	//ICfor1dRM(mesh,RiemannProblem1D_Leblanc());
	ICfor1dRM(mesh,RiemannProblem1D_ShuOsher());

	GKSSubcellFrameworkDiag1D diag;
	int final_step = 0;
	const bool ok = AdvanceCase(mesh, 0.02, 1.8, gksfr_transmissive_strict, config, diag, true, &final_step);
	if (!ok)
	{
		std::cout << "GKS-subcell run stopped because of an invalid state." << std::endl;
		return;
	}

	WriteSodTecplot(mesh, StepTaggedPath("gks_subcell", final_step, ".plt").c_str());
}
