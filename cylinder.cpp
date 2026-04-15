#include"cylinder.h"

void cylinder()
{
	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = false;
	block.ghost = 3;

	double tstop = 1;
	block.CFL = 0.5;

	K = 3;
	Gamma = 1.4;
	Pr = 1.0;

	double renum = 1e3;
	double den_ref = 1.0;
	double u_ref = 5.0;
	double L = 1;  //cylinder diameter 1

	Fluid2d icvalue;
	icvalue.primvar[0] = den_ref;
	icvalue.primvar[1] = u_ref;
	icvalue.primvar[2] = 0.0;
	icvalue.primvar[3] = den_ref / Gamma;

	Mu = den_ref * u_ref * L / renum;
	cout << Mu << endl;

	gausspoint = 1;
	SetGuassPoint();
	reconstruction_variable = conservative;
	wenotype = wenoz;


	cellreconstruction_2D_normal = Vanleer_normal;
	cellreconstruction_2D_tangent = Vanleer_tangent;
	g0reconstruction_2D_normal = Center_do_nothing_normal;
	g0reconstruction_2D_tangent = Center_all_collision_multi;

	is_reduce_order_warning = false;
	flux_function_2d = GKS2D;
	gks2dsolver = gks2nd_2d;
	tau_type = NS;
	c1_euler = 0.05;
	c2_euler = 5; //强间断 c2建议取大于等于5

	timecoe_list_2d = S1O2_2D;

	Initial_stages(block);

	Fluid2d* fluids = NULL;
	Interface2d* xinterfaces = NULL;
	Interface2d* yinterfaces = NULL;
	Flux2d_gauss** xfluxes = NULL;
	Flux2d_gauss** yfluxes = NULL;

	string grid = add_mesh_directory_modify_for_linux()
		+ "C:/code/gks2d-str/structured-mesh/half-cylinder-str-1.plt";

	Read_structured_mesh
	(grid, &fluids, &xinterfaces, &yinterfaces, &xfluxes, &yfluxes, block);

	ICfor_uniform_2d(fluids, icvalue.primvar, block);

	runtime.finish_initial = clock();
	block.t = 0;
	block.step = 0;

	int inputstep = 1;
	while (block.t < tstop)
	{
		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output2d_center(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}

		CopyFluid_new_to_old(fluids, block);
		Convar_to_Primvar(fluids, block);

		block.dt = Get_CFL(block, fluids, tstop);

		for (int i = 0; i < block.stages; i++)
		{

			boundary_for_cylinder(fluids, xinterfaces, block, icvalue);

			Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

			Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block);

			Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i);

			Update(fluids, xfluxes, yfluxes, block, i);

		}


		++block.step;
		block.t = block.t + block.dt;

		if (block.step % 100 == 0)
		{
			cout << "step 10 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		Residual2d(fluids, block, 10);

	}
	runtime.finish_compute = clock();
	cout << "\n the total running time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << " second!" << endl;
	cout << "\n the time for initializing is " << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << " second!" << endl;
	cout << "\n the time for computing is " << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << " second!" << endl;

	//output_prim_variable_at_final_time
	output2d_center(fluids, block);
}


void boundary_for_cylinder(Fluid2d* fluids, Interface2d* faces, Block2d block, Fluid2d bcvalue)
{
	Fluid2d wall_vel;
	wall_vel.primvar[0] = bcvalue.primvar[0];
	wall_vel.primvar[1] = 0.0;
	wall_vel.primvar[2] = 0.0;
	wall_vel.primvar[3] = bcvalue.primvar[3];

	noslip_adiabatic_boundary_right(fluids, faces, block, wall_vel);

	inflow_boundary_left(fluids, block, bcvalue);

	free_boundary_down(fluids, block, bcvalue);
	free_boundary_up(fluids, block, bcvalue);


}
