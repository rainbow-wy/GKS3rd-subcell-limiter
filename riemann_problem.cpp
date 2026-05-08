#include"riemann_problem.h"

namespace
{
	const double kPi = 3.14159265358979323846;

	double CellCenterX(const GKSFRMesh1D& mesh, int e)
	{
		return mesh.x_left + (e + 0.5) * mesh.dx;
	}

	double SolutionPointX(const GKSFRMesh1D& mesh, int e, int i)
	{
		return CellCenterX(mesh, e) + 0.5 * mesh.dx * GKSFR_GL_Point(i);
	}

	void SetFRPointPrimitive(double* Q, double rho, double u, double p)
	{
		double prim[3] = { rho, u, p };
		Primvar_to_convar_1D(Q, prim);
	}

	void SetFRPointPrimitive2D(double* Q, const double* prim)
	{
		double local_prim[4] = { prim[0], prim[1], prim[2], prim[3] };
		Primvar_to_convar_2D(Q, local_prim);
	}

	void PrimitiveForProblemAtX(double* prim, const RiemannProblem1D& problem, double x)
	{
		if (problem.profile_type == RiemannProblem1D::three_constant_states)
		{
			if (x < problem.x_discontinuity)
			{
				Copy_Array(prim, const_cast<double*>(problem.left_prim), 3);
				return;
			}
			if (x < problem.x_discontinuity_2)
			{
				Copy_Array(prim, const_cast<double*>(problem.middle_prim), 3);
				return;
			}
			Copy_Array(prim, const_cast<double*>(problem.right_prim), 3);
			return;
		}

		if (x < problem.x_discontinuity)
		{
			Copy_Array(prim, const_cast<double*>(problem.left_prim), 3);
			return;
		}

		Copy_Array(prim, const_cast<double*>(problem.right_prim), 3);
		if (problem.profile_type == RiemannProblem1D::shu_osher)
		{
			prim[0] = problem.right_prim[0]
				+ problem.right_density_amplitude * sin(problem.right_density_wavenumber * x);
		}
	}
}

RiemannProblem1D RiemannProblem1D_Sod()
{
	RiemannProblem1D problem{};
	problem.x_discontinuity = 0.5;
	problem.left_prim[0] = 1.0;
	problem.left_prim[1] = 0.0;
	problem.left_prim[2] = 1.0;
	problem.right_prim[0] = 0.125;
	problem.right_prim[1] = 0.0;
	problem.right_prim[2] = 0.1;
	problem.profile_type = RiemannProblem1D::constant_states;
	return problem;
}

RiemannProblem1D RiemannProblem1D_DoubleRarefaction()
{
	RiemannProblem1D problem{};
	problem.x_discontinuity = 0;
	problem.left_prim[0] = 7.0;
	problem.left_prim[1] = -1.0;
	problem.left_prim[2] = 0.2;
	problem.right_prim[0] = 7.0;
	problem.right_prim[1] = 1.0;
	problem.right_prim[2] = 0.2;
	problem.profile_type = RiemannProblem1D::constant_states;
	return problem;
}

RiemannProblem1D RiemannProblem1D_Leblanc()
{
	RiemannProblem1D problem{};
	problem.x_discontinuity = 0;
	problem.left_prim[0] = 2.0;
	problem.left_prim[1] = 0.0;
	problem.left_prim[2] = 1.0e9;
	problem.right_prim[0] = 0.001;
	problem.right_prim[1] = 0.0;
	problem.right_prim[2] = 1.0;
	problem.profile_type = RiemannProblem1D::constant_states;
	return problem;
}

RiemannProblem1D RiemannProblem1D_ShuOsher()
{
	RiemannProblem1D problem{};
	// Standard 1D Shu-Osher shock-entropy wave interaction:
	// left state is the post-shock constant state, and the right density is
	// 1 + 0.2 sin(5 x) with u = 0, p = 1 for x > -4.
	problem.x_discontinuity = -4.0;
	problem.left_prim[0] = 3.857143;
	problem.left_prim[1] = 2.629369;
	problem.left_prim[2] = 10.333333;
	problem.right_prim[0] = 1.0;
	problem.right_prim[1] = 0.0;
	problem.right_prim[2] = 1.0;
	problem.profile_type = RiemannProblem1D::shu_osher;
	problem.right_density_amplitude = 0.2;
	problem.right_density_wavenumber = 5.0;
	return problem;
}

RiemannProblem1D RiemannProblem1D_BlastWave()
{
	RiemannProblem1D problem{};
	problem.x_discontinuity = 0.1;
	problem.x_discontinuity_2 = 0.9;
	problem.left_prim[0] = 1.0;
	problem.left_prim[1] = 0.0;
	problem.left_prim[2] = 1000.0;
	problem.middle_prim[0] = 1.0;
	problem.middle_prim[1] = 0.0;
	problem.middle_prim[2] = 0.01;
	problem.right_prim[0] = 1.0;
	problem.right_prim[1] = 0.0;
	problem.right_prim[2] = 100.0;
	problem.profile_type = RiemannProblem1D::three_constant_states;
	return problem;
}

RiemannProblem1D RiemannProblem1D_SedovBlastWave(double x_left,double dx)
{
	RiemannProblem1D problem{};
	(void)x_left;
	const double ambient_energy = 1.0e-12;
	const double blast_energy = 3.2e6 / dx;
	const double ambient_pressure = (Gamma - 1.0) * ambient_energy;
	const double blast_pressure = (Gamma - 1.0) * blast_energy;

	problem.x_discontinuity = -0.5 * dx;
	problem.x_discontinuity_2 = 0.5 * dx;
	problem.left_prim[0] = 1.0;
	problem.left_prim[1] = 0.0;
	problem.left_prim[2] = ambient_pressure;
	problem.middle_prim[0] = 1.0;
	problem.middle_prim[1] = 0.0;
	problem.middle_prim[2] = blast_pressure;
	problem.right_prim[0] = 1.0;
	problem.right_prim[1] = 0.0;
	problem.right_prim[2] = ambient_pressure;
	problem.profile_type = RiemannProblem1D::three_constant_states;
	return problem;
}

RiemannProblem2D RiemannProblem2D_SubcellLimiterReference()
{
	RiemannProblem2D problem{};
	problem.x_discontinuity = 0.5;
	problem.y_discontinuity = 0.5;

	problem.upper_right_prim[0] = 0.5313;
	problem.upper_right_prim[1] = 0.0;
	problem.upper_right_prim[2] = 0.0;
	problem.upper_right_prim[3] = 0.4;

	problem.upper_left_prim[0] = 1.0;
	problem.upper_left_prim[1] = 0.7276;
	problem.upper_left_prim[2] = 0.0;
	problem.upper_left_prim[3] = 1.0;

	problem.lower_left_prim[0] = 0.8;
	problem.lower_left_prim[1] = 0.0;
	problem.lower_left_prim[2] = 0.0;
	problem.lower_left_prim[3] = 1.0;

	problem.lower_right_prim[0] = 1.0;
	problem.lower_right_prim[1] = 0.0;
	problem.lower_right_prim[2] = 0.7276;
	problem.lower_right_prim[3] = 1.0;
	return problem;
}

void ICfor2dRM(GKSFRMesh2D& mesh, const RiemannProblem2D& problem)
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
					const double* prim = nullptr;
					if (x >= problem.x_discontinuity && y >= problem.y_discontinuity)
					{
						prim = problem.upper_right_prim;
					}
					else if (x < problem.x_discontinuity && y >= problem.y_discontinuity)
					{
						prim = problem.upper_left_prim;
					}
					else if (x < problem.x_discontinuity && y < problem.y_discontinuity)
					{
						prim = problem.lower_left_prim;
					}
					else
					{
						prim = problem.lower_right_prim;
					}
					SetFRPointPrimitive2D(cell.Q[p][q], prim);
				}
			}
		}
	}
}


void riemann_problem_1d()
{
	Runtime runtime;//statement for recording the running time
	runtime.start_initial = clock();

	Block1d block;
	block.nodex = 400;
	block.ghost = 4;

	double tstop = 1.8;
	block.CFL = 0.5;
	Fluid1d* bcvalue = new Fluid1d[2];
	K = 4;
	Gamma = 1.4;
	R_gas = 1.0;

	gks1dsolver = gks2nd;
	//g0type = collisionn;
	tau_type = Euler;
	//Smooth = false;
	c1_euler = 0.05;
	c2_euler = 1;

	//prepare the boundary condtion function
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	//prepare the reconstruction function

	cellreconstruction = WENO5_AO;
	wenotype = wenoz;
	reconstruction_variable = characteristic;
	g0reconstruction = Center_3rd;
	is_reduce_order_warning = true;
	//prepare the flux function
	flux_function = GKS;
	//prepare time marching stratedgy
	timecoe_list = S1O2;
	Initial_stages(block);

	// allocate memory for 1-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N
	block.nx = block.nodex + 2 * block.ghost;
	block.nodex_begin = block.ghost;
	block.nodex_end = block.nodex + block.ghost - 1;
	Fluid1d* fluids = new Fluid1d[block.nx];
	// then the interfaces reconstructioned value, N+1
	Interface1d* interfaces = new Interface1d[block.nx + 1];
	// then the flux, which the number is identical to interfaces
	Flux1d** fluxes = Setflux_array(block);
	//end the allocate memory part

	//bulid or read mesh,
	//since the mesh is all structured from left and right,
	//there is no need for buliding the complex topology between cells and interfaces
	//just using the index for address searching

	block.left = -5.0;
	block.right = 5.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);

	//ended mesh part

	// then its about initializing, lets first initialize a sod test case
	//you can initialize whatever kind of 1d test case as you like
	//const RiemannProblem1D problem = RiemannProblem1D_Sod();
	//const RiemannProblem1D problem = RiemannProblem1D_DoubleRarefaction();
	//const RiemannProblem1D problem = RiemannProblem1D_Leblanc();
	const RiemannProblem1D problem = RiemannProblem1D_ShuOsher();
	ICfor1dRM(fluids, problem, block);
	//initializing part end

	// === 验证 WENO5_AO_poly 与 weno_5th_ao_left/right 的一致性 ===
	leftboundary(fluids, block, bcvalue[0]);
	rightboundary(fluids, block, bcvalue[1]);
	verify_WENO5_AO_poly(fluids, block);
	// === 验证结束，确认通过后可删除上面三行 ===

	//then we are about to do the loop
	block.t = 0;//the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result
	runtime.finish_initial = clock();
	while (block.t < tstop)
	{
		// assume you are using command window,
		// you can specify a running step conveniently
		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output1d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}

		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);


		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			// here the boudary type, you shall go above the search the key words"BoundaryCondition leftboundary;"
			// to see the pointer to the corresponding function

			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(interfaces, fluids, block);
			Reconstruction_forg0(interfaces, fluids, block);

			//then is solver part
			Calculate_flux(fluxes, interfaces, block, i);

			//then is update flux part
			Update(fluids, fluxes, block, i);
			//output1d_checking(fluids, interfaces, fluxes, block);
		}

		block.step++;
		block.t = block.t + block.dt;
		if (block.step % 10 == 0)
		{
			cout << "step 10 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		if ((block.t - tstop) > 0)
		{
			runtime.finish_compute = clock();
			cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
			cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
			output1d(fluids, block);
			//output1d_checking(fluids, interfaces, fluxes, block);
		}
	}

}

void ICfor1dRM(Fluid1d* fluids, Fluid1d zone1, Fluid1d zone2, Block1d block)
{
	RiemannProblem1D problem{};
	problem.x_discontinuity = 0.5 * (block.left + block.right);
	Copy_Array(problem.left_prim, zone1.primvar, 3);
	Copy_Array(problem.right_prim, zone2.primvar, 3);
	ICfor1dRM(fluids, problem, block);
}

void ICfor1dRM(Fluid1d* fluids, const RiemannProblem1D& problem, Block1d block)
{
	double prim[3];
	for (int i = 0; i < block.nx; i++)
	{
		PrimitiveForProblemAtX(prim, problem, fluids[i].cx);
		Copy_Array(fluids[i].primvar, prim, 3);
	}

	for (int i = 0; i < block.nx; i++)
	{
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);
	}
}

void ICfor1dRM(GKSFRMesh1D& mesh, const RiemannProblem1D& problem)
{
	double prim[3];
	for (int e = 0; e < mesh.cells; ++e)
	{
		for (int i = 0; i < 3; ++i)
		{
			const double x = SolutionPointX(mesh, e, i);
			PrimitiveForProblemAtX(prim, problem, x);
			SetFRPointPrimitive(
				mesh.cell[e].Q[i],
				prim[0],
				prim[1],
				prim[2]);
		}
	}
}

void riemann_problem_2d()
{

	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = true;
	block.nodex = 200;
	block.nodey = 200;
	block.ghost = 3;



	block.CFL = 0.5;
	Fluid2d* bcvalue = new Fluid2d[4];

	K = 3;
	Gamma = 1.4;




	//prepare the boundary condtion function
	BoundaryCondition2dnew leftboundary(0);
	BoundaryCondition2dnew rightboundary(0);
	BoundaryCondition2dnew downboundary(0);
	BoundaryCondition2dnew upboundary(0);

	leftboundary = free_boundary_left;
	rightboundary = free_boundary_right;
	downboundary = free_boundary_down;
	upboundary = free_boundary_up;

	//prepare the reconstruction function

	gausspoint = 1;
	SetGuassPoint();

	reconstruction_variable = characteristic;
	wenotype = wenoz;

	cellreconstruction_2D_normal = WENO5_normal;
	cellreconstruction_2D_tangent = WENO5_tangent;
	g0reconstruction_2D_normal = Center_do_nothing_normal;
	g0reconstruction_2D_tangent = Center_all_collision_multi;

	is_reduce_order_warning = true; //重构出负后 输出出负的单元

	//prepare the flux function
	gks2dsolver = gks2nd_2d;
	tau_type = Euler;
	c1_euler = 0.05;
	c2_euler = 1;
	flux_function_2d = GKS2D;

	//prepare time marching stratedgy


	//time coe list must be 2d
	timecoe_list_2d = S2O4_2D;
	Initial_stages(block);


	// allocate memory for 2-D fluid field
	// in a standard finite element method, we have 
	// first the cell average value, N*N
	block.nx = block.nodex + 2 * block.ghost;
	block.ny = block.nodey + 2 * block.ghost;

	Fluid2d* fluids = Setfluid(block);
	// then the interfaces reconstructioned value, (N+1)*(N+1)
	Interface2d* xinterfaces = Setinterface_array(block);
	Interface2d* yinterfaces = Setinterface_array(block);
	// then the flux, which the number is identical to interfaces
	Flux2d_gauss** xfluxes = Setflux_gauss_array(block);
	Flux2d_gauss** yfluxes = Setflux_gauss_array(block);
	//end the allocate memory part

	block.left = 0.0;
	block.right = 2.0;
	block.down = 0.0;
	block.up = 2.0;
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes);
	//ended mesh part

	//RM 1 T=0.2
	// double tstop[] = { 0.2 };
	// double zone1[] = { 0.8, 0, 0, 1.0 };
	// double zone2[]{ 1.0, 0, 0.7276, 1.0 };
	// double zone3[]{ 0.5313, 0, 0, 0.4 };
	// double zone4[]{ 1.0, 0.7276, 0, 1.0 };

	////RM config 1 four rarefaction wave T=0.2
	// double tstop[] = { 0.2 };
	// double zone1[]={ 0.1072,-0.7259,-1.4045,0.0439 };
	// double zone2[]={ 0.2579,0,-1.4045,0.15};
	// double zone3[]={ 1,0,0,1};
	// double zone4[]={ 0.5197,-0.7259,0,0.4};

	// double tstop[] = { 0.2 };
	// double zone1[]={ 1.0, 0, 0, 1.0 };
	// double zone4[]={ 1.0, 0, 0, 1.0 };
	// double zone3[]={ 0.125, 0, 0, 0.1 };
	// double zone2[]={ 0.125, 0, 0, 0.1 };

	//double tstop[] = { 0.2 };
	//Initial zone1{ 1.0, 0, 0, 1.0 };
	//Initial zone2{ 1.0, 0, 0, 1.0 };
	//Initial zone3{ 0.125, 0, 0, 0.1 };
	//Initial zone4{ 0.125, 0, 0, 0.1 };

	////RM 3 T=0.8 with x,y=0.5
	double tstop[] = { 0.6 };
	double p0 = 1.0; // this p0 is for adjust the speed of shear layer
	double zone1[]={ 1, -0.75, 0.5, p0 };
	double zone2[]={ 3.0, -0.75, -0.5, p0 };
	double zone3[]={ 1.0, 0.75, -0.5, p0 };
	double zone4[]={ 2.0, 0.75, 0.5, p0 };

	////RM 2 T=0.6 with x,y = 0.7
	//double tstop[] = { 0.3, 0.4, 0.5 };
	//Initial zone1{ 0.138, 1.206, 1.206, 0.029 };
	//Initial zone2{ 0.5323, 0, 1.206, 0.3 };
	//Initial zone3{ 1.5, 0, 0, 1.5 };
	//Initial zone4{ 0.5323, 1.206, 0, 0.3 };

	////RM 4 T=0.25 with x,y=0.5
	//double tstop[] = { 0.25,0.3};
	//double p0 = 1.0; // this p0 is for adjust the speed of shear layer
	//Initial zone3{ 1, -0.75, -0.5, p0 };
	//Initial zone4{ 2.0, -0.75, 0.5, p0 };
	//Initial zone1{ 1.0, 0.75, 0.5, p0 };
	//Initial zone2{ 3.0, 0.75, -0.5, p0 };

	IC_for_riemann_2d(fluids, zone1, zone2, zone3, zone4, block);

	runtime.finish_initial = clock();
	block.t = 0;//the current simulation time
	block.step = 0; //the current step
	int tstop_dim = sizeof(tstop) / sizeof(double);

	int inputstep = 1000000;
	//int inputstep = 1;//input a certain step,
					  //initialize inputstep=1, to avoid a 0 result

	//在若干个模拟时间输出
	for (int instant = 0; instant < tstop_dim; ++instant)
	{

		if (inputstep == 0)
		{
			break;
		}
		while (block.t < tstop[instant])
		{

			if (block.step % inputstep == 0)
			{
				cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
				cin >> inputstep;
				if (inputstep == 0)
				{
					output2d_binary(fluids, block);
					break;
				}
			}
			if (runtime.start_compute == 0.0)
			{
				runtime.start_compute = clock();
				cout << "runtime-start " << endl;
			}

			CopyFluid_new_to_old(fluids, block);

			block.dt = Get_CFL(block, fluids, tstop[instant]);

			for (int i = 0; i < block.stages; i++)
			{
				leftboundary(fluids, block, bcvalue[0]);
				rightboundary(fluids, block, bcvalue[1]);
				downboundary(fluids, block, bcvalue[2]);
				upboundary(fluids, block, bcvalue[3]);

				Convar_to_Primvar(fluids, block);

				Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block);

				Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block);

				Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i);

				Update(fluids, xfluxes, yfluxes, block, i);

			}

			block.step++;
			block.t = block.t + block.dt;

			if (block.step % 10 == 0)
			{
				cout << "step 10 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
			}

			if ((block.t - tstop[instant]) > 0)
			{
				output2d(fluids, block);
			}
			//if (block.step % 100 == 0)
			//{
			//	output2d(fluids, block);
			//}
			if (fluids[block.nx / 2 * (block.ny) + block.ny / 2].convar[0] != fluids[block.nx / 2 * (block.ny) + block.ny / 2].convar[0])
			{
				cout << "the program blows up!" << endl;
				output2d(fluids, block);
				break;
			}
		}
	}
	runtime.finish_compute = clock();
	cout << "the total run time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "initializing time is" << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "the pure computational time is" << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << " second !" << endl;

}

void IC_for_riemann_2d(Fluid2d* fluid, double* zone1, double* zone2, double* zone3, double* zone4, Block2d block)
{
	double discon[2];
	discon[0] = 0.5;
	discon[1] = 0.5;
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			if (i < discon[0] * block.nx && j < discon[1] * block.ny)
			{
				Copy_Array(fluid[i * block.ny + j].primvar, zone1, 4);
			}
			else
			{
				if (i >= discon[0] * block.nx && j < discon[1] * block.ny)
				{
					Copy_Array(fluid[i * block.ny + j].primvar, zone2, 4);
				}
				else
				{
					if (i >= discon[0] * block.nx && j >= discon[1] * block.ny)
					{
						Copy_Array(fluid[i * block.ny + j].primvar, zone3, 4);
					}
					else
					{
						Copy_Array(fluid[i * block.ny + j].primvar, zone4, 4);
					}
				}
			}

		}

	}

	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			Primvar_to_convar_2D(fluid[i * block.ny + j].convar, fluid[i * block.ny + j].primvar);

		}
	}
}
