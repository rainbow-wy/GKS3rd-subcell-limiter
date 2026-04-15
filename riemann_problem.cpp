#include"riemann_problem.h"

void riemann_problem_1d()
{
	Runtime runtime;//statement for recording the running time
	runtime.start_initial = clock();

	Block1d block;
	block.nodex = 100;
	block.ghost = 4;

	double tstop = 0.2;
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

	cellreconstruction = Vanleer;
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

	block.left = 0.0;
	block.right = 1.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);

	//ended mesh part

	// then its about initializing, lets first initialize a sod test case
	//you can initialize whatever kind of 1d test case as you like
	Fluid1d zone1; zone1.primvar[0] = 1.0; zone1.primvar[1] = 0.0; zone1.primvar[2] = 1.0;
	Fluid1d zone2; zone2.primvar[0] = 0.125; zone2.primvar[1] = 0.0; zone2.primvar[2] = 0.1;

	//Fluid1d zone1; zone1.primvar[0] = 1.0; zone1.primvar[1] = -2.0; zone1.primvar[2] = 0.4;
	//Fluid1d zone2; zone2.primvar[0] = 1.0; zone2.primvar[1] = 2.0; zone2.primvar[2] = 0.4;
	ICfor1dRM(fluids, zone1, zone2, block);
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
	for (int i = 0; i < block.nx; i++)
	{
		if (i < 0.5 * block.nx)
		{
			Copy_Array(fluids[i].primvar, zone1.primvar, 3);
		}
		else
		{
			Copy_Array(fluids[i].primvar, zone2.primvar, 3);
		}

	}

	for (int i = 0; i < block.nx; i++)
	{
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);
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