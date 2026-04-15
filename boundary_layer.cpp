#include"boundary_layer.h"


void boundary_layer()
{

	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block;
	block.uniform = false;
	block.ghost = 3;

	double tstop = 20000;
	block.CFL = 0.3;

	K = 3;
	Gamma = 1.4;
	Pr = 1.0;

	double renum = 1e5;
	double den_ref = 1.0;
	double u_ref = 0.15;
	double L = 100;

	//prepare the boundary condtion function
	Fluid2d* bcvalue = new Fluid2d[1];
	bcvalue[0].primvar[0] = den_ref;
	bcvalue[0].primvar[1] = u_ref;
	bcvalue[0].primvar[2] = 0.0;
	bcvalue[0].primvar[3] = den_ref / Gamma;

	Mu = den_ref * u_ref * L / renum;
	cout << Mu << endl;



	//prepare the reconstruction function

	gausspoint = 1;
	SetGuassPoint();
	reconstruction_variable = conservative;
	wenotype = wenoz;//只有WENO才会起作用
	cellreconstruction_2D_normal = Vanleer_normal;
	cellreconstruction_2D_tangent = Vanleer_tangent;
	g0reconstruction_2D_normal = Center_3rd_normal;
	g0reconstruction_2D_tangent = Center_3rd_tangent;

	is_reduce_order_warning = true;

	flux_function_2d = GKS2D;
	gks2dsolver = gks2nd_2d;
	tau_type = NS;
	c1_euler = 0.0;
	c2_euler = 0;

	timecoe_list_2d = S1O2_2D;


	Initial_stages(block);


	Fluid2d* fluids = NULL;
	Interface2d* xinterfaces = NULL;
	Interface2d* yinterfaces = NULL;
	Flux2d_gauss** xfluxes = NULL;
	Flux2d_gauss** yfluxes = NULL;

	//读入2维结构化网格，先找到网格名字
	string grid = add_mesh_directory_modify_for_linux()
		+ "C:/code/gks2d-str/structured-mesh/flat-plate-str.plt";

	//读入2维结构化网格，格式是tecplot的3维binary，z方向只有一层。可以通过pointwise输出tecplot的binary格式实现
	Read_structured_mesh
	(grid, &fluids, &xinterfaces, &yinterfaces, &xfluxes, &yfluxes, block);

	ICfor_uniform_2d(fluids, bcvalue[0].primvar, block);

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
				output2d(fluids, block);
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}

		CopyFluid_new_to_old(fluids, block);
		block.dt = Get_CFL(block, fluids, tstop);

		for (int i = 0; i < block.stages; i++)
		{

			boundaryforBoundary_layer(fluids, block, bcvalue[0]);

			Convar_to_Primvar(fluids, block);

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
		output_wall_info_boundary_layer(xinterfaces, yinterfaces, xfluxes, yfluxes, block, L);
		output_skin_friction_boundary_layer(fluids, yinterfaces, yfluxes, block, u_ref, renum, L);
		if (block.step % 5000 == 0 || abs(block.t - tstop) < 1e-14)
		{
			output2d(fluids, block);
			char loc1[] = "0050";
			char loc2[] = "0105";
			char loc3[] = "0165";
			char loc4[] = "0265";
			char loc5[] = "0305";
			char loc6[] = "0365";
			char loc7[] = "0465";
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.050 * L, loc1);
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.105 * L, loc2);
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.165 * L, loc3);
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.265 * L, loc4);
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.305 * L, loc5);
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.365 * L, loc6);
			output_boundary_layer(fluids, block, bcvalue[0].primvar, 0.465 * L, loc7);

		}

	}
	runtime.finish_compute = clock();
	;
	cout << "\n the total running time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << "秒！" << endl;
	cout << "\n the time for initializing is " << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "秒！" << endl;
	cout << "\n the time for computing is " << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "秒！" << endl;

}

void boundaryforBoundary_layer(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{

	inflow_boundary_left(fluids, block, bcvalue);

	inflow_boundary_right(fluids, block, bcvalue);
	//free_boundary_up(fluids, block, bcvalue);
	//outflow_boundary_up(fluids, block, bcvalue);
	inflow_boundary_up(fluids, block, bcvalue);

	//the boundary at the bottom
	int order = block.ghost;

	for (int j = order - 1; j >= 0; j--)
	{
		for (int i = order; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			int ref = i * block.ny + 2 * order - 1 - j;
			if (fluids[index].coordx < 0.0)
			{

				fluids[index].convar[0] = fluids[ref].convar[0];
				fluids[index].convar[1] = fluids[ref].convar[1];
				fluids[index].convar[2] = -fluids[ref].convar[2];
				fluids[index].convar[3] = fluids[ref].convar[3];



			}
			else
			{
				fluids[index].convar[0] = fluids[ref].convar[0];
				fluids[index].convar[1] = -fluids[ref].convar[1];
				fluids[index].convar[2] = -fluids[ref].convar[2];
				fluids[index].convar[3] = fluids[ref].convar[3];


			}


		}

	}
}

void output_boundary_layer(Fluid2d* fluids, Block2d block, double* bcvalue, double coordx, char* label)
{
	int order = block.ghost;
	int i;
	for (int k = 0; k < block.nx; k++)
	{
		int index = k * block.ny + order;
		if (fluids[index].coordx > coordx)
		{
			i = k;
			break;
		}
	}

	stringstream name;
	name << "C:/code/gks2d-str/code/result/boundarylayer_line" << label << setfill('0') << setw(5) << "step" << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "variables =ys,density,us,vs,pressure" << endl;
	out << "zone i = " << 1 << ",j = " << block.nodey << ", F=POINT" << endl;

	double vel = bcvalue[1];
	double nu = Mu / bcvalue[0];
	for (int j = order; j < block.ny - order; j++)
	{
		int index = i * block.ny + j;
		Convar_to_primvar_2D(fluids[index].primvar, fluids[index].convar);

		out << fluids[index].coordy / (sqrt(nu * fluids[index].coordx / vel)) << " "
			<< fluids[index].primvar[0] << " "
			<< fluids[index].primvar[1] / vel << " "
			<< fluids[index].primvar[2] / (sqrt(nu * vel / fluids[index].coordx)) << " "
			<< fluids[index].primvar[3] << " "
			<< endl;
	}

	out.close();

}

void output_wall_info_boundary_layer(Interface2d* xInterfaces, Interface2d* yInterfaces,
	Flux2d_gauss** xFluxes, Flux2d_gauss** yFluxes, Block2d block, double L_ref)
{
	//there are four blocks
	int outstep = 20;
	if (block.step % outstep != 0)
	{
		return;
	}
	double t = block.t;
	double dt = block.timecoefficient[block.stages - 1][0][0] * block.dt;
	double var[8];
	Array_zero(var, 8);


	//the contribute from the bottom
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{

		int ref = i * (block.ny + 1) + block.ghost;

		if (yInterfaces[ref].x > 0.0)
		{
			double f_local[4];
			double derf_local[4];
			for (int m = 0; m < gausspoint; m++)
			{
				Local_to_Global(f_local, yFluxes[ref][0].gauss[m].f, yFluxes[ref][0].gauss[m].normal);
				Local_to_Global(derf_local, yFluxes[ref][0].gauss[m].derf, yFluxes[ref][0].gauss[m].normal);
			}
			double l = yInterfaces[ref].length;
			double sign = 1;
			for (int k = 0; k < 4; k++)
			{
				for (int m = 0; m < gausspoint; m++)
				{
					var[k] += sign * l * yFluxes[ref][block.stages - 1].gauss[m].x[k] / dt;
					var[k + 4] += sign * l * gauss_weight[m] * (f_local[k]) / block.dt;
				}
			}
		}

	}

	ofstream out;
	if (block.step == outstep)
	{
		out.open("result/boundary-layer-force-history.plt", ios::out);
	}
	else
	{
		out.open("result/boundary-layer-force-history.plt", ios::ate | ios::out | ios::in);
	}

	if (!out.is_open())
	{
		cout << "cannot find square-cylinder-force-history.plt" << endl;
		cout << "a new case will start" << endl;
	}

	if (block.step == outstep)
	{
		out << "# CFL number is " << block.CFL << endl;
		out << "variables = t,mass_transfer,cd,cl,heatflux,";
		out << "mass_transfer_instant, cd_instant, cl_instant, heatflux_instant" << endl;
	}

	Array_scale(var, 8, L_ref);

	out << block.t << " "
		<< var[0] << " " << var[1] << " "
		<< var[2] << " " << var[3] << " ";
	out << var[4] << " " << var[5] << " "
		<< var[6] << " " << var[7] << " ";
	out << endl;
	out.close();
}


void output_skin_friction_boundary_layer
(Fluid2d* fluids, Interface2d* yInterfaces, Flux2d_gauss** yFluxes, Block2d block, double u_ref, double re_ref, double L_ref)
{
	int outstep = 1000;
	if (block.step % outstep != 0)
	{
		return;
	}
	double t = block.t;
	double dt = block.timecoefficient[0][0][0] * block.dt;
	double var[3];
	Array_zero(var, 3);

	ofstream out;
	out.open("result/boundary-layer-skin-friction.plt", ios::out);
	out << "variables = log<sub>10</sub>Re<sub>x</sub>, log<sub>10</sub>C<sub>f</sub>, C<sub>p</sub>";
	out << endl;
	//output the data
	//the contribute from the bottom
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		int ref = i * (block.ny + 1) + block.ghost;
		int ref_cell = i * (block.ny) + block.ghost;
		if (yInterfaces[ref].x > 0.0)
		{
			Array_zero(var, 3);
			for (int k = 1; k < 3; k++)
			{
				for (int m = 0; m < gausspoint; m++)
				{
					var[k] += yFluxes[ref][0].gauss[m].x[k] / dt;
				}
			}
			var[1] = -var[1];
			var[2] = (var[2] - 1.0 / Gamma);
			double rex = yInterfaces[ref].x / L_ref * re_ref;

			for (int k = 0; k < 3; ++k)
			{
				var[k] /= (0.5 * u_ref * u_ref);
			}
			out << log10(rex) << " " << log10(var[1]) << " " << var[2] << endl;

		}
	}
	out.close();
}