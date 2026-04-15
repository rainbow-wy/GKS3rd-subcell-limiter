#include"accuracy_test.h"

void accuracy_sinwave_1d()
{
	int mesh_set = 3; //用一系列网格尺寸等比例减小的网格来测试格式的精度
	int mesh_number_start = 10; //最粗网格数
	double length = 2.0; //计算域总长度为2，[0,2]
	double CFL = 0.5; //CFL数

	// 注意，为了防止时间推进和空间离散的精度不一致，dt = CFL * （dx^dt_ratio），这一点的原因也可以参看
	// X. Ji et al. Journal of Computational Physics 356 (2018) 150–173 的 subsection 4.1
	double dt_ratio = 1.0;
	//end

	double mesh_size_start = length / mesh_number_start; //均匀网格大小

	//各套网格数量-尺寸-误差 开辟数组存一下
	int* mesh_number = new int[mesh_set];
	double* mesh_size = new double[mesh_set];
	double** error = new double* [mesh_set];
	//end

	for (int i = 0; i < mesh_set; i++)
	{
		error[i] = new double[3]; //每一套网格都计算 L1 L2 L_inf 误差
		mesh_size[i] = mesh_size_start / pow(2, i); //等比变化
		mesh_number[i] = mesh_number_start * pow(2, i); //等比变化

		//计算每一套网格的工况,得到对应的误差
		accuracy_sinwave_1d(CFL, dt_ratio, mesh_number[i], error[i]);
		//end
	}

	output_error_form(CFL, dt_ratio, mesh_set, mesh_number, error);
}
void accuracy_sinwave_1d(double& CFL, double& dt_ratio, int& mesh_number, double* error)
{
	Runtime runtime;//记录一下计算时间
	runtime.start_initial = clock();

	Block1d block; //存储一些一维算例相关的基本信息

	block.nodex = mesh_number;
	block.ghost = 3; // ghost cell 数。一般来讲 二阶-三阶 需要两个，四阶-五阶 需要三个。具体算法则需具体讨论

	double tstop = 2.0; //计算一个sinwave的周期

	block.CFL = CFL; // CFL数



	K = 4; //内部自由度
	Gamma = 1.4; //比热容是由内部自由度决定的，这里为了简单，直接手动给出

	tau_type = Euler; //无粘还是粘性问题

	//准备通量函数
	//函数指针 flux_function 可以方便的切换 使用的通量求解器
	//如 L-F， HLLC, GKS
	flux_function = GKS; //GKS
	//通常 GKS是指 2阶gks，但通过一些近似，可以变形为
	//一阶GKS（gks1st） 一阶KFVS(kfvs1st) 等，在gks1dsolver中枚举出来
	gks1dsolver = gks3rd; //现在选择三阶gks
	c1_euler = 0.0; //仅对于无粘流动，gks里面的人工碰撞时间系数c1
	c2_euler = 0.0; //对于无粘或者粘性流动，gks里面的人工碰撞时间系数c2
	//end

	//边界条件函数设定
	Fluid1d* bcvalue = new Fluid1d[2]; //一维左右dirichlet边界的值
	BoundaryCondition leftboundary(0);
	BoundaryCondition rightboundary(0);
	leftboundary = periodic_boundary_left;
	rightboundary = periodic_boundary_right;
	//end

	//设置重构函数
	cellreconstruction = WENO5_AO; //以单位为中心的重构，重构出界面左值和右值
	wenotype = wenoz; //非线性权的种类,linear就意味着采用线性权
	reconstruction_variable = characteristic; //重构变量的类型，分为特征量和守恒量
	g0reconstruction = Center_collision;//平衡态函数g0的重构方式（传统黎曼求解器用此函数得到粘性项所需要的一阶导数）
	//end

	//显式时间推进选取，支持N步M阶时间精度，N目前写死了小于等于5
	timecoe_list = S1O3;
	//timecoe_list = S2O4;//这叫做2步4阶 two-stage fourth-order
	Initial_stages(block); //赋值一下
	//end

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
	block.right = 2.0;
	block.dx = (block.right - block.left) / block.nodex;
	//set the uniform geometry information
	SetUniformMesh(block, fluids, interfaces, fluxes);


	//初始化一下
	ICfor_sinwave(fluids, block);

	//初始化时间和时间步长
	block.t = 0;//the current simulation time
	block.step = 0; //the current step
	int inputstep = 1;//输入一个运行的时间步数

	runtime.finish_initial = clock(); //记录一下初始化需要的cpu时间....

	while (block.t < tstop) //显式时间步循环，from t^n to t^n+1
	{
		//输入一个运行时间步数的主要好处是，你可以在开发阶段，先让程序跑几步停一下看看结果
		//避免跑了一万步，发现从第一部就爆掉了
		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output1d(fluids, block); //这个是输出一维的流场数据
				break;
			}
		}

		//记录运算时间
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}
		//end

		//Copy the fluid vales to fluid old
		CopyFluid_new_to_old(fluids, block);

		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop);


		for (int i = 0; i < block.stages; ++i)
		{
			//边界条件
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);

			//重构
			Reconstruction_within_cell(interfaces, fluids, block);
			Reconstruction_forg0(interfaces, fluids, block);

			//计算通量
			Calculate_flux(fluxes, interfaces, block, i);

			//更新
			Update(fluids, fluxes, block, i);
		}

		block.step++;
		block.t = block.t + block.dt;
		if ((block.t - tstop) > 0)
		{
			error_for_sinwave(fluids, block, tstop, error);
		}
	}
	runtime.finish_compute = clock();
	cout << "initializiation time is " << (float)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << "seconds" << endl;
	cout << "computational time is " << (float)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << "seconds" << endl;
	output1d(fluids, block);
}

void ICfor_sinwave(Fluid1d* fluids, Block1d block)
{
#pragma omp parallel for 
	for (int i = block.nodex_begin; i <= block.nodex_end; i++)
	{
		double pi = 3.14159265358979323846;
		fluids[i].primvar[0]
			= 1.0 - 0.2 / pi / block.dx *
			(cos(pi * (i + 1 - block.ghost) * block.dx) - cos(pi * (i - block.ghost) * block.dx));
		fluids[i].primvar[1] = 1;
		fluids[i].primvar[2] = 1;
		Primvar_to_convar_1D(fluids[i].convar, fluids[i].primvar);

	}
}

void error_for_sinwave(Fluid1d* fluids, Block1d block, double tstop, double* error)
{

	if (abs(block.t - tstop) <= (1e-10))
	{

		double error1 = 0;
		double error2 = 0;
		double error3 = 0;

		for (int i = block.ghost; i < block.nx - block.ghost; i++)
		{
			double pi = 3.14159265358979323846;
			int index = i;
			double primvar0 = 1 - 0.2 / pi / block.dx * (cos(pi * (i + 1 - block.ghost) * block.dx) - cos(pi * (i - block.ghost) * block.dx));
			error1 = error1 + abs(fluids[index].convar[0] - primvar0);
			error2 = error2 + pow(fluids[index].convar[0] - primvar0, 2);
			error3 = (error3 > abs(fluids[index].convar[0] - primvar0)) ? error3 : abs(fluids[index].convar[0] - primvar0);
		}
		error1 /= block.nodex;
		error2 = sqrt(error2 / block.nodex);

		error[0] = error1; error[1] = error2; error[2] = error3;

		cout << scientific << "L1 error=" << error1 << " L2 error=" << error2 << " Linf error=" << error3 << endl;

	}


}

void accuracy_sinwave_2d()
{
	int mesh_set = 3;
	int mesh_number_start = 10;
	double length = 2.0;
	double CFL = 0.5;
	double dt_ratio = 1.0;

	double mesh_size_start = length / mesh_number_start;
	int* mesh_number = new int[mesh_set];
	double* mesh_size = new double[mesh_set];
	double** error = new double* [mesh_set];

	for (int i = 0; i < mesh_set; ++i)
	{
		error[i] = new double[3];
		mesh_size[i] = mesh_size_start / pow(2, i);
		mesh_number[i] = mesh_number_start * pow(2, i);
		sinwave_2d(CFL, dt_ratio, mesh_number[i], error[i]);
	}

	output_error_form(CFL, dt_ratio, mesh_set, mesh_number, error);
}

void sinwave_2d(double& CFL, double& dt_ratio, int& mesh_number, double* error)
{

	Runtime runtime;
	runtime.start_initial = clock();
	Block2d block; // Class, geometry variables
	block.uniform = false;
	block.nodex = mesh_number;
	block.nodey = mesh_number;
	block.ghost = 3; // 5th-order reconstruction, should 3 ghost cell

	double tstop = 2;
	block.CFL = CFL;


	K = 3; // 1d K=4, 2d K=3, 3d K=2;
	Gamma = 1.4; // diatomic gas, r=1.4

	//this part should rewritten ad gks2dsolver blabla
	gks2dsolver = gks2nd_2d; // Emumeration, choose solver type
	// 2nd means spatical 2nd-order, the traditional GKS; 2d means two dimensions

	tau_type = Euler; // Emumeration, choose collision time tau type
	// for accuracy test case, the tau type should be Euler and more accurately ZERO
	// for inviscid case, use Euler type tau; for viscous case, use NS type tau

	c1_euler = 0.0; // first coefficient in Euler type tau calculation
	c2_euler = 0.0; // second coefficient in Euler type tau calculation
	// when Euler type, tau is always zero;
	// if the c1 c2 were given, (c1, c2, should be given together), the tau_num is determined by c1 c2, (c1*dt + c2*deltaP)
	// if c1, c2 are both zero, (or not given value either), tau_num is also zero, which might be useful for accuracy test.
	// when NS type, tau is always physical tau, relating to Mu or Nu;
	// besides, tau_num would be determined by tau and c2 part, (tau + c2*deltaP),
	// so, c2 should be given (c1 is useless now), and Mu or Nu should be given ONLY ONE;
	// for NS type, specially, if Smooth is true, tau_num is determined ONLY by tau, (tau_num = tau), which means c2 is zero

	Mu = 0.0;
	Pr = 0.73;
	R_gas = 1;

	//prepare the boundary condtion function
	Fluid2d* bcvalue = new Fluid2d[4]; // Class, primitive variables only for boundary
	BoundaryCondition2dnew leftboundary(0); // 声明 边界条件的 函数指针
	BoundaryCondition2dnew rightboundary(0); // 声明 边界条件的 函数指针
	BoundaryCondition2dnew downboundary(0); // 声明 边界条件的 函数指针
	BoundaryCondition2dnew upboundary(0); // 声明 边界条件的 函数指针

	leftboundary = periodic_boundary_left;  // 给函数指针 赋值
	rightboundary = periodic_boundary_right;  // 给函数指针 赋值
	downboundary = periodic_boundary_down;  // 给函数指针 赋值
	upboundary = periodic_boundary_up;  // 给函数指针 赋值

	//prepare the reconstruction
	gausspoint = 2; // fifth-order or sixth-order use THREE gauss points
	// WENO5 has the function relating to arbitrary gausspoints
	// WENO5_AO supports 2 gausspoint now, so fourth-order at most for spacial reconstruction (enough for two step fourth-order GKS)
	SetGuassPoint(); // Function, set Gauss points coordinates and weight factor

	reconstruction_variable = conservative; // Emumeration, choose the variables used for reconstruction type
	wenotype = linear; // Emumeration, choose reconstruction type

	cellreconstruction_2D_normal = WENO5_normal;  // reconstruction in normal directon
	cellreconstruction_2D_tangent = WENO5_tangent;  // reconstruction in tangential directon
	g0reconstruction_2D_normal = Center_5th_normal;  // reconstruction for g0 in normal directon
	g0reconstruction_2D_tangent = Center_5th_tangent;  // reconstruction for g0 in tangential directon


	//prepare the flux function
	flux_function_2d = GKS2D_smooth; // 给函数指针 赋值 flux calculation type
	// solver is GKS
	// GKS2D means full solver (used for shock), GKS2D_smooth means smooth solver (used for non_shock)
	// for accuracy test case, tau is ZERO, so GKS2D is equal to GKS2D_smooth for results, but the later is faster
	// for inviscid flow, Euler type tau equals artifical viscosity;
	// for viscous flow, NS type tau equals the modified collision time;
	//end

	//prepare time marching stratedgy
	//time coe list must be 2d, end by _2D
	timecoe_list_2d = S2O4_2D; // 两步四阶
	Initial_stages(block);


	// allocate memory for 2-D fluid field
	// in a standard finite element method, we have
	// first the cell average value, N*N
	block.nx = block.nodex + 2 * block.ghost;
	block.ny = block.nodey + 2 * block.ghost;

	Fluid2d* fluids = Setfluid(block); // Function, input a class (geometry), output the pointer of one class (conservative variables)
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
	SetUniformMesh(block, fluids, xinterfaces, yinterfaces, xfluxes, yfluxes); // Function, set mesh

	//end

	ICfor_sinwave_2d(fluids, block);

	runtime.finish_initial = clock();
	block.t = 0; //the current simulation time
	block.step = 0; //the current step

	int inputstep = 1;//input a certain step

	while (block.t < tstop)
	{

		if (block.step % inputstep == 0)
		{
			cout << "pls cin interation step, if input is 0, then the program will exit " << endl;
			cin >> inputstep;
			if (inputstep == 0)
			{
				output2d_binary(fluids, block); // 二进制输出，用单元中心点代替单元
				break;
			}
		}
		if (runtime.start_compute == 0.0)
		{
			runtime.start_compute = clock();
			cout << "runtime-start " << endl;
		}

		CopyFluid_new_to_old(fluids, block); // Function, copy variables
		//determine the cfl condtion
		block.dt = Get_CFL(block, fluids, tstop); // Function, get real CFL number

		if ((block.t + block.dt - tstop) > 0)
		{
			block.dt = tstop - block.t + 1e-16;
		}

		for (int i = 0; i < block.stages; i++)
		{
			//after determine the cfl condition, let's implement boundary condtion
			leftboundary(fluids, block, bcvalue[0]);
			rightboundary(fluids, block, bcvalue[1]);
			downboundary(fluids, block, bcvalue[2]);
			upboundary(fluids, block, bcvalue[3]);
			//cout << "after bc " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			Convar_to_Primvar(fluids, block); // Function
			//then is reconstruction part, which we separate the left or right reconstrction
			//and the center reconstruction
			Reconstruction_within_cell(xinterfaces, yinterfaces, fluids, block); // Function
			//cout << "after cell recon " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			Reconstruction_forg0(xinterfaces, yinterfaces, fluids, block); // Function
			//cout << "after g0 recon " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			//then is solver part
			Calculate_flux(xfluxes, yfluxes, xinterfaces, yinterfaces, block, i); // Function
			//cout << "after flux calcu " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

			//then is update flux part
			Update(fluids, xfluxes, yfluxes, block, i); // Function
			//cout << "after updating " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;

		}

		block.step++;
		block.t = block.t + block.dt;

		if (block.step % 1000 == 0)
		{
			cout << "step 1000 time is " << (double)(clock() - runtime.start_compute) / CLOCKS_PER_SEC << endl;
		}
		if ((block.t - tstop) >= 0)
		{
			output2d_binary(fluids, block);
		}

	}

	error_for_sinwave_2d(fluids, block, tstop, error); // Function, get error
	runtime.finish_compute = clock();

	cout << "the total run time is " << (double)(runtime.finish_compute - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "initializing time is" << (double)(runtime.finish_initial - runtime.start_initial) / CLOCKS_PER_SEC << " second !" << endl;
	cout << "the pure computational time is" << (double)(runtime.finish_compute - runtime.start_compute) / CLOCKS_PER_SEC << " second !" << endl;

}

void ICfor_sinwave_2d(Fluid2d* fluids, Block2d block)
{
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			double pi = 3.14159265358979323846;
			int index = i * (block.ny) + j;
			double xleft = (i - block.ghost) * block.dx;
			double xright = (i + 1 - block.ghost) * block.dx;
			double yleft = (j - block.ghost) * block.dy;
			double yright = (j + 1 - block.ghost) * block.dy;

			//case in two dimensional x-y-plane
			double k1 = sin(pi * (xright + yright));
			double k2 = sin(pi * (xright + yleft));
			double k3 = sin(pi * (xleft + yright));
			double k4 = sin(pi * (xleft + yleft));
			fluids[index].primvar[0] = 1.0 - 0.2 / pi / pi / block.dx / block.dy * ((k1 - k2) - (k3 - k4));

			fluids[index].primvar[1] = 1;
			fluids[index].primvar[2] = 1;
			fluids[index].primvar[3] = 1;
		}
	}
#pragma omp parallel for
	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			int index = i * (block.ny) + j;
			Primvar_to_convar_2D(fluids[index].convar, fluids[index].primvar);
		}
	}
}

void error_for_sinwave_2d(Fluid2d* fluids, Block2d block, double tstop, double* error)
{
	if (abs(block.t - tstop) <= (1e-10))
	{
		cout << "Accuracy-residual-file-output" << endl;
		double error1 = 0;
		double error2 = 0;
		double error3 = 0;
		for (int i = block.ghost; i < block.nx - block.ghost; i++)
		{
			for (int j = block.ghost; j < block.ny - block.ghost; j++)
			{
				double pi = 3.14159265358979323846;
				int index = i * block.ny + j;
				double xleft = (i - block.ghost) * block.dx;
				double xright = (i + 1 - block.ghost) * block.dx;
				double yleft = (j - block.ghost) * block.dy;
				double yright = (j + 1 - block.ghost) * block.dy;

				//case in two dimensional x-y-plane
				double k1 = sin(pi * (xright + yright));
				double k2 = sin(pi * (xright + yleft));
				double k3 = sin(pi * (xleft + yright));
				double k4 = sin(pi * (xleft + yleft));
				double primvar0 = 1.0 - 0.2 / pi / pi / block.dx / block.dy * ((k1 - k2) - (k3 - k4));


				error1 = error1 + abs(fluids[index].convar[0] - primvar0);
				error2 = error2 + pow(fluids[index].convar[0] - primvar0, 2);
				error3 = (error3 > abs(fluids[index].convar[0] - primvar0)) ? error3 : abs(fluids[index].convar[0] - primvar0);
			}
		}
		error1 /= (block.nodex * block.nodey);
		error2 = sqrt(error2 / (block.nodex * block.nodey));
		error[0] = error1; error[1] = error2; error[2] = error3;
	}
}
