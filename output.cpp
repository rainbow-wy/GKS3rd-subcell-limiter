#include"output.h"



void make_directory_for_result()
{
	string newfolder = "result";

#if defined(_WIN32)
	//一个简单的宏定义判断，如果是win系统 则运行以下代码
	if (_mkdir(newfolder.c_str()) != 0)
	{
		_mkdir(newfolder.c_str());
	}
	else
	{
		cout << "fail to create new folder for result file" << endl;
	}
	//end
#else 
	//如果是linux系统 则运行以下代码
	mkdir(newfolder.c_str(), 0777);
	if (mkdir(newfolder.c_str(), 0777) != 0)
	{
		mkdir(newfolder.c_str(), 0777);
	}
	else
	{
		cout << "fail to create new folder for result file" << endl;
	}
	//end
#endif

}

void couthere()
{
	// to test that the program can run up to here
	cout << "here" << endl;
}
void cinhere()
{
	// to test taht the program can run up here
	//and only show once
	double a;
	cout << "the program run here" << endl;
	cin >> a;
}

void output_error_form(double CFL, double dt_ratio, int mesh_set, int* mesh_number, double** error)
{
	cout << "Accuracy-residual-file-output" << endl;

	ofstream error_out("result/error.txt");
	error_out << "the CFL number is " << CFL << endl;
	error_out << "the dt ratio over dx is " << dt_ratio << endl;


	double** order = new double* [mesh_set - 1];
	for (int i = 0; i < mesh_set - 1; i++)
	{
		order[i] = new double[3];
		for (int j = 0; j < 3; j++)
		{
			order[i][j] = log(error[i][j] / error[i + 1][j]) / log(2);
		}
	}


	for (int i = 0; i < mesh_set; i++)
	{
		if (i == 0)
		{
			error_out << "1/" << mesh_number[i]
				<< scientific << setprecision(6)
				<< " & " << error[i][0] << " & ~"
				<< " & " << error[i][1] << " & ~"
				<< " & " << error[i][2] << " & ~"
				<< " \\\\ " << endl;
		}
		else
		{
			error_out << "1/" << mesh_number[i]
				<< " & " << scientific << setprecision(6) << error[i][0];
			error_out << " & " << fixed << setprecision(2) << order[i - 1][0];
			error_out << " & " << scientific << setprecision(6) << error[i][1];
			error_out << " & " << fixed << setprecision(2) << order[i - 1][1];
			error_out << " & " << scientific << setprecision(6) << error[i][2];
			error_out << " & " << fixed << setprecision(2) << order[i - 1][2];
			error_out << " \\\\ " << endl;
		}

	}
	error_out.close();
}


void output1d(Fluid1d* fluids, Block1d block)
{
	stringstream name;
	name << "C:/code/gks2d-str/code/result/Result1D-" << setfill('0') << setw(5) << block.step << ".plt"<< endl ;
	string s;
	name >> s;
	cout << "Output file: " << s << endl;  // 在 `name >> s;` 后添加
	ofstream out(s);
	out << "variables = x,density,u,pressure,temperature,entropy,Ma" << endl;
	out << "zone i = " << block.nodex << ", F=POINT" << endl;

	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);
	//output the data
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		out << fluids[i].cx << " "
			<< setprecision(10)
			<< fluids[i].primvar[0] << " "
			<< fluids[i].primvar[1] << " "
			<< fluids[i].primvar[2] << " "
			<< Temperature(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
			<< entropy(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
			<< sqrt(pow(fluids[i].primvar[1], 2)) / Soundspeed(fluids[i].primvar[0], fluids[i].primvar[2]) << " "
			<< endl;
	}
	out.close();
}

void continue_compute(Fluid2d* fluids, Block2d& block)
{
	fstream in("result/continue-in.txt", ios::in);
	if (!in.is_open())
	{
		cout << "cannot find continue.txt" << endl;
		cout << "a new case will start" << endl;
		return;
	}
	int in_mesh[2];
	in >> in_mesh[0] >> in_mesh[1];

	if (in_mesh[0] == block.nodex && in_mesh[1] == block.nodey)
	{

		fstream dat("result/continue-in.plt", ios::in);
		if (!dat.is_open())
		{
			cout << "cannot find continue.plt, which contains fluid data" << endl;
			cout << "a new case will start" << endl;
			return;
		}
		string s;
		getline(dat, s);
		dat.ignore(256, '\n');

		if (block.ghost <= 0 || block.ghost >= 5)
		{
			cout << "a wrong ghost set up initializing" << endl;
		}

		//cell avg information
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int i = block.ghost; i < block.ghost + block.nodex; i++)
			{
				int index = i * block.ny + j;
				dat >> fluids[index].coordx
					>> fluids[index].coordy;

				for (int var = 0; var < 4; var++)
				{
					dat >> fluids[index].convar[var];
				}



			}

		}


		in >> block.t;
		cout << "successfully read in continue-compute-file in t= " << block.t << endl;

		Convar_to_Primvar(fluids, block);
		cout << "change convar to primvar" << endl;
	}
	else
	{
		cout << "mesh size is not same" << endl;
		cout << "a new case will start" << endl;
		return;
	}
	in.close();
	return;
}
void output2d(Fluid2d* fluids, Block2d block)
{
	//for create a output file, pls include the follow standard lib
	//#include <fstream>
	//#include <sstream>
	//#include<iomanip>
	output_record(block);
	stringstream name;
	name << "C:/code/gks2d-str/code/result/Result2D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "# CFL number is " << block.CFL << endl;
	out << "# tau type (0 refer euler, 1 refer NS) is " << tau_type << endl;
	if (tau_type == Euler)
	{
		out << "# c1_euler " << c1_euler << " c2_euler " << c2_euler << endl;
	}
	if (tau_type == NS)
	{
		out << "# Mu " << Mu << " Nu " << Nu << endl;
	}
	out << "# time marching strategy is " << block.stages << " stage method" << endl;
	out << "variables = x,y,density,u,v,pressure,temperature,entropy,Ma" << endl;
	out << "zone i = " << block.nodex << " " << ", j = " << block.nodey << ", F=POINT" << endl;

	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);
	//output the data
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{

			out << fluids[i * block.ny + j].coordx << " "
				<< fluids[i * block.ny + j].coordy << " "
				<< setprecision(10)
				<< fluids[i * block.ny + j].primvar[0] << " "
				<< fluids[i * block.ny + j].primvar[1] << " "
				<< fluids[i * block.ny + j].primvar[2] << " "
				<< fluids[i * block.ny + j].primvar[3] << " "
				<< Temperature(fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[3]) << " "
				<< entropy(fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[3]) << " "
				<< sqrt(pow(fluids[i * block.ny + j].primvar[1], 2) + pow(fluids[i * block.ny + j].primvar[2], 2)) /
				Soundspeed(fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[3]) << " "
				<< endl;
		}
	}

	//fclose(fp);//关闭文件
	out.close();

}
void output2d_center(Fluid2d* fluids, Block2d block)
{
	output_record(block);
	stringstream name;
	name << "result/result2D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);


	out << "# CFL number is " << block.CFL << endl;
	out << "# tau type (0 refer euler, 1 refer NS) is " << tau_type << endl;
	if (tau_type == Euler)
	{
		out << "# c1_euler " << c1_euler << " c2_euler " << c2_euler << endl;
	}
	if (tau_type == NS)
	{
		out << "# Mu " << Mu << " Nu " << Nu << endl;
	}
	out << "# time marching strategy is " << block.stages << " stage method" << endl;
	out << "variables = x,y,density,u,v,pressure,temperature,entropy,Ma" << endl;
	out << "zone i = " << block.nodex + 1 << " " << ", j = " << block.nodey + 1 << ", DATAPACKING=BLOCK, VARLOCATION=([3-9]=CELLCENTERED)" << endl;

	//output locations
	for (int j = block.ghost; j < block.ghost + block.nodey + 1; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex + 1; i++)
		{
			out << fluids[i * block.ny + j].node[0] << " ";
		}
		out << endl;
	}

	for (int j = block.ghost; j < block.ghost + block.nodey + 1; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex + 1; i++)
		{
			out << fluids[i * block.ny + j].node[1] << " ";
		}
		out << endl;
	}


	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);
	//output the data

	//output four primary variables
	for (int var = 0; var < 4; var++)
	{
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int i = block.ghost; i < block.ghost + block.nodex; i++)
			{
				out << setprecision(10)
					<< fluids[i * block.ny + j].primvar[var] << " ";
			}
			out << endl;
		}
	}
	//output temperature
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			out << setprecision(10)
				<< Temperature(fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[3]) << " ";
		}
		out << endl;
	}
	//output entropy
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			out << setprecision(10)
				<< entropy(fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[3]) << " ";
		}
		out << endl;
	}

	//output Ma number
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			out << setprecision(10)
				<< sqrt(pow(fluids[i * block.ny + j].primvar[1], 2))
				/ Soundspeed(fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[3]) << " ";
		}
		out << endl;
	}

	//fclose(fp);//关闭文件
	out.close();
}

void output2d_blow_up(Fluid2d* fluids, Block2d block)
{
	output_record(block);
	stringstream name;
	name << "result/BlowUp-2D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "# CFL number is " << block.CFL << endl;
	out << "# tau type (0 refer euler, 1 refer NS) is " << tau_type << endl;
	if (tau_type == Euler)
	{
		out << "# c1_euler " << c1_euler << " c2_euler " << c2_euler << endl;
	}
	if (tau_type == NS)
	{
		out << "# Mu " << Mu << " Nu " << Nu << endl;
	}
	out << "# time marching strategy is " << block.stages << " stage method" << endl;
	out << "variables = x,y,blow_point,density,u,v,pressure" << endl;
	out << "zone i = " << block.nodex + 1 << " " << ", j = " << block.nodey + 1 << ", DATAPACKING=BLOCK, VARLOCATION=([3-7]=CELLCENTERED)" << endl;

	//output locations
	for (int j = block.ghost; j < block.ghost + block.nodey + 1; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex + 1; i++)
		{
			out << fluids[i * block.ny + j].node[0] << " ";
		}
		out << endl;
	}

	for (int j = block.ghost; j < block.ghost + block.nodey + 1; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex + 1; i++)
		{
			out << fluids[i * block.ny + j].node[1] << " ";
		}
		out << endl;
	}


	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);
	//output the data

	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			if (negative_density_or_pressure(fluids[i * block.ny + j].primvar))
			{
				out << 1.0 << " ";
			}
			else
			{
				out << 0.0 << " ";
			}
		}
		out << endl;
	}
	//output four primary variables
	for (int var = 0; var < 4; var++)
	{
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int i = block.ghost; i < block.ghost + block.nodex; i++)
			{
				if (negative_density_or_pressure(fluids[i * block.ny + j].primvar))
				{
					out << -1.0 << " ";
				}
				else
				{
					out << setprecision(10)
						<< fluids[i * block.ny + j].primvar[var] << " ";
				}
			}
			out << endl;
		}
	}


	//fclose(fp);//关闭文件
	out.close();
}
void output_record(Block2d block)
{
	//输出一下二维的模拟工况，为了方便记录之用
	stringstream name;
	name << "result/Record2D-" << setfill('0') << setw(5) << block.step << ".txt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "# This is a two dimensional result" << endl;
	if (block.uniform == true)
	{
		out << "# uniform mesh is used for calculation" << endl;
		out << "# left region coordinate " << endl;
		out << block.left << endl;
		out << "# right region coordinate " << endl;
		out << block.right << endl;
		out << "# down region coordinate " << endl;
		out << block.down << endl;
		out << "# up region coordinate " << endl;
		out << block.up << endl;
	}
	else
	{
		out << "# non-uniform mesh is used for calculation" << endl;
	}
	out << "# ghost cell number is " << endl;
	out << block.ghost << endl;

	out << "# Mesh size nodex and nodey" << endl;
	out << block.nodex << " " << block.nodey << endl;
	out << "# CFL number is " << endl;
	out << block.CFL << endl;
	out << "# simulation time is " << endl;
	out << block.t << endl;
	out << "# simulation step is " << endl;
	out << block.step << endl;
	out << "# tau type (0 refer euler, 1 refer NS) is " << endl;
	out << tau_type << endl;

	if (tau_type == Euler)
	{
		out << "# c1_euler " << " c2_euler " << endl;
		out << c1_euler << " " << c2_euler << endl;
	}

	if (tau_type == NS)
	{
		out << "# Mu " << " Nu " << endl;
		out << Mu << " " << Nu << endl;
	}
	out << "# the number of gaussian point it use: " << endl;
	out << gausspoint << endl;

	out << "# normal_reconstruction_for_left_right_states is ";
	if (cellreconstruction_2D_normal == First_order_normal)
	{
		out << "FirstOrder_splitting " << endl;
	}
	else if (cellreconstruction_2D_normal == Vanleer_normal)
	{
		out << "Vanleer_splitting " << endl;
	}
	else
	{
		out << "not_specified" << endl;
	}


	out << "# normal_reconstruction_for_g0_states is ";
	if (g0reconstruction_2D_normal == Center_collision_normal)
	{
		out << "Center_collision_normal " << endl;
	}
	else if (g0reconstruction_2D_normal == Center_3rd_normal)
	{
		out << "Center_3rd_normal " << endl;
	}
	else
	{
		out << "not_specified" << endl;
	}

	out << "# tangential_reconstruction_for_left_right_states is ";
	if (cellreconstruction_2D_tangent == Vanleer_tangent)
	{
		out << "Vanleer_tangent " << endl;
	}
	else if (cellreconstruction_2D_tangent == First_order_tangent)
	{
		out << "First_order_tangent " << endl;
	}
	else
	{
		out << "not_specified" << endl;
	}

	out << "# g0reconstruction_2D_multi is ";
	if (g0reconstruction_2D_tangent == Center_all_collision_multi)
	{
		out << "Center_all_collision_tangent " << endl;
	}
	else if (g0reconstruction_2D_tangent == Center_3rd_tangent)
	{
		out << "Center_all_collision_tangent " << endl;
	}
	else
	{
		out << "not_specified" << endl;
	}

	out << "# time marching strategy is ";
	if (timecoe_list_2d == S1O1_2D)
	{
		out << "S1O1 " << endl;
	}
	if (timecoe_list_2d == S2O4_2D)
	{
		out << "S2O4 " << endl;
	}
	if (timecoe_list_2d == RK4_2D)
	{
		out << "RK4 " << endl;
	}

	out << "# the flux function is ";
	if (flux_function_2d == GKS2D_smooth)
	{
		out << "gks2nd_smooth" << endl;
	}
	if (flux_function_2d == GKS2D)
	{
		if (gks2dsolver == kfvs1st_2d)
		{
			out << "kfvs1st " << endl;
		}
		if (gks2dsolver == kfvs2nd_2d)
		{
			out << "kfvs2nd " << endl;
		}
		if (gks2dsolver == gks1st_2d)
		{
			out << "gks1st " << endl;
		}
		if (gks2dsolver == gks2nd_2d)
		{
			out << "gks2nd " << endl;
		}
	}

	out << "# the limitation for current scheme is " << endl;
	out.close();
}
void output2d_binary(Fluid2d* fluids, Block2d block)
{
	cout << "Output binary result..." << endl;
	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);

	output_record(block);
	stringstream name;
	name << "result/Result2D-" << setfill('0') << setw(5) << block.step << ".plt" << endl;
	string s;
	name >> s;
	ofstream out;
	out.open(s, ios::out | ios::binary);

	int IMax = block.nodex;
	int JMax = block.nodey;
	int KMax = 1;
	int NumVar = 7;

	char Title[] = "structure2d";


	char Varname1[] = "X";
	char Varname2[] = "Y";
	char Varname3[] = "Z";
	char Varname4[] = "density";
	char Varname5[] = "u";
	char Varname6[] = "v";
	char Varname7[] = "pressure";

	char Zonename1[] = "Zone 001";

	float ZONEMARKER = 299.0;

	float EOHMARKER = 357.0;

	//==============Header Secontion =================//

	//------1.1 Magic number, Version number

	char MagicNumber[] = "#!TDV101";

	// the magic number must be output with 8 lengths without NULL end
	out.write(MagicNumber, 8);
	//---- - 1.2.Integer value of 1.----------------------------------------------------------

	int IntegerValue = 1;
	out.write((char*)&IntegerValue, sizeof(IntegerValue));

	//---- - 1.3.Title and variable names.------------------------------------------------ -

	//---- - 1.3.1.The TITLE.
	write_character_into_plt(Title, out);
	//---- - 1.3.2 Number of variables(NumVar) in the c_strfile.



	out.write((char*)&NumVar, sizeof(NumVar));


	//------1.3.3 Variable names.N = L[1] + L[2] + ....L[NumVar]

	write_character_into_plt(Varname1, out);
	write_character_into_plt(Varname2, out);
	write_character_into_plt(Varname3, out);
	write_character_into_plt(Varname4, out);
	write_character_into_plt(Varname5, out);
	write_character_into_plt(Varname6, out);
	write_character_into_plt(Varname7, out);

	//---- - 1.4.Zones------------------------------------------------------------------ -

	//--------Zone marker.Value = 299.0

	out.write((char*)&ZONEMARKER, sizeof(ZONEMARKER));

	//--------Zone name.
	write_character_into_plt(Zonename1, out);
	//--------Zone color

	int ZoneColor = -1;

	out.write((char*)&ZoneColor, sizeof(ZoneColor));

	//--------ZoneType

	int ZoneType = 0;

	out.write((char*)&ZoneType, sizeof(ZoneType));

	//--------DaraPacking 0=Block, 1=Point

	int DaraPacking = 1;

	out.write((char*)&DaraPacking, sizeof(DaraPacking));
	//--------Specify Var Location. 0 = Don't specify, all c_str is located at the nodes. 1 = Specify

	int SpecifyVarLocation = 0;

	out.write((char*)&SpecifyVarLocation, sizeof(SpecifyVarLocation));
	//--------Number of user defined face neighbor connections(value >= 0)

	int NumOfNeighbor = 0;

	out.write((char*)&NumOfNeighbor, sizeof(NumOfNeighbor));

	//-------- - IMax, JMax, KMax

	out.write((char*)&IMax, sizeof(IMax));
	out.write((char*)&JMax, sizeof(JMax));
	out.write((char*)&KMax, sizeof(KMax));

	//----------// -1 = Auxiliary name / value pair to follow 0 = No more Auxiliar name / value pairs.

	int AuxiliaryName = 0;

	out.write((char*)&AuxiliaryName, sizeof(AuxiliaryName));

	//----I HEADER OVER--------------------------------------------------------------------------------------------

	//=============================Geometries section=======================

	//=============================Text section======================

	// EOHMARKER, value = 357.0

	out.write((char*)&EOHMARKER, sizeof(EOHMARKER));

	//================II.Data section===============//

	//---- - 2.1 zone---------------------------------------------------------------------- -

	out.write((char*)&ZONEMARKER, sizeof(ZONEMARKER));

	//--------variable c_str format, 1 = Float, 2 = Double, 3 = LongInt, 4 = ShortInt, 5 = Byte, 6 = Bit

	int fomat[20];
	for (int ifomat = 0; ifomat < 20; ifomat++)
	{
		fomat[ifomat] = 2;
	}
	for (int ifomat = 0; ifomat < NumVar; ifomat++)
	{
		out.write((char*)&fomat[ifomat], sizeof(fomat[ifomat]));
	}


	//--------Has variable sharing 0 = no, 1 = yes.

	int HasVarSharing = 0;

	out.write((char*)&HasVarSharing, sizeof(HasVarSharing));

	//----------Zone number to share connectivity list with(-1 = no sharing).

	int ZoneNumToShareConnectivity = -1;

	out.write((char*)&ZoneNumToShareConnectivity, sizeof(ZoneNumToShareConnectivity));
	//----------Zone c_str.Each variable is in c_str format asspecified above.
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			double coordz = 0.0;
			out.write((char*)&fluids[i * block.ny + j].coordx, sizeof(double));
			out.write((char*)&fluids[i * block.ny + j].coordy, sizeof(double));
			out.write((char*)&coordz, sizeof(double));
			out.write((char*)&fluids[i * block.ny + j].primvar[0], sizeof(double));
			out.write((char*)&fluids[i * block.ny + j].primvar[1], sizeof(double));
			out.write((char*)&fluids[i * block.ny + j].primvar[2], sizeof(double));
			out.write((char*)&fluids[i * block.ny + j].primvar[3], sizeof(double));
		}

	}
	out.close();

}
void write_character_into_plt(char* str, ofstream& out)
{
	int value = 0;
	while ((*str) != '\0')
	{
		value = (int)*str;
		out.write((char*)&value, sizeof(value));
		str++;
	}

	char null_char[] = "";
	value = (int)*null_char;
	out.write((char*)&value, sizeof(value));
}
void continue_compute_output(Fluid2d* fluids, Block2d block)
{
	//for create a output file, pls include the follow standard lib
	//#include <fstream>
	//#include <sstream>
	//#include<iomanip>
	continue_compute_output_param(block);
	stringstream name;
	name << "result/continue-out" << ".plt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << "variables = x,y,convar0,convar1,convar2,convar3";
	out << endl;

	out << "zone i = " << block.nodex << " " << ", j = " << block.nodey << ", F=POINT" << endl;

	//before output data, converting conservative varabile to primitive variable
	Convar_to_primvar(fluids, block);

	//output the data
	for (int j = block.ghost; j < block.ghost + block.nodey; j++)
	{
		for (int i = block.ghost; i < block.ghost + block.nodex; i++)
		{
			out << fluids[i * block.ny + j].coordx << " "
				<< fluids[i * block.ny + j].coordy << " "
				<< setprecision(10);
			for (int var = 0; var < 4; var++)
			{
				out << fluids[i * block.ny + j].convar[var] << " ";
			}

			out << endl;
		}
	}

	out.close();
}
void continue_compute_output_param(Block2d block)
{
	stringstream name;
	name << "result/continue-out" << ".txt" << endl;
	string s;
	name >> s;
	ofstream out(s);
	out << block.nodex << " " << block.nodey << endl; //mesh point
	out << setprecision(14);
	out << block.t << endl;//simulation time
}

