#pragma once
#include"fvm_gks2d.h"
#include"time.h"
#include<string.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#if defined(_WIN32)
#include<direct.h>
#else 
#include<sys/stat.h>
//#include<sys/types.h>
#endif
struct Runtime {
	//for record wall clock time
	clock_t start_initial;
	clock_t finish_initial;
	clock_t start_compute;
	clock_t finish_compute;
	Runtime()
	{
		memset(this, 0, sizeof(Runtime));
	}
};

void make_directory_for_result();

void couthere();
void cinhere();

void output_error_form(double CFL, double dt_ratio, int mesh_set, int* mesh_number, double** error);

//tecplot 一维输出
void output1d(Fluid1d* fluids, Block1d block);
void continue_compute(Fluid2d* fluids, Block2d& block);
//tecplot 一维输出，单元中心点值
void output2d(Fluid2d* fluids, Block2d block);
//tecplot 一维输出，单元值
void output2d_center(Fluid2d* fluids, Block2d block);
//输出哪里爆掉了
void output2d_blow_up(Fluid2d* fluids, Block2d block);
void output_record(Block2d block);
void output2d_binary(Fluid2d* fluids, Block2d block);
void write_character_into_plt(char* str, ofstream& out);
void continue_compute_output(Fluid2d* fluids, Block2d block);
void continue_compute_output_param(Block2d block);


