#include"fvm_gks2d.h"

//默认参数
GKS1d_type gks1dsolver = nothing; //initialization
Reconstruction_within_Cell cellreconstruction = Vanleer; //initialization
Reconstruction_forG0 g0reconstruction = Center_collision; //initialization
//G0_construct_type g0type = collisionn; //initialization
Flux_function flux_function = GKS; //initialization
Reconstruction_variable reconstruction_variable = conservative; //initialization
WENOtype wenotype = wenojs; //initialization
TimeMarchingCoefficient timecoe_list = S2O4; //initialization
bool is_reduce_order_warning = false; //initialization
bool is_show_1d_timestep = true; //initialization


void Convar_to_primvar(Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		Convar_to_primvar_1D(fluids[i].primvar, fluids[i].convar);
	}
}




void free_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	//for this case, the input variable can also use Fluid1d fluids
	//the call is the same with the current method
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		fluids[i].convar[0] = fluids[i + 1].convar[0];
		fluids[i].convar[1] = fluids[i + 1].convar[1];
		fluids[i].convar[2] = fluids[i + 1].convar[2];
	}
}
void free_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		fluids[i].convar[0] = fluids[i - 1].convar[0];
		fluids[i].convar[1] = fluids[i - 1].convar[1];
		fluids[i].convar[2] = fluids[i - 1].convar[2];
	}
}
void reflection_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		int ref = 2 * block.ghost - 1 - i;
		fluids[i].convar[0] = fluids[ref].convar[0];
		fluids[i].convar[1] = -fluids[ref].convar[1];
		fluids[i].convar[2] = fluids[ref].convar[2];
	}
}
void reflection_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		int ref = 2 * (block.nx - block.ghost) - 1 - i;
		fluids[i].convar[0] = fluids[ref].convar[0];
		fluids[i].convar[1] = -fluids[ref].convar[1];
		fluids[i].convar[2] = fluids[ref].convar[2];
	}
}
void periodic_boundary_left(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nodex_begin - 1; i >= 0; i--)
	{
		int ref = i + block.nodex;
		Copy_Array(fluids[i].convar, fluids[ref].convar, 3);
	}
}
void periodic_boundary_right(Fluid1d* fluids, Block1d block, Fluid1d bcvalue)
{
	for (int i = block.nodex_end + 1; i < block.nx; i++)
	{
		int ref = i - block.nodex;
		Copy_Array(fluids[i].convar, fluids[ref].convar, 3);
	}
}


void Check_Order_Reduce(Point1d& left, Point1d& right, Fluid1d& fluid)
{
	bool order_reduce[2];
	Check_Order_Reduce_by_Lambda_1D(order_reduce[0], right.convar);
	Check_Order_Reduce_by_Lambda_1D(order_reduce[1], left.convar);
	//if lambda <0, then reduce to the first order
	if (order_reduce[0] == true || order_reduce[1] == true)
	{
		if (is_reduce_order_warning == true)
			cout << "reconstruction reduce order in location  " << fluid.cx << endl;
		for (int m = 0; m < 3; m++)
		{
			left.convar[m] = fluid.convar[m];
			right.convar[m] = fluid.convar[m];
			left.der1[m] = 0.0;
			right.der1[m] = 0.0;
			left.der2[m] = 0.0;
			right.der2[m] = 0.0;
		}
	}
}
void Check_Order_Reduce_by_Lambda_1D(bool& order_reduce, double* convar)
{
	order_reduce = false;
	double lambda;
	lambda = Lambda(convar[0], convar[1] / convar[0], convar[2]);
	//if lambda <0, then reduce to the first order
	if (lambda <= 0.0 || (lambda == lambda) == false)
	{
		order_reduce = true;
	}
}


void Reconstruction_within_cell(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
		cellreconstruction(interfaces[i].right, interfaces[i + 1].left, &fluids[i]);
	}
}
void Vanleer(Point1d& left, Point1d& right, Fluid1d* fluids)
{
	//Note: function by vanleer reconstruction
	//the final updated variables is conservative variables
	Fluid1d wn1 = fluids[-1];
	Fluid1d w = fluids[0];
	Fluid1d wp1 = fluids[1];
	double splus[3], sminus[3];

	for (int i = 0; i < 3; i++)
	{
		splus[i] = (wp1.convar[i] - w.convar[i]) / ((wp1.dx + w.dx) / 2.0);
		sminus[i] = (w.convar[i] - wn1.convar[i]) / ((wn1.dx + w.dx) / 2.0);

		if ((splus[i] * sminus[i]) > 0)
		{
			left.der1[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
			right.der1[i] = left.der1[i];
		}
		else
		{
			left.der1[i] = 0.0;
			right.der1[i] = 0.0;
		}
		left.der2[i] = 0.0;
		right.der2[i] = 0.0;
		left.convar[i] = w.convar[i] - 0.5 * w.dx * left.der1[i];
		right.convar[i] = w.convar[i] + 0.5 * w.dx * right.der1[i];
	}

	Check_Order_Reduce(left, right, fluids[0]);
}

void weno_3rd_left(double& var, double& der1, double& der2,
	double wn1, double w0, double wp1, double h)
{
	//-- - parameter of WENO-- -
	double beta[2], d[2], ww[2], alpha[2];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[2], px[2], pxx[2], tempvar;
	double sum_alpha;

	//two small stencil
	d[0] = 2.0 / 3.0;
	d[1] = 1.0 / 3.0;


	if (wenotype == linear)
	{
		for (int k = 0; k < 2; k++)
		{
			ww[k] = d[k];
		}
	}
	else
	{

		beta[0] = (w0 - wn1) * (w0 - wn1);
		beta[1] = (w0 - wp1) * (w0 - wp1);

		double tau3 = (abs(beta[0] - beta[1]));


		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 2; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}

		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 2; i++)
			{
				double global_div = pow(tau3 / (beta[i] + epsilonW), 2);
				alpha[i] = d[i] * (1 + global_div);
				sum_alpha += alpha[i];
			}
		}


		for (int k = 0; k < 2; k++)
		{
			ww[k] = alpha[k] / sum_alpha;
		}
	}

	//-- - candidate polynomial-- -
	p[0] = (w0 + wn1) / 2.0;
	p[1] = (3 * w0 - wp1) / 2.0;

	px[0] = (w0 - wn1) / h;
	px[1] = (wp1 - w0) / h;

	pxx[0] = 0.0; pxx[1] = 0.0;

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	der2 = 0.0;
	double final_weight[2];
	for (int k = 0; k < 2; k++)
	{
		final_weight[k] = ww[k];
	}

	for (int k = 0; k < 2; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
		der2 += final_weight[k] * pxx[k];
	}
}

void weno_3rd_right(double& var, double& der1, double& der2,
	double wn1, double w0, double wp1,
	double h)
{
	//-- - parameter of WENO-- -
	double beta[2], d[2], ww[2], alpha[2];
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double p[2], px[2], pxx[2], tempvar;
	double sum_alpha;

	//two small stencil
	d[0] = 1.0 / 3.0;
	d[1] = 2.0 / 3.0;


	if (wenotype == linear)
	{
		for (int k = 0; k < 3; k++)
		{
			ww[k] = d[k];
		}
	}
	else
	{
		beta[0] = (w0 - wn1) * (w0 - wn1);
		beta[1] = (w0 - wp1) * (w0 - wp1);

		double tau3 = (abs(beta[0] - beta[1]));


		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 2; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}

		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 2; i++)
			{
				double global_div = tau3 / (beta[i] + epsilonW);
				alpha[i] = d[i] * (1 + global_div * global_div);
				sum_alpha += alpha[i];
			}
		}


		for (int k = 0; k < 2; k++)
		{
			ww[k] = alpha[k] / sum_alpha;
		}
	}

	//-- - candidate polynomial-- -
	p[0] = (3 * w0 - wn1) / 2.0;
	p[1] = (w0 + wp1) / 2.0;

	px[0] = (w0 - wn1) / h;
	px[1] = (wp1 - w0) / h;

	pxx[0] = 0.0; pxx[1] = 0.0;

	////-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	der2 = 0.0;
	double final_weight[2];
	for (int k = 0; k < 2; k++)
	{
		final_weight[k] = ww[k];
	}

	for (int k = 0; k < 2; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
		der2 += final_weight[k] * pxx[k];
	}
}

void WENO3_left(double& left, double w_back, double w, double w_forward, double h)
{
	//left poistion according to the center of w.
	//-- - parameter of WENO-- -
	double beta[2], d[2], ww[2], alpha[2];
	double p = 2;
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double h0, h1, h2, tempvar;
	d[0] = 2.0 / 3.0;
	d[1] = 1.0 / 3.0;

	beta[0] = (w - w_back) * (w - w_back);
	beta[1] = (w_forward - w) * (w_forward - w);

	tempvar = abs(beta[0] - beta[1]);
	h0 = 0.0;
	for (int k = 0; k < 2; k++)
	{
		alpha[k] = d[k] * (1 + tempvar / (beta[k] + epsilonW));
		h0 = h0 + alpha[k];
	}
	for (int k = 0; k < 2; k++)
	{
		ww[k] = alpha[k] / h0;
	}
	//-- - candidate polynomial-- -
	h0 = (w_back + w) / 2.0;
	h1 = (3.0 * w - w_forward) / 2.0;
	//-- - combination-- -
	left = ww[0] * h0 + ww[1] * h1;
}

void WENO3_right(double& right, double w_back, double w, double w_forward, double h)
{
	//left poistion according to the center of w.
	//-- - parameter of WENO-- -
	double beta[2], d[2], ww[2], alpha[2];
	double p = 2;
	double epsilonW = 1e-6;
	//-- - intermediate parameter-- -
	double h0, h1, h2, tempvar;
	d[0] = 1.0 / 3.0;
	d[1] = 2.0 / 3.0;


	beta[0] = (w - w_back) * (w - w_back);
	beta[1] = (w_forward - w) * (w_forward - w);


	tempvar = abs(beta[0] - beta[1]);

	h0 = 0.0;
	for (int k = 0; k < 2; k++)
	{
		alpha[k] = d[k] * (1 + tempvar / (beta[k] + epsilonW));
		h0 = h0 + alpha[k];
	}


	for (int k = 0; k < 2; k++)
	{
		ww[k] = alpha[k] / h0;
	}

	//-- - candidate polynomial-
	h0 = (3.0 * w - w_back) / 2.0;
	h1 = (w + w_forward) / 2.0;
	//-- - combination-- -
	right = ww[0] * h0 + ww[1] * h1;
}

void WENO5_AO(Point1d& left, Point1d& right, Fluid1d* fluids)
{
	//Note: function by WENO5_AO reconstruction

	double wn2[3]; Copy_Array(wn2, fluids[-2].convar, 3);
	double  wn1[3]; Copy_Array(wn1, fluids[-1].convar, 3);
	double  w0[3];  Copy_Array(w0, fluids[0].convar, 3);
	double  wp1[3]; Copy_Array(wp1, fluids[1].convar, 3);
	double wp2[3]; Copy_Array(wp2, fluids[2].convar, 3);
	double tmp;
	//non-uniform grid was treated as uniform grid
	double  h = fluids[0].dx;

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_right(right.convar[i], right.der1[i], right.der2[i], wn2[i], wn1[i], w0[i], wp1[i], wp2[i], h);
			weno_5th_ao_left(left.convar[i], left.der1[i], left.der2[i], wn2[i], wn1[i], w0[i], wp1[i], wp2[i], h);
		}
	}

	if (reconstruction_variable == characteristic)
	{

		double ren3[3], ren2[3], ren1[3], re0[3], rep1[3], rep2[3], rep3[3];
		double var[3], der1[3], der2[3];

		double base_left[3], base_right[3];

		double wn1_primvar[3], w_primvar[3], wp1_primvar[3];
		Convar_to_primvar_1D(wn1_primvar, wn1);
		Convar_to_primvar_1D(w_primvar, w0);
		Convar_to_primvar_1D(wp1_primvar, wp1);

		for (int i = 0; i < 3; i++)
		{
			base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
			base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
		}

		Convar_to_char1D(ren2, base_left, wn2);
		Convar_to_char1D(ren1, base_left, wn1);
		Convar_to_char1D(re0, base_left, w0);
		Convar_to_char1D(rep1, base_left, wp1);
		Convar_to_char1D(rep2, base_left, wp2);

		// left_reconstruction
		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
		}
		Char_to_convar1D(left.convar, base_left, var);
		Char_to_convar1D(left.der1, base_left, der1);
		Char_to_convar1D(left.der2, base_left, der2);

		// right reconstruction

		Convar_to_char1D(ren2, base_right, wn2);
		Convar_to_char1D(ren1, base_right, wn1);
		Convar_to_char1D(re0, base_right, w0);
		Convar_to_char1D(rep1, base_right, wp1);
		Convar_to_char1D(rep2, base_right, wp2);

		for (int i = 0; i < 3; i++)
		{
			weno_5th_ao_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
		}
		Char_to_convar1D(right.convar, base_right, var);
		Char_to_convar1D(right.der1, base_right, der1);
		Char_to_convar1D(right.der2, base_right, der2);

	}

	Check_Order_Reduce(left, right, fluids[0]);

}

//=============================================================================
// WENO5_AO 重构多项式计算
// 在系数层面做 WENO-AO 加权，输出完整的 Cellpoly1d 多项式
// 不修改原有 WENO5_AO / weno_5th_ao_left / weno_5th_ao_right 函数
//=============================================================================

// 辅助函数：对单个标量分量，计算 WENO5_AO 加权后的多项式系数
// 输出 (q_out, A_out, B_out, C_out, D_out)，对应 Cellpoly1d 的正交基
//
// 候选多项式定义在局部坐标 ξ = (x - x_i) / h ∈ [-0.5, 0.5]
// Cellpoly1d 基函数：
//   φ_0 = 1,  φ_1 = ξ,  φ_2 = ξ²-1/12,  φ_3 = ξ³-3ξ/20,  φ_4 = ξ⁴-3ξ²/14+3/560
// 
// 每个候选多项式 P_k(ξ) 通过保均值条件唯一确定：
//   ∫_{-1/2}^{1/2} P_k(ξ - j) dξ = w_{i+j}  对模板中所有单元 j
//
// 小模板 (二次多项式，C=D=0):
//   S0: {i-2, i-1, i}    → P0(ξ) = q + A0*ξ + B0*(ξ²-1/12)
//   S1: {i-1, i,   i+1}  → P1(ξ) = q + A1*ξ + B1*(ξ²-1/12)
//   S2: {i,   i+1, i+2}  → P2(ξ) = q + A2*ξ + B2*(ξ²-1/12)
// 大模板 (四次多项式):
//   S3: {i-2, i-1, i, i+1, i+2} → P3(ξ) = q + A3*ξ + B3*(ξ²-1/12) + C3*(ξ³-3ξ/20) + D3*(ξ⁴-3ξ²/14+3/560)
//
static void weno_5th_ao_cellpoly(
	double& q_out, double& A_out, double& B_out, double& C_out, double& D_out,
	double wn2, double wn1, double w0, double wp1, double wp2)
{
	double dhi = 0.85;
	double dlo = 0.85;

	// WENO 参数
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-8;
	double sum_alpha;

	// 线性权重
	d[0] = (1 - dhi) * (1 - dlo) / 2.0;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) / 2.0;
	d[3] = dhi;

	// --- 光滑指标（与原函数完全一致）---
	if (wenotype == linear)
	{
		for (int k = 0; k < 4; k++)
			ww[k] = d[k];
	}
	else
	{
		beta[0] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w0), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((w0 - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w0 - 4 * wp1 + wp2), 2);

		beta[3] = (1.0 / 5040.0) * (231153 * w0 * w0 + 104963 * wn1 * wn1 + 6908 * wn2 * wn2 -
			38947 * wn2 * wp1 + 104963 * wp1 * wp1 +
			wn1 * (-51001 * wn2 + 179098 * wp1 - 38947 * wp2) -
			3 * w0 * (99692 * wn1 - 22641 * wn2 + 99692 * wp1 - 22641 * wp2) +
			8209 * wn2 * wp2 - 51001 * wp1 * wp2 + 6908 * wp2 * wp2);

		double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 4; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}
		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 4; i++)
			{
				double global_div = tau5 / (beta[i] + epsilonW);
				alpha[i] = d[i] * (1 + global_div * global_div);
				sum_alpha += alpha[i];
			}
		}

		for (int k = 0; k < 4; k++)
			ww[k] = alpha[k] / sum_alpha;
	}

	// --- 各候选多项式的 Cellpoly1d 系数 ---
	// q_k = w0 对所有模板（保均值条件: q 恒等于当前单元均值）
	// 
	// 小模板 S0: {i-2, i-1, i}
	//   A0 = (wn2 - 4*wn1 + 3*w0) / 2    （注意这是 ξ 坐标下的系数，不含 1/h）
	//   B0 = (wn2 - 2*wn1 + w0) / 2
	//   C0 = D0 = 0
	//
	// 小模板 S1: {i-1, i, i+1}
	//   A1 = (wp1 - wn1) / 2
	//   B1 = (wn1 - 2*w0 + wp1) / 2
	//   C1 = D1 = 0
	//
	// 小模板 S2: {i, i+1, i+2}
	//   A2 = (-3*w0 + 4*wp1 - wp2) / 2
	//   B2 = (w0 - 2*wp1 + wp2) / 2
	//   C2 = D2 = 0
	//
	// S3: {i-2, i-1, i, i+1, i+2}
    //   A3 = (-11*wp2 + 82*wp1 - 82*wn1 + 11*wn2) / 120
    //   B3 = (-3*wp2 + 40*wp1 - 74*w0 + 40*wn1 - 3*wn2) / 56
    //   C3 = (wp2 - 2*wp1 + 2*wn1 - wn2) / 12
    //   D3 = (wp2 - 4*wp1 + 6*w0 - 4*wn1 + wn2) / 24

	double q_k[4], A_k[4], B_k[4], C_k[4], D_k[4];

	// 所有模板的 q = w0
	for (int k = 0; k < 4; k++)
		q_k[k] = w0;

	// S0: {i-2, i-1, i}
	A_k[0] = (wn2 - 4.0 * wn1 + 3.0 * w0) / 2.0;
	B_k[0] = (wn2 - 2.0 * wn1 + w0) / 2.0;
	C_k[0] = 0.0;
	D_k[0] = 0.0;

	// S1: {i-1, i, i+1}
	A_k[1] = (wp1 - wn1) / 2.0;
	B_k[1] = (wn1 - 2.0 * w0 + wp1) / 2.0;
	C_k[1] = 0.0;
	D_k[1] = 0.0;

	// S2: {i, i+1, i+2}
	A_k[2] = (-3.0 * w0 + 4.0 * wp1 - wp2) / 2.0;
	B_k[2] = (w0 - 2.0 * wp1 + wp2) / 2.0;
	C_k[2] = 0.0;
	D_k[2] = 0.0;

	// S3: {i-2, i-1, i, i+1, i+2}
	A_k[3] = (-11.0 * wp2 + 82.0 * wp1 - 82.0 * wn1 + 11.0 * wn2) / 120.0;
    B_k[3] = (-3.0 * wp2 + 40.0 * wp1 - 74.0 * w0 + 40.0 * wn1 - 3.0 * wn2) / 56.0;
    C_k[3] = (wp2 - 2.0 * wp1 + 2.0 * wn1 - wn2) / 12.0;
    D_k[3] = (wp2 - 4.0 * wp1 + 6.0 * w0 - 4.0 * wn1 + wn2) / 24.0;

	// --- AO 权重修正（与原函数一致）---
	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	// --- 在系数层面做加权组合 ---
	q_out = 0.0;
	A_out = 0.0;
	B_out = 0.0;
	C_out = 0.0;
	D_out = 0.0;
	for (int k = 0; k < 4; k++)
	{
		q_out += final_weight[k] * q_k[k];
		A_out += final_weight[k] * A_k[k];
		B_out += final_weight[k] * B_k[k];
		C_out += final_weight[k] * C_k[k];
		D_out += final_weight[k] * D_k[k];
	}
}

// 计算单元 i 的 WENO5_AO 重构多项式（对外接口）
// 输入: fluids 指针指向当前单元 fluids[0]，需要 fluids[-2..2]
// 输出: poly 填入完整的 Cellpoly1d 多项式系数
void WENO5_AO_poly(Cellpoly1d& poly, Fluid1d* fluids)
{
	double wn2[3], wn1[3], w0[3], wp1[3], wp2[3];
	Copy_Array(wn2, fluids[-2].convar, 3);
	Copy_Array(wn1, fluids[-1].convar, 3);
	Copy_Array(w0,  fluids[0].convar,  3);
	Copy_Array(wp1, fluids[1].convar,  3);
	Copy_Array(wp2, fluids[2].convar,  3);

	poly.reset();

	for (int m = 0; m < 3; m++)
	{
		weno_5th_ao_cellpoly(
			poly.q[m], poly.A[m], poly.B[m], poly.C[m], poly.D[m],
			wn2[m], wn1[m], w0[m], wp1[m], wp2[m]);
	}
}

// 验证函数：对比 WENO5_AO_poly 与 weno_5th_ao_left/right 的边界值
// 在初始化/第一步之后调用一次即可，例如：
//   verify_WENO5_AO_poly(fluids, block);
void verify_WENO5_AO_poly(Fluid1d* fluids, Block1d block)
{
	cout << "====== verify_WENO5_AO_poly ======" << endl;
	double max_err_val = 0.0, max_err_der = 0.0;
	int test_count = 0;

	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		double h = fluids[i].dx;

		// --- 方法1: 原函数 weno_5th_ao_left/right ---
		Point1d pt_left, pt_right;
		double tmp;
		for (int m = 0; m < 3; m++)
		{
			weno_5th_ao_left(pt_left.convar[m], pt_left.der1[m], tmp,
				fluids[i - 2].convar[m], fluids[i - 1].convar[m],
				fluids[i].convar[m], fluids[i + 1].convar[m],
				fluids[i + 2].convar[m], h);

			weno_5th_ao_right(pt_right.convar[m], pt_right.der1[m], tmp,
				fluids[i - 2].convar[m], fluids[i - 1].convar[m],
				fluids[i].convar[m], fluids[i + 1].convar[m],
				fluids[i + 2].convar[m], h);
		}

		// --- 方法2: 多项式系数 ---
		Cellpoly1d poly;
		WENO5_AO_poly(poly, &fluids[i]);

		// --- 对比 ---
		for (int m = 0; m < 3; m++)
		{
			double val_left_poly  = poly.value(m, -0.5);
			double val_right_poly = poly.value(m,  0.5);
			double der_left_poly  = poly.dxValue(m, -0.5, h);
			double der_right_poly = poly.dxValue(m,  0.5, h);

			double err_vl = abs(val_left_poly  - pt_left.convar[m]);
			double err_vr = abs(val_right_poly - pt_right.convar[m]);
			double err_dl = abs(der_left_poly  - pt_left.der1[m]);
			double err_dr = abs(der_right_poly - pt_right.der1[m]);

			max_err_val = max(max_err_val, max(err_vl, err_vr));
			max_err_der = max(max_err_der, max(err_dl, err_dr));

			if (err_vl > 1e-10 || err_vr > 1e-10 || err_dl > 1e-10 || err_dr > 1e-10)
			{
				cout << "  MISMATCH at cell " << i << " comp " << m << ":"
					<< " err_vl=" << err_vl << " err_vr=" << err_vr
					<< " err_dl=" << err_dl << " err_dr=" << err_dr << endl;
			}
		}
		test_count++;
	}

	cout << "  Tested " << test_count << " cells." << endl;
	cout << "  Max point-value error: " << max_err_val << endl;
	cout << "  Max derivative  error: " << max_err_der << endl;
	if (max_err_val < 1e-10 && max_err_der < 1e-10)
		cout << "  Result: ALL PASSED" << endl;
	else
		cout << "  Result: FAILED - see MISMATCH above" << endl;
	cout << "==================================" << endl;
}

void weno_5th_ao_left(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-8;
	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1 - dhi) * (1 - dlo) / 2.0;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	if (wenotype == linear)
	{
		for (int k = 0; k < 4; k++)
		{
			ww[k] = d[k];
		}
	}
	else
	{
		//cout << "here" << endl;
		beta[0] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w0), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((w0 - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w0 - 4 * wp1 + wp2), 2);

		beta[3] = (1.0 / 5040.0) * (231153 * w0 * w0 + 104963 * wn1 * wn1 + 6908 * wn2 * wn2 -
			38947 * wn2 * wp1 + 104963 * wp1 * wp1 +
			wn1 * (-51001 * wn2 + 179098 * wp1 - 38947 * wp2) -
			3 * w0 * (99692 * wn1 - 22641 * wn2 + 99692 * wp1 - 22641 * wp2) +
			8209 * wn2 * wp2 - 51001 * wp1 * wp2 + 6908 * wp2 * wp2);

		double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 4; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}
		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 4; i++)
			{
				double global_div = tau5 / (beta[i] + epsilonW);
				alpha[i] = d[i] * (1 + global_div * global_div);
				sum_alpha += alpha[i];
			}
		}

		for (int k = 0; k < 4; k++)
		{
			ww[k] = alpha[k] / sum_alpha;
		}
	}
	//-- - candidate polynomial-- -
	p[0] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w0;
	p[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w0 - 1.0 / 6.0 * wp1;
	p[2] = 11.0 / 6.0 * w0 - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	p[3] = (1.0 / 60.0) * (47 * w0 + 27 * wn1 - 3 * wn2 - 13 * wp1 + 2 * wp2);

	px[0] = (w0 - wn1) / h;
	px[1] = (w0 - wn1) / h;
	px[2] = -((2 * w0 - 3 * wp1 + wp2) / h);
	px[3] = (15 * w0 - 15 * wn1 + wn2 - wp1) / (12 * h);

	pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
	pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
	pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
	pxx[3] = ((-8 * w0 + 2 * wn1 + wn2 + 6 * wp1 - wp2) / (4 * h * h));

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	der2 = 0.0;
	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	for (int k = 0; k < 4; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
		der2 += final_weight[k] * pxx[k];
	}
}
void weno_5th_ao_right(double& var, double& der1, double& der2, double wn2, double wn1, double w0, double wp1, double wp2, double h)
{
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-8;

	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1 - dhi) * (1 - dlo) / 2.0;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) / 2.0;
	//one big stencil
	d[3] = dhi;

	if (wenotype == linear)
	{
		for (int k = 0; k < 4; k++)
		{
			ww[k] = d[k];
		}
	}
	else
	{

		beta[0] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w0), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w0), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w0 + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((w0 - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w0 - 4 * wp1 + wp2), 2);

		beta[3] = (1.0 / 5040.0) * (231153 * w0 * w0 + 104963 * wn1 * wn1 + 6908 * wn2 * wn2 -
			38947 * wn2 * wp1 + 104963 * wp1 * wp1 +
			wn1 * (-51001 * wn2 + 179098 * wp1 - 38947 * wp2) -
			3 * w0 * (99692 * wn1 - 22641 * wn2 + 99692 * wp1 - 22641 * wp2) +
			8209 * wn2 * wp2 - 51001 * wp1 * wp2 + 6908 * wp2 * wp2);

		double tau5 = 1.0 / 3.0 * (abs(beta[3] - beta[0]) + abs(beta[3] - beta[1]) + abs(beta[3] - beta[2]));

		if (wenotype == wenojs)
		{
			sum_alpha = 0.0;
			for (int k = 0; k < 4; k++)
			{
				alpha[k] = d[k] / ((epsilonW + beta[k]) * (epsilonW + beta[k]));
				sum_alpha += alpha[k];
			}
		}
		else if (wenotype == wenoz)
		{
			sum_alpha = 0.0;
			for (int i = 0; i < 4; i++)
			{
				double global_div = tau5 / (beta[i] + epsilonW);
				alpha[i] = d[i] * (1 + global_div * global_div);
				sum_alpha += alpha[i];
			}
		}

		for (int k = 0; k < 4; k++)
		{
			ww[k] = alpha[k] / sum_alpha;
		}

	}
	//-- - candidate polynomial-- -

	p[0] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w0;
	p[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w0 + 1.0 / 3.0 * wp1;
	p[2] = 1.0 / 3.0 * w0 + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	p[3] = (1.0 / 60.0) * (47 * w0 - 13 * wn1 + 2 * wn2 + 27 * wp1 - 3 * wp2);

	px[0] = (2 * w0 - 3 * wn1 + wn2) / h;
	px[1] = (-w0 + wp1) / h;
	px[2] = (-w0 + wp1) / h;
	px[3] = (-15 * w0 + wn1 + 15 * wp1 - wp2) / (12 * h);

	pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
	pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
	pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
	pxx[3] = (-8 * w0 + 6 * wn1 - wn2 + 2 * wp1 + wp2) / (4 * h * h);

	//-- - combination-- -
	var = 0.0;
	der1 = 0.0;
	der2 = 0.0;
	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	for (int k = 0; k < 4; k++)
	{
		var += final_weight[k] * p[k];
		der1 += final_weight[k] * px[k];
		der2 += final_weight[k] * pxx[k];
	}
}

namespace
{
	void Evaluate_Cellpoly_Interface_1D(double* convar, const Cellpoly1d& poly, double xi)
	{
		for (int m = 0; m < 3; ++m)
		{
			convar[m] = poly.value(m, xi);
		}
	}

	void Collision_Interface_1D(double* wbar, const double* convar_left, const double* convar_right)
	{
		double left_copy[3], right_copy[3];
		for (int m = 0; m < 3; ++m)
		{
			left_copy[m] = convar_left[m];
			right_copy[m] = convar_right[m];
		}

		double prim_left[3], prim_right[3];
		Convar_to_ULambda_1d(prim_left, left_copy);
		Convar_to_ULambda_1d(prim_right, right_copy);

		MMDF1d ml(prim_left[1], prim_left[2]);
		MMDF1d mr(prim_right[1], prim_right[2]);
		double unit[3]{ 1.0, 0.0, 0.0 };
		double gl[3], gr[3];
		GL(0, 0, gl, unit, ml);
		GR(0, 0, gr, unit, mr);
		for (int m = 0; m < 3; ++m)
		{
			wbar[m] = convar_left[0] * gl[m] + convar_right[0] * gr[m];
		}
	}
}


void Reconstruction_forg0(Interface1d* interfaces, Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; ++i)
	{
		g0reconstruction(interfaces[i], &fluids[i - 1]);
	}
}

void Center_2nd_collisionless(Interface1d& interfaces, Fluid1d* fluids)
{
	double w[3], wp1[3];

	//assume this is uniform mesh,
	//when mesh is no uniform, accuracy decrease
	double h = (fluids[0].dx + fluids[1].dx) / 2.0;
	for (int i = 0; i < 3; i++)
	{
		w[i] = fluids[0].convar[i];
		wp1[i] = fluids[1].convar[i];
	}
	//collisionless
	//at this time, only second order accuracy obtained
	for (int i = 0; i < 3; i++)
	{
		interfaces.center.convar[i] = (wp1[i] + w[i]) / 2.0;
		interfaces.center.der1[i] = (wp1[i] - w[i]) / h;
		interfaces.center.der2[i] = 0.0;

	}
}

void Center_3rd(Interface1d& interfaces, Fluid1d* fluids)
{
	double w[3], wp1[3];

	//assume this is uniform mesh,
	//when mesh is no uniform, accuracy decrease
	double h = (fluids[0].dx + fluids[1].dx) / 2.0;
	for (int i = 0; i < 3; i++)
	{
		w[i] = fluids[0].convar[i];
		wp1[i] = fluids[1].convar[i];
	}
	// to separate the reconstruction part from the flux part,
	// here we do a collision to get g0
	double convar_left[3], convar_right[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = interfaces.left.convar[i];
		convar_right[i] = interfaces.right.convar[i];
	}

	double prim_left[3], prim_right[3];
	Convar_to_ULambda_1d(prim_left, convar_left);
	Convar_to_ULambda_1d(prim_right, convar_right);

	MMDF1d ml(prim_left[1], prim_left[2]);
	MMDF1d mr(prim_right[1], prim_right[2]);

	double unit[3]{ 1.0, 0.0, 0.0 };

	double gl[3], gr[3];
	GL(0, 0, gl, unit, ml);
	GR(0, 0, gr, unit, mr);
	for (int i = 0; i < 3; i++)
	{
		interfaces.center.convar[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
	}
	for (int i = 0; i < 3; i++)
	{
		interfaces.center.der1[i] = (wp1[i] - w[i]) / h;
		interfaces.center.der2[i] = 0.0;

	}

}

void Center_GKS3rd(Interface1d& interfaces, Fluid1d* fluids)
{
	// Paper Sec. 3.4 discusses a dedicated equilibrium reconstruction at the
	// interface. In this 1D minimal implementation, we first construct the
	// equilibrium state Wbar at three neighboring interfaces and then recover
	// Wbar_x and Wbar_xx by centered finite differences on the uniform mesh.
	const double h = fluids[0].dx;

	Cellpoly1d poly_m1, poly_0, poly_p1, poly_p2;
	WENO5_AO_poly(poly_m1, &fluids[-1]);
	WENO5_AO_poly(poly_0, &fluids[0]);
	WENO5_AO_poly(poly_p1, &fluids[1]);
	WENO5_AO_poly(poly_p2, &fluids[2]);

	double wl_imh[3], wr_imh[3];
	double wl_iph[3], wr_iph[3];
	double wl_iph3[3], wr_iph3[3];
	Evaluate_Cellpoly_Interface_1D(wl_imh, poly_m1, 0.5);
	Evaluate_Cellpoly_Interface_1D(wr_imh, poly_0, -0.5);
	Evaluate_Cellpoly_Interface_1D(wl_iph, poly_0, 0.5);
	Evaluate_Cellpoly_Interface_1D(wr_iph, poly_p1, -0.5);
	Evaluate_Cellpoly_Interface_1D(wl_iph3, poly_p1, 0.5);
	Evaluate_Cellpoly_Interface_1D(wr_iph3, poly_p2, -0.5);

	double wbar_imh[3], wbar_iph[3], wbar_iph3[3];
	Collision_Interface_1D(wbar_imh, wl_imh, wr_imh);
	Collision_Interface_1D(wbar_iph, wl_iph, wr_iph);
	Collision_Interface_1D(wbar_iph3, wl_iph3, wr_iph3);

	for (int m = 0; m < 3; ++m)
	{
		interfaces.center.convar[m] = wbar_iph[m];
		interfaces.center.der1[m] = (wbar_iph3[m] - wbar_imh[m]) / (2.0 * h);
		interfaces.center.der2[m] = (wbar_iph3[m] - 2.0 * wbar_iph[m] + wbar_imh[m]) / (h * h);
	}
}

void Center_simple_avg(Interface1d& interfaces, Fluid1d* fluids)
{
	for (int i = 0; i < 3; i++)
	{
		interfaces.center.convar[i] = (interfaces.left.convar[i] + interfaces.right.convar[i]) / 2.0;
		interfaces.center.der1[i] = (interfaces.left.der1[i] + interfaces.right.der1[i]) / 2.0;
		interfaces.center.der2[i] = 0.0;

	}
}

void Center_collision(Interface1d& interfaces, Fluid1d* fluids)
{
	double convar_left[3], convar_right[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = interfaces.left.convar[i];
		convar_right[i] = interfaces.right.convar[i];
	}

	double prim_left[3], prim_right[3]; //rho, U, lambda
	Convar_to_ULambda_1d(prim_left, convar_left);
	Convar_to_ULambda_1d(prim_right, convar_right);

	MMDF1d ml(prim_left[1], prim_left[2]);
	MMDF1d mr(prim_right[1], prim_right[2]);

	double unit[3]{ 1.0, 0.0, 0.0 };

	double gl[3], gr[3];
	GL(0, 0, gl, unit, ml); // gl, means Wl(u>0), by input uint
	GR(0, 0, gr, unit, mr); // gr, means Wr(u<0), by input uint
	double axl[3], axr[3];
	Microslope(axl, interfaces.left.der1, prim_left); // axl, means a coefficient indicating slope
	Microslope(axr, interfaces.right.der1, prim_right); // axr, means a coefficient indicating slope
	double ax0l[3], ax0r[3];
	GL(0, 0, ax0l, axl, ml);  // ax0l, means Wlx(u>0), by input axl
	GR(0, 0, ax0r, axr, mr); // ax0r, means Wrx(u<0), by input axr
	for (int i = 0; i < 3; i++)
	{
		interfaces.center.convar[i] = convar_left[0] * gl[i] + convar_right[0] * gr[i];
		interfaces.center.der1[i] = convar_left[0] * ax0l[i] + convar_right[0] * ax0r[i];
		interfaces.center.der2[i] = 0.0;

	}
}


void Calculate_flux(Flux1d** fluxes, Interface1d* interfaces, Block1d& block, int stage)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
	{
		flux_function(fluxes[i][stage], interfaces[i], block.dt);

	}
}

namespace
{
	struct GKS3rdTimeCoeff1D
	{
		double neq_g;
		double neq_au;
		double neq_axx;
		double eq_g;
		double eq_au;
		double eq_A;
		double eq_axx;
		double eq_axt;
		double eq_att;
	};

	double TauEulerPaper1D(const double* convar_left, const double* convar_right, double dt)
	{
		double prim_left[3], prim_right[3];
		double left_copy[3], right_copy[3];
		for (int m = 0; m < 3; ++m)
		{
			left_copy[m] = convar_left[m];
			right_copy[m] = convar_right[m];
		}
		Convar_to_primvar_1D(prim_left, left_copy);
		Convar_to_primvar_1D(prim_right, right_copy);
		const double pl = prim_left[2];
		const double pr = prim_right[2];
		const double denom = fabs(pl + pr);
		const double jump = (denom > 1e-14) ? fabs(pl - pr) / denom : 0.0;
		return c1_euler * dt + c2_euler * jump * dt;
	}

	void GetGKS3rdTimeCoeff1D(GKS3rdTimeCoeff1D& tc, double T, double tau_n)
	{
		if (tau_n <= 1e-14)
		{
			tc.neq_g = 0.0;
			tc.neq_au = 0.0;
			tc.neq_axx = 0.0;
			tc.eq_g = T;
			tc.eq_au = 0.0;
			tc.eq_A = 0.5 * T * T;
			tc.eq_axx = 0.0;
			tc.eq_axt = 0.0;
			tc.eq_att = T * T * T / 6.0;
			return;
		}

		const double tau = tau_n;
		const double eta = exp(-T / tau);
		tc.neq_g = tau * (1.0 - eta);
		tc.neq_au = tau * ((T + tau) * eta - tau);
		tc.neq_axx = tau * tau * tau - 0.5 * tau * eta * (T * T + 2.0 * tau * T + 2.0 * tau * tau);

		tc.eq_g = T - tau + tau * eta;
		tc.eq_au = 2.0 * tau * tau * (1.0 - eta) - tau * T * (1.0 + eta);
		tc.eq_A = 0.5 * T * T - tau * T + tau * tau * (1.0 - eta);
		tc.eq_axx = tau * tau * T - 3.0 * tau * tau * tau
			+ 0.5 * tau * eta * (T * T + 4.0 * tau * T + 6.0 * tau * tau);
		tc.eq_axt = -0.5 * tau * T * T + 2.0 * tau * tau * T - 3.0 * tau * tau * tau
			+ tau * tau * (T + 3.0 * tau) * eta;
		tc.eq_att = T * T * T / 6.0 - 0.5 * tau * T * T + tau * tau * T - tau * tau * tau * (1.0 - eta);
	}
}

void GKS(Flux1d& flux, Interface1d& interface, double dt)
{
	if (gks1dsolver == nothing)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	if (gks1dsolver == gks3rd)
	{
		double convar_left[3], convar_right[3], convar_bar[3];
		double der1_left[3], der1_right[3], der1_bar[3];
		double der2_left[3], der2_right[3], der2_bar[3];
		for (int m = 0; m < 3; ++m)
		{
			convar_left[m] = interface.left.convar[m];
			convar_right[m] = interface.right.convar[m];
			convar_bar[m] = interface.center.convar[m];
			der1_left[m] = interface.left.der1[m];
			der1_right[m] = interface.right.der1[m];
			der1_bar[m] = interface.center.der1[m];
			der2_left[m] = interface.left.der2[m];
			der2_right[m] = interface.right.der2[m];
			der2_bar[m] = interface.center.der2[m];
			flux.derf[m] = 0.0;
			flux.der2f[m] = 0.0;
		}

		double prim_left[3], prim_right[3], prim_bar[3];
		Convar_to_ULambda_1d(prim_left, convar_left);
		Convar_to_ULambda_1d(prim_right, convar_right);
		Convar_to_ULambda_1d(prim_bar, convar_bar);

		const double p_left = Pressure(convar_left[0], convar_left[1], convar_left[2]);
		const double p_right = Pressure(convar_right[0], convar_right[1], convar_right[2]);
		const double pressure_jump = fabs(p_left - p_right) / (fabs(p_left + p_right) + 1e-14);

		const double tau_n = TauEulerPaper1D(convar_left, convar_right, dt);
		GKS3rdTimeCoeff1D tc;
		GetGKS3rdTimeCoeff1D(tc, dt, tau_n);

		double a_left[3], a_right[3], a_bar[3], A_bar[3];
		double axx_left[3], axx_right[3], axx_bar[3], axt_bar[3], att_bar[3];
		SolveGKS_a1D(a_left, der1_left, prim_left);
		SolveGKS_a1D(a_right, der1_right, prim_right);
		SolveGKS_a1D(a_bar, der1_bar, prim_bar);
		SolveGKS_A1D(A_bar, a_bar, prim_bar);
		SolveGKS_axx1D(axx_left, der2_left, prim_left);
		SolveGKS_axx1D(axx_right, der2_right, prim_right);
		SolveGKS_axx1D(axx_bar, der2_bar, prim_bar);
		SolveGKS_axt1D(axt_bar, axx_bar, prim_bar);
		SolveGKS_att1D(att_bar, axt_bar, prim_bar);

		if (pressure_jump > 0.05)
		{
			// The paper remarks that the one-step 3rd-order GKS is less robust in
			// discontinuous regions. For the current minimal 1D implementation, we
			// locally switch off the genuinely 3rd-order terms and keep the 2nd-order
			// GKS core near strong jumps.
			Array_zero(axx_left, 3);
			Array_zero(axx_right, 3);
			Array_zero(axx_bar, 3);
			Array_zero(axt_bar, 3);
			Array_zero(att_bar, 3);
			tc.neq_axx = 0.0;
			tc.eq_axx = 0.0;
			tc.eq_axt = 0.0;
			tc.eq_att = 0.0;
		}

		MMDF1d ml(prim_left[1], prim_left[2]);
		MMDF1d mr(prim_right[1], prim_right[2]);
		MMDF1d mb(prim_bar[1], prim_bar[2]);
		double unit[3]{ 1.0, 0.0, 0.0 };
		double moment[3];
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] = 0.0;
		}

		GL(1, 0, moment, unit, ml);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_left[0] * tc.neq_g * moment[m];
		}
		GR(1, 0, moment, unit, mr);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_right[0] * tc.neq_g * moment[m];
		}
		GL(2, 0, moment, a_left, ml);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_left[0] * tc.neq_au * moment[m];
		}
		GR(2, 0, moment, a_right, mr);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_right[0] * tc.neq_au * moment[m];
		}
		GL(3, 0, moment, axx_left, ml);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_left[0] * tc.neq_axx * moment[m];
		}
		GR(3, 0, moment, axx_right, mr);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_right[0] * tc.neq_axx * moment[m];
		}

		G(1, 0, moment, unit, mb);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_bar[0] * tc.eq_g * moment[m];
		}
		G(2, 0, moment, a_bar, mb);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_bar[0] * tc.eq_au * moment[m];
		}
		G(1, 0, moment, A_bar, mb);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_bar[0] * tc.eq_A * moment[m];
		}
		G(3, 0, moment, axx_bar, mb);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_bar[0] * tc.eq_axx * moment[m];
		}
		G(2, 0, moment, axt_bar, mb);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_bar[0] * tc.eq_axt * moment[m];
		}
		G(1, 0, moment, att_bar, mb);
		for (int m = 0; m < 3; ++m)
		{
			flux.f[m] += prim_bar[0] * tc.eq_att * moment[m];
		}
		return;
	}
	double Flux[2][3];
	//change conservative variables to rho u lambda
	double convar_left[3], convar_right[3], convar0[3];
	for (int i = 0; i < 3; i++)
	{
		convar_left[i] = interface.left.convar[i];
		convar_right[i] = interface.right.convar[i];
		convar0[i] = interface.center.convar[i];
	}

	double prim_left[3], prim_right[3], prim0[3];
	Convar_to_ULambda_1d(prim_left, convar_left);
	Convar_to_ULambda_1d(prim_right, convar_right);
	Convar_to_ULambda_1d(prim0, convar0);

	//then lets get the coefficient of time intergation factors

	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[2]);
	tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[2], prim_right[2], prim0[2], dt);
	double eta = exp(-dt / tau_num);
	double t[10];
	// non equ part time coefficient for gks_2nd algorithm (f0)
	t[0] = tau_num * (1 - eta); // this refers glu, gru part
	t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
	t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
	// then, equ part time coefficient for gks 2nd (g0)
	t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
	t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
	t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt; //this refers A0u part

	if (gks1dsolver == kfvs1st)
	{
		t[0] = dt;
		for (int i = 1; i < 6; i++)
		{
			t[i] = 0.0;
		}
		//do nothing, kfvs1st only use t[0]=dt part;
	}
	else if (gks1dsolver == kfvs2nd)
	{
		t[0] = dt;
		t[1] = -dt * dt / 2.0;
		for (int i = 2; i < 6; i++)
		{
			t[i] = 0.0;
		}
	}
	MMDF1d ml(prim_left[1], prim_left[2]);
	MMDF1d mr(prim_right[1], prim_right[2]);

	double unit[3] = { 1, 0.0, 0.0 };

	double glu[3], gru[3];
	GL(1, 0, glu, unit, ml);
	GR(1, 0, gru, unit, mr);

	//only one part, the kfvs1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
	}

	if (gks1dsolver == kfvs1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// kfvs1st part ended

	//now the equ part added, m0 term added, gks1st part begin
	MMDF1d m0(prim0[1], prim0[2]);

	double g0u[3];
	G(1, 0, g0u, unit, m0);

	//the equ g0u part, the gks1st part
	for (int i = 0; i < 3; i++)
	{
		Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
	}

	if (gks1dsolver == gks1st)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// gks1d solver ended

	//for kfvs2nd part
	double der1left[3], der1right[3];
	for (int i = 0; i < 3; i++)
	{
		der1left[i] = interface.left.der1[i];
		der1right[i] = interface.right.der1[i];
	}

	double alx[3];
	Microslope(alx, der1left, prim_left);

	double alxuul[3];
	GL(2, 0, alxuul, alx, ml);

	double arx[3];
	Microslope(arx, der1right, prim_right);
	double arxuur[3];
	GR(2, 0, arxuur, arx, mr);

	for (int i = 0; i < 3; i++)
	{	// t1 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
	}
	if (gks1dsolver == kfvs2nd)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// the kfvs2nd part ended

	// then we still need t[2], t[4] t[5] part for gks 2nd
	//for t[2] Aru,Alu part
	double alxu[3];
	double arxu[3];

	//take <u> moment for al, ar
	G(1, 0, alxu, alx, ml);
	G(1, 0, arxu, arx, mr);

	double Al[3], Ar[3];
	double der_AL[3], der_AR[3];

	//using compatability condition to get the time derivative
	for (int i = 0; i < 3; i++)
	{
		der_AL[i] = -prim_left[0] * (alxu[i]);
		der_AR[i] = -prim_right[0] * (arxu[i]);
	}
	// solve the coefficient martix b=ma
	Microslope(Al, der_AL, prim_left);
	Microslope(Ar, der_AR, prim_right);

	//to obtain the Alu and Aru
	double Alul[3];
	double Arur[3];
	GL(1, 0, Alul, Al, ml);
	GR(1, 0, Arur, Ar, mr);

	for (int i = 0; i < 3; i++)
	{	// t2 part
		Flux[0][i] = Flux[0][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
	}

	// for t[4] a0xuu part

	double a0x[3];
	double der1[3];

	for (int i = 0; i < 3; i++)
	{
		der1[i] = interface.center.der1[i];
	}

	//solve the microslope
	Microslope(a0x, der1, prim0);
	//a0x <u> moment
	double a0xu[3];
	G(1, 0, a0xu, a0x, m0);
	//a0x <u^2> moment
	double a0xuu[3];
	G(2, 0, a0xuu, a0x, m0);

	for (int i = 0; i < 3; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i]);
	}


	// for t[5] A0u part
	double derA0[3];

	for (int i = 0; i < 3; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i]);
	}
	double A0[3];
	Microslope(A0, derA0, prim0);
	double A0u[3];
	G(1, 0, A0u, A0, m0);
	for (int i = 0; i < 3; i++)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}
	if (gks1dsolver == gks2nd && timecoe_list == S1O1)
	{
		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	if (gks1dsolver == gks2nd)
	{
		double dt2 = 0.5 * dt;
		eta = exp(-dt2 / tau_num);
		// non equ part time coefficient for gks_2nd algorithm
		t[0] = tau_num * (1 - eta); // this refers glu, gru part
		t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
		t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
		// then, equ part time coefficient for gks 2nd
		t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
		t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
		t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

		for (int i = 0; i < 3; i++)
		{
			// t0 part
			Flux[1][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
			// t1 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[1] * (alxuul[i]) + prim_right[0] * t[1] * (arxuur[i]);
			// t2 part
			Flux[1][i] = Flux[1][i] + prim_left[0] * t[2] * (Alul[i]) + prim_right[0] * t[2] * (Arur[i]);
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}

		for (int i = 0; i < 3; i++)
		{
			flux.f[i] = (4.0 * Flux[1][i] - Flux[0][i]);
			flux.derf[i] = 4.0 * (Flux[0][i] - 2.0 * Flux[1][i]);
		}

		return;
	}
	else
	{
		cout << "no valid solver specify" << endl;
		exit(0);
	}
}
void HLLC(Flux1d& flux, Interface1d& interface, double dt)
{
	double pl[3], pr[3];
	Convar_to_primvar_1D(pl, interface.left.convar);
	Convar_to_primvar_1D(pr, interface.right.convar);

	double al, ar, pvars, pstar, tmp1, tmp2, tmp3, qk, sl, sr, star;
	al = sqrt(Gamma * pl[2] / pl[0]); //sound speed
	ar = sqrt(Gamma * pr[2] / pr[0]);
	tmp1 = 0.5 * (al + ar);         //avg of sound and density
	tmp2 = 0.5 * (pl[0] + pr[0]);

	pvars = 0.5 * (pl[2] + pr[2]) - 0.5 * (pr[1] - pl[1]) * tmp1 * tmp2;
	pstar = fmax(0.0, pvars);

	double flxtmp[3], qstar[3];


	ESTIME(sl, star, sr, pl[0], pl[1], pl[2], al, pr[0], pr[1], pr[2], ar);

	tmp1 = pr[2] - pl[2] + pl[0] * pl[1] * (sl - pl[1]) - pr[0] * pr[1] * (sr - pr[1]);
	tmp2 = pl[0] * (sl - pl[1]) - pr[0] * (sr - pr[1]);
	star = tmp1 / tmp2;

	if (sl >= 0.0) {
		get_Euler_flux(pl, flux.f);
	}
	else if (sr <= 0.0)
	{
		get_Euler_flux(pr, flux.f);
	}
	else if ((star >= 0.0) && (sl <= 0.0))
	{
		get_Euler_flux(pl, flxtmp);
		ustarforHLLC(pl[0], pl[1], pl[2], sl, star, qstar);

		for (int m = 0; m < 3; m++)
		{
			flux.f[m] = flxtmp[m] + sl * (qstar[m] - interface.left.convar[m]);
		}
	}
	else if ((star <= 0.0) && (sr >= 0.0))
	{
		get_Euler_flux(pr, flxtmp);
		ustarforHLLC(pr[0], pr[1], pr[2], sr, star, qstar);
		for (int m = 0; m < 3; m++)
		{
			flux.f[m] = flxtmp[m] + sr * (qstar[m] - interface.right.convar[m]);
		}
	}
	for (int m = 0; m < 3; m++)
	{
		flux.f[m] *= dt;
	}
}
void get_Euler_flux(double p[3], double* flux)
{
	//the flux of Euler equation, density*u, density*u*u+p, (p+0.5*density*u*u+p/(r-1))*u
	//p[0 1 2], means, density, u, p
	flux[0] = p[0] * p[1];
	flux[1] = p[0] * p[1] * p[1] + p[2];
	double ENERGS = 0.5 * (p[1] * p[1]) * p[0] + p[2] / (Gamma - 1.0);
	flux[2] = p[1] * (ENERGS + p[2]);
}
void ustarforHLLC(double d1, double u1, double p1, double s1, double star1, double* ustar)
{
	double tmp1, tmp2, tmp3;
	tmp1 = d1 * (s1 - u1) / (s1 - star1);
	tmp2 = 0.5 * (u1 * u1) + p1 / ((Gamma - 1.0) * d1);
	tmp3 = star1 + p1 / (d1 * (s1 - u1));

	ustar[0] = tmp1;
	ustar[1] = tmp1 * star1;
	ustar[2] = tmp1 * (tmp2 + (star1 - u1) * tmp3);
}
void LF(Flux1d& flux, Interface1d& interface, double dt)
{
	double pl[3], pr[3];
	Convar_to_primvar_1D(pl, interface.left.convar);
	Convar_to_primvar_1D(pr, interface.right.convar);

	//Reimann invariants, u +/- 2c/(r-1)
	//Sound speed, c = sqrt(r*p/density)
	double k[2];
	k[0] = abs(pl[1]) + sqrt(Gamma * pl[2] / pl[0]); //abs(u)+c, left
	k[1] = abs(pr[1]) + sqrt(Gamma * pr[2] / pr[0]); //abs(u)+c, right

	double beta = k[0]; // beta, means, abs(partialF/parialX)
	if (k[1] > k[0]) { beta = k[1]; }
	double flux_l[3], flux_r[3];
	get_Euler_flux(pl, flux_l);
	get_Euler_flux(pr, flux_r);

	for (int m = 0; m < 3; m++)
	{
		flux.f[m] = 0.5 * ((flux_l[m] + flux_r[m]) - beta * (interface.right.convar[m] - interface.left.convar[m]));
		flux.f[m] *= dt;
	}
}


void S1O1(Block1d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;
}
void S1O2(Block1d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;
}
void S1O3(Block1d& block)
{
	// Paper Sec. 2.3: the one-stage 3rd-order accuracy comes from the
	// time-accurate GKS flux itself. Here Update() still consumes a single
	// time-integrated flux over the whole step, so the stage coefficient is one.
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;
	block.timecoefficient[0][0][2] = 0.0;
}
void S2O4(Block1d& block)
{
	// two stages, so the extra stages coefficients are zero.
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0 / 8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0 / 6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;
}
void RK4(Block1d& block)
{
	block.stages = 4;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;
	block.timecoefficient[2][2][0] = 1.0;
	block.timecoefficient[3][0][0] = 1.0 / 6.0;
	block.timecoefficient[3][1][0] = 1.0 / 3.0;
	block.timecoefficient[3][2][0] = 1.0 / 3.0;
	block.timecoefficient[3][3][0] = 1.0 / 6.0;
}

void Initial_stages(Block1d& block)
{
	// initialize five stages (5 steps) and set all them as ZERO
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage
		{
			for (int k = 0; k < 3; k++) //refers f, derf
			{
				block.timecoefficient[i][j][k] = 0.0;
			}
		}
	}
	// by timecoe_list set the correct stages (others still be ZERO)
	timecoe_list(block);
}

Flux1d** Setflux_array(Block1d block)
{
	Flux1d** fluxes = new Flux1d * [block.nx + 1];  // dynamic variable (since block.nx is not determined)

	for (int i = 0; i <= block.nx; i++)
	{
		// for m th step time marching schemes, m subflux needed
		fluxes[i] = new Flux1d[block.stages];
		for (int j = 0; j < block.stages; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				fluxes[i][j].f[k] = 0.0;
				fluxes[i][j].derf[k] = 0.0;
				fluxes[i][j].der2f[k] = 0.0;
			}
		}
	}
	if (fluxes == 0)
	{
		cout << "memory allocation failed for muli-stage flux";
		return NULL;
	}
	cout << "the memory for muli-stage flux has been allocated..." << endl;
	return fluxes;
}

void SetUniformMesh(Block1d block, Fluid1d* fluids, Interface1d* interfaces, Flux1d** fluxes)
{
	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		fluids[i].dx = block.dx; //cell size
		fluids[i].cx = block.left + (i + 0.5 - block.ghost) * block.dx; //cell center location
	}
	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		interfaces[i].x = block.left + (i - block.ghost) * block.dx;
		interfaces[i].left.x = interfaces[i].x;
		interfaces[i].right.x = interfaces[i].x;
		interfaces[i].center.x = interfaces[i].x;
		interfaces[i].flux = fluxes[i];
	}
}


void CopyFluid_new_to_old(Fluid1d* fluids, Block1d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fluids[i].convar_old[j] = fluids[i].convar[j];
		}
	}
}


double Get_CFL(Block1d& block, Fluid1d* fluids, double tstop)
{
	double dt = block.dx;
	for (int i = 0; i < block.nodex; i++)
	{
		dt = Dtx(dt, block.dx, block.CFL, fluids[i + block.ghost].convar);
	}
	if (block.t + dt > tstop)
	{
		dt = tstop - block.t + 1e-15;
	}
	//print time step information
	if (is_show_1d_timestep)
	{
		cout << "step= " << block.step
			<< "time size is " << dt
			<< " time= " << block.t << endl;
	}
	return dt;
}
double Dtx(double dtx, double dx, double CFL, double convar[3])
{
	double tmp;
	double prim[3];
	Convar_to_primvar_1D(prim, convar);
	tmp = abs(prim[1]) + sqrt(Gamma * prim[2] / prim[0]);
	if (tmp > CFL * dx / dtx) // if abs(u)+c (tmp) > abs(u)+c (old)
	{
		dtx = CFL * dx / tmp; // abs(u)+c determine one smaller time step
	}
	//consider viscous problem (tau_type == NS)
	if (dtx > 0.25 * CFL * dx * dx / Mu && tau_type == NS && Mu > 0)
	{
		// Note: if dtx (above) > CFL * dx * dx / 4
		// replace dxt by the smaller time step determined by viscous term
		dtx = 0.25 * CFL * dx * dx / Mu;
	}
	return dtx;
}


void Update(Fluid1d* fluids, Flux1d** fluxes, Block1d block, int stage)
{
	if (stage > block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			double Flux = 0.0;
			for (int k = 0; k < stage + 1; ++k)
			{
				Flux = Flux
					+ block.timecoefficient[stage][k][0] * fluxes[i][k].f[j]
					+ block.timecoefficient[stage][k][1] * fluxes[i][k].derf[j];
			}
			fluxes[i][stage].F[j] = Flux;
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			fluids[i].convar[j] = fluids[i].convar_old[j] + 1.0 / fluids[i].dx * (fluxes[i][stage].F[j] - fluxes[i + 1][stage].F[j]);
		}
	}

}






void ESTIME(double& SL, double& SM, double& SR, double DL, double UL, double PL, double CL, double DR, double UR, double PR, double CR)
{
	//this program is copied from toro
	//to compute wave speed estimates for the HLLC Riemann
	//solver usingand adaptive approximate - state Riemann
	//solver including the PVRS, TRRSand TSRS solvers
	double QUSER = 2.0;
	double G1 = (Gamma - 1.0) / (2.0 * Gamma);
	double G2 = (Gamma + 1.0) / (2.0 * Gamma);
	double G3 = 2.0 * Gamma / (Gamma - 1.0);
	double G4 = 2.0 / (Gamma - 1.0);
	double G5 = 2.0 / (Gamma + 1.0);
	double G6 = (Gamma - 1.0) / (Gamma + 1.0);
	double G7 = (Gamma - 1.0) / 2.0;
	double G8 = Gamma - 1.0;
	//     Compute guess pressure from PVRS Riemann solver
	double	CUP = 0.25 * (DL + DR) * (CL + CR);
	double	PPV = 0.5 * (PL + PR) + 0.5 * (UL - UR) * CUP;
	PPV = fmax(0.0, PPV);
	double	PMIN = fmin(PL, PR);
	double	PMAX = fmax(PL, PR);
	double	QMAX = PMAX / PMIN;
	double PM;
	double UM;


	double PQ = pow((PL / PR), G1);
	UM = (PQ * UL / CL + UR / CR + G4 * (PQ - 1.0)) / (PQ / CL + 1.0 / CR);
	double	PTL = 1.0 + G7 * (UL - UM) / CL;
	double	PTR = 1.0 + G7 * (UM - UR) / CR;
	PM = 0.5 * (PL * pow(PTL, G3) + PR * pow(PTR, G3));


	PM = pow((CL + CR - G7 * (UR - UL)) / (CL / pow(PL, G1) + CR / pow(PR, G1)), G3);
	//Find speeds
	if (PM <= PL)
	{
		SL = UL - CL;
	}
	else
	{
		SL = UL - CL * sqrt(1.0 + G2 * (PM / PL - 1.0));
	}

	SM = UM;

	if (PM <= PR)
	{
		SR = UR + CR;
	}
	else
	{
		SR = UR + CR * sqrt(1.0 + G2 * (PM / PR - 1.0));
	}
}

void get_flux(double p[4], double* flux)
{
	flux[0] = p[0] * p[1];
	flux[1] = p[0] * p[1] * p[1] + p[3];
	flux[2] = p[0] * p[1] * p[2];
	double ENERGS = 0.5 * (p[1] * p[1] + p[2] * p[2]) * p[0] + p[3] / (Gamma - 1.0);
	flux[3] = p[1] * (ENERGS + p[3]);
}

void ustarforHLLC(double d1, double u1, double v1, double p1, double s1, double star1, double* ustar)
{

	double tmp1, tmp2, tmp3;

	tmp1 = d1 * (s1 - u1) / (s1 - star1);
	tmp2 = 0.5 * (u1 * u1 + v1 * v1) + p1 / ((Gamma - 1.0) * d1);
	tmp3 = star1 + p1 / (d1 * (s1 - u1));

	ustar[0] = tmp1;
	ustar[1] = tmp1 * star1;
	ustar[2] = tmp1 * v1;
	ustar[3] = tmp1 * (tmp2 + (star1 - u1) * tmp3);
}
