
//说明指南:
//这是一个气体动理学格式的求解器，支持二维单块结构化网格模拟，目的是方便有需求的科研工作者入门之用。
//程序基本采用面向过程编写，没有什么弯弯绕，
//故建议从main函数 顺序 看起，搭配注释学习。注释采用土洋结合，以简洁方便的传递作者意图为原则。
//
//如用作学术用途，请引用任意下列文献：
//（1）X.JI, F.ZHAO, W.SHYY, & K.XU(2018).
//A family of high - order gas - kinetic schemes and its comparison with Riemann solver based high - order methods.
//Journal of Computational Physics, 356, 150 - 173.
//（2）X.JI, & K.XU(2020).Performance Enhancement for High - order Gas - kinetic Scheme Based on WENO - adaptive - order Reconstruction.
//Communication in Computational Physics, 28, 2, 539 - 590


#include<omp.h>  //omp并行，尖括号<>代表外部库文件
#include"output.h" //双引号代表该项目中的头文件
#include "gks3rd_demo.h"
#include "gks_fr_driver.h"
#include "gks_subcell_blend_driver.h"
//以下是几个测试算例的头文件
#include"accuracy_test.h"
#include"boundary_layer.h"
#include"cylinder.h"
#include"riemann_problem.h"
//end

using namespace std; //默认std的命名规则

//程序运行开始，手动输入omp并行的线程数
//如果在函数名上面注释，那么函数在别处调用时，可以看到注释
void Set_omp_thread() 
{
	//在小于等于cpu核数的时候，n线程可以理解为n核并行
	int num_thread;
	cout << "please_specify threads number for omp parallel:  ";
	cin >> num_thread;
	omp_set_num_threads(num_thread);
}

int main()
{
	Set_omp_thread(); //设置omp并行线程
	//make_directory_for_result(); //兼容linux的结果文件夹的临时方案
	//accuracy_sinwave_1d(); //一维精度测试，通过周期线性的正弦波传播算例，测试格式在光滑（smooth flow）无粘工况下的精度
	//accuracy_sinwave_1d_gks3rd();//一维FV-单步3阶GKS精度测试
	//accuracy_sinwave_1d_gksfr();//一维FR-单步3阶GKS精度测试
	//accuracy_sinwave_1d_gks_subcell();//一维子单元限制器精度测试
	//riemann_problem_1d_gksfr();//一维FR-单步3阶GKS黎曼问题测试（由于未做子单元限制器，产生振荡基本跑不通）
	//riemann_problem_1d_gks_subcell();//一维子单元限制器单步3阶GKS黎曼问题测试
	//riemann_problem_1d_gks3rd();//一维FV-单步3阶GKS黎曼问题测试
	//riemann_problem_1d();  //一维黎曼问题测试，测试格式对于可压缩间断问题的鲁棒性和分辨率

	//accuracy_sinwave_2d_gks3rd();//二维FV-单步3阶GKS精度测试
	//accuracy_sinwave_2d_gksfr();//二维FR-单步3阶GKS光滑精度测试
	//riemann_problem_2d_gks3rd();//二维FV-单步3阶GKS黎曼问题测试riemann
    //accuracy_sinwave_2d(); //二维精度测试，通过二维周期线性的正弦波传播算例，测试格式在光滑（smooth flow）无粘工况下的精度
	//riemann_problem_2d();  //二维黎曼问题测试，测试格式对于多维可压缩间断问题的鲁棒性和分辨率
	//riemann_problem_2d_gks_subcell();//二维子单元限制器Riemann问题测试
	//double_mach_reflection_2d_gks_subcell();//二维子单元限制器双马赫反射问题测试
	detonation_shock_diffraction_2d_gks_subcell();//二维缺失矩形区域 Detonation Shock Diffraction 测试


	//boundary_layer();      //二维边界层问题测试，测试格式对于光滑（smooth flow）粘性边界的分辨率，非均匀直角网格
	//cylinder();              //二维超音速圆柱绕流问题测试，测试格式对于间断问题和非均匀网格的鲁棒性
    return 0;
}
