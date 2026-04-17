#include "gks_subcell_blend_driver.h"

#include <omp.h>

int main()
{
	omp_set_num_threads(1);
	riemann_problem_1d_gks_subcell();
	return 0;
}
