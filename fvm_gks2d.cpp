#include"fvm_gks2d.h"
#include"output.h"

int gausspoint = 0; //initilization
double* gauss_loc = new double[gausspoint]; //initilization
double* gauss_weight = new double[gausspoint]; //initilization
GKS2d_type gks2dsolver = kfvs1st_2d; //initilization

Reconstruction_within_Cell_2D_normal cellreconstruction_2D_normal = Vanleer_normal; //initilization
Reconstruction_within_Cell_2D_tangent cellreconstruction_2D_tangent = Vanleer_tangent; //initilization
Reconstruction_forG0_2D_normal g0reconstruction_2D_normal = Center_do_nothing_normal; //initilization
Reconstruction_forG0_2D_tangent g0reconstruction_2D_tangent = Center_all_collision_multi; //initilization

Flux_function_2d flux_function_2d = GKS2D_smooth; //initilization
TimeMarchingCoefficient_2d timecoe_list_2d = S1O1_2D; //initilization




void SetGuassPoint()
{
	if (gausspoint == 0)
	{
		cout << "no gausspoint specify" << endl;
		exit(0);
	}
	if (gausspoint == 1)
	{
		gauss_loc = new double[1]; gauss_loc[0] = 0.0;
		gauss_weight = new double[1]; gauss_weight[0] = 1.0;
	}
	if (gausspoint == 2)
	{
		gauss_loc = new double[2]; gauss_loc[0] = -sqrt(1.0 / 3.0); gauss_loc[1] = -gauss_loc[0];
		gauss_weight = new double[2]; gauss_weight[0] = 0.5; gauss_weight[1] = 0.5;
	}
	if (gausspoint == 3)
	{
		cout << "the usage of 3 gausspoint may be imcompatiable for some reconstructions" << endl;
		//exit(0);
		gauss_loc = new double[3]; gauss_loc[0] = 0.0; gauss_loc[1] = sqrt(3.0 / 5.0); gauss_loc[2] = -gauss_loc[1];
		gauss_weight = new double[3];
		gauss_weight[0] = 4.0 / 9.0; gauss_weight[1] = 0.5 - 0.5 * gauss_weight[0]; gauss_weight[2] = gauss_weight[1];
	}
	if (gausspoint == 4)
	{
		cout << "the usage of 4 gausspoint may be imcompatiable for some reconstructions" << endl;
		gauss_loc = new double[4]; gauss_weight = new double[4];

		gauss_loc[0] = -1.0; gauss_loc[1] = -gauss_loc[0];
		gauss_loc[2] = -sqrt(5) / 5.0; gauss_loc[3] = -gauss_loc[2];
		gauss_weight[0] = 1.0 / 12.0;		gauss_weight[1] = 1.0 / 12.0;
		gauss_weight[2] = 5.0 / 12.0;		gauss_weight[3] = 5.0 / 12.0;

	}
	cout << gausspoint << " gausspoint(s) specify" << endl;
}
double Get_CFL(Block2d& block, Fluid2d* fluids, double tstop)
{
	double dt = fluids[block.ghost * block.ny + block.ghost].dx;

	for (int i = block.ghost; i < block.nx - block.ghost; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost; j++)
		{
			dt = Dtx(dt, fluids[i * block.ny + j].dx, block.CFL, fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[1],
				fluids[i * block.ny + j].primvar[2], fluids[i * block.ny + j].primvar[3]);
			dt = Dtx(dt, fluids[i * block.ny + j].dy, block.CFL, fluids[i * block.ny + j].primvar[0], fluids[i * block.ny + j].primvar[1],
				fluids[i * block.ny + j].primvar[2], fluids[i * block.ny + j].primvar[3]);
		}
	}

	if (block.t + dt > tstop)
	{
		dt = tstop - block.t + 1e-16;
		cout << "last step here" << endl;
	}
	//print time step information
	if (block.step % 100 == 0)
	{
		cout << "step= " << block.step
			<< "time size is " << dt
			<< " time= " << block.t
			<< endl;
	}
	return dt;
}



void A_point(double* a, double* der, double* prim)
{
	double R4, R3, R2;
	double overden = 1.0 / prim[0];
	R4 = der[3] * overden - 0.5 * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]) * der[0] * overden;
	R3 = (der[2] - prim[2] * der[0]) * overden;
	R2 = (der[1] - prim[1] * der[0]) * overden;
	a[3] = (4.0 / (K + 2)) * prim[3] * prim[3] * (2 * R4 - 2 * prim[1] * R2 - 2 * prim[2] * R3);
	a[2] = 2 * prim[3] * R3 - prim[2] * a[3];
	a[1] = 2 * prim[3] * R2 - prim[1] * a[3];
	a[0] = der[0] * overden - prim[1] * a[1] - prim[2] * a[2] - 0.5 * a[3] * (prim[1] * prim[1] + prim[2] * prim[2] + 0.5 * (K + 2) / prim[3]);
}
void G_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{

	psi[0] = a[0] * m.uvxi[no_u][no_v][no_xi] + a[1] * m.uvxi[no_u + 1][no_v][no_xi] + a[2] * m.uvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 2][no_v][no_xi] + m.uvxi[no_u][no_v + 2][no_xi] + m.uvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.uvxi[no_u + 1][no_v][no_xi] + a[1] * m.uvxi[no_u + 2][no_v][no_xi] + a[2] * m.uvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 3][no_v][no_xi] + m.uvxi[no_u + 1][no_v + 2][no_xi] + m.uvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.uvxi[no_u][no_v + 1][no_xi] + a[1] * m.uvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.uvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.uvxi[no_u + 2][no_v + 1][no_xi] + m.uvxi[no_u][no_v + 3][no_xi] + m.uvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.uvxi[no_u + 2][no_v][no_xi] + m.uvxi[no_u][no_v + 2][no_xi] + m.uvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.uvxi[no_u + 3][no_v][no_xi] + m.uvxi[no_u + 1][no_v + 2][no_xi] + m.uvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.uvxi[no_u + 2][no_v + 1][no_xi] + m.uvxi[no_u][no_v + 3][no_xi] + m.uvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.uvxi[no_u + 4][no_v][no_xi] + m.uvxi[no_u][no_v + 4][no_xi] + m.uvxi[no_u][no_v][no_xi + 2] + 2 * m.uvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.uvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.uvxi[no_u][no_v + 2][no_xi + 1]));
}
void GL_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{

	psi[0] = a[0] * m.upvxi[no_u][no_v][no_xi] + a[1] * m.upvxi[no_u + 1][no_v][no_xi] + a[2] * m.upvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 2][no_v][no_xi] + m.upvxi[no_u][no_v + 2][no_xi] + m.upvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.upvxi[no_u + 1][no_v][no_xi] + a[1] * m.upvxi[no_u + 2][no_v][no_xi] + a[2] * m.upvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 3][no_v][no_xi] + m.upvxi[no_u + 1][no_v + 2][no_xi] + m.upvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.upvxi[no_u][no_v + 1][no_xi] + a[1] * m.upvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.upvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.upvxi[no_u + 2][no_v + 1][no_xi] + m.upvxi[no_u][no_v + 3][no_xi] + m.upvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.upvxi[no_u + 2][no_v][no_xi] + m.upvxi[no_u][no_v + 2][no_xi] + m.upvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.upvxi[no_u + 3][no_v][no_xi] + m.upvxi[no_u + 1][no_v + 2][no_xi] + m.upvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.upvxi[no_u + 2][no_v + 1][no_xi] + m.upvxi[no_u][no_v + 3][no_xi] + m.upvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.upvxi[no_u + 4][no_v][no_xi] + m.upvxi[no_u][no_v + 4][no_xi] + m.upvxi[no_u][no_v][no_xi + 2] + 2 * m.upvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.upvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.upvxi[no_u][no_v + 2][no_xi + 1]));
}
void GR_address(int no_u, int no_v, int no_xi, double* psi, double a[4], MMDF& m)
{
	// Similar to the GL_address
	psi[0] = a[0] * m.unvxi[no_u][no_v][no_xi] + a[1] * m.unvxi[no_u + 1][no_v][no_xi] + a[2] * m.unvxi[no_u][no_v + 1][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 2][no_v][no_xi] + m.unvxi[no_u][no_v + 2][no_xi] + m.unvxi[no_u][no_v][no_xi + 1]);
	psi[1] = a[0] * m.unvxi[no_u + 1][no_v][no_xi] + a[1] * m.unvxi[no_u + 2][no_v][no_xi] + a[2] * m.unvxi[no_u + 1][no_v + 1][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 3][no_v][no_xi] + m.unvxi[no_u + 1][no_v + 2][no_xi] + m.unvxi[no_u + 1][no_v][no_xi + 1]);
	psi[2] = a[0] * m.unvxi[no_u][no_v + 1][no_xi] + a[1] * m.unvxi[no_u + 1][no_v + 1][no_xi] + a[2] * m.unvxi[no_u][no_v + 2][no_xi] + a[3] * 0.5 * (m.unvxi[no_u + 2][no_v + 1][no_xi] + m.unvxi[no_u][no_v + 3][no_xi] + m.unvxi[no_u][no_v + 1][no_xi + 1]);
	psi[3] = 0.5 * (a[0] * (m.unvxi[no_u + 2][no_v][no_xi] + m.unvxi[no_u][no_v + 2][no_xi] + m.unvxi[no_u][no_v][no_xi + 1]) +
		a[1] * (m.unvxi[no_u + 3][no_v][no_xi] + m.unvxi[no_u + 1][no_v + 2][no_xi] + m.unvxi[no_u + 1][no_v][no_xi + 1]) +
		a[2] * (m.unvxi[no_u + 2][no_v + 1][no_xi] + m.unvxi[no_u][no_v + 3][no_xi] + m.unvxi[no_u][no_v + 1][no_xi + 1]) +
		a[3] * 0.5 * (m.unvxi[no_u + 4][no_v][no_xi] + m.unvxi[no_u][no_v + 4][no_xi] + m.unvxi[no_u][no_v][no_xi + 2] + 2 * m.unvxi[no_u + 2][no_v + 2][no_xi] + 2 * m.unvxi[no_u + 2][no_v][no_xi + 1] + 2 * m.unvxi[no_u][no_v + 2][no_xi + 1]));
}


void free_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[(i + 1) * block.ny + j].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[(i + 1) * block.ny + j].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[(i + 1) * block.ny + j].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[(i + 1) * block.ny + j].convar[3];
		}
	}
}
void free_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[(i - 1) * block.ny + j].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[(i - 1) * block.ny + j].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[(i - 1) * block.ny + j].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[(i - 1) * block.ny + j].convar[3];
		}

	}

}
void free_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ghost - 1; j >= 0; j--)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[i * block.ny + j + 1].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[i * block.ny + j + 1].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[i * block.ny + j + 1].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[i * block.ny + j + 1].convar[3];
		}
	}

}
void free_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			fluids[i * block.ny + j].convar[0] = fluids[i * block.ny + j - 1].convar[0];
			fluids[i * block.ny + j].convar[1] = fluids[i * block.ny + j - 1].convar[1];
			fluids[i * block.ny + j].convar[2] = fluids[i * block.ny + j - 1].convar[2];
			fluids[i * block.ny + j].convar[3] = fluids[i * block.ny + j - 1].convar[3];
		}
	}

}

void periodic_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		int tar = i + block.nodex;
		for (int j = 0; j < block.ny; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[tar * block.ny + j].convar[var];
			}
		}
	}

}
void periodic_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		int tar = i - block.nodex;
		for (int j = 0; j < block.ny; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[tar * block.ny + j].convar[var];
			}
		}
	}

}
void periodic_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ghost - 1; j >= 0; j--)
	{
		int tar = j + block.nodey;
		for (int i = 0; i < block.nx; i++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[i * block.ny + tar].convar[var];
			}
		}
	}

}
void periodic_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		int tar = j - block.nodey;
		for (int i = 0; i < block.nx; i++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar[var] = fluids[i * block.ny + tar].convar[var];
			}
		}
	}

}
void noslip_adiabatic_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;
	for (int i = order - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;
			int ref = (2 * order - 1 - i) * block.ny + j;
			fluids[index].convar[0] = fluids[ref].convar[0];
			fluids[index].convar[1] = 2 * fluids[index].convar[0] * bcvalue.primvar[1] - fluids[ref].convar[1];
			fluids[index].convar[2] = 2 * fluids[index].convar[0] * bcvalue.primvar[2] - fluids[ref].convar[2];
			fluids[index].convar[3] = fluids[ref].convar[3];

		}
	}

}

void noslip_adiabatic_boundary_right(Fluid2d* fluids, Interface2d* face, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;

	for (int i = block.nx - order; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;
			int ref = (2 * block.nx - 2 * order - 1 - i) * block.ny + j;
			int face_index = (2 * block.nx - 2 * order - 1 - i) * (block.ny + 1) + j;
			noslip_adiabatic_boundary(fluids[index].convar, fluids[ref].convar, face[face_index].normal, &bcvalue.primvar[1]);

		}
	}
}

void noslip_adiabatic_boundary(double* w_ghost, double* w_inner, double* normal, double* local_wall_velocity)
{
	double inner[4];
	Copy_Array(inner, w_inner, 4);
	Global_to_Local(inner, normal);
	inner[1] = 2 * local_wall_velocity[0] - inner[1];
	inner[2] = 2 * local_wall_velocity[1] - inner[2];
	Local_to_Global(inner, normal);
	Copy_Array(w_ghost, inner, 4);
}

void noslip_adiabatic_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;

	for (int j = order - 1; j >= 0; j--)
	{
		for (int i = 0; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			int ref = (i)*block.ny + 2 * order - 1 - j;
			fluids[index].convar[0] = fluids[ref].convar[0];
			fluids[index].convar[1] = 2 * fluids[index].convar[0] * bcvalue.primvar[1] - fluids[ref].convar[1];
			fluids[index].convar[2] = 2 * fluids[index].convar[0] * bcvalue.primvar[2] - fluids[ref].convar[2];
			fluids[index].convar[3] = fluids[ref].convar[3];

		}
	}

}
void noslip_adiabatic_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;

	for (int j = block.ny - order; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			int ref = (i)*block.ny + 2 * block.ny - 2 * order - 1 - j;
			fluids[index].convar[0] = fluids[ref].convar[0];
			fluids[index].convar[1] = 2 * fluids[index].convar[0] * bcvalue.primvar[1] - fluids[ref].convar[1];
			fluids[index].convar[2] = 2 * fluids[index].convar[0] * bcvalue.primvar[2] - fluids[ref].convar[2];
			fluids[index].convar[3] = fluids[ref].convar[3];

		}
	}

}
void reflection_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;
	for (int j = block.ny - order; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			int ref = (i)*block.ny + 2 * block.ny - 2 * order - 1 - j;

			fluids[index].convar[0] = fluids[ref].convar[0];
			fluids[index].convar[1] = fluids[ref].convar[1];
			fluids[index].convar[2] = -fluids[ref].convar[2];
			fluids[index].convar[3] = fluids[ref].convar[3];

		}

	}

}
void reflection_boundary_down(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;

	for (int j = order - 1; j >= 0; j--)
	{
		for (int i = 0; i < block.nx; i++)
		{
			int index = i * block.ny + j;
			int ref = (i)*block.ny + 2 * order - 1 - j;

			fluids[index].convar[0] = fluids[ref].convar[0];
			fluids[index].convar[1] = fluids[ref].convar[1];
			fluids[index].convar[2] = -fluids[ref].convar[2];
			fluids[index].convar[3] = fluids[ref].convar[3];

		}
	}
}
void reflection_boundary_right(Fluid2d* fluids, Interface2d* face, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;

	for (int i = block.nx - order; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;
			int ref = (2 * block.nx - 2 * order - 1 - i) * block.ny + j;
			int face_index = (2 * block.nx - 2 * order - 1 - i) * (block.ny + 1) + j;
			reflection_boundary(fluids[index].convar, fluids[ref].convar, face[face_index].normal);

		}
	}
}
void reflection_boundary_left(Fluid2d* fluids, Interface2d* face, Block2d block, Fluid2d bcvalue)
{
	int order = block.ghost;
	for (int i = order - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int index = i * block.ny + j;
			int ref = (2 * order - 1 - i) * block.ny + j;
			int face_index = order * (block.ny + 1) + j;

			reflection_boundary(fluids[index].convar, fluids[ref].convar, face[face_index].normal);

		}
	}
}

void reflection_boundary(double* w_ghost, double* w_inner, double* normal)
{
	double inner[4];
	Copy_Array(inner, w_inner, 4);
	Global_to_Local(inner, normal);
	inner[1] = -inner[1];
	Local_to_Global(inner, normal);
	Copy_Array(w_ghost, inner, 4);
}



void inflow_boundary_left(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.ghost - 1; i >= 0; i--)
	{
		for (int j = 0; j < block.ny; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].primvar[var] = bcvalue.primvar[var];
			}

			Primvar_to_convar_2D(fluids[i * block.ny + j].convar, fluids[i * block.ny + j].primvar);
		}
	}

}

void inflow_boundary_right(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int i = block.nx - block.ghost; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].primvar[var] = bcvalue.primvar[var];
			}
			Primvar_to_convar_2D(fluids[i * block.ny + j].convar, fluids[i * block.ny + j].primvar);
		}

	}

}


void inflow_boundary_up(Fluid2d* fluids, Block2d block, Fluid2d bcvalue)
{
	for (int j = block.ny - block.ghost; j < block.ny; j++)
	{
		for (int i = 0; i < block.nx; i++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].primvar[var] = bcvalue.primvar[var];
			}
			Primvar_to_convar_2D(fluids[i * block.ny + j].convar, fluids[i * block.ny + j].primvar);
		}
	}

}

void ICfor_uniform_2d(Fluid2d* fluid, double* prim, Block2d block)
{

	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			fluid[i * block.ny + j].primvar[0] = prim[0];
			fluid[i * block.ny + j].primvar[1] = prim[1];
			fluid[i * block.ny + j].primvar[2] = prim[2];
			fluid[i * block.ny + j].primvar[3] = prim[3];
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

bool negative_density_or_pressure(double* primvar)
{
	//detect whether density or pressure is negative
	if (primvar[0] < 0 || primvar[3] < 0 ||
		isnan(primvar[0])
		|| isnan(primvar[3]))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Convar_to_Primvar(Fluid2d* fluids, Block2d block)
{
	bool all_positive = true;
#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			Convar_to_primvar_2D(fluids[i * block.ny + j].primvar, fluids[i * block.ny + j].convar);
			if (negative_density_or_pressure(fluids[i * block.ny + j].primvar))
			{
				all_positive = false;
			}
		}
	}

	if (all_positive == false)
	{
		cout << "the program blows up at t=" << block.t << "!" << endl;
		output2d_blow_up(fluids, block);
		exit(0);
	}
}
void Convar_to_primvar(Fluid2d* fluid, Block2d block)
{
#pragma omp parallel  for
	for (int i = 0; i < block.nx * block.ny; i++)
	{
		Convar_to_primvar_2D(fluid[i].primvar, fluid[i].convar);
	}
}
void YchangetoX(double* fluidtmp, double* fluid)
{
	fluidtmp[0] = fluid[0];
	fluidtmp[1] = fluid[2];
	fluidtmp[2] = -fluid[1];
	fluidtmp[3] = fluid[3];
}
void XchangetoY(double* fluidtmp, double* fluid)
{
	fluidtmp[0] = fluid[0];
	fluidtmp[1] = -fluid[2];
	fluidtmp[2] = fluid[1];
	fluidtmp[3] = fluid[3];
}
void CopyFluid_new_to_old(Fluid2d* fluids, Block2d block)
{
#pragma omp parallel  for
	for (int i = block.ghost; i < block.ghost + block.nodex; i++)
	{
		for (int j = block.ghost; j < block.ghost + block.nodey; j++)
		{
			for (int var = 0; var < 4; var++)
			{
				fluids[i * block.ny + j].convar_old[var] = fluids[i * block.ny + j].convar[var];
			}
		}
	}
}




void Reconstruction_within_cell(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			cellreconstruction_2D_normal
			(xinterfaces[i * (block.ny + 1) + j], xinterfaces[(i + 1) * (block.ny + 1) + j],
				yinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost - 1; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost - 1; j < block.ny - block.ghost + 1; j++)
		{
			cellreconstruction_2D_tangent
			(&xinterfaces[i * (block.ny + 1) + j], &xinterfaces[(i + 1) * (block.ny + 1) + j],
				&yinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j + 1], &fluids[i * block.ny + j], block);
		}
	}
}

void First_order_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	first_order(left.line.right, right.line.left,
		left.normal, right.normal, fluids[0].convar);

	first_order(down.line.right, up.line.left,
		down.normal, up.normal, fluids[0].convar);
}

void first_order(Point2d& left, Point2d& right, double* normal_l, double* normal_r, double* w)
{
	double splus[4], sminus[4];

	double w_l[4];
	Global_to_Local(w_l, w, normal_l);
	Copy_Array(left.convar, w_l, 4);
	Array_zero(left.der1x, 4);

	double w_r[4];
	Global_to_Local(w_r, w, normal_r);
	Copy_Array(right.convar, w_r, 4);
	Array_zero(right.der1x, 4);

}

void Vanleer_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			Vanleer(left.line.right, right.line.left, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, block.dx);
		}
		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double wn1tmp[4], wtmp[4], wp1tmp[4];
			YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);

			Vanleer(down.line.right, up.line.left, wn1tmp, wtmp, wp1tmp, block.dy);
		}
	}
	else
	{

		double h = fluids[0].dx;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			Vanleer_non_uniform(left.line.right, right.line.left,
				left.normal, right.normal, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, h);
		}
		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			h = fluids[0].dy;
			Vanleer_non_uniform(down.line.right, up.line.left,
				down.normal, up.normal, fluids[-1].convar, fluids[0].convar, fluids[1].convar, h);
		}

	}
}
void Vanleer(Point2d& left, Point2d& right, double* wn1, double* w, double* wp1, double h)
{
	double splus[4], sminus[4];

	for (int i = 0; i < 4; i++)
	{
		splus[i] = (wp1[i] - w[i]) / h;
		sminus[i] = (w[i] - wn1[i]) / h;
	}

	for (int i = 0; i < 4; i++)
	{
		if ((splus[i] * sminus[i]) > 0)
		{
			left.der1x[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
			right.der1x[i] = left.der1x[i];
		}

		else
		{
			left.der1x[i] = 0.0;
			right.der1x[i] = 0.0;
		}
		left.convar[i] = w[i] - 0.5 * h * left.der1x[i];
		right.convar[i] = w[i] + 0.5 * h * right.der1x[i];
	}

	//if lambda <0, then reduce to the first order
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);
	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);

	if (left.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
		{
			cout << "order reduce on left interface x= " << left.x << " y= " << left.y << endl;
		}
		for (int m = 0; m < 4; m++)
		{
			left.convar[m] = wn1[m];
		}
		left.is_reduce_order = true;
	}
	if (right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
		{
			cout << "order reduce on right interface x= " << right.x << " y= " << right.y << endl;
		}
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = wp1[m];
		}
		right.is_reduce_order = true;
	}

}

void WENO3_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{

	//for non-uniform mesh.
	Point2d voidpoint;
	if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
	{
		double dx = fluids[0].dx;
		double normal[2];
		double  wn1tmp[4], wtmp[4], wp1tmp[4];

		//cell left reconstruction
		Copy_Array(normal, left.normal, 2);

		Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
		Global_to_Local(wtmp, fluids[0].convar, normal);
		Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);

		WENO3(left.line.right, voidpoint, wn1tmp, wtmp, wp1tmp, dx);
		//cell right reconstruction
		Copy_Array(normal, right.normal, 2);

		Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
		Global_to_Local(wtmp, fluids[0].convar, normal);
		Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);

		WENO3(voidpoint, right.line.left, wn1tmp, wtmp, wp1tmp, dx);

	}


	if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
	{
		double dy = fluids[0].dy;
		double normal[2];
		double  wn1tmp[4], wtmp[4], wp1tmp[4];

		//down interface reconstruction
		Copy_Array(normal, down.normal, 2);

		Global_to_Local(wn1tmp, fluids[-1].convar, normal);
		Global_to_Local(wtmp, fluids[0].convar, normal);
		Global_to_Local(wp1tmp, fluids[1].convar, normal);

		WENO3(down.line.right, voidpoint, wn1tmp, wtmp, wp1tmp, dy);

		//up interface reconstruction
		Copy_Array(normal, up.normal, 2);

		Global_to_Local(wn1tmp, fluids[-1].convar, normal);
		Global_to_Local(wtmp, fluids[0].convar, normal);
		Global_to_Local(wp1tmp, fluids[1].convar, normal);

		WENO3(voidpoint, up.line.left, wn1tmp, wtmp, wp1tmp, dy);
	}

}

void WENO3(Point2d& left, Point2d& right, double* wn1, double* w, double* wp1, double h)
{
	//we denote that   |left...cell-center...right|
	double ren1[4], re0[4], rep1[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{

			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];

		}
	}
	else if (reconstruction_variable == characteristic)
	{

		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);

	}


	for (int i = 0; i < 4; i++)
	{
		weno_3rd_left(var[i], der1[i], der2[i], ren1[i], re0[i], rep1[i], h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];

		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
		}
	}
	else
	{
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_3rd_right(var[i], der1[i], der2[i], ren1[i], re0[i], rep1[i], h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];

		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " WENO5-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
		}
	}

}

void Vanleer_non_uniform(Point2d& left, Point2d& right, double* normal_l, double* normal_r, double* wn1, double* w, double* wp1, double h)
{
	double splus[4], sminus[4];

	double wn1_l[4], w_l[4], wp1_l[4];
	Global_to_Local(wn1_l, wn1, normal_l);
	Global_to_Local(w_l, w, normal_l);
	Global_to_Local(wp1_l, wp1, normal_l);

	for (int i = 0; i < 4; i++)
	{
		splus[i] = (wp1_l[i] - w_l[i]) / h;
		sminus[i] = (w_l[i] - wn1_l[i]) / h;
	}

	for (int i = 0; i < 4; i++)
	{
		if ((splus[i] * sminus[i]) > 0)
		{
			left.der1x[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
		}
		else
		{
			left.der1x[i] = 0.0;
		}
		left.convar[i] = w_l[i] - 0.5 * h * left.der1x[i];

	}

	double wn1_r[4], w_r[4], wp1_r[4];
	Global_to_Local(wn1_r, wn1, normal_r);
	Global_to_Local(w_r, w, normal_r);
	Global_to_Local(wp1_r, wp1, normal_r);


	for (int i = 0; i < 4; i++)
	{
		splus[i] = (wp1_r[i] - w_r[i]) / h;
		sminus[i] = (w_r[i] - wn1_r[i]) / h;
	}

	for (int i = 0; i < 4; i++)
	{
		if ((splus[i] * sminus[i]) > 0)
		{
			right.der1x[i] = 2 * splus[i] * sminus[i] / (splus[i] + sminus[i]);
		}
		else
		{
			right.der1x[i] = 0.0;
		}
		right.convar[i] = w_r[i] + 0.5 * h * right.der1x[i];

	}

	//if lambda <0, then reduce to the first order
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);
	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);

	if (left.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
		{
			cout << "order reduce on left interface x= " << left.x << " y= " << left.y << endl;
		}
		for (int m = 0; m < 4; m++)
		{
			left.convar[m] = w_l[m];
		}
		left.is_reduce_order = true;
	}
	if (right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
		{
			cout << "order reduce on right interface x= " << right.x << " y= " << right.y << endl;
		}
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w_r[m];
		}
		right.is_reduce_order = true;
	}

}

void WENO5_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			WENO(left.line.right, right.line.left,
				fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[2 * block.ny].convar
				, fluids[0].dx);
		}


		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];
			YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wn2tmp, fluids[-2].convar); YchangetoX(wp2tmp, fluids[2].convar);

			WENO(down.line.right, up.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, fluids[0].dy);
		}
	}
	else
	{


		//for non-uniform mesh.
		Point2d voidpoint;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			double dx = fluids[0].dx;
			double normal[2];
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];

			//cell left reconstruction
			Copy_Array(normal, left.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			WENO(left.line.right, voidpoint, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dx);
			//cell right reconstruction
			Copy_Array(normal, right.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			WENO(voidpoint, right.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dx);

		}


		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double dy = fluids[0].dy;
			double normal[2];
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];

			//down interface reconstruction
			Copy_Array(normal, down.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			WENO(down.line.right, voidpoint, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dy);

			//up interface reconstruction
			Copy_Array(normal, up.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			WENO(voidpoint, up.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dy);
		}
	}

}
void WENO(Point2d& left, Point2d& right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	//we denote that   |left...cell-center...right|
	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		WENO5_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
			left.der2xx[i] = der2[i];
		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
		Char_to_convar(left.der2xx, base_left, der2);
	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		WENO5_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
			right.der2xx[i] = der2[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
		Char_to_convar(right.der2xx, base_right, der2);
	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " WENO5-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
			right.der2xx[m] = 0.0;
			left.der2xx[m] = 0.0;
		}
	}
}
void WENO5_left(double& var, double& der1, double& der2, double wn2, double wn1, double w, double wp1, double wp2, double h)
{
	double qleft[3];
	double qright[3];
	double dleft[3];
	double dright[3];

	qleft[0] = 11.0 / 6.0 * w - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	qleft[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w - 1.0 / 6.0 * wp1;
	qleft[2] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w;


	dleft[0] = 0.1;
	dleft[1] = 0.6;
	dleft[2] = 0.3;

	qright[0] = 1.0 / 3.0 * w + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	qright[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w + 1.0 / 3.0 * wp1;
	qright[2] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w;

	dright[0] = 0.3;
	dright[1] = 0.6;
	dright[2] = 0.1;

	double 	beta[3];

	beta[0] = 13.0 / 12.0 * pow((w - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w - 4 * wp1 + wp2), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w), 2);

	double epsilon = 1e-8;
	double alphaleft[3];	double alpharight[3];
	if (wenotype == wenojs)
	{
		for (int i = 0; i < 3; i++)
		{
			alphaleft[i] = dleft[i] / ((beta[i] + epsilon) * (beta[i] + epsilon));
			alpharight[i] = dright[i] / ((beta[i] + epsilon) * (beta[i] + epsilon));
		}
	}
	if (wenotype == wenoz)
	{
		double tau5 = abs(beta[0] - beta[2]);
		for (int i = 0; i < 3; i++)
		{
			alphaleft[i] = dleft[i] * (1 + tau5 / (beta[i] + epsilon));
			alpharight[i] = dright[i] * (1 + tau5 / (beta[i] + epsilon));
		}
	}

	if (wenotype == linear)
	{
		for (int k = 0; k < 3; k++)
		{
			alphaleft[k] = dleft[k];
			alpharight[k] = dright[k];
		}
	}

	double alphal = alphaleft[0] + alphaleft[1] + alphaleft[2];

	double omegaleft[3];
	omegaleft[0] = alphaleft[0] / alphal;
	omegaleft[1] = alphaleft[1] / alphal;
	omegaleft[2] = alphaleft[2] / alphal;

	double left = omegaleft[0] * qleft[0] + omegaleft[1] * qleft[1] + omegaleft[2] * qleft[2];



	double alphar = alpharight[0] + alpharight[1] + alpharight[2];

	double omegaright[3];

	omegaright[0] = alpharight[0] / alphar;
	omegaright[1] = alpharight[1] / alphar;
	omegaright[2] = alpharight[2] / alphar;

	double right = omegaright[0] * qright[0] + omegaright[1] * qright[1] + omegaright[2] * qright[2];
	var = left;
	der2 = 6 * (left + right - 2 * w) / (h * h);
	der1 = (right - left) / h - 0.5 * h * der2;

}
void WENO5_right(double& var, double& der1, double& der2, double wn2, double wn1, double w, double wp1, double wp2, double h)
{
	double qleft[3];
	double qright[3];
	double dleft[3];
	double dright[3];

	qleft[0] = 11.0 / 6.0 * w - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	qleft[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w - 1.0 / 6.0 * wp1;
	qleft[2] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w;


	dleft[0] = 0.1;
	dleft[1] = 0.6;
	dleft[2] = 0.3;

	qright[0] = 1.0 / 3.0 * w + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	qright[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w + 1.0 / 3.0 * wp1;
	qright[2] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w;

	dright[0] = 0.3;
	dright[1] = 0.6;
	dright[2] = 0.1;

	double 	beta[3];

	beta[0] = 13.0 / 12.0 * pow((w - 2 * wp1 + wp2), 2) + 0.25 * pow((3 * w - 4 * wp1 + wp2), 2);
	beta[1] = 13.0 / 12.0 * pow((wn1 - 2 * w + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
	beta[2] = 13.0 / 12.0 * pow((wn2 - 2 * wn1 + w), 2) + 0.25 * pow((wn2 - 4 * wn1 + 3 * w), 2);

	double epsilon = 1e-8;
	double alphaleft[3];	double alpharight[3];
	if (wenotype == wenojs)
	{
		for (int i = 0; i < 3; i++)
		{
			alphaleft[i] = dleft[i] / ((beta[i] + epsilon) * (beta[i] + epsilon));
			alpharight[i] = dright[i] / ((beta[i] + epsilon) * (beta[i] + epsilon));
		}
	}
	if (wenotype == wenoz)
	{
		double tau5 = abs(beta[0] - beta[2]);
		for (int i = 0; i < 3; i++)
		{
			alphaleft[i] = dleft[i] * (1 + tau5 / (beta[i] + epsilon));
			alpharight[i] = dright[i] * (1 + tau5 / (beta[i] + epsilon));
		}
	}

	if (wenotype == linear)
	{
		for (int k = 0; k < 3; k++)
		{
			alphaleft[k] = dleft[k];
			alpharight[k] = dright[k];
		}
	}

	double alphal = alphaleft[0] + alphaleft[1] + alphaleft[2];

	double omegaleft[3];
	omegaleft[0] = alphaleft[0] / alphal;
	omegaleft[1] = alphaleft[1] / alphal;
	omegaleft[2] = alphaleft[2] / alphal;

	double left = omegaleft[0] * qleft[0] + omegaleft[1] * qleft[1] + omegaleft[2] * qleft[2];



	double alphar = alpharight[0] + alpharight[1] + alpharight[2];

	double omegaright[3];

	omegaright[0] = alpharight[0] / alphar;
	omegaright[1] = alpharight[1] / alphar;
	omegaright[2] = alpharight[2] / alphar;

	double right = omegaright[0] * qright[0] + omegaright[1] * qright[1] + omegaright[2] * qright[2];
	var = right;
	der2 = 6 * (left + right - 2 * w) / (h * h);
	der1 = (right - left) / h + 0.5 * h * der2;

}

void WENO5_AO_normal(Interface2d& left, Interface2d& right, Interface2d& down, Interface2d& up, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			WENO5_AO(left.line.right, right.line.left, fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[2 * block.ny].convar, fluids[0].dx);
		}


		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];
			YchangetoX(wn1tmp, fluids[-1].convar); YchangetoX(wtmp, fluids[0].convar); YchangetoX(wp1tmp, fluids[1].convar);
			YchangetoX(wn2tmp, fluids[-2].convar); YchangetoX(wp2tmp, fluids[2].convar);

			WENO5_AO(down.line.right, up.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, fluids[0].dy);
		}
	}
	else
	{
		//for non-uniform mesh.
		Point2d voidpoint;
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			double dx = fluids[0].dx;
			double normal[2];
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];

			//cell left reconstruction
			Copy_Array(normal, left.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			WENO5_AO(left.line.right, voidpoint, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dx);
			//cell right reconstruction
			Copy_Array(normal, right.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2 * block.ny].convar, normal);
			Global_to_Local(wn1tmp, fluids[-block.ny].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[block.ny].convar, normal);
			Global_to_Local(wp2tmp, fluids[2 * block.ny].convar, normal);
			WENO5_AO(voidpoint, right.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dx);

		}


		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double dy = fluids[0].dy;
			double normal[2];
			double wn2tmp[4], wn1tmp[4], wtmp[4], wp1tmp[4], wp2tmp[4];

			//down interface reconstruction
			Copy_Array(normal, down.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			WENO5_AO(down.line.right, voidpoint, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dy);

			//up interface reconstruction
			Copy_Array(normal, up.normal, 2);
			Global_to_Local(wn2tmp, fluids[-2].convar, normal);
			Global_to_Local(wn1tmp, fluids[-1].convar, normal);
			Global_to_Local(wtmp, fluids[0].convar, normal);
			Global_to_Local(wp1tmp, fluids[1].convar, normal);
			Global_to_Local(wp2tmp, fluids[2].convar, normal);
			WENO5_AO(voidpoint, up.line.left, wn2tmp, wn1tmp, wtmp, wp1tmp, wp2tmp, dy);
		}
	}
}
void WENO5_AO(Point2d& left, Point2d& right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	//we denote that   |left...cell-center...right|
	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_5th_ao_left(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);

	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left.convar[i] = var[i];
			left.der1x[i] = der1[i];
			left.der2xx[i] = der2[i];

		}
	}
	else
	{
		Char_to_convar(left.convar, base_left, var);
		Char_to_convar(left.der1x, base_left, der1);
		Char_to_convar(left.der2xx, base_left, der2);

	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		weno_5th_ao_right(var[i], der1[i], der2[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right.convar[i] = var[i];
			right.der1x[i] = der1[i];
			right.der2xx[i] = der2[i];
		}
	}
	else
	{
		Char_to_convar(right.convar, base_right, var);
		Char_to_convar(right.der1x, base_right, der1);
		Char_to_convar(right.der2xx, base_right, der2);

	}

	Check_Order_Reduce_by_Lambda_2D(right.is_reduce_order, right.convar);
	Check_Order_Reduce_by_Lambda_2D(left.is_reduce_order, left.convar);

	if (left.is_reduce_order == true || right.is_reduce_order == true)
	{
		if (is_reduce_order_warning == true)
			cout << " WENO5-cell-splitting order reduce" << endl;
		for (int m = 0; m < 4; m++)
		{
			right.convar[m] = w[m];
			left.convar[m] = w[m];
			right.der1x[m] = 0.0;
			left.der1x[m] = 0.0;
			right.der2xx[m] = 0.0;
			left.der2xx[m] = 0.0;
		}
	}

}


void First_order_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
	{
		first_order_tangent(left[0].gauss[num_gauss].right, left[0].line.right);
		first_order_tangent(right[0].gauss[num_gauss].left, right[0].line.left);
		first_order_tangent(down[0].gauss[num_gauss].right, down[0].line.right);
		first_order_tangent(up[0].gauss[num_gauss].left, up[0].line.left);
	}
}

void first_order_tangent(Point2d& gauss, Point2d& w0)
{
	//tangential
	Copy_Array(gauss.convar, w0.convar, 4);
	Copy_Array(gauss.der1x, w0.der1x, 4);
	Array_zero(gauss.der1y, 4);
}

void Vanleer_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		// along x direction tangitial recontruction,
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			Vanleer(left[0].gauss[num_gauss].right, left[-1].line.right,
				left[0].line.right, left[1].line.right, left[0].length);
			Vanleer(right[0].gauss[num_gauss].left, right[-1].line.left,
				right[0].line.left, right[1].line.left, right[0].length);
		}
		//since we already do the coordinate transform, along y, no transform needed.
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			Vanleer(down[0].gauss[num_gauss].right, down[block.ny + 1].line.right,
				down[0].line.right, down[-block.ny - 1].line.right, down[0].length);
			Vanleer(up[0].gauss[num_gauss].left, up[block.ny + 1].line.left,
				up[0].line.left, up[-block.ny - 1].line.left, up[0].length);
		}
	}
	else
	{
		Interface2d re[2];
		int index[2]; index[0] = -1; index[1] = 1;

		for (int iface = 0; iface < 2; ++iface)
		{
			Local_to_Global(re[iface].line.left.convar, right[index[iface]].line.left.convar,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.left.convar,
				right[0].normal);
			Local_to_Global(re[iface].line.left.der1x, right[index[iface]].line.left.der1x,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.left.der1x,
				right[0].normal);
			Local_to_Global(re[iface].line.right.convar, right[index[iface]].line.right.convar,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.right.convar,
				right[0].normal);
			Local_to_Global(re[iface].line.right.der1x, right[index[iface]].line.right.der1x,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.right.der1x,
				right[0].normal);
		}

		for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
		{
			Vanleer(right[0].gauss[num_gauss].left, re[0].line.left,
				right[0].line.left, re[1].line.left, right[0].length);
			Vanleer(right[0].gauss[num_gauss].right, re[0].line.right,
				right[0].line.right, re[1].line.right, right[0].length);
		}

		index[0] = (block.ny + 1); index[1] = -(block.ny + 1);

		for (int iface = 0; iface < 2; ++iface)
		{
			Local_to_Global(re[iface].line.left.convar, up[index[iface]].line.left.convar,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.left.convar,
				up[0].normal);
			Local_to_Global(re[iface].line.left.der1x, up[index[iface]].line.left.der1x,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.left.der1x,
				up[0].normal);
			Local_to_Global(re[iface].line.right.convar, up[index[iface]].line.right.convar,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.right.convar,
				up[0].normal);
			Local_to_Global(re[iface].line.right.der1x, up[index[iface]].line.right.der1x,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.right.der1x,
				up[0].normal);
		}

		for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
		{
			Vanleer(up[0].gauss[num_gauss].left, re[0].line.left,
				up[0].line.left, re[1].line.left, up[0].length);
			Vanleer(up[0].gauss[num_gauss].right, re[0].line.right,
				up[0].line.right, re[1].line.right, up[0].length);
		}
	}

}
void Vanleer(Point2d& gauss, Point2d& wn1, Point2d& w0, Point2d& wp1, double h)
{
	//tangential
	for (int var = 0; var < 4; var++)
	{
		double coe[2];
		Vanleer_tangential(coe, wn1.convar[var], w0.convar[var], wp1.convar[var], h);
		gauss.convar[var] = coe[0] + gauss.x * coe[1];
		gauss.der1y[var] = coe[1];

		Vanleer_tangential(coe, wn1.der1x[var], w0.der1x[var], wp1.der1x[var], h);
		gauss.der1x[var] = coe[0] + gauss.x * coe[1];
	}

	Check_Order_Reduce_by_Lambda_2D(gauss.is_reduce_order, gauss.convar);
	//if lambda <0, then reduce to the first order
	if (gauss.is_reduce_order == true)
	{

		for (int m = 0; m < 4; ++m)
		{
			gauss.convar[m] = w0.convar[m];
			gauss.der1x[m] = w0.der1x[m];
			gauss.der1y[m] = 0.0;
		}
	}
}
void Vanleer_tangential(double* coe, double wn1, double w0, double wp1, double h)
{
	if ((w0 - wn1) * (wp1 - w0) > 0)
	{
		double slope_left = (w0 - wn1) / h;
		double slope_right = (wp1 - w0) / h;
		coe[1] = 2 * slope_left * slope_right / (slope_left + slope_right);
		coe[0] = w0;
	}
	else
	{
		coe[0] = w0;
		coe[1] = 0.0;
	}
}

void WENO5_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		// along x direction tangitial recontruction,
		WENO5_tangential(right[0].gauss, right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, right[0].length);

		//since we already do the coordinate transform, along y, no transform needed.
		WENO5_tangential(up[0].gauss, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
			up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, up[0].length);
	}
	else
	{
		Interface2d re[4];
		int index[4]; index[0] = -2; index[1] = -1; index[2] = 1; index[3] = 2;

		for (int iface = 0; iface < 4; iface++)
		{
			Local_to_Global(re[iface].line.left.convar, right[index[iface]].line.left.convar,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.left.convar,
				right[0].normal);
			Local_to_Global(re[iface].line.left.der1x, right[index[iface]].line.left.der1x,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.left.der1x,
				right[0].normal);
			Local_to_Global(re[iface].line.right.convar, right[index[iface]].line.right.convar,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.right.convar,
				right[0].normal);
			Local_to_Global(re[iface].line.right.der1x, right[index[iface]].line.right.der1x,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.right.der1x,
				right[0].normal);
		}

		WENO5_tangential(right[0].gauss, re[0].line,
			re[1].line, right[0].line, re[2].line, re[3].line,
			right[0].length);


		index[0] = 2 * (block.ny + 1); index[1] = (block.ny + 1);
		index[2] = -(block.ny + 1); index[3] = -2 * (block.ny + 1);

		for (int iface = 0; iface < 4; iface++)
		{
			Local_to_Global(re[iface].line.left.convar, up[index[iface]].line.left.convar,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.left.convar,
				up[0].normal);
			Local_to_Global(re[iface].line.left.der1x, up[index[iface]].line.left.der1x,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.left.der1x,
				up[0].normal);
			Local_to_Global(re[iface].line.right.convar, up[index[iface]].line.right.convar,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.right.convar,
				up[0].normal);
			Local_to_Global(re[iface].line.right.der1x, up[index[iface]].line.right.der1x,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.right.der1x,
				up[0].normal);
		}

		WENO5_tangential(up[0].gauss, re[0].line, re[1].line,
			up[0].line, re[2].line, re[3].line,
			up[0].length);
	}


}

void WENO5_tangential(double* left, double* right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	//we denote that   |left...cell-center...right|
	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2);
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
		Convar_to_char(rep2, base_left, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		WENO5_left(var[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left[i] = var[i];
		}
	}
	else
	{
		Char_to_convar(left, base_left, var);

	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren2[i] = wn2[i];
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
			rep2[i] = wp2[i];
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2);
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
		Convar_to_char(rep2, base_right, wp2);
	}

	for (int i = 0; i < 4; i++)
	{
		WENO5_right(var[i], ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right[i] = var[i];

		}
	}
	else
	{
		Char_to_convar(right, base_right, var);
	}



}

void WENO5_tangential_for_slope(double* left, double* right, double* wn2, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	//we denote that   |left...cell-center...right|

	for (int i = 0; i < 4; i++)
	{
		WENO5_left(left[i], wn2[i], wn1[i], w[i], wp1[i], wp2[i], h);
	}
	for (int i = 0; i < 4; i++)
	{
		WENO5_right(right[i], wn2[i], wn1[i], w[i], wp1[i], wp2[i], h);
	}
}

void Polynomial_3rd(double* coefficent, double pn1, double w0, double pp1, double h)
{
	coefficent[0] = (1.0 / 4.0) * (-pn1 - pp1 + 6.0 * w0);
	coefficent[1] = -(pn1 - pp1) / h;
	coefficent[2] = (3 * (pn1 + pp1 - 2 * w0)) / (h * h);
}


void Polynomial_5th(double* coefficent, double wn2, double wn1, double w0, double wp1, double wp2, double h)
{
	coefficent[0] = 1.0 / 1920.0 * (2134 * w0 - 116 * (wn1 + wp1) + 9 * (wn2 + wp2));
	coefficent[1] = 1.0 / h * (-5.0 / 48.0 * (wp2 - wn2) + 17.0 / 24.0 * (wp1 - wn1));
	coefficent[2] = 1.0 / 16.0 / h / h * (12 * (wp1 + wn1) - (wp2 + wn2) - 22 * w0);
	coefficent[3] = 1.0 / 12.0 / h / h / h * ((wp2 - wn2) - 2 * (wp1 - wn1));
	coefficent[4] = 1.0 / 24.0 / h / h / h / h * ((wp2 + wn2) - 4 * (wp1 + wn1) + 6 * w0);
}
double Value_Polynomial(int order, double x0, double* coefficient)
{
	if (order <= 0) { cout << "wrong input for the order of polynomial" << endl; return 0; }
	if (order == 1) { return coefficient[0]; }
	if (order == 2) { return coefficient[0] + coefficient[1] * x0; };
	if (order == 3) { return coefficient[0] + x0 * (coefficient[1] + coefficient[2] * x0); };
	if (order == 4) { return coefficient[0] + x0 * (coefficient[1] + x0 * (coefficient[2] + coefficient[3] * x0)); };
	if (order == 5) {
		return coefficient[0] + x0 * (coefficient[1] +
			x0 * (coefficient[2] + x0 * (coefficient[3] + coefficient[4] * x0)));
	};
	if (order == 6) {
		return coefficient[0] + x0 * (coefficient[1] +
			x0 * (coefficient[2] + x0 * (coefficient[3]
				+ x0 * (coefficient[4] + coefficient[5] * x0))));
	};
	if (order == 7) {
		return coefficient[0] + x0 * (coefficient[1] +
			x0 * (coefficient[2] + x0 * (coefficient[3]
				+ x0 * (coefficient[4] + x0 * (coefficient[5] + coefficient[6] * x0)))));
	};
	return 0.0;
}

double Der1_Polynomial(int order, double x0, double* coefficient)
{
	if (order <= 0) { cout << "wrong input for the order of polynomial" << endl; return 0; }
	if (order == 1) { return 0.0; }
	if (order == 2) { return coefficient[1]; };
	if (order == 3) { return coefficient[1] + 2 * coefficient[2] * x0; };
	if (order == 4) { return coefficient[1] + x0 * (2.0 * coefficient[2] + 3.0 * coefficient[3] * x0); };
	if (order == 5)
	{
		return  coefficient[1] +
			x0 * (2.0 * coefficient[2] + x0 * (3.0 * coefficient[3] + 4.0 * coefficient[4] * x0));
	};
	if (order == 6)
	{
		return coefficient[1] +
			x0 * (2.0 * coefficient[2] + x0 * (3.0 * coefficient[3]
				+ x0 * (4.0 * coefficient[4] + 5.0 * coefficient[5] * x0)));
	};
	if (order == 7)
	{
		return coefficient[1] +
			x0 * (2.0 * coefficient[2] + x0 * (3.0 * coefficient[3]
				+ x0 * (4.0 * coefficient[4] + x0 * (5.0 * coefficient[5] + 6.0 * coefficient[6] * x0))));
	};
	return 0.0;
}

double Der2_Polynomial(int order, double x0, double* coefficient)
{
	if (order <= 0) { cout << "wrong input for the order of polynomial" << endl; return 0; }
	if (order <= 2) { return 0.0; }
	if (order == 3) { return 2.0 * coefficient[2]; }
	if (order == 4) { return 2.0 * coefficient[2] + 6.0 * coefficient[3] * x0; }
	if (order == 5)
	{
		return 2.0 * coefficient[2] + x0 * (6.0 * coefficient[3] + 12.0 * coefficient[4] * x0);
	}
	if (order == 6)
	{
		return 2.0 * coefficient[2] + x0 * (6.0 * coefficient[3]
			+ x0 * (12.0 * coefficient[4] + 20.0 * coefficient[5] * x0));
	}
	if (order == 7)
	{
		return 2.0 * coefficient[2] + x0 * (6.0 * coefficient[3]
			+ x0 * (12.0 * coefficient[4] + x0 * (20.0 * coefficient[5] + 30.0 * coefficient[6] * x0)));
	}
	return 0.0;
}



void Polynomial_3rd_avg(double* coefficent, double wn1, double w0, double wp1, double h)
{
	coefficent[0] = w0;
	coefficent[1] = 1.0 / 2.0 / h * (wp1 - wn1);
	coefficent[2] = 1.0 / 2.0 / h / h * ((wp1 + wn1) - 2 * w0);
}
double Value_Polynomial_3rd(double x, double coefficient[3])
{
	return coefficient[0] + coefficient[1] * x + coefficient[2] * x * x;
}
double Der1_Polynomial_3rd(double x, double coefficient[3])
{
	return coefficient[1] + 2 * coefficient[2] * x;
}

double Der2_Polynomial_3rd(double coefficient[3])
{
	return 2.0 * coefficient[2];
}


void WENO5_left(double& var, double wn2, double wn1, double w, double wp1, double wp2, double h)
{
	double qleft[3];
	double dleft[3];

	qleft[0] = 11.0 / 6.0 * w - 7.0 / 6.0 * wp1 + 1.0 / 3.0 * wp2;
	qleft[1] = 1.0 / 3.0 * wn1 + 5.0 / 6.0 * w - 1.0 / 6.0 * wp1;
	qleft[2] = -1.0 / 6.0 * wn2 + 5.0 / 6.0 * wn1 + 1.0 / 3.0 * w;

	dleft[0] = 0.1;
	dleft[1] = 0.6;
	dleft[2] = 0.3;

	double omegaleft[3];

	if (wenotype == linear)
	{
		omegaleft[0] = dleft[0];
		omegaleft[1] = dleft[1];
		omegaleft[2] = dleft[2];
	}
	else
	{
		double 	beta[3];

		beta[0] = 13.0 / 12.0 * pow((w - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w - 4.0 * wp1 + wp2), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w), 2);

		double epsilon = 1e-6;

		double alphaleft[3];
		if (wenotype == wenojs)
		{
			for (int i = 0; i < 3; i++)
			{
				alphaleft[i] = dleft[i] / ((beta[i] + epsilon) * (beta[i] + epsilon));
			}
		}
		if (wenotype == wenoz)
		{
			double tau5 = abs(beta[0] - beta[2]);
			for (int i = 0; i < 3; i++)
			{
				alphaleft[i] = dleft[i] * (1 + tau5 / (beta[i] + epsilon));
			}
		}

		double alphal = alphaleft[0] + alphaleft[1] + alphaleft[2];


		omegaleft[0] = alphaleft[0] / alphal;
		omegaleft[1] = alphaleft[1] / alphal;
		omegaleft[2] = alphaleft[2] / alphal;
	}


	double left = omegaleft[0] * qleft[0] + omegaleft[1] * qleft[1] + omegaleft[2] * qleft[2];

	var = left;
}
void WENO5_right(double& var, double wn2, double wn1, double w, double wp1, double wp2, double h)
{
	double qright[3];
	double dright[3];

	qright[0] = 1.0 / 3.0 * w + 5.0 / 6.0 * wp1 - 1.0 / 6.0 * wp2;
	qright[1] = -1.0 / 6.0 * wn1 + 5.0 / 6.0 * w + 1.0 / 3.0 * wp1;
	qright[2] = 1.0 / 3.0 * wn2 - 7.0 / 6.0 * wn1 + 11.0 / 6.0 * w;

	dright[0] = 0.3;
	dright[1] = 0.6;
	dright[2] = 0.1;

	double omegaright[3];

	if (wenotype == linear)
	{
		omegaright[0] = dright[0];
		omegaright[1] = dright[1];
		omegaright[2] = dright[2];
	}
	else
	{

		double 	beta[3];

		beta[0] = 13.0 / 12.0 * pow((w - 2.0 * wp1 + wp2), 2) + 0.25 * pow((3.0 * w - 4.0 * wp1 + wp2), 2);
		beta[1] = 13.0 / 12.0 * pow((wn1 - 2.0 * w + wp1), 2) + 0.25 * pow((wn1 - wp1), 2);
		beta[2] = 13.0 / 12.0 * pow((wn2 - 2.0 * wn1 + w), 2) + 0.25 * pow((wn2 - 4.0 * wn1 + 3.0 * w), 2);

		double epsilon = 1e-6;
		double alphaleft[3];	double alpharight[3];
		if (wenotype == wenojs)
		{
			for (int i = 0; i < 3; i++)
			{
				alpharight[i] = dright[i] / ((beta[i] + epsilon) * (beta[i] + epsilon));
			}
		}
		if (wenotype == wenoz)
		{
			double tau5 = abs(beta[0] - beta[2]);
			for (int i = 0; i < 3; i++)
			{
				alpharight[i] = dright[i] * (1 + tau5 / (beta[i] + epsilon));
			}
		}

		double alphar = alpharight[0] + alpharight[1] + alpharight[2];

		omegaright[0] = alpharight[0] / alphar;
		omegaright[1] = alpharight[1] / alphar;
		omegaright[2] = alpharight[2] / alphar;
	}


	double right = omegaright[0] * qright[0] + omegaright[1] * qright[1] + omegaright[2] * qright[2];
	var = right;

}
void WENO5_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double h)
{
	double left[4], right[4];
	double coe[3];

	// let's reconstruct the left first
	WENO5_tangential(left, right, wn2.left.convar, wn1.left.convar, w0.left.convar, wp1.left.convar, wp2.left.convar, h);

	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.left.convar[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].left.convar[var] = Value_Polynomial_3rd(re[num_gauss].left.x, coe);
			re[num_gauss].left.der1y[var] = Der1_Polynomial_3rd(re[num_gauss].left.x, coe);
		}
	}

	WENO5_tangential_for_slope(left, right, wn2.left.der1x, wn1.left.der1x, w0.left.der1x, wp1.left.der1x, wp2.left.der1x, h);
	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.left.der1x[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].left.der1x[var] = Value_Polynomial_3rd(re[num_gauss].left.x, coe);

		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
			}
		}
	}

	// then is the right part
	WENO5_tangential(left, right, wn2.right.convar, wn1.right.convar, w0.right.convar, wp1.right.convar, wp2.right.convar, h);

	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.right.convar[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].right.convar[var] = Value_Polynomial_3rd(re[num_gauss].right.x, coe);
			re[num_gauss].right.der1y[var] = Der1_Polynomial_3rd(re[num_gauss].right.x, coe);
		}
	}

	WENO5_tangential_for_slope(left, right, wn2.right.der1x, wn1.right.der1x, w0.right.der1x, wp1.right.der1x, wp2.right.der1x, h);

	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.right.der1x[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].right.der1x[var] = Value_Polynomial_3rd(re[num_gauss].right.x, coe);

		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
			}
		}
	}

	//right part compelete
}

void WENO5_AO_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		// along x direction tangitial recontruction,
		WENO5_AO_tangential(right[0].gauss,
			right[-2].line, right[-1].line, right[0].line, right[1].line, right[2].line, right[0].length);

		//since we already do the coordinate transform, along y, no transform needed.
		WENO5_AO_tangential(up[0].gauss, up[2 * (block.ny + 1)].line, up[block.ny + 1].line,
			up[0].line, up[-(block.ny + 1)].line, up[-2 * (block.ny + 1)].line, up[0].length);
	}
	else
	{
		Interface2d re[4];
		int index[4]; index[0] = -2; index[1] = -1; index[2] = 1; index[3] = 2;

		for (int iface = 0; iface < 4; iface++)
		{
			Local_to_Global(re[iface].line.left.convar, right[index[iface]].line.left.convar,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.left.convar,
				right[0].normal);
			Local_to_Global(re[iface].line.left.der1x, right[index[iface]].line.left.der1x,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.left.der1x,
				right[0].normal);
			Local_to_Global(re[iface].line.right.convar, right[index[iface]].line.right.convar,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.right.convar,
				right[0].normal);
			Local_to_Global(re[iface].line.right.der1x, right[index[iface]].line.right.der1x,
				right[index[iface]].normal);
			Global_to_Local(re[iface].line.right.der1x,
				right[0].normal);
		}

		WENO5_AO_tangential(right[0].gauss, re[0].line,
			re[1].line, right[0].line, re[2].line, re[3].line,
			right[0].length);


		index[0] = 2 * (block.ny + 1); index[1] = (block.ny + 1);
		index[2] = -(block.ny + 1); index[3] = -2 * (block.ny + 1);

		for (int iface = 0; iface < 4; iface++)
		{
			Local_to_Global(re[iface].line.left.convar, up[index[iface]].line.left.convar,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.left.convar,
				up[0].normal);
			Local_to_Global(re[iface].line.left.der1x, up[index[iface]].line.left.der1x,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.left.der1x,
				up[0].normal);
			Local_to_Global(re[iface].line.right.convar, up[index[iface]].line.right.convar,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.right.convar,
				up[0].normal);
			Local_to_Global(re[iface].line.right.der1x, up[index[iface]].line.right.der1x,
				up[index[iface]].normal);
			Global_to_Local(re[iface].line.right.der1x,
				up[0].normal);
		}

		WENO5_AO_tangential(up[0].gauss, re[0].line, re[1].line,
			up[0].line, re[2].line, re[3].line,
			up[0].length);

	}
}
void WENO5_AO_tangential(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double h)
{
	//lets first reconstruction the left value

	double ren2[4], ren1[4], re0[4], rep1[4], rep2[4];
	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(w_primvar, w0.left.convar);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			if (gausspoint == 2)
			{
				weno_5th_ao_2gauss(re[0].left.convar[i], re[0].left.der1y[i], re[0].left.der2yy[i],
					re[1].left.convar[i], re[1].left.der1y[i], re[1].left.der2yy[i],
					wn2.left.convar[i], wn1.left.convar[i], w0.left.convar[i], wp1.left.convar[i], wp2.left.convar[i], h, 2);
			}
		}
	}
	else
	{
		Convar_to_char(ren2, base_left, wn2.left.convar);
		Convar_to_char(ren1, base_left, wn1.left.convar);
		Convar_to_char(re0, base_left, w0.left.convar);
		Convar_to_char(rep1, base_left, wp1.left.convar);
		Convar_to_char(rep2, base_left, wp2.left.convar);
		if (gausspoint == 2)
		{
			double var[2][4], der1[2][4], der2[2][4];
			for (int i = 0; i < 4; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], der2[0][i],
					var[1][i], der1[1][i], der2[1][i],
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h, 2);
			}
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				Char_to_convar(re[igauss].left.convar, base_left, var[igauss]);
				Char_to_convar(re[igauss].left.der1y, base_left, der1[igauss]);
				Char_to_convar(re[igauss].left.der2yy, base_left, der2[igauss]);
			}
		}

	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		if (gausspoint == 2)
		{
			weno_5th_ao_2gauss(re[0].left.der1x[i], tmp[0], tmp[1],
				re[1].left.der1x[i], tmp[2], tmp[3],
				wn2.left.der1x[i], wn1.left.der1x[i], w0.left.der1x[i], wp1.left.der1x[i], wp2.left.der1x[i], h, 1);
			re[0].left.der2xy[i] = tmp[0];
			re[1].left.der2xy[i] = tmp[2];
		}
	}

	//then let's construct the right part....
	Convar_to_primvar_2D(w_primvar, w0.right.convar);
	for (int i = 0; i < 4; i++)
	{
		base_right[i] = (w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		double tmp[2];
		for (int i = 0; i < 4; i++)
		{
			if (gausspoint == 2)
			{
				weno_5th_ao_2gauss(re[0].right.convar[i], re[0].right.der1y[i], re[0].right.der2yy[i],
					re[1].right.convar[i], re[1].right.der1y[i], re[1].right.der2yy[i],
					wn2.right.convar[i], wn1.right.convar[i], w0.right.convar[i], wp1.right.convar[i], wp2.right.convar[i], h, 2);
			}
		}
	}
	else
	{
		Convar_to_char(ren2, base_right, wn2.right.convar);
		Convar_to_char(ren1, base_right, wn1.right.convar);
		Convar_to_char(re0, base_right, w0.right.convar);
		Convar_to_char(rep1, base_right, wp1.right.convar);
		Convar_to_char(rep2, base_right, wp2.right.convar);
		if (gausspoint == 2)
		{
			double var[2][4], der1[2][4], der2[2][4];
			for (int i = 0; i < 4; i++)
			{
				weno_5th_ao_2gauss(var[0][i], der1[0][i], der2[0][i],
					var[1][i], der1[1][i], der2[1][i],
					ren2[i], ren1[i], re0[i], rep1[i], rep2[i], h, 2);
			}
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				Char_to_convar(re[igauss].right.convar, base_right, var[igauss]);
				Char_to_convar(re[igauss].right.der1y, base_right, der1[igauss]);
				Char_to_convar(re[igauss].right.der2yy, base_right, der2[igauss]);
			}
		}
	}

	for (int i = 0; i < 4; i++)
	{
		double tmp[4];
		if (gausspoint == 2)
		{
			weno_5th_ao_2gauss(re[0].right.der1x[i], tmp[0], tmp[1],
				re[1].right.der1x[i], tmp[2], tmp[3],
				wn2.right.der1x[i], wn1.right.der1x[i], w0.right.der1x[i], wp1.right.der1x[i], wp2.right.der1x[i], h, 1);
			re[0].right.der2xy[i] = tmp[0];
			re[1].right.der2xy[i] = tmp[2];
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{

			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
				re[num_gauss].left.der2xy[var] = 0.0;
				re[num_gauss].left.der2yy[var] = 0.0;
			}
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{

			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
				re[num_gauss].right.der2xy[var] = 0.0;
				re[num_gauss].right.der2yy[var] = 0.0;
			}
		}
	}

}
void weno_5th_ao_2gauss(double& g1, double& g1x, double& g1xx, double& g2, double& g2x, double& g2xx, double wn2, double wn1, double w0, double wp1, double wp2, double h, int order)
{
	//the parameter order constrols up to which order you want construct
	// order from 0, 1, 2
	if (order > 2 || order < 0)
	{
		cout << "invalid order input for the function " << __FUNCTION__ << endl;
		exit(0);
	}
	double dhi = 0.85;
	double dlo = 0.85;
	//-- - parameter of WENO-- -
	double beta[4], d[4], ww[4], alpha[4];
	double epsilonW = 1e-8;
	//-- - intermediate parameter-- -
	double p[4], px[4], pxx[4], tempvar;
	double sum_alpha;

	//three small stencil
	d[0] = (1 - dhi) * (1 - dlo) * 0.5;
	d[1] = (1 - dhi) * dlo;
	d[2] = (1 - dhi) * (1 - dlo) * 0.5;
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
	////-- - combination-- -

	double final_weight[4];
	final_weight[3] = ww[3] / d[3];
	for (int k = 0; k < 3; k++)
	{
		final_weight[k] = ww[k] - ww[3] / d[3] * d[k];
	}

	g1 = 0; g1x = 0; g1xx = 0;
	g2 = 0; g2x = 0; g2xx = 0;

	double sqrt3 = sqrt(3);
	//-- - candidate polynomial-- for gauss 1 point
	p[0] = w0 - (sqrt3 * w0) / 4 + (4 * wn1 - wn2) / (4 * sqrt3);
	p[1] = w0 + (wn1 - wp1) / (4 * sqrt3);
	p[2] = (1.0 / 12.0) * (3 * (4 + sqrt3) * w0 + sqrt3 * (-4 * wp1 + wp2));
	p[3] = (4314 * w0 + (4 + 500 * sqrt3) * wn1 - wn2 - 70 * sqrt3 * wn2 +
		4 * wp1 - 500 * sqrt3 * wp1 - wp2 + 70 * sqrt3 * wp2) / 4320;

	for (int k = 0; k < 4; k++)
	{
		g1 += final_weight[k] * p[k];
	}

	if (order > 0)
	{
		px[0] = -((-9 + sqrt3) * w0 - 2 * (-6 + sqrt3) * wn1 + (-3 + sqrt3) * wn2) / (6 * h);
		px[1] = -(-2 * sqrt3 * w0 + (3 + sqrt3) * wn1 + (-3 + sqrt3) * wp1) / (6 * h);
		px[2] = -((9 + sqrt3) * w0 - 2 * (6 + sqrt3) * wp1 + (3 + sqrt3) * wp2) / (6 * h);
		px[3] = (48 * sqrt3 * w0 - 72 * wn1 - 26 * sqrt3 * wn1 + 9 * wn2 +
			2 * sqrt3 * wn2 + 72 * wp1 - 26 * sqrt3 * wp1 - 9 * wp2 +
			2 * sqrt3 * wp2) / (108 * h);

		for (int k = 0; k < 4; k++)
		{
			g1x += final_weight[k] * px[k];
		}
		if (order == 2)
		{
			pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
			pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
			pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
			pxx[3] = -((30 * w0 + 2 * (-8 + sqrt3) * wn1 + wn2 - sqrt3 * wn2 - 16 * wp1 -
				2 * sqrt3 * wp1 + wp2 + sqrt3 * wp2) / (12 * h * h));

			for (int k = 0; k < 4; k++)
			{
				g1xx += final_weight[k] * pxx[k];
			}
		}
	}

	//-- - candidate polynomial-- for gauss 2 point
	p[0] = (1.0 / 12.0) * (3 * (4 + sqrt3) * w0 + sqrt3 * (-4 * wn1 + wn2));
	p[1] = w0 + (-wn1 + wp1) / (4 * sqrt3);
	p[2] = w0 - (sqrt3 * w0) / 4 + (4 * wp1 - wp2) / (4 * sqrt3);
	p[3] = (4314 * w0 + (4 - 500 * sqrt3) * wn1 - wn2 + 70 * sqrt3 * wn2 +
		4 * wp1 + 500 * sqrt3 * wp1 - wp2 - 70 * sqrt3 * wp2) / 4320;
	for (int k = 0; k < 4; k++)
	{
		g2 += final_weight[k] * p[k];
	}
	if (order > 0)
	{
		px[0] = ((9 + sqrt3) * w0 - 2 * (6 + sqrt3) * wn1 + (3 + sqrt3) * wn2) /
			(6 * h);
		px[1] = (-2 * sqrt3 * w0 + (-3 + sqrt3) * wn1 + (3 + sqrt3) * wp1) / (6 * h);
		px[2] = ((-9 + sqrt3) * w0 - 2 * (-6 + sqrt3) * wp1 + (-3 + sqrt3) * wp2) / (6 * h);
		px[3] = -((48 * sqrt3 * w0 + 72 * wn1 - 26 * sqrt3 * wn1 - 9 * wn2 +
			2 * sqrt3 * wn2 - 72 * wp1 - 26 * sqrt3 * wp1 + 9 * wp2 +
			2 * sqrt3 * wp2) / (108 * h));
		for (int k = 0; k < 4; k++)
		{
			g2x += final_weight[k] * px[k];
		}
		if (order == 2)
		{
			pxx[0] = (w0 - 2 * wn1 + wn2) / h / h;
			pxx[1] = (-2 * w0 + wn1 + wp1) / h / h;
			pxx[2] = (w0 - 2 * wp1 + wp2) / h / h;
			pxx[3] = -((30 * w0 - 2 * (8 + sqrt3) * wn1 + wn2 + sqrt3 * wn2 - 16 * wp1 +
				2 * sqrt3 * wp1 + wp2 - sqrt3 * wp2) / (12 * h * h));
		}
		for (int k = 0; k < 4; k++)
		{
			g2xx += final_weight[k] * pxx[k];
		}
	}
}

void WENO3_tangent(Interface2d* left, Interface2d* right, Interface2d* down, Interface2d* up, Fluid2d* fluids, Block2d block)
{

	Interface2d re[4];
	int index[4]; index[0] = -2; index[1] = -1; index[2] = 1; index[3] = 2;

	for (int iface = 0; iface < 4; iface++)
	{
		Local_to_Global(re[iface].line.left.convar, right[index[iface]].line.left.convar,
			right[index[iface]].normal);
		Global_to_Local(re[iface].line.left.convar,
			right[0].normal);
		Local_to_Global(re[iface].line.left.der1x, right[index[iface]].line.left.der1x,
			right[index[iface]].normal);
		Global_to_Local(re[iface].line.left.der1x,
			right[0].normal);
		Local_to_Global(re[iface].line.right.convar, right[index[iface]].line.right.convar,
			right[index[iface]].normal);
		Global_to_Local(re[iface].line.right.convar,
			right[0].normal);
		Local_to_Global(re[iface].line.right.der1x, right[index[iface]].line.right.der1x,
			right[index[iface]].normal);
		Global_to_Local(re[iface].line.right.der1x,
			right[0].normal);
	}

	WENO3_tangential(right[0].gauss,
		re[1].line, right[0].line, re[2].line,
		right[0].length);


	index[0] = 2 * (block.ny + 1); index[1] = (block.ny + 1);
	index[2] = -(block.ny + 1); index[3] = -2 * (block.ny + 1);

	for (int iface = 0; iface < 4; iface++)
	{
		Local_to_Global(re[iface].line.left.convar, up[index[iface]].line.left.convar,
			up[index[iface]].normal);
		Global_to_Local(re[iface].line.left.convar,
			up[0].normal);
		Local_to_Global(re[iface].line.left.der1x, up[index[iface]].line.left.der1x,
			up[index[iface]].normal);
		Global_to_Local(re[iface].line.left.der1x,
			up[0].normal);
		Local_to_Global(re[iface].line.right.convar, up[index[iface]].line.right.convar,
			up[index[iface]].normal);
		Global_to_Local(re[iface].line.right.convar,
			up[0].normal);
		Local_to_Global(re[iface].line.right.der1x, up[index[iface]].line.right.der1x,
			up[index[iface]].normal);
		Global_to_Local(re[iface].line.right.der1x,
			up[0].normal);
	}

	WENO3_tangential(up[0].gauss, re[1].line,
		up[0].line, re[2].line,
		up[0].length);

}

void WENO3_tangential(Recon2d* re, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, double h)
{
	double left[4], right[4];
	double coe[3];

	// let's reconstruct the left first
	WENO3_tangential(left, right, wn1.left.convar, w0.left.convar, wp1.left.convar, h);

	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.left.convar[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].left.convar[var] = Value_Polynomial_3rd(re[num_gauss].left.x, coe);
			re[num_gauss].left.der1y[var] = Der1_Polynomial_3rd(re[num_gauss].left.x, coe);
		}
	}

	WENO3_tangential_for_slope(left, right, wn1.left.der1x, w0.left.der1x, wp1.left.der1x, h);
	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.left.der1x[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].left.der1x[var] = Value_Polynomial_3rd(re[num_gauss].left.x, coe);

		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].left.is_reduce_order, re[num_gauss].left.convar);
		if (re[num_gauss].left.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].left.convar[var] = w0.left.convar[var];
				re[num_gauss].left.der1x[var] = w0.left.der1x[var];
				re[num_gauss].left.der1y[var] = 0.0;
			}
		}
	}

	// then is the right part
	WENO3_tangential(left, right, wn1.right.convar, w0.right.convar, wp1.right.convar, h);

	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.right.convar[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].right.convar[var] = Value_Polynomial_3rd(re[num_gauss].right.x, coe);
			re[num_gauss].right.der1y[var] = Der1_Polynomial_3rd(re[num_gauss].right.x, coe);
		}
	}

	WENO3_tangential_for_slope(left, right, wn1.right.der1x, w0.right.der1x, wp1.right.der1x, h);

	for (int var = 0; var < 4; var++)
	{
		Polynomial_3rd(coe, left[var], w0.right.der1x[var], right[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].right.der1x[var] = Value_Polynomial_3rd(re[num_gauss].right.x, coe);
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].right.is_reduce_order, re[num_gauss].right.convar);
		if (re[num_gauss].right.is_reduce_order == true)
		{
			for (int var = 0; var < 4; var++)
			{
				re[num_gauss].right.convar[var] = w0.right.convar[var];
				re[num_gauss].right.der1x[var] = w0.right.der1x[var];
				re[num_gauss].right.der1y[var] = 0.0;
			}
		}
	}
	//right part compelete
}

void WENO3_tangential(double* left, double* right, double* wn1, double* w, double* wp1, double h)
{
	//we denote that   |left...cell-center...right|
	double  ren1[4], re0[4], rep1[4];
	double var[4], der1[4], der2[4];

	double base_left[4];
	double base_right[4];
	double wn1_primvar[4], w_primvar[4], wp1_primvar[4];
	Convar_to_primvar_2D(wn1_primvar, wn1);
	Convar_to_primvar_2D(w_primvar, w);
	Convar_to_primvar_2D(wp1_primvar, wp1);

	for (int i = 0; i < 4; i++)
	{
		base_left[i] = 0.5 * (wn1_primvar[i] + w_primvar[i]);
		base_right[i] = 0.5 * (wp1_primvar[i] + w_primvar[i]);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
		}
	}
	else
	{
		Convar_to_char(ren1, base_left, wn1);
		Convar_to_char(re0, base_left, w);
		Convar_to_char(rep1, base_left, wp1);
	}

	for (int i = 0; i < 4; i++)
	{
		WENO3_left(var[i], ren1[i], re0[i], rep1[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			left[i] = var[i];
		}
	}
	else
	{
		Char_to_convar(left, base_left, var);
	}

	// cell right
	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			ren1[i] = wn1[i];
			re0[i] = w[i];
			rep1[i] = wp1[i];
		}
	}
	else
	{
		Convar_to_char(ren1, base_right, wn1);
		Convar_to_char(re0, base_right, w);
		Convar_to_char(rep1, base_right, wp1);
	}

	for (int i = 0; i < 4; i++)
	{
		WENO3_right(var[i], ren1[i], re0[i], rep1[i], h);
	}

	if (reconstruction_variable == conservative)
	{
		for (int i = 0; i < 4; i++)
		{
			right[i] = var[i];
		}
	}
	else
	{
		Char_to_convar(right, base_right, var);
	}



}

void WENO3_tangential_for_slope(double* left, double* right, double* wn1, double* w, double* wp1, double h)
{
	//we denote that   |left...cell-center...right|

	for (int i = 0; i < 4; i++)
	{
		WENO3_left(left[i], wn1[i], w[i], wp1[i], h);
	}
	for (int i = 0; i < 4; i++)
	{
		WENO3_right(right[i], wn1[i], w[i], wp1[i], h);
	}
}

void Reconstruction_forg0(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{

#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			g0reconstruction_2D_normal(&xinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j], &fluids[i * (block.ny) + j], block);
		}
	}

	// then get the guass point value. That is, so called multi-dimensional property
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nx - block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.ny - block.ghost + 1; j++)
		{
			g0reconstruction_2D_tangent(&xinterfaces[i * (block.ny + 1) + j], &yinterfaces[i * (block.ny + 1) + j], &fluids[i * (block.ny) + j], block);
		}
	}
}
void Center_collision_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
	{
		Center_all_collision_2d_multi(xinterfaces[0].line);
	}
	if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
	{
		Center_all_collision_2d_multi(yinterfaces[0].line);
	}
}
void Center_3rd_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	//assume this is uniform mesh,
	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			Center_3rd_Splitting(xinterfaces[0], fluids[-block.ny].convar, fluids[0].convar, block.dx);
		}
		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double tmp_w[4], tmp_wp1[4];
			YchangetoX(tmp_wp1, fluids[0].convar);
			YchangetoX(tmp_w, fluids[-1].convar);
			Center_3rd_Splitting(yinterfaces[0], tmp_w, tmp_wp1, block.dy);
		}
	}
	else
	{
		double tmp_w[4], tmp_wp1[4];
		double normal[2];
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{

			double dx = 0.5 * (fluids[-block.ny].dx + fluids[0].dx);

			Copy_Array(normal, xinterfaces[0].normal, 2);
			Global_to_Local(tmp_w, fluids[-block.ny].convar, normal);
			Global_to_Local(tmp_wp1, fluids[0].convar, normal);
			double over_cos = -dx / ((fluids[-block.ny].coordx - fluids[0].coordx) * normal[0] + (fluids[-block.ny].coordy - fluids[0].coordy) * normal[1]);

			Center_3rd_Splitting(xinterfaces[0], tmp_w, tmp_wp1, over_cos * dx);
		}
		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double dy = 0.5 * (fluids[-1].dy + fluids[0].dy);

			Copy_Array(normal, yinterfaces[0].normal, 2);
			Global_to_Local(tmp_w, fluids[-1].convar, normal);
			Global_to_Local(tmp_wp1, fluids[0].convar, normal);
			double over_cos = -dy / ((fluids[-1].coordx - fluids[0].coordx) * normal[0] + (fluids[-1].coordy - fluids[0].coordy) * normal[1]);

			Center_3rd_Splitting(yinterfaces[0], tmp_w, tmp_wp1, over_cos * dy);
		}

	}
}
void Center_3rd_Splitting(Interface2d& interface, double* w, double* wp1, double h)
{
	double prim_left[4], prim_right[4];
	Convar_to_ULambda_2d(prim_left, interface.line.left.convar);
	Convar_to_ULambda_2d(prim_right, interface.line.right.convar);

	MMDF1st ml(prim_left[1], prim_left[2], prim_left[3]);
	MMDF1st mr(prim_right[1], prim_right[2], prim_right[3]);

	Collision(interface.line.center.convar, prim_left[0], prim_right[0], ml, mr);

	for (int var = 0; var < 4; var++)
	{
		interface.line.center.der1x[var] = (wp1[var] - w[var]) / h;
		interface.line.center.der2xx[var] = 0.0;
	}
}
void Center_do_nothing_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	// do nothing;
}

void Center_4th_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{

	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			Center_4th_normal(xinterfaces[0], fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[0].dx);
		}
		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double tmp_wn1[4], tmp_w[4], tmp_wp1[4], tmp_wp2[4];
			YchangetoX(tmp_wn1, fluids[-2].convar);
			YchangetoX(tmp_w, fluids[-1].convar);
			YchangetoX(tmp_wp1, fluids[0].convar);
			YchangetoX(tmp_wp2, fluids[1].convar);
			Center_4th_normal(yinterfaces[0], tmp_wn1, tmp_w, tmp_wp1, tmp_wp2, fluids[0].dy);
		}
	}
	else
	{
		double tmp_wn1[4], tmp_w[4], tmp_wp1[4], tmp_wp2[4];
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			double dx = 0.5 * (fluids[-block.ny].dx + fluids[0].dx);

			Global_to_Local(tmp_wn1, fluids[-2 * block.ny].convar, xinterfaces[0].normal);
			Global_to_Local(tmp_w, fluids[-1 * block.ny].convar, xinterfaces[0].normal);
			Global_to_Local(tmp_wp1, fluids[0].convar, xinterfaces[0].normal);
			Global_to_Local(tmp_wp2, fluids[block.ny].convar, xinterfaces[0].normal);
			Center_4th_normal(xinterfaces[0], tmp_wn1, tmp_w, tmp_wp1, tmp_wp2, dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double dy = 0.5 * (fluids[-1].dy + fluids[0].dy);
			Global_to_Local(tmp_wn1, fluids[-2].convar, yinterfaces[0].normal);
			Global_to_Local(tmp_w, fluids[-1].convar, yinterfaces[0].normal);
			Global_to_Local(tmp_wp1, fluids[0].convar, yinterfaces[0].normal);
			Global_to_Local(tmp_wp2, fluids[1].convar, yinterfaces[0].normal);
			Center_4th_normal(yinterfaces[0], tmp_wn1, tmp_w, tmp_wp1, tmp_wp2, dy);
		}
	}
}

void Center_4th_normal(Interface2d& interface, double* wn1, double* w, double* wp1, double* wp2, double h)
{
	for (int var = 0; var < 4; var++)
	{
		interface.line.center.convar[var] = (-1.0 / 12.0 * (wp2[var] + wn1[var]) + 7.0 / 12.0 * (wp1[var] + w[var]));
		interface.line.center.der1x[var] = (-1.0 / 12.0 * (wp2[var] - wn1[var]) + 1.25 * (wp1[var] - w[var])) / h;
	}
}

void Center_5th_normal(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{

	if (block.uniform == true)
	{
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			Center_5th_normal(xinterfaces[0], fluids[-2 * block.ny].convar, fluids[-block.ny].convar, fluids[0].convar, fluids[block.ny].convar, fluids[0].dx);
		}
		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double tmp_wn1[4], tmp_w[4], tmp_wp1[4], tmp_wp2[4];
			YchangetoX(tmp_wn1, fluids[-2].convar);
			YchangetoX(tmp_w, fluids[-1].convar);
			YchangetoX(tmp_wp1, fluids[0].convar);
			YchangetoX(tmp_wp2, fluids[1].convar);
			Center_5th_normal(yinterfaces[0], tmp_wn1, tmp_w, tmp_wp1, tmp_wp2, fluids[0].dy);
		}
	}
	else
	{
		double tmp_wn1[4], tmp_w[4], tmp_wp1[4], tmp_wp2[4];
		if ((fluids[0].xindex > block.ghost - 2) && (fluids[0].xindex < block.nx - block.ghost + 1))
		{
			double dx = 0.5 * (fluids[-block.ny].dx + fluids[0].dx);

			Global_to_Local(tmp_wn1, fluids[-2 * block.ny].convar, xinterfaces[0].normal);
			Global_to_Local(tmp_w, fluids[-1 * block.ny].convar, xinterfaces[0].normal);
			Global_to_Local(tmp_wp1, fluids[0].convar, xinterfaces[0].normal);
			Global_to_Local(tmp_wp2, fluids[block.ny].convar, xinterfaces[0].normal);
			Center_5th_normal(xinterfaces[0], tmp_wn1, tmp_w, tmp_wp1, tmp_wp2, dx);
		}

		if ((fluids[0].yindex > block.ghost - 2) && (fluids[0].yindex < block.ny - block.ghost + 1))
		{
			double dy = 0.5 * (fluids[-1].dy + fluids[0].dy);
			Global_to_Local(tmp_wn1, fluids[-2].convar, yinterfaces[0].normal);
			Global_to_Local(tmp_w, fluids[-1].convar, yinterfaces[0].normal);
			Global_to_Local(tmp_wp1, fluids[0].convar, yinterfaces[0].normal);
			Global_to_Local(tmp_wp2, fluids[1].convar, yinterfaces[0].normal);
			Center_5th_normal(yinterfaces[0], tmp_wn1, tmp_w, tmp_wp1, tmp_wp2, dy);
		}
	}
}


void Center_5th_normal(Interface2d& interface, double* wn1, double* w, double* wp1, double* wp2, double h)
{




	double prim_left[4], prim_right[4];
	Convar_to_ULambda_2d(prim_left, interface.line.left.convar);
	Convar_to_ULambda_2d(prim_right, interface.line.right.convar);

	MMDF1st ml(prim_left[1], prim_left[2], prim_left[3]);
	MMDF1st mr(prim_right[1], prim_right[2], prim_right[3]);

	Collision(interface.line.center.convar, prim_left[0], prim_right[0], ml, mr);
	for (int var = 0; var < 4; var++)
	{
		interface.line.center.der1x[var] = (-1.0 / 12.0 * (wp2[var] - wn1[var]) + 1.25 * (wp1[var] - w[var])) / h;
		interface.line.center.der2xx[var] = (wn1[var] - w[var] - wp1[var] + wp2[var]) / (2.0 * h * h);
	}


	bool order_reduce = false;

	if (interface.line.left.is_reduce_order == true
		|| interface.line.right.is_reduce_order == true)
	{
		order_reduce = true;
	}

	if (order_reduce == true)
	{
		for (int var = 0; var < 4; var++)
		{
			interface.line.center.der1x[var] = (wp1[var] - w[var]) / h;
			interface.line.center.der2xx[var] = 0.0;
		}
	}

}


void Center_all_collision_multi(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	if (flux_function_2d == GKS2D_smooth)
	{
		for (int m = 0; m < gausspoint; m++)
		{
			Center_all_collision_2d_multi(xinterfaces[0].gauss[m]);
			Center_all_collision_2d_multi(yinterfaces[0].gauss[m]);
		}
	}
}
void Center_all_collision_2d_multi(Recon2d& gauss)
{
	double prim_left[4], prim_right[4];
	Convar_to_ULambda_2d(prim_left, gauss.left.convar);
	Convar_to_ULambda_2d(prim_right, gauss.right.convar);

	MMDF ml(prim_left[1], prim_left[2], prim_left[3]);
	MMDF mr(prim_right[1], prim_right[2], prim_right[3]);
	Collision(gauss.center.convar, prim_left[0], prim_right[0], ml, mr);

	double a0[4] = { 1.0, 0.0, 0.0, 0.0 };
	double alx[4], arx[4];
	//w_x
	A_point(alx, gauss.left.der1x, prim_left);
	A_point(arx, gauss.right.der1x, prim_right);
	double al0x[4];
	double ar0x[4];
	GL_address(0, 0, 0, al0x, alx, ml);
	GR_address(0, 0, 0, ar0x, arx, mr);


	double aly[4], ary[4];
	//w_y
	A_point(aly, gauss.left.der1y, prim_left);
	A_point(ary, gauss.right.der1y, prim_right);
	double al0y[4];
	double ar0y[4];
	GL_address(0, 0, 0, al0y, aly, ml);
	GR_address(0, 0, 0, ar0y, ary, mr);





	for (int var = 0; var < 4; var++)
	{
		gauss.center.der1x[var] = prim_left[0] * al0x[var] + prim_right[0] * ar0x[var];
		gauss.center.der1y[var] = prim_left[0] * al0y[var] + prim_right[0] * ar0y[var];
	}

}

void Center_first_order_tangent(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
	{
		first_order_tangent(xinterfaces[0].gauss[num_gauss].center, xinterfaces[0].line.center);
		first_order_tangent(yinterfaces[0].gauss[num_gauss].center, yinterfaces[0].line.center);
	}
}

void Center_3rd_tangent(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			//assume this is uniform mesh,
			// along x direction tangitial recontruction,
			Center_3rd_Multi(xinterfaces[0].gauss[num_gauss].center, xinterfaces[-1].line.center,
				xinterfaces[0].line.center, xinterfaces[1].line.center, xinterfaces[0].length);
			//since we already do the coordinate transform, along y, no transform needed.
			Center_3rd_Multi(yinterfaces[0].gauss[num_gauss].center, yinterfaces[block.ny + 1].line.center,
				yinterfaces[0].line.center, yinterfaces[-block.ny - 1].line.center, yinterfaces[0].length);
		}
	}
	else
	{
		Interface2d re[2];
		int index[2]; index[0] = -1; index[1] = 1;

		for (int iface = 0; iface < 2; iface++)
		{
			Local_to_Global(re[iface].line.center.convar, xinterfaces[index[iface]].line.center.convar,
				xinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.convar,
				xinterfaces[0].normal);
			Local_to_Global(re[iface].line.center.der1x, xinterfaces[index[iface]].line.center.der1x,
				xinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.der1x,
				xinterfaces[0].normal);

		}
		for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
		{
			//assume this is uniform mesh,
			Center_3rd_Multi(xinterfaces[0].gauss[num_gauss].center, re[0].line.center,
				xinterfaces[0].line.center, re[1].line.center, xinterfaces[0].length);
		}

		index[0] = (block.ny + 1); index[1] = -(block.ny + 1);
		for (int iface = 0; iface < 2; iface++)
		{
			Local_to_Global(re[iface].line.center.convar, yinterfaces[index[iface]].line.center.convar,
				yinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.convar,
				yinterfaces[0].normal);
			Local_to_Global(re[iface].line.center.der1x, yinterfaces[index[iface]].line.center.der1x,
				yinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.der1x,
				yinterfaces[0].normal);

		}
		for (int num_gauss = 0; num_gauss < gausspoint; ++num_gauss)
		{
			//assume this is uniform mesh,
			Center_3rd_Multi(yinterfaces[0].gauss[num_gauss].center, re[0].line.center,
				yinterfaces[0].line.center, re[1].line.center, yinterfaces[0].length);

		}
	}
}
void Center_5th_tangent(Interface2d* xinterfaces, Interface2d* yinterfaces, Fluid2d* fluids, Block2d block)
{
	if (block.uniform == true)
	{
		Center_5th_Multi(xinterfaces[0].gauss, xinterfaces[-2].line, xinterfaces[-1].line, xinterfaces[0].line, xinterfaces[1].line, xinterfaces[2].line, xinterfaces[0].length);
		Center_5th_Multi(yinterfaces[0].gauss, yinterfaces[2 * (block.ny + 1)].line, yinterfaces[block.ny + 1].line, yinterfaces[0].line, yinterfaces[-block.ny - 1].line, yinterfaces[-2 * (block.ny + 1)].line, yinterfaces[0].length);
		// x direction
		bool order_reduce = false;
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			if (xinterfaces[0].gauss[num_gauss].center.is_reduce_order == true)
			{
				order_reduce = true;
			}
		}
		if (order_reduce == true)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Center_3rd_Multi(xinterfaces[0].gauss[num_gauss].center, xinterfaces[-1].line.center, xinterfaces[0].line.center, xinterfaces[1].line.center, xinterfaces[0].length);
			}
		}

		// y  direction
		order_reduce = false;
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			if (yinterfaces[0].gauss[num_gauss].center.is_reduce_order == true)
			{
				order_reduce = true;
			}
		}
		if (order_reduce == true)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Center_3rd_Multi(yinterfaces[0].gauss[num_gauss].center, yinterfaces[block.ny + 1].line.center, yinterfaces[0].line.center, yinterfaces[-block.ny - 1].line.center, yinterfaces[0].length);
			}
		}

	}
	else
	{
		Interface2d re[4];
		int index[4]; index[0] = -2; index[1] = -1; index[2] = 1; index[3] = 2;

		for (int iface = 0; iface < 4; ++iface)
		{
			Local_to_Global(re[iface].line.center.convar, xinterfaces[index[iface]].line.center.convar,
				xinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.convar,
				xinterfaces[0].normal);
			Local_to_Global(re[iface].line.center.der1x, xinterfaces[index[iface]].line.center.der1x,
				xinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.der1x,
				xinterfaces[0].normal);

		}

		//assume this is uniform mesh,
		Center_5th_Multi(xinterfaces[0].gauss, re[0].line, re[1].line,
			xinterfaces[0].line, re[2].line, re[3].line, xinterfaces[0].length);

		index[0] = 2 * (block.ny + 1); index[1] = (block.ny + 1); index[2] = -(block.ny + 1); index[3] = -2 * (block.ny + 1);
		for (int iface = 0; iface < 4; ++iface)
		{
			Local_to_Global(re[iface].line.center.convar, yinterfaces[index[iface]].line.center.convar,
				yinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.convar,
				yinterfaces[0].normal);
			Local_to_Global(re[iface].line.center.der1x, yinterfaces[index[iface]].line.center.der1x,
				yinterfaces[index[iface]].normal);
			Global_to_Local(re[iface].line.center.der1x,
				yinterfaces[0].normal);

		}
		Center_5th_Multi(yinterfaces[0].gauss, re[0].line, re[1].line,
			yinterfaces[0].line, re[2].line, re[3].line, yinterfaces[0].length);

	}
}
void Center_5th_Multi(Recon2d* re, Recon2d& wn2, Recon2d& wn1, Recon2d& w0, Recon2d& wp1, Recon2d& wp2, double h)
{

	for (int var = 0; var < 4; var++)
	{
		double coe[5];
		Polynomial_5th(coe, wn2.center.convar[var], wn1.center.convar[var],
			w0.center.convar[var], wp1.center.convar[var], wp2.center.convar[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].center.convar[var] = Value_Polynomial(5, re[num_gauss].center.x, coe);
			re[num_gauss].center.der1y[var] = Der1_Polynomial(5, re[num_gauss].center.x, coe);
			re[num_gauss].center.der2yy[var] = Der2_Polynomial(5, re[num_gauss].center.x, coe);
		}
		Polynomial_5th(coe, wn2.center.der1x[var], wn1.center.der1x[var],
			w0.center.der1x[var], wp1.center.der1x[var], wp2.center.der1x[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].center.der1x[var] = Value_Polynomial(5, re[num_gauss].center.x, coe);
			re[num_gauss].center.der2xy[var] = Der1_Polynomial(5, re[num_gauss].center.x, coe);
		}
		Polynomial_5th(coe, wn2.center.der2xx[var], wn1.center.der2xx[var],
			w0.center.der2xx[var], wp1.center.der2xx[var], wp2.center.der2xx[var], h);
		for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
		{
			re[num_gauss].center.der2xx[var] = Value_Polynomial(5, re[num_gauss].center.x, coe);
		}
	}

	for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
	{
		Check_Order_Reduce_by_Lambda_2D(re[num_gauss].center.is_reduce_order, re[num_gauss].center.convar);
		if (re[num_gauss].center.is_reduce_order == true)
		{
			for (int m = 0; m < 4; m++)
			{
				re[num_gauss].center.convar[m] = w0.center.convar[m];
				re[num_gauss].center.der1x[m] = w0.center.der1x[m];
				re[num_gauss].center.der1y[m] = 0.0;
				re[num_gauss].center.der2xx[m] = w0.center.der2xx[m];
				re[num_gauss].center.der2xy[m] = 0.0;
				re[num_gauss].center.der2yy[m] = 0.0;
			}
		}
	}
}
void Center_3rd_Multi(Point2d& gauss, Point2d& wn1, Point2d& w0, Point2d& wp1, double h)
{
	for (int var = 0; var < 4; var++)
	{
		double coe[3];
		Polynomial_3rd_avg(coe, wn1.convar[var], w0.convar[var], wp1.convar[var], h);
		gauss.convar[var] = Value_Polynomial_3rd(gauss.x, coe);
		gauss.der1y[var] = Der1_Polynomial_3rd(gauss.x, coe);
		gauss.der2yy[var] = Der2_Polynomial_3rd(coe);
		Polynomial_3rd_avg(coe, wn1.der1x[var], w0.der1x[var], wp1.der1x[var], h);
		gauss.der1x[var] = Value_Polynomial_3rd(gauss.x, coe);
		gauss.der2xy[var] = Der1_Polynomial_3rd(gauss.x, coe);
		Polynomial_3rd_avg(coe, wn1.der2xx[var], w0.der2xx[var], wp1.der2xx[var], h);
		gauss.der2xx[var] = Value_Polynomial_3rd(gauss.x, coe);
	}

	Check_Order_Reduce_by_Lambda_2D(gauss.is_reduce_order, gauss.convar);
	//if lambda <0, then reduce to the first order
	if (gauss.is_reduce_order == true)
	{
		for (int m = 0; m < 4; m++)
		{
			gauss.convar[m] = w0.convar[m];
			gauss.der1x[m] = w0.der1x[m];
			gauss.der1y[m] = 0.0;
			gauss.der2xx[m] = w0.der2xx[m];
			gauss.der2xy[m] = 0.0;
			gauss.der2yy[m] = 0.0;
		}
	}
}




void Calculate_flux(Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Interface2d* xinterfaces, Interface2d* yinterfaces, Block2d block, int stage)
{
	// Note : calculate the final flux of gauss points, in xfluxes,  by the reconstructed variables in xinterfaces; the same for y
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				flux_function_2d(xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss], xinterfaces[i * (block.ny + 1) + j].gauss[num_gauss], block.dt);
				// calculate the final flux in xfluxes,  by the reconstructed variables in xinterfaces; the same for y
			}
		}
	}
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				flux_function_2d(yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss], yinterfaces[i * (block.ny + 1) + j].gauss[num_gauss], block.dt);
			}
		}
	}
}

namespace
{
	const double gks2d_tau_eps = 1.0e-14;

	enum GKS2DMomentSpace
	{
		gks2d_full_moment,
		gks2d_left_half_moment,
		gks2d_right_half_moment
	};

	struct GKS2DLocalThirdOrderState
	{
		double prim[4];
		double ax[4];
		double ay[4];
		double at[4];
		double axx[4];
		double axy[4];
		double ayy[4];
		double axt[4];
		double ayt[4];
		double att[4];
		double tau_phys;
	};

	struct GKS3rdTimeCoeff2D
	{
		double eq_g;
		double eq_grad;
		double eq_A;
		double eq_quad;
		double eq_grad_t;
		double eq_tt;
	};

	double GetTauNumericalDissipation2D(const double* prim_left, const double* prim_right, double dt)
	{
		double denom = fabs(prim_left[0] / prim_left[3] + prim_right[0] / prim_right[3]);
		double jump = 0.0;
		if (denom > 1.0e-14)
		{
			jump = fabs(prim_left[0] / prim_left[3] - prim_right[0] / prim_right[3]) / denom;
		}
		if (tau_type == Euler)
		{
			if (c1_euler <= 0.0 && c2_euler <= 0.0)
			{
				return 0.0;
			}
			return c1_euler * dt + c2_euler * jump * dt;
		}
		if (tau_type == NS || tau_type == Sutherland || tau_type == Power_law)
		{
			double tau_n = c2_euler * jump * dt;
			if (tau_n != tau_n)
			{
				tau_n = 0.0;
			}
			return tau_n;
		}
		return 0.0;
	}

	double GetTauPhysFromPrim2D(const double* prim)
	{
		if ((tau_type != NS && tau_type != Sutherland && tau_type != Power_law) ||
			((Mu <= 0.0) && (Nu <= 0.0)))
		{
			return 0.0;
		}
		if ((Mu > 0.0 && Nu > 0.0) || (Mu < 0.0 && Nu < 0.0))
		{
			return 0.0;
		}
		return Get_Tau_NS(prim[0], prim[3]);
	}

	void GetExpTimeIntegrals2D(double T, double tau_n, double& e0, double& e1, double& e2)
	{
		if (tau_n <= gks2d_tau_eps)
		{
			e0 = 0.0;
			e1 = 0.0;
			e2 = 0.0;
			return;
		}
		const double eta = exp(-T / tau_n);
		e0 = tau_n * (1.0 - eta);
		e1 = tau_n * tau_n - tau_n * (T + tau_n) * eta;
		e2 = tau_n * tau_n * tau_n - 0.5 * tau_n * eta *
			(T * T + 2.0 * tau_n * T + 2.0 * tau_n * tau_n);
	}

	void GetGKS3rdEquilibriumTimeCoeff2D(GKS3rdTimeCoeff2D& tc, double T, double tau_n)
	{
		if (tau_n <= gks2d_tau_eps)
		{
			tc.eq_g = T;
			tc.eq_grad = 0.0;
			tc.eq_A = 0.5 * T * T;
			tc.eq_quad = 0.0;
			tc.eq_grad_t = 0.0;
			tc.eq_tt = T * T * T / 6.0;
			return;
		}

		const double tau = tau_n;
		const double eta = exp(-T / tau);
		tc.eq_g = T - tau + tau * eta;
		tc.eq_grad = 2.0 * tau * tau * (1.0 - eta) - tau * T * (1.0 + eta);
		tc.eq_A = 0.5 * T * T - tau * T + tau * tau * (1.0 - eta);
		tc.eq_quad = tau * tau * T - 3.0 * tau * tau * tau
			+ 0.5 * tau * eta * (T * T + 4.0 * tau * T + 6.0 * tau * tau);
		tc.eq_grad_t = -0.5 * tau * T * T + 2.0 * tau * tau * T - 3.0 * tau * tau * tau
			+ tau * tau * (T + 3.0 * tau) * eta;
		tc.eq_tt = T * T * T / 6.0 - 0.5 * tau * T * T + tau * tau * T - tau * tau * tau * (1.0 - eta);
	}

	void SolveCoeffFromMoment2D(double* coeff, const double* moment, const double* prim)
	{
		double rhs[4];
		double prim_local[4];
		Copy_Array(rhs, const_cast<double*>(moment), 4);
		Copy_Array(prim_local, const_cast<double*>(prim), 4);
		A(coeff, rhs, prim_local);
	}

	void AddMoment2D(
		double* accum,
		GKS2DMomentSpace space,
		double rho,
		MMDF& moment,
		int no_u,
		int no_v,
		const double* coeff,
		double scale)
	{
		if (fabs(scale) <= 1.0e-30)
		{
			return;
		}
		double coeff_local[4];
		double psi[4];
		Copy_Array(coeff_local, const_cast<double*>(coeff), 4);
		if (space == gks2d_full_moment)
		{
			G_address(no_u, no_v, 0, psi, coeff_local, moment);
		}
		else if (space == gks2d_left_half_moment)
		{
			GL_address(no_u, no_v, 0, psi, coeff_local, moment);
		}
		else
		{
			GR_address(no_u, no_v, 0, psi, coeff_local, moment);
		}
		for (int m = 0; m < 4; ++m)
		{
			accum[m] += scale * rho * psi[m];
		}
	}

	void SolveGKSAt2D(double* at, const double* ax, const double* ay, const double* prim)
	{
		double rhs[4] = { 0.0, 0.0, 0.0, 0.0 };
		MMDF moment(prim[1], prim[2], prim[3]);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 1, 0, ax, -1.0);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 0, 1, ay, -1.0);
		SolveCoeffFromMoment2D(at, rhs, prim);
	}

	void SolveGKSAxt2D(double* axt, const double* axx, const double* axy, const double* prim)
	{
		double rhs[4] = { 0.0, 0.0, 0.0, 0.0 };
		MMDF moment(prim[1], prim[2], prim[3]);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 1, 0, axx, -1.0);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 0, 1, axy, -1.0);
		SolveCoeffFromMoment2D(axt, rhs, prim);
	}

	void SolveGKSAyt2D(double* ayt, const double* axy, const double* ayy, const double* prim)
	{
		double rhs[4] = { 0.0, 0.0, 0.0, 0.0 };
		MMDF moment(prim[1], prim[2], prim[3]);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 1, 0, axy, -1.0);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 0, 1, ayy, -1.0);
		SolveCoeffFromMoment2D(ayt, rhs, prim);
	}

	void SolveGKSAtt2D(double* att, const double* axt, const double* ayt, const double* prim)
	{
		double rhs[4] = { 0.0, 0.0, 0.0, 0.0 };
		MMDF moment(prim[1], prim[2], prim[3]);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 1, 0, axt, -1.0);
		AddMoment2D(rhs, gks2d_full_moment, prim[0], moment, 0, 1, ayt, -1.0);
		SolveCoeffFromMoment2D(att, rhs, prim);
	}

	void BuildGKS2DThirdOrderState(
		GKS2DLocalThirdOrderState& state,
		const Point2d& point,
		double tau_phys,
		bool need_att)
	{
		double convar[4];
		Copy_Array(convar, const_cast<double*>(point.convar), 4);
		Convar_to_ULambda_2d(state.prim, convar);
		SolveCoeffFromMoment2D(state.ax, point.der1x, state.prim);
		SolveCoeffFromMoment2D(state.ay, point.der1y, state.prim);
		SolveGKSAt2D(state.at, state.ax, state.ay, state.prim);
		SolveCoeffFromMoment2D(state.axx, point.der2xx, state.prim);
		SolveCoeffFromMoment2D(state.axy, point.der2xy, state.prim);
		SolveCoeffFromMoment2D(state.ayy, point.der2yy, state.prim);
		SolveGKSAxt2D(state.axt, state.axx, state.axy, state.prim);
		SolveGKSAyt2D(state.ayt, state.axy, state.ayy, state.prim);
		if (need_att)
		{
			SolveGKSAtt2D(state.att, state.axt, state.ayt, state.prim);
		}
		else
		{
			Array_zero(state.att, 4);
		}
		state.tau_phys = tau_phys;
	}

	void AddInitialTransportContribution2D(
		double* flux,
		const GKS2DLocalThirdOrderState& state,
		GKS2DMomentSpace space,
		double e0,
		double e1,
		double e2)
	{
		double unit[4] = { 1.0, 0.0, 0.0, 0.0 };
		MMDF moment(state.prim[1], state.prim[2], state.prim[3]);
		const double second_scale = state.tau_phys * e1 + e2;

		AddMoment2D(flux, space, state.prim[0], moment, 1, 0, unit, e0);
		AddMoment2D(flux, space, state.prim[0], moment, 1, 0, state.at, -state.tau_phys * e0);
		AddMoment2D(flux, space, state.prim[0], moment, 2, 0, state.ax, -state.tau_phys * e0 - e1);
		AddMoment2D(flux, space, state.prim[0], moment, 1, 1, state.ay, -state.tau_phys * e0 - e1);
		AddMoment2D(flux, space, state.prim[0], moment, 2, 0, state.axt, state.tau_phys * e1);
		AddMoment2D(flux, space, state.prim[0], moment, 1, 1, state.ayt, state.tau_phys * e1);
		AddMoment2D(flux, space, state.prim[0], moment, 3, 0, state.axx, second_scale);
		AddMoment2D(flux, space, state.prim[0], moment, 2, 1, state.axy, 2.0 * second_scale);
		AddMoment2D(flux, space, state.prim[0], moment, 1, 2, state.ayy, second_scale);
	}

	void ComputeGKS3rdFluxIntegral2DLocal(
		double* flux,
		const GKS2DLocalThirdOrderState& gbar_state,
		const GKS2DLocalThirdOrderState& gl_state,
		const GKS2DLocalThirdOrderState& gr_state,
		double T,
		double tau_num)
	{
		Array_zero(flux, 4);
		double e0, e1, e2;
		GetExpTimeIntegrals2D(T, tau_num, e0, e1, e2);
		GKS3rdTimeCoeff2D tc_eq;
		GetGKS3rdEquilibriumTimeCoeff2D(tc_eq, T, tau_num);

		const double T2 = T * T;
		double unit[4] = { 1.0, 0.0, 0.0, 0.0 };
		MMDF moment_bar(gbar_state.prim[1], gbar_state.prim[2], gbar_state.prim[3]);

		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 1, 0, unit, tc_eq.eq_g);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 2, 0, gbar_state.ax,
			tc_eq.eq_grad - gbar_state.tau_phys * T);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 1, 1, gbar_state.ay,
			tc_eq.eq_grad - gbar_state.tau_phys * T);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 1, 0, gbar_state.at,
			tc_eq.eq_A - gbar_state.tau_phys * T);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 3, 0, gbar_state.axx,
			tc_eq.eq_quad);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 2, 1, gbar_state.axy,
			2.0 * tc_eq.eq_quad);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 1, 2, gbar_state.ayy,
			tc_eq.eq_quad);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 2, 0, gbar_state.axt,
			tc_eq.eq_grad_t - 0.5 * gbar_state.tau_phys * T2);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 1, 1, gbar_state.ayt,
			tc_eq.eq_grad_t - 0.5 * gbar_state.tau_phys * T2);
		AddMoment2D(flux, gks2d_full_moment, gbar_state.prim[0], moment_bar, 1, 0, gbar_state.att,
			tc_eq.eq_tt - 0.5 * gbar_state.tau_phys * T2);

		AddInitialTransportContribution2D(flux, gl_state, gks2d_left_half_moment, e0, e1, e2);
		AddInitialTransportContribution2D(flux, gr_state, gks2d_right_half_moment, e0, e1, e2);
	}

	void RecoverGKS3rdFluxTimeDerivatives2D(
		double* F,
		double* Ft,
		double* Ftt,
		const double* I_dt,
		const double* I_2dt3,
		const double* I_dt3,
		double dt)
	{
		for (int m = 0; m < 4; ++m)
		{
			F[m] = (I_dt[m] - 4.5 * I_2dt3[m] + 9.0 * I_dt3[m]) / dt;
			Ft[m] = -9.0 * (I_dt[m] - 4.0 * I_2dt3[m] + 5.0 * I_dt3[m]) / (dt * dt);
			Ftt[m] = 9.0 * (3.0 * I_dt[m] - 9.0 * I_2dt3[m] + 9.0 * I_dt3[m]) / (dt * dt * dt);
		}
	}
}

void GKS2D_smooth(Flux2d& flux, Recon2d& interface, double dt)
{
	if (gks2dsolver == nothing_2d)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	for (int m = 0; m < 4; ++m)
	{
		flux.der2f[m] = 0.0;
	}
	double Flux[2][4];
	//change conservative variables to rho u lambda
	double convar0[4];
	for (int i = 0; i < 4; i++)
	{

		convar0[i] = interface.center.convar[i];
	}
	//cout << endl;
	double prim0[4];
	Convar_to_ULambda_2d(prim0, convar0);

	//then lets get the coefficient of time intergation factors

	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[3]);
	tau_num = tau;
	double eta = 0.0;
	double t[10];
	// non equ part time coefficient for gks_2nd algorithm
	t[0] = tau_num * (1 - eta); // this refers glu, gru part
	t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part

	t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
								  // then, equ part time coefficient for gks 2nd
	t[3] = dt; //this refers g0u part
	t[4] = -dt * tau; //this refers a0uu part
	t[5] = 0.5 * dt * dt - tau * dt; //this refers A0u part


	double unit[4] = { 1, 0.0, 0.0, 0.0 };



	//only one part, the kfvs1st part
	for (int i = 0; i < 4; i++)
	{
		Flux[0][i] = 0.0;
	}

	//now the equ part added, m0 term added, gks1st part begin
	MMDF m0(prim0[1], prim0[2], prim0[3]);

	double g0u[4];
	G_address(1, 0, 0, g0u, unit, m0);

	//the equ g0u part, the gks1st part
	for (int i = 0; i < 4; i++)
	{
		Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
	}

	// then we still need t[2], t[4] t[5] part for gks 2nd

	// for t[4] a0xuu part

	double a0x[4];
	double der1[4];
	for (int i = 0; i < 4; i++)
	{
		der1[i] = interface.center.der1x[i];
	}
	//solve the microslope
	A(a0x, der1, prim0);
	//a0x <u> moment
	double a0xu[4];
	G_address(1, 0, 0, a0xu, a0x, m0);
	//a0x <u^2> moment
	double a0xuu[4];
	G_address(2, 0, 0, a0xuu, a0x, m0);

	double a0y[4];
	double dery[4];
	for (int i = 0; i < 4; i++)
	{
		dery[i] = interface.center.der1y[i];
	}
	A(a0y, dery, prim0);
	double a0yv[4];
	G_address(0, 1, 0, a0yv, a0y, m0);
	//a0x <u^2> moment
	double a0yuv[4];
	G_address(1, 1, 0, a0yuv, a0y, m0);

	for (int i = 0; i < 4; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i] + a0yuv[i]);
	}


	// for t[5] A0u part
	double derA0[4];

	for (int i = 0; i < 4; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i] + a0yv[i]);
	}
	double A0[4];
	A(A0, derA0, prim0);
	double A0u[4];
	G_address(1, 0, 0, A0u, A0, m0);
	for (int i = 0; i < 4; ++i)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}
	if (gks2dsolver == gks2nd_2d && timecoe_list_2d == S1O1_2D)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
			flux.derf[i] = 0.0;
		}

		return;
	}


	if (gks2dsolver == gks2nd_2d && (timecoe_list_2d != S1O1_2D))
	{
		double dt2 = 0.5 * dt;
		t[3] = dt2; //this refers g0u part
		t[4] = dt2 * tau; //this refers a0uu part
		t[5] = 0.5 * dt2 * dt2 - tau * dt2; //this refers A0u part

		for (int i = 0; i < 4; i++)
		{
			// t0 part
			Flux[1][i] = 0.0;
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i] + a0yuv[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}



		for (int i = 0; i < 4; i++)
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
void GKS2D(Flux2d& flux, Recon2d& interface, double dt)
{
	// Note : // calculate the final flux in flux,  by the reconstructed variables in interface
	if (gks2dsolver == nothing_2d)
	{
		cout << "no gks solver specify" << endl;
		exit(0);
	}
	for (int m = 0; m < 4; ++m)
	{
		flux.der2f[m] = 0.0;
	}
	if (gks2dsolver == gks3rd_2d)
	{
		double convar_left[4], convar_right[4], convar_bar[4];
		double prim_left[4], prim_right[4], prim_bar[4];
		for (int m = 0; m < 4; ++m)
		{
			convar_left[m] = interface.left.convar[m];
			convar_right[m] = interface.right.convar[m];
			convar_bar[m] = interface.center.convar[m];
		}
		Convar_to_ULambda_2d(prim_left, convar_left);
		Convar_to_ULambda_2d(prim_right, convar_right);
		Convar_to_ULambda_2d(prim_bar, convar_bar);

		const double tau_phys_left = GetTauPhysFromPrim2D(prim_left);
		const double tau_phys_right = GetTauPhysFromPrim2D(prim_right);
		const double tau_phys_bar = GetTauPhysFromPrim2D(prim_bar);
		const double tau_num = tau_phys_bar + GetTauNumericalDissipation2D(prim_left, prim_right, dt);

		GKS2DLocalThirdOrderState left_state, right_state, bar_state;
		BuildGKS2DThirdOrderState(left_state, interface.left, tau_phys_left, false);
		BuildGKS2DThirdOrderState(right_state, interface.right, tau_phys_right, false);
		BuildGKS2DThirdOrderState(bar_state, interface.center, tau_phys_bar, true);

		double I_dt[4];
		double I_2dt3[4];
		double I_dt3[4];
		double F0[4];
		double Ft[4];
		double Ftt[4];
		ComputeGKS3rdFluxIntegral2DLocal(I_dt, bar_state, left_state, right_state, dt, tau_num);
		ComputeGKS3rdFluxIntegral2DLocal(I_2dt3, bar_state, left_state, right_state, 2.0 * dt / 3.0, tau_num);
		ComputeGKS3rdFluxIntegral2DLocal(I_dt3, bar_state, left_state, right_state, dt / 3.0, tau_num);
		RecoverGKS3rdFluxTimeDerivatives2D(F0, Ft, Ftt, I_dt, I_2dt3, I_dt3, dt);

		for (int m = 0; m < 4; ++m)
		{
			flux.f[m] = I_dt[m];
			flux.derf[m] = Ft[m];
			flux.der2f[m] = Ftt[m];
		}
		return;
	}
	double Flux[2][4];
	double convar_left[4], convar_right[4], convar0[4];
	for (int i = 0; i < 4; ++i)
	{
		convar_left[i] = interface.left.convar[i];
		convar_right[i] = interface.right.convar[i];
	}
	//cout << endl;
	//change conservative variables to rho u lambda
	double prim_left[4], prim_right[4], prim0[4];
	Convar_to_ULambda_2d(prim_left, convar_left);
	Convar_to_ULambda_2d(prim_right, convar_right);

	//prim_left[1], prim_left[2], prim_left[3], means U, V, Lambda
	MMDF ml(prim_left[1], prim_left[2], prim_left[3]);
	MMDF mr(prim_right[1], prim_right[2], prim_right[3]);

	if (g0reconstruction_2D_tangent == Center_all_collision_multi)
	{
		Collision(interface.center.convar, prim_left[0], prim_right[0], ml, mr); // ml, mr, pass the whole class into the function
		double a0[4] = { 1.0, 0.0, 0.0, 0.0 };
		double alx_t[4], arx_t[4];
		//w_x
		A_point(alx_t, interface.left.der1x, prim_left); // input the slope of macroscopic variables, output the a coefficient
		A_point(arx_t, interface.right.der1x, prim_right);
		double al0x_t[4];
		double ar0x_t[4];
		GL_address(0, 0, 0, al0x_t, alx_t, ml); // ml, mr, pass the whole class into the function
		GR_address(0, 0, 0, ar0x_t, arx_t, mr); // ml, mr, pass the whole class into the function
		//GL_address, GR_address, G_address, are calculating the macroscopic variables by the microscopic variables based on moment calculation
		//Note: when the input of moment calculation is a coefficient, the result is the derivative of W; when the input is (0 0 0 0), the result is W
		//when the input is (1 0 0 0), the result is the relative flux, without considering the integral of time

		//w_y
		double aly_t[4], ary_t[4];
		A_point(aly_t, interface.left.der1y, prim_left);
		A_point(ary_t, interface.right.der1y, prim_right);
		//The difference of A_point with A is just the input form, matrix or pointer
		//The content, input variables, output variables are the same

		double al0y_t[4];
		double ar0y_t[4];
		GL_address(0, 0, 0, al0y_t, aly_t, ml);
		GR_address(0, 0, 0, ar0y_t, ary_t, mr);
		for (int var = 0; var < 4; ++var)
		{
			interface.center.der1x[var] = prim_left[0] * al0x_t[var] + prim_right[0] * ar0x_t[var];
			interface.center.der1y[var] = prim_left[0] * al0y_t[var] + prim_right[0] * ar0y_t[var];
		}
	}

	for (int i = 0; i < 4; ++i)
	{
		convar0[i] = interface.center.convar[i];
	}

	Convar_to_ULambda_2d(prim0, convar0);
	//then lets get the coefficient of time intergation factors

	double tau;
	double tau_num;
	tau = Get_Tau_NS(prim0[0], prim0[3]);
	tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[3], prim_right[3], prim0[3], dt);
	double eta = exp(-dt / tau_num);
	double t[10];
	// non equ part time coefficient for gks_2nd algorithm
	t[0] = tau_num * (1 - eta); // this refers glu, gru part
	t[1] = tau_num * (eta * (dt + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part

	t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
	// then, equ part time coefficient for gks 2nd
	t[3] = tau_num * eta + dt - tau_num; //this refers g0u part
	t[4] = tau_num * (tau_num - eta * (dt + tau_num) - tau * (eta - 1)) - dt * tau; //this refers a0uu part
	t[5] = 0.5 * dt * dt - tau * tau_num * (eta - 1) - tau * dt; //this refers A0u part

	//reset the t[0~5], used ONLY for kfvs method
	//in kfvs, only use time step deltat, but NOT collision time
	if (gks2dsolver == kfvs1st_2d)
	{
		t[0] = dt;
		for (int i = 1; i < 6; i++)
		{
			t[i] = 0.0;
		}
		//do nothing, kfvs1st only use t[0]=dt part;
	}
	else if (gks2dsolver == kfvs2nd_2d)
	{
		t[0] = dt;
		t[1] = -0.5 * dt * dt;
		for (int i = 2; i < 6; i++)
		{
			t[i] = 0.0;
		}
	}

	double unit[4] = { 1, 0.0, 0.0, 0.0 };
	// alx_t = a(0)*1 + a(1)*u + a(2)*v + a(3)*e; arx_t, aly_t, ary_t are the same
	// the slope of the (macroscopic) conservative variables value (/density) can be obtained by (in al0x_t), GL_address(0, 0, 0, al0x_t, alx_t, ml);
	// thus, when replace a by unit (=[1 0 0 0]), GL_address can get the relative (macroscopic) conservative variables (/density); but pay attention to the (1 0 0) but not (0 0 0)
	// G_address, GR_address, are the same with GL_address

	double glu[4], gru[4];

	GL_address(1, 0, 0, glu, unit, ml); // (1 0 0) get GLu
	GR_address(1, 0, 0, gru, unit, mr); // (1 0 0) get GRu

	//only one part, the kfvs1st part
	for (int i = 0; i < 4; i++)
	{
		Flux[0][i] = prim_left[0] * t[0] * glu[i] + prim_right[0] * t[0] * gru[i];
	}
	//cout<< endl;
	if (gks2dsolver == kfvs1st_2d)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
		}

		return;
	}
	// kfvs1st part ended

	//now the equ part added, m0 term added, gks1st part begin
	MMDF m0(prim0[1], prim0[2], prim0[3]);

	double g0u[4];
	G_address(1, 0, 0, g0u, unit, m0);

	if (gks2dsolver != kfvs2nd_2d)
	{
		//the equ g0u part, the gks1st part
		for (int i = 0; i < 4; i++)
		{
			Flux[0][i] = Flux[0][i] + prim0[0] * t[3] * g0u[i];
		}
	}

	if (gks2dsolver == gks1st_2d)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// gks1d solver ended

	//for kfvs2nd part
	double der1xleft[4], der1xright[4];
	double der1yleft[4], der1yright[4];
	for (int i = 0; i < 4; i++)
	{
		der1xleft[i] = interface.left.der1x[i];
		der1xright[i] = interface.right.der1x[i];
		der1yleft[i] = interface.left.der1y[i];
		der1yright[i] = interface.right.der1y[i];
	}

	double alx[4];
	A(alx, der1xleft, prim_left);
	double alxuul[4];
	GL_address(2, 0, 0, alxuul, alx, ml);

	double arx[4];
	A(arx, der1xright, prim_right);
	double arxuur[4];
	GR_address(2, 0, 0, arxuur, arx, mr);

	double aly[4];
	A(aly, der1yleft, prim_left);
	double alyuvl[4];
	GL_address(1, 1, 0, alyuvl, aly, ml);

	double ary[4];
	A(ary, der1yright, prim_right);
	double aryuvr[4];
	GR_address(1, 1, 0, aryuvr, ary, mr);

	for (int i = 0; i < 4; i++)
	{	// t1 part
		Flux[0][i] = Flux[0][i] + t[1] * (prim_left[0] * (alxuul[i] + alyuvl[i]) + prim_right[0] * (arxuur[i] + aryuvr[i]));
	}
	if (gks2dsolver == kfvs2nd_2d)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
		}
		return;
	}
	// the kfvs2nd part ended

	// then we still need t[2], t[4] t[5] part for gks 2nd
	//for t[2] Aru,Alu part
	double alxu[4];		double alyv[4];
	double arxu[4];		double aryv[4];

	//take <u> moment for al, ar
	G_address(1, 0, 0, alxu, alx, ml);
	G_address(1, 0, 0, arxu, arx, mr);
	G_address(0, 1, 0, alyv, aly, ml);
	G_address(0, 1, 0, aryv, ary, mr);

	double Al[4], Ar[4];
	double der_AL[4], der_AR[4];

	//using compatability condition to get the time derivative
	for (int i = 0; i < 4; i++)
	{
		der_AL[i] = -prim_left[0] * (alxu[i] + alyv[i]);
		der_AR[i] = -prim_right[0] * (arxu[i] + aryv[i]);
	}
	// solve the coefficient martix b=ma
	A(Al, der_AL, prim_left);
	A(Ar, der_AR, prim_right);

	//to obtain the Alu and Aru
	double Alul[4];
	double Arur[4];
	GL_address(1, 0, 0, Alul, Al, ml);
	GR_address(1, 0, 0, Arur, Ar, mr);

	for (int i = 0; i < 4; i++)
	{	// t2 part
		Flux[0][i] = Flux[0][i] + t[2] * (prim_left[0] * Alul[i] + prim_right[0] * Arur[i]);
	}

	// for t[4] a0xuu part
	double a0x[4];
	double der1[4];
	for (int i = 0; i < 4; i++)
	{
		der1[i] = interface.center.der1x[i]; //Only the averaged value for g0
	}
	//solve the microslope
	A(a0x, der1, prim0);
	//a0x <u> moment
	double a0xu[4];
	G_address(1, 0, 0, a0xu, a0x, m0); //get a0xu, used for the following determination of derA0, and then A0
	//a0x <u^2> moment
	double a0xuu[4];
	G_address(2, 0, 0, a0xuu, a0x, m0);

	double a0y[4];
	double dery[4];
	for (int i = 0; i < 4; i++)
	{
		dery[i] = interface.center.der1y[i];
	}
	A(a0y, dery, prim0);
	double a0yv[4];
	G_address(0, 1, 0, a0yv, a0y, m0); //get a0yv, used for the following determination of derA0, and then A0
	//a0x <u^2> moment
	double a0yuv[4];
	G_address(1, 1, 0, a0yuv, a0y, m0);

	for (int i = 0; i < 4; i++)
	{	// t4 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[4] * (a0xuu[i] + a0yuv[i]);
	}

	// for t[5] A0u part
	double derA0[4];
	for (int i = 0; i < 4; i++)
	{
		derA0[i] = -prim0[0] * (a0xu[i] + a0yv[i]);
	}
	double A0[4];
	A(A0, derA0, prim0);
	double A0u[4];
	G_address(1, 0, 0, A0u, A0, m0);
	for (int i = 0; i < 4; i++)
	{	// t5 part
		Flux[0][i] = Flux[0][i] + prim0[0] * t[5] * (A0u[i]);
	}

	if (is_Prandtl_fix == true)
	{
		double qs;
		qs = a0xuu[3] + A0u[3] - prim0[1] * (a0xuu[1] + A0u[1]) - prim0[2] * (a0xuu[2] + A0u[2]);
		qs *= -prim0[0] * tau * dt;
		Flux[0][3] += (1.0 / Pr - 1.0) * qs;
	}
	if (gks2dsolver == gks2nd_2d && timecoe_list_2d == S1O1_2D)
	{
		for (int i = 0; i < 4; i++)
		{
			flux.f[i] = Flux[0][i];
			flux.derf[i] = 0.0;
		}
		return;
	}

	if (gks2dsolver == gks2nd_2d && timecoe_list_2d != S1O1_2D)
	{
		double dt2 = 0.5 * dt; // the following is dt2
		tau_num = Get_Tau(prim_left[0], prim_right[0], prim0[0], prim_left[3], prim_right[3], prim0[3], dt2);
		eta = exp(-dt2 / tau_num);
		// non equ part time coefficient for gks_2nd algorithm
		t[0] = tau_num * (1 - eta); // this refers glu, gru part
		t[1] = tau_num * (eta * (dt2 + tau_num) - tau_num) + tau * tau_num * (eta - 1); //this refers aluu, aruu part
		t[2] = tau * tau_num * (eta - 1); //this refers Alu, Aru part
		// then, equ part time coefficient for gks 2nd
		t[3] = tau_num * eta + dt2 - tau_num; //this refers g0u part
		t[4] = tau_num * (tau_num - eta * (dt2 + tau_num) - tau * (eta - 1)) - dt2 * tau; //this refers a0uu part
		t[5] = 0.5 * dt2 * dt2 - tau * tau_num * (eta - 1) - tau * dt2; //this refers A0u part

		for (int i = 0; i < 4; i++)
		{
			// t0 part
			Flux[1][i] = t[0] * (prim_left[0] * glu[i] + prim_right[0] * gru[i]);
			// t1 part
			Flux[1][i] = Flux[1][i] + t[1] * (prim_left[0] * (alxuul[i] + alyuvl[i]) + prim_right[0] * (arxuur[i] + aryuvr[i]));
			// t2 part
			Flux[1][i] = Flux[1][i] + t[2] * (prim_left[0] * Alul[i] + prim_right[0] * Arur[i]);
			// t3 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[3] * g0u[i];
			// t4 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[4] * (a0xuu[i] + a0yuv[i]);
			// t5 part
			Flux[1][i] = Flux[1][i] + prim0[0] * t[5] * (A0u[i]);
		}

		if (is_Prandtl_fix == true)
		{
			double qs;
			qs = a0xuu[3] + A0u[3] - prim0[1] * (a0xuu[1] + A0u[1]) - prim0[2] * (a0xuu[2] + A0u[2]);
			qs *= -prim0[0] * tau * 0.5 * dt;
			Flux[1][3] += (1.0 / Pr - 1.0) * qs;
		}

		for (int i = 0; i < 4; i++)
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
void LF2D(Flux2d& flux, Recon2d& interface, double dt)
{
	double pl[4], pr[4];
	Convar_to_primvar_2D(pl, interface.left.convar);
	Convar_to_primvar_2D(pr, interface.right.convar);

	double k[2];
	k[0] = abs(pl[1]) + sqrt(Gamma * pl[3] / pl[0]);
	k[1] = abs(pr[1]) + sqrt(Gamma * pr[3] / pr[0]);
	double beta = k[0];
	if (k[1] > k[0]) { beta = k[1]; }
	double flux_l[4], flux_r[4];
	get_flux(pl, flux_l);
	get_flux(pr, flux_r);

	for (int m = 0; m < 4; m++)
	{
		flux.f[m] = 0.5 * ((flux_l[m] + flux_r[m]) - beta * (interface.right.convar[m] - interface.left.convar[m]));
		flux.f[m] *= dt;
	}
	if (tau_type == NS)
	{
		NS_by_central_difference_prim_2D(flux, interface, dt);
	}

}
void HLLC2D(Flux2d& flux, Recon2d& interface, double dt)
{
	double pl[4], pr[4];
	Convar_to_primvar_2D(pl, interface.left.convar);
	Convar_to_primvar_2D(pr, interface.right.convar);

	double al, ar, pvars, pstar, tmp1, tmp2, tmp3, qk, sl, sr, star;
	al = sqrt(Gamma * pl[3] / pl[0]); //sound speed
	ar = sqrt(Gamma * pr[3] / pr[0]);
	tmp1 = 0.5 * (al + ar);         //avg of sound and density
	tmp2 = 0.5 * (pl[0] + pr[0]);

	pvars = 0.5 * (pl[3] + pr[3]) - 0.5 * (pr[1] - pl[1]) * tmp1 * tmp2;
	pstar = fmax(0.0, pvars);

	double flxtmp[4], qstar[4];

	ESTIME(sl, star, sr, pl[0], pl[1], pl[3], al, pr[0], pr[1], pr[3], ar);

	tmp1 = pr[3] - pl[3] + pl[0] * pl[1] * (sl - pl[1]) - pr[0] * pr[1] * (sr - pr[1]);
	tmp2 = pl[0] * (sl - pl[1]) - pr[0] * (sr - pr[1]);
	star = tmp1 / tmp2;

	if (sl >= 0.0)
	{
		get_flux(pl, flux.f);
	}
	else if (sr <= 0.0)
	{
		get_flux(pr, flux.f);
	}
	else if ((star >= 0.0) && (sl <= 0.0))
	{
		get_flux(pl, flxtmp);
		ustarforHLLC(pl[0], pl[1], pl[2], pl[3], sl, star, qstar);

		for (int m = 0; m < 4; m++)
		{
			flux.f[m] = flxtmp[m] + sl * (qstar[m] - interface.left.convar[m]);
		}
	}
	else if ((star <= 0.0) && (sr >= 0.0))
	{
		get_flux(pr, flxtmp);
		ustarforHLLC(pr[0], pr[1], pr[2], pr[3], sr, star, qstar);
		for (int m = 0; m < 4; m++)
		{
			flux.f[m] = flxtmp[m] + sr * (qstar[m] - interface.right.convar[m]);
		}
	}
	else
	{
		cout << "couldnt be possible that hllc cannot give any result" << endl;
	}

	for (int m = 0; m < 4; m++)
	{
		flux.f[m] *= dt;
	}
	if (tau_type == NS)
	{
		NS_by_central_difference_convar_2D(flux, interface, dt);
	}
}
void NS_by_central_difference_prim_2D(Flux2d& flux, Recon2d& interface, double dt)
{
	// Note : Sixth-order centrial differential scheme for viscous term
	// here convar[i] represents density, u, v, and temperature
	double mu = Mu;
	if (Nu > 0) { mu = Nu * interface.center.convar[0]; }
	double u = interface.center.convar[1];
	double v = interface.center.convar[2];
	double ux, uy, vx, vy;
	ux = interface.center.der1x[1];
	vx = interface.center.der1x[2];
	uy = interface.center.der1y[1];
	vy = interface.center.der1y[2];
	double tau_xx = 2 * mu * ux - 2.0 / 3.0 * mu * (ux + vy);
	double tau_xy = mu * (uy + vx);
	double q = u * tau_xx + v * tau_xy + (K + 4) / (2 * Pr) * mu * interface.center.der1x[3];
	flux.f[1] += tau_xx * dt;
	flux.f[2] += tau_xy * dt;
	flux.f[3] += q * dt;
}
void NS_by_central_difference_convar_2D(Flux2d& flux, Recon2d& interface, double dt)
{
	// Note: Write a function calculating the flux of viscous term in NS, used for HLLC Riemann solver
	// If g0type = collisionnless, means the conservative variables used for viscous flux are those of Gauss points, which were interpolated by the line-averaged variables in the tangential direction
	// If g0type = collisionn, means, the conservative variables used for viscous flux are averaged by the left and right variables (line-ageraged)
	// This is because of the first method might cause negative variables for density or pressure, which is not allowed.
	double convar[4];
	for (int m = 0; m < 4; m++)
	{
		convar[m] = 0.5 * (interface.left.convar[m] + interface.right.convar[m]);
	}

	// here convar[i] represents the conservative variables
	double den = convar[0];
	double mu = Mu;
	if (Nu > 0) { mu = Nu * den; }
	double u = convar[1] / den;
	double v = convar[2] / den;
	double ux, uy, vx, vy;
	ux = (interface.center.der1x[1] - interface.center.der1x[0] * u) / den;
	vx = (interface.center.der1x[2] - interface.center.der1x[0] * v) / den;
	uy = (interface.center.der1y[1] - interface.center.der1y[0] * u) / den;
	vy = (interface.center.der1y[2] - interface.center.der1y[0] * v) / den;
	double tau_xx = 2 * mu * ux - 2.0 / 3.0 * mu * (ux + vy);
	double tau_xy = mu * (uy + vx);
	double Tx = interface.center.der1x[3] / den - u * ux - u * vx
		- convar[3] * interface.center.der1x[0] / den / den;
	double q = u * tau_xx + v * tau_xy + (K + 4) / (2.0 * Pr) * mu * Tx;
	flux.f[1] -= tau_xx * dt;
	flux.f[2] -= tau_xy * dt;
	flux.f[3] -= q * dt;
}

//forward_euler
void S1O1_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;

}
void S1O2_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.5;

}
void S1O3_2D(Block2d& block)
{
	block.stages = 1;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[0][0][1] = 0.0;
	block.timecoefficient[0][0][2] = 0.0;
}
void RK2_2D(Block2d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 1.0;
	block.timecoefficient[1][0][0] = 0.5;
	block.timecoefficient[1][1][0] = 0.5;

}
void S2O4_2D(Block2d& block)
{
	block.stages = 2;
	block.timecoefficient[0][0][0] = 0.5;
	block.timecoefficient[0][0][1] = 1.0 / 8.0;
	block.timecoefficient[1][0][0] = 1.0;
	block.timecoefficient[1][1][0] = 0.0;
	block.timecoefficient[1][0][1] = 1.0 / 6.0;
	block.timecoefficient[1][1][1] = 1.0 / 3.0;

}
void RK4_2D(Block2d& block)
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
void Initial_stages(Block2d& block)
{
	for (int i = 0; i < 5; i++) //refers the n stage
	{
		for (int j = 0; j < 5; j++) //refers the nth coefficient at n stage
		{
			for (int k = 0; k < 3; k++) //refers f, derf, der2f
			{
				block.timecoefficient[i][j][k] = 0.0;
			}
		}
	}
	timecoe_list_2d(block);
}



Fluid2d* Setfluid(Block2d& block)
{
	Fluid2d* var = new Fluid2d[block.nx * block.ny]; // dynamic variable (since block.nx is not determined)
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			var[i * block.ny + j].xindex = i;
			var[i * block.ny + j].yindex = j;
		}
	}

	cout << "fluid variable allocate done..." << endl;
	return var;
}
Interface2d* Setinterface_array(Block2d block)
{
	Interface2d* var = new Interface2d[(block.nx + 1) * (block.ny + 1)];  // dynamic variable (since block.nx is not determined)
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				var[i * (block.ny + 1) + j].line.left.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.left.der1y[k] = 0.0;
				var[i * (block.ny + 1) + j].line.left.der2xx[k] = 0.0;
				var[i * (block.ny + 1) + j].line.left.der2xy[k] = 0.0;
				var[i * (block.ny + 1) + j].line.left.der2yy[k] = 0.0;

				var[i * (block.ny + 1) + j].line.right.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.right.der1y[k] = 0.0;
				var[i * (block.ny + 1) + j].line.right.der2xx[k] = 0.0;
				var[i * (block.ny + 1) + j].line.right.der2xy[k] = 0.0;
				var[i * (block.ny + 1) + j].line.right.der2yy[k] = 0.0;

				var[i * (block.ny + 1) + j].line.center.der1x[k] = 0.0;
				var[i * (block.ny + 1) + j].line.center.der1y[k] = 0.0;
				var[i * (block.ny + 1) + j].line.center.der2xx[k] = 0.0;
				var[i * (block.ny + 1) + j].line.center.der2xy[k] = 0.0;
				var[i * (block.ny + 1) + j].line.center.der2yy[k] = 0.0;
			}
		}
	}
	cout << "interface variable allocate done..." << endl;
	return var;
}
Flux2d_gauss** Setflux_gauss_array(Block2d block)
{
	Flux2d_gauss** var = new Flux2d_gauss * [(block.nx + 1) * (block.ny + 1)];  // dynamic variable (since block.nx is not determined)

	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			// for m th step time marching schemes, m subflux needed
			var[i * (block.ny + 1) + j] = new Flux2d_gauss[block.stages];
		}
	}

	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			for (int k = 0; k < block.stages; k++)
			{
				if (gausspoint == 0)
				{
					var[i * (block.ny + 1) + j][k].gauss = new Flux2d[1];
					for (int m = 0; m < 4; m++)
					{
						var[i * (block.ny + 1) + j][k].gauss[0].f[m] = 0.0;
						var[i * (block.ny + 1) + j][k].gauss[0].derf[m] = 0.0;
						var[i * (block.ny + 1) + j][k].gauss[0].der2f[m] = 0.0;

					}
				}
				else
				{
					var[i * (block.ny + 1) + j][k].gauss = new Flux2d[gausspoint];
					for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
					{
						for (int m = 0; m < 4; m++)
						{
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].f[m] = 0.0;
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[m] = 0.0;
							var[i * (block.ny + 1) + j][k].gauss[num_gauss].der2f[m] = 0.0;
						}
					}
				}
			}
		}
	}
	if (var == 0)
	{
		cout << "fluid variable allocate fail...";
		return NULL;
	}
	cout << "flux with gausspoint variable allocate done..." << endl;
	return var;

}




void SetUniformMesh(Block2d& block, Fluid2d* fluids, Interface2d* xinterfaces, Interface2d* yinterfaces, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes)
{
	//nodex, nodey are the real node
	//interface number = cell number + 1
	block.dx = (block.right - block.left) / block.nodex;
	block.dy = (block.up - block.down) / block.nodey;
	block.overdx = 1 / block.dx;
	block.overdy = 1 / block.dy;

	block.xcell_begin = block.ghost;
	block.xcell_end = block.ghost + block.nodex - 1;
	block.ycell_begin = block.ghost;
	block.ycell_end = block.ghost + block.nodey - 1;

	block.xinterface_begin_n = block.ghost;
	block.xinterface_end_n = block.ghost + block.nodex;
	block.xinterface_begin_t = block.ghost;
	block.xinterface_end_t = block.ghost + block.nodex - 1;

	block.yinterface_begin_n = block.ghost;
	block.yinterface_end_n = block.ghost + block.nodey;
	block.yinterface_begin_t = block.ghost;
	block.yinterface_end_t = block.ghost + block.nodey - 1;

	//cell avg information
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			//two dimension geometry to one dimension store matrix, y direciton first and x direction second
			fluids[i * block.ny + j].dx = block.dx; //cell size
			fluids[i * block.ny + j].dy = block.dy; //cell size
			fluids[i * block.ny + j].coordx = block.left + (i + 0.5 - block.ghost) * block.dx; //cell center location
			fluids[i * block.ny + j].coordy = block.down + (j + 0.5 - block.ghost) * block.dy; //cell center location
			fluids[i * block.ny + j].area = block.dx * block.dy;
			fluids[i * block.ny + j].node[0] = block.left + (i - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[1] = block.down + (j - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[2] = block.left + (i + 1 - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[3] = block.down + (j - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[4] = block.left + (i + 1 - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[5] = block.down + (j + 1 - block.ghost) * block.dy;
			fluids[i * block.ny + j].node[6] = block.left + (i - block.ghost) * block.dx;
			fluids[i * block.ny + j].node[7] = block.down + (j + 1 - block.ghost) * block.dy;
		}
	}

	// interface information
	for (int i = 0; i <= block.nx; i++)
	{
		for (int j = 0; j <= block.ny; j++)
		{
			xinterfaces[i * (block.ny + 1) + j].x = block.left + (i - block.ghost) * block.dx;
			xinterfaces[i * (block.ny + 1) + j].y = block.down + (j - block.ghost + 0.5) * block.dy;
			xinterfaces[i * (block.ny + 1) + j].length = block.dy;
			xinterfaces[i * (block.ny + 1) + j].normal[0] = 1.0;
			xinterfaces[i * (block.ny + 1) + j].normal[1] = 0.0;

			Copy_geo_from_interface_to_line(xinterfaces[i * (block.ny + 1) + j]);
			xinterfaces[i * (block.ny + 1) + j].gauss = new Recon2d[gausspoint];

			Copy_geo_from_interface_to_flux
			(xinterfaces[i * (block.ny + 1) + j], xfluxes[i * (block.ny + 1) + j], block.stages);

			yinterfaces[i * (block.ny + 1) + j].y = block.down + (j - block.ghost) * block.dy;
			yinterfaces[i * (block.ny + 1) + j].x = block.left + (i - block.ghost + 0.5) * block.dx;
			yinterfaces[i * (block.ny + 1) + j].length = block.dx;
			yinterfaces[i * (block.ny + 1) + j].normal[0] = 0.0;
			yinterfaces[i * (block.ny + 1) + j].normal[1] = 1.0;

			Copy_geo_from_interface_to_line(yinterfaces[i * (block.ny + 1) + j]);

			yinterfaces[i * (block.ny + 1) + j].gauss = new Recon2d[gausspoint];

			Copy_geo_from_interface_to_flux
			(yinterfaces[i * (block.ny + 1) + j], yfluxes[i * (block.ny + 1) + j], block.stages);

			Set_Gauss_Coordinate(xinterfaces[i * (block.ny + 1) + j], yinterfaces[i * (block.ny + 1) + j]);
		}
	}
	cout << "set uniform information done..." << endl;
}
void Copy_geo_from_interface_to_line(Interface2d& interface)
{
	interface.line.x = interface.x;
	interface.line.y = interface.y;
	interface.line.normal[0] = interface.normal[0];
	interface.line.normal[1] = interface.normal[1];
}
void Copy_geo_from_interface_to_flux(Interface2d& interface, Flux2d_gauss* flux, int stages)
{
	for (int istage = 0; istage < stages; istage++)
	{
		if (gausspoint == 0)
		{
			int igauss = 0;
			flux[istage].gauss[igauss].normal[0] = interface.normal[0];
			flux[istage].gauss[igauss].normal[1] = interface.normal[1];
			flux[istage].gauss[igauss].length = interface.length;
		}
		else
		{
			for (int igauss = 0; igauss < gausspoint; igauss++)
			{
				flux[istage].gauss[igauss].normal[0] = interface.normal[0];
				flux[istage].gauss[igauss].normal[1] = interface.normal[1];
				flux[istage].gauss[igauss].length = interface.length;
			}
		}
	}
}
void Set_Gauss_Coordinate(Interface2d& xinterface, Interface2d& yinterface)
{
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		// first is gauss parameter
		xinterface.gauss[num_guass].x = xinterface.x;
		xinterface.gauss[num_guass].y = xinterface.y + gauss_loc[num_guass] * 0.5 * xinterface.length;
		xinterface.gauss[num_guass].normal[0] = xinterface.normal[0];
		xinterface.gauss[num_guass].normal[1] = xinterface.normal[1];
		// each gauss point contain left center right point
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].left, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].right, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].center, num_guass, xinterface.length);

		yinterface.gauss[num_guass].x = yinterface.x + gauss_loc[num_guass] * 0.5 * yinterface.length;
		yinterface.gauss[num_guass].y = yinterface.y;
		yinterface.gauss[num_guass].normal[0] = yinterface.normal[0];
		yinterface.gauss[num_guass].normal[1] = yinterface.normal[1];
		// each gauss point contain left center right point
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].left, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].right, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].center, num_guass, yinterface.length);

	}
}
void Set_Gauss_Intergation_Location_x(Point2d& xgauss, int index, double h)
{
	xgauss.x = gauss_loc[index] * h / 2.0;
}
void Set_Gauss_Intergation_Location_y(Point2d& ygauss, int index, double h)
{
	ygauss.x = gauss_loc[index] * h / 2.0;
}



void Update(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	// Note : write the function to update the conservative variables of each cell
	if (stage > block.stages)
	{
		cout << "wrong middle stage,pls check the time marching setting" << endl;
		exit(0);
	}

	double dt = block.dt;

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double Flux = 0.0;
					for (int k = 0; k < stage + 1; k++)
					{
	
							Flux = Flux
							+ gauss_weight[num_gauss] *
							(block.timecoefficient[stage][k][0] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
								+ block.timecoefficient[stage][k][1] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]
								+ block.timecoefficient[stage][k][2] * xfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].der2f[var]);

					}
					xfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux;
					// calculate the final flux of the interface, in x in xfluxes, by the obtained the flux and its derivative (f, def, der2f) at guass points, in xfluxes, and the corresponding weight factors
					// calculate by several stages according to the time marching method. same for yfluxes
				}
			}
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				for (int var = 0; var < 4; var++)
				{
					double Flux = 0.0;
					for (int k = 0; k < stage + 1; k++)
					{
							Flux = Flux
							+ gauss_weight[num_gauss] *
							(block.timecoefficient[stage][k][0] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].f[var]
								+ block.timecoefficient[stage][k][1] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].derf[var]
								+ block.timecoefficient[stage][k][2] * yfluxes[i * (block.ny + 1) + j][k].gauss[num_gauss].der2f[var]);
					}
					yfluxes[i * (block.ny + 1) + j][stage].gauss[num_gauss].x[var] = Flux;
				}
			}
		}
	}

	Update_with_gauss(fluids, xfluxes, yfluxes, block, stage);
}
void Update_with_gauss(Fluid2d* fluids, Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, Block2d block, int stage)
{
	//Note : calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost + 1; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost + 1; j++)
		{
			int face = i * (block.ny + 1) + j;
			for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
			{
				Local_to_Global(yfluxes[face][stage].gauss[num_gauss].x, yfluxes[face][stage].gauss[num_gauss].normal);
				Local_to_Global(xfluxes[face][stage].gauss[num_gauss].x, xfluxes[face][stage].gauss[num_gauss].normal);
			}
		}
	}

#pragma omp parallel  for
	for (int i = block.ghost; i < block.nodex + block.ghost; i++)
	{
		for (int j = block.ghost; j < block.nodey + block.ghost; j++)
		{
			int cell = i * (block.ny) + j;
			int face = i * (block.ny + 1) + j;

			for (int var = 0; var < 4; var++)
			{
				fluids[cell].convar[var] = fluids[cell].convar_old[var]; //get the Wn from convar_old
				double total_flux = 0.0;
				for (int num_gauss = 0; num_gauss < gausspoint; num_gauss++)
				{
					total_flux += yfluxes[face][stage].gauss[num_gauss].length * yfluxes[face][stage].gauss[num_gauss].x[var];
					total_flux += -yfluxes[face + 1][stage].gauss[num_gauss].length * yfluxes[face + 1][stage].gauss[num_gauss].x[var];
					total_flux += xfluxes[face][stage].gauss[num_gauss].length * xfluxes[face][stage].gauss[num_gauss].x[var];
					total_flux += -xfluxes[face + block.ny + 1][stage].gauss[num_gauss].length * xfluxes[face + block.ny + 1][stage].gauss[num_gauss].x[var];
				}
				fluids[cell].convar[var] += total_flux / fluids[cell].area;
				// calculate the final flux of the cell, in fluids, by the obtained final flux of interface, in xfluxes and yfluxes
				// in 2d, the total flux be updated by the line-averaged flux of four interface
			}
		}
	}
}

void Check_Order_Reduce_by_Lambda_2D(bool& order_reduce, double* convar)
{
	order_reduce = false;
	double lambda;
	lambda = Lambda(convar[0], convar[1] / convar[0], convar[2] / convar[0], convar[3]);
	//if lambda <0, then reduce to the first order
	if (lambda <= 0.0 || (lambda == lambda) == false)
	{
		order_reduce = true;
	}
}
void Recompute_KFVS_1st(Fluid2d& fluid, double* center, double* left, double* right, double* down, double* up, Block2d& block)
{
	double flux_left[4];
	double flux_right[4];
	double flux_down[4];
	double flux_up[4];
	KFVS_1st(flux_left, left, center, block.dt);
	KFVS_1st(flux_right, center, right, block.dt);

	double center_tmp[4], up_tmp[4], down_tmp[4];
	YchangetoX(center_tmp, center);
	YchangetoX(up_tmp, up);
	YchangetoX(down_tmp, down);

	KFVS_1st(flux_down, down_tmp, center_tmp, block.dt);
	KFVS_1st(flux_up, center_tmp, up_tmp, block.dt);

	fluid.convar[0] = fluid.convar_old[0] + 1.0 / block.dy * (flux_down[0] - flux_up[0]);
	fluid.convar[1] = fluid.convar_old[1] + 1.0 / block.dy * (-flux_down[2] + flux_up[2]);
	fluid.convar[2] = fluid.convar_old[2] + 1.0 / block.dy * (flux_down[1] - flux_up[1]);
	fluid.convar[3] = fluid.convar_old[3] + 1.0 / block.dy * (flux_down[3] - flux_up[3]);

	fluid.convar[0] = fluid.convar[0] + 1.0 / block.dx * (flux_left[0] - flux_right[0]);
	fluid.convar[1] = fluid.convar[1] + 1.0 / block.dx * (flux_left[1] - flux_right[1]);
	fluid.convar[2] = fluid.convar[2] + 1.0 / block.dx * (flux_left[2] - flux_right[2]);
	fluid.convar[3] = fluid.convar[3] + 1.0 / block.dx * (flux_left[3] - flux_right[3]);
}
void KFVS_1st(double* flux, double* left, double* right, double dt)
{
	double density_left, density_right;

	double u_left, u_right, v_left, v_right;
	//λ
	double lambda_left, lambda_right;

	density_left = left[0];
	density_right = right[0];

	u_left = U(density_left, left[1]);
	v_left = V(density_left, left[2]);
	lambda_left = Lambda(density_left, u_left, v_left, left[3]);

	u_right = U(density_right, right[1]);
	v_right = V(density_right, right[2]);
	lambda_right = Lambda(density_right, u_right, v_right, right[3]);

	MMDF1st m2(u_left, v_left, lambda_left);
	MMDF1st m3(u_right, v_right, lambda_right);

	//t4part
	flux[0] =
		density_left * dt * m2.uplus[1] +
		density_right * dt * m3.uminus[1];

	flux[1] =
		density_left * dt * m2.uplus[2] +
		density_right * dt * m3.uminus[2];

	flux[2] =
		density_left * dt * m2.uplus[1] * m2.vwhole[1] +
		density_right * dt * m3.uminus[1] * m3.vwhole[1];

	flux[3] =
		density_left * 0.5 * dt * (m2.uplus[3] + m2.uplus[1] * m2.vwhole[2] + m2.uplus[1] * m2.xi2) +
		density_right * 0.5 * dt * (m3.uminus[3] + m3.uminus[1] * m3.vwhole[2] + m3.uminus[1] * m3.xi2);
}





void SetNonUniformMesh(Block2d& block, Fluid2d* fluids,
	Interface2d* xinterfaces, Interface2d* yinterfaces,
	Flux2d_gauss** xfluxes, Flux2d_gauss** yfluxes, double** node2d)
{
	block.xcell_begin = block.ghost;
	block.xcell_end = block.ghost + block.nodex - 1;
	block.ycell_begin = block.ghost;
	block.ycell_end = block.ghost + block.nodey - 1;

	block.xinterface_begin_n = block.ghost;
	block.xinterface_end_n = block.ghost + block.nodex;
	block.xinterface_begin_t = block.ghost;
	block.xinterface_end_t = block.ghost + block.nodex - 1;

	block.yinterface_begin_n = block.ghost;
	block.yinterface_end_n = block.ghost + block.nodey;
	block.yinterface_begin_t = block.ghost;
	block.yinterface_end_t = block.ghost + block.nodey - 1;

#pragma omp parallel  for
	for (int i = 0; i < block.nx; i++)
	{
		for (int j = 0; j < block.ny; j++)
		{
			int cindex = i * block.ny + j;
			int node_index[4];
			node_index[0] = i * (block.ny + 1) + j;
			node_index[1] = (i + 1) * (block.ny + 1) + j;
			node_index[2] = (i + 1) * (block.ny + 1) + j + 1;
			node_index[3] = (i) * (block.ny + 1) + j + 1;
			Set_cell_geo_from_quad_node(fluids[cindex],
				node2d[node_index[0]], node2d[node_index[1]],
				node2d[node_index[2]], node2d[node_index[3]]);
		}
	}

#pragma omp parallel  for
	for (int i = 0; i < block.nx + 1; i++)
	{
		for (int j = 0; j < block.ny + 1; j++)
		{
			int interface_index = i * (block.ny + 1) + j;
			int node_index[3];
			node_index[0] = i * (block.ny + 1) + j;
			node_index[1] = (i + 1) * (block.ny + 1) + j;
			node_index[2] = (i) * (block.ny + 1) + j + 1;
			if (j < block.ny)
			{
				Set_interface_geo_from_two_node
				(xinterfaces[interface_index], node2d[node_index[0]], node2d[node_index[2]], 0);
			}
			if (i < block.nx)
			{
				Set_interface_geo_from_two_node
				(yinterfaces[interface_index], node2d[node_index[0]], node2d[node_index[1]], 1);
			}

			Copy_geo_from_interface_to_line(xinterfaces[interface_index]);
			Copy_geo_from_interface_to_line(yinterfaces[interface_index]);

			xinterfaces[interface_index].gauss = new Recon2d[gausspoint];
			yinterfaces[interface_index].gauss = new Recon2d[gausspoint];
			if (j < block.ny)
			{
				Set_Gauss_Coordinate_general_mesh_x(xinterfaces[interface_index],
					node2d[node_index[0]], node2d[node_index[2]]);
			}
			if (i < block.nx)
			{
				Set_Gauss_Coordinate_general_mesh_y(yinterfaces[interface_index],
					node2d[node_index[0]], node2d[node_index[1]]);
			}


			Copy_geo_from_interface_to_flux
			(xinterfaces[interface_index], xfluxes[interface_index], block.stages);
			Copy_geo_from_interface_to_flux
			(yinterfaces[interface_index], yfluxes[interface_index], block.stages);

		}
	}

}


void Set_cell_geo_from_quad_node
(Fluid2d& fluid, double* n1, double* n2, double* n3, double* n4)
{
	//4            //3
	/////////////////
	//             //
	//             //
	//             //
	//             //
	//             //
	/////////////////
	//1            //2
	Copy_Array(&fluid.node[0], n1, 2);
	Copy_Array(&fluid.node[2], n2, 2);
	Copy_Array(&fluid.node[4], n3, 2);
	Copy_Array(&fluid.node[6], n4, 2);
	fluid.coordx = (n1[0] + n2[0] + n4[0] + n3[0]) / 4.0;
	fluid.coordy = (n1[1] + n2[1] + n4[1] + n3[1]) / 4.0;
	fluid.area = 0.5 * sqrt((pow((n2[0] - n1[0]), 2) + pow((n2[1] - n1[1]), 2)) *
		(pow((n4[0] - n1[0]), 2) + pow((n4[1] - n1[1]), 2)) -
		pow(((n2[0] - n1[0]) * (n4[0] - n1[0]) +
			(n2[1] - n1[1]) * (n4[1] - n1[1])), 2)) +
		0.5 * sqrt((pow((n2[0] - n3[0]), 2) + pow((n2[1] - n3[1]), 2)) *
			(pow((n4[0] - n3[0]), 2) + pow((n4[1] - n3[1]), 2)) -
			pow(((n2[0] - n3[0]) * (n4[0] - n3[0]) +
				(n2[1] - n3[1]) * (n4[1] - n3[1])), 2));
	fluid.dx = sqrt(pow((0.5 * (n2[0] + n3[0]) - 0.5 * (n1[0] + n4[0])), 2)
		+ pow((0.5 * (n3[1] - n4[1]) + 0.5 * (n2[1] - n1[1])), 2));
	fluid.dy = sqrt(pow((0.5 * (n4[0] + n3[0]) - 0.5 * (n1[0] + n2[0])), 2) +
		pow((0.5 * (n4[1] + n3[1]) - 0.5 * (n1[1] + n2[1])), 2));

}


void Set_interface_geo_from_two_node
(Interface2d& interface, double* n1, double* n2, int direction)
{
	interface.x = (n1[0] + n2[0]) / 2.0;
	interface.y = (n1[1] + n2[1]) / 2.0;
	interface.length =
		sqrt(pow((n2[0] - n1[0]), 2) + pow((n2[1] - n1[1]), 2));
	if (direction == 0) //we see it as x direction
	{
		interface.normal[0] =
			(n2[1] - n1[1]) / interface.length;
		interface.normal[1] =
			-(n2[0] - n1[0]) / interface.length;
	}
	else
	{
		interface.normal[0] =
			-(n2[1] - n1[1]) / interface.length;
		interface.normal[1] =
			(n2[0] - n1[0]) / interface.length;
	}

}


void Set_Gauss_Coordinate_general_mesh_x
(Interface2d& xinterface, double* node0, double* nodex1)
{
	double xoff[2];
	xoff[0] = nodex1[0] - node0[0]; xoff[1] = nodex1[1] - node0[1];
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		// first is gauss parameter
		xinterface.gauss[num_guass].x = xinterface.x + gauss_loc[num_guass] * 0.5 * xoff[0];
		xinterface.gauss[num_guass].y = xinterface.y + gauss_loc[num_guass] * 0.5 * xoff[1];
		xinterface.gauss[num_guass].normal[0] = xinterface.normal[0];
		xinterface.gauss[num_guass].normal[1] = xinterface.normal[1];
		// each gauss point contain left center right point 
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].left, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].right, num_guass, xinterface.length);
		Set_Gauss_Intergation_Location_x(xinterface.gauss[num_guass].center, num_guass, xinterface.length);
	}
}

void Set_Gauss_Coordinate_general_mesh_y
(Interface2d& yinterface, double* node0, double* nodey1)
{
	double yoff[2];
	yoff[0] = nodey1[0] - node0[0]; yoff[1] = nodey1[1] - node0[1];
	for (int num_guass = 0; num_guass < gausspoint; num_guass++)
	{
		yinterface.gauss[num_guass].x = yinterface.x + gauss_loc[num_guass] * 0.5 * yoff[0];
		yinterface.gauss[num_guass].y = yinterface.y + gauss_loc[num_guass] * 0.5 * yoff[1];
		yinterface.gauss[num_guass].normal[0] = yinterface.normal[0];
		yinterface.gauss[num_guass].normal[1] = yinterface.normal[1];
		// each gauss point contain left center right point 
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].left, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].right, num_guass, yinterface.length);
		Set_Gauss_Intergation_Location_y(yinterface.gauss[num_guass].center, num_guass, yinterface.length);
	}
}


void Residual2d(Fluid2d* fluids, Block2d block, int outputstep)
{
	if (block.step % outputstep == 0)
	{
		int order = block.ghost;
		double residual[4];
		double sum_old[4];
		for (int k = 0; k < 4; k++)
		{
			residual[k] = 0.0;
			sum_old[k] = 0.0;
		}
		for (int i = order; i < block.nx - order; i++)
		{
			for (int j = order; j < block.ny - order; j++)
			{
				int index = i * block.ny + j;
				for (int k = 0; k < 4; k++)
				{
					residual[k] = residual[k] + abs(fluids[index].convar[k] - fluids[index].convar_old[k]);
					sum_old[k] = sum_old[k] + abs(fluids[index].convar_old[k]);
				}
			}
		}
		for (int k = 0; k < 4; k++)
		{
			residual[k] = residual[k] / (sum_old[k] + 1e-15) / block.dt;
		}

		cout << "step= " << block.step
			<< " density " << residual[0] << " u " << residual[1]
			<< " v " << residual[2] << " densityE " << residual[3] << endl;


		ofstream out;
		if (block.step == outputstep)
		{
			out.open("result/residual-2D.plt", ios::out);
		}
		else
		{
			out.open("result/residual-2D.plt", ios::ate | ios::out | ios::in);
		}

		if (!out.is_open())
		{
			cout << "cannot find residual-2D.plt" << endl;
			cout << "a new case will start" << endl;
		}

		if (block.step == outputstep)
		{
			out << "# CFL number is " << block.CFL << endl;
			out << "# tau type (0 refer euler, 1 refer NS) is " << tau_type << endl;
			out << "# time marching strategy is " << block.stages << " stage method" << endl;
			out << "variables = step,density_residual,m_residual,n_residual,E_residual" << endl;

		}

		//output the data
		out << block.step << " "
			<< log10(residual[0]) << " " << log10(residual[1]) << " "
			<< log10(residual[2]) << " " << log10(residual[3]) << endl;
		out.close();
	}
}
