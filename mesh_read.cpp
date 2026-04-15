#include"mesh_read.h"

//兼容linux
string add_mesh_directory_modify_for_linux()
{
#if defined(_WIN32)
	return "";
#else
	return "../";
#endif
}

void Read_structured_mesh(string& filename,
	Fluid2d** fluids, Interface2d** xinterfaces, Interface2d** yinterfaces,
	Flux2d_gauss*** xfluxes, Flux2d_gauss*** yfluxes, Block2d& block)
{
	//参考 因特网无名前辈
	int block_num = 1;

	ifstream in;
	in.open(filename, ifstream::binary);
	if (!in.is_open())
	{
		cout << "cannot find input structured grid file (.plt) " << endl;
		while (1);
		exit(0);
	}
	cout << "reading mesh..." << endl;

	//prepare node information
	double*** node2d = new double** [block_num];


	// Header section
	// section i
	char buffer[9];
	in.read(buffer, 8);
	buffer[8] = '\0';
	string version = string(buffer);
	if (version != "#!TDV112")
		cout << "Parser based on version #!TDV112, current version is: " << version << endl;

	// section ii
	int IntegerValue;
	in.read((char*)&IntegerValue, sizeof(IntegerValue));
	if (IntegerValue != 1)
		cout << "the interger value is not equal to 1 but " << IntegerValue << endl;

	// section iii
	//file type  0=FULL 1=GRID 2=SOLUTION
	int file_type;
	in.read((char*)&file_type, sizeof(file_type));
	//file title ---> defined name of the plt
	string title;
	read_ascii_to_string(title, in);
	cout << "the mesh title is " << title;
	//the number of total variables (x, y, z) as defult
	int number_variables;
	in.read((char*)&number_variables, sizeof(number_variables));
	if (number_variables != 3)
	{
		cout << "warning, the input grid file has the number of variables " << number_variables << endl;
	}
	string variable_names[3];
	for (int i = 0; i < number_variables; i++)
		read_ascii_to_string(variable_names[i], in);

	// section iv
	// we try to read multi-zone; the marker for a start of new zone is a float number 
	// validation_marker==299.0
	float validation_marker;
	in.read((char*)&validation_marker, sizeof(validation_marker));
	assert(validation_marker == 299.0);
	int zone_num = 0;
	int imax[100];
	int jmax[100];
	int kmax[100];
	size_t number_nodes[100];
	while (validation_marker == 299.0)
	{
		zone_num++;
		string zone_name;
		int parent_zone;
		int strand_id;
		double solution_time;
		int not_used;
		int zone_type;  // 0=ORDERED 1=FELINESEG 2=FETRIANGLE 3=FEQUADRILATERAL 
						// 4=FETETRAHEDRON 5=FEBRICK 6=FEPOLYGON 7=FEPOLYHEDRON
		int data_packing; // 0=Block, 1=Point
		int var_location;

		read_ascii_to_string(zone_name, in);
		in.read((char*)&parent_zone, sizeof(parent_zone));
		in.read((char*)&strand_id, sizeof(strand_id));
		in.read((char*)&solution_time, sizeof(solution_time));
		in.read((char*)&not_used, sizeof(not_used));
		in.read((char*)&zone_type, sizeof(zone_type));
		in.read((char*)&data_packing, sizeof(data_packing));
		in.read((char*)&var_location, sizeof(var_location));


		if (var_location == 0)
		{
		}
		else
		{
			cout << "we donnt support this function" << endl;
			exit(0);
		}

		// face neighbors
		int face_neighbors;
		in.read((char*)&face_neighbors, sizeof(face_neighbors));

		if (face_neighbors == 0)
		{
		}
		else
		{
			cout << "we donnt support this function" << endl;
			exit(0);
		}

		// ordered zone
		if (zone_type == 0)
		{
			in.read((char*)&imax[zone_num - 1], sizeof(imax[zone_num - 1]));
			in.read((char*)&jmax[zone_num - 1], sizeof(jmax[zone_num - 1]));
			in.read((char*)&kmax[zone_num - 1], sizeof(kmax[zone_num - 1]));
			number_nodes[zone_num - 1]
				= imax[zone_num - 1] * jmax[zone_num - 1] * kmax[zone_num - 1];
			cout << zone_name << " has total "
				<< imax[zone_num - 1] << " * " << jmax[zone_num - 1] << " * " << kmax[zone_num - 1] << "points " << endl;

		}
		else
		{
			cout << "we donnt support unstructure mesh currently" << endl;
			exit(0);
		}
		// to do create an auxiliary structure
		// add all the structures to a vector
		in.read((char*)&validation_marker, sizeof(validation_marker));


		// read the next marker to check if 
		// there is more than one zone
		in.read((char*)&validation_marker, sizeof(validation_marker));
	} // end of zone section

	int sum_mesh_number = 0;

	for (int iblock = 0; iblock < block_num; iblock++)
	{
		block.nodex = imax[iblock] - 1;
		block.nodey = jmax[iblock] - 1;
		for (int iblock = 0; iblock < block_num; iblock++)
		{
			block.nx = block.nodex + 2 * block.ghost;
			block.ny = block.nodey + 2 * block.ghost;
			sum_mesh_number += block.nodex * block.nodey;
		}
		node2d[iblock] = new double* [(block.nx + 1) * (block.ny + 1)];
		for (int icell = 0; icell < (block.nx + 1) * (block.ny + 1); icell++)
		{
			node2d[iblock][icell] = new double[2];
		}
	}



	// section v
	if (validation_marker == 399)
	{
		// read geometries
		cout << "Geometry section not implemented. Stopping." << endl;
		exit(-1);
	}
	// end of geometrie

	////////////////////////
	// end of header section
	////////////////////////
	cout << "total " << zone_num << " zones are read in" << endl;
	cout << "the total cell for all zones is " << sum_mesh_number << endl;
	// Data section
	assert(validation_marker == 357.0);
	for (size_t index = 0; index < zone_num; index++)
	{
		in.read((char*)&validation_marker, sizeof(validation_marker));
		assert(validation_marker == 299.0);
		int format[3];
		for (int i_var = 0; i_var < number_variables; i_var++)
		{
			in.read((char*)&format[i_var], sizeof(format[i_var]));
			if (format[i_var] != 2)
			{
				cout << "we only support double pricision coordinates now" << endl;
				exit(0);
			}
		}
		int has_passive_variables;
		in.read((char*)&has_passive_variables, sizeof(has_passive_variables));

		if (has_passive_variables)
		{
			cout << "we donnt support has_passive_variables != 0 now" << endl;
			exit(0);
		}
		int has_variable_sharing;
		in.read((char*)&has_variable_sharing, sizeof(has_variable_sharing));

		if (has_variable_sharing)
		{
			cout << "we donnt support has_variable_sharing != 0 now" << endl;
			exit(0);
		}
		int zone_share_connectivity;
		in.read((char*)&zone_share_connectivity, sizeof(zone_share_connectivity));

		double  min_value[3], max_value[3];

		for (int i = 0; i < number_variables; i++)
		{
			in.read((char*)&min_value[i], sizeof(min_value[i]));
			in.read((char*)&max_value[i], sizeof(max_value[i]));
		}

		for (int i_var = 0; i_var < number_variables; i_var++)
		{
			for (int k = 0; k < kmax[index]; k++)
			{
				for (int j = 0; j < jmax[index]; j++)
				{
					for (int i = 0; i < imax[index]; i++)
					{
						double coord;
						in.read((char*)&coord, sizeof(coord));
						int icell = block.ghost + i;
						int jcell = block.ghost + j;
						int inode = icell * (block.ny + 1) + jcell;
						if (i_var < 2)
						{
							node2d[index][inode][i_var] = coord;
						}

					}
				}
			}
		}

	}

	create_mirror_ghost_cell(block, node2d[0]);

	//once we have the dimension infomation, we can use it 
	// to allocate memory
	*fluids = Setfluid(block);
	// then the interfaces reconstructioned value, (N+1)*(N+1)
	*xinterfaces = Setinterface_array(block);
	*yinterfaces = Setinterface_array(block);
	// then the flux, which the number is identical to interfaces
	*xfluxes = Setflux_gauss_array(block);
	*yfluxes = Setflux_gauss_array(block);

	//we shall use the node coordinate to compute the geometric infos for cells/interfaces


	SetNonUniformMesh(block, fluids[0],
		xinterfaces[0], yinterfaces[0],
		xfluxes[0], yfluxes[0], node2d[0]);

}

void read_ascii_to_string(string& value, ifstream& in)
{
	int ascii(-1);

	in.read((char*)&ascii, sizeof(ascii));
	while (ascii != 0)
	{
		char temp = (char)ascii;
		value.append(sizeof(char), temp);
		in.read(reinterpret_cast<char*>(&ascii), sizeof(ascii));
	}
}

void create_mirror_ghost_cell(Block2d& block, double** node2d)
{

	//create the coordinate for ghost nodes
	for (int j = 0; j < block.ny + 1; j++)
	{
		for (int i = 0; i < block.nx + 1; i++)
		{
			int inode = i * (block.ny + 1) + j;
			int base_node = 0, bc_node = 0;
			int base_index_x = i, base_index_y = j;
			int bc_index_x = i, bc_index_y = j;
			if (i < block.ghost)
			{
				base_index_x = 2 * block.ghost - i;
				bc_index_x = block.ghost;
			}
			else if (i > block.nx - block.ghost)
			{
				base_index_x = 2 * (block.nx - block.ghost) - i;
				bc_index_x = block.nx - block.ghost;
			}
			if (j < block.ghost)
			{
				base_index_y = 2 * block.ghost - j;
				bc_index_y = block.ghost;
			}
			else if (j > block.ny - block.ghost)
			{
				base_index_y = 2 * (block.ny - block.ghost) - j;
				bc_index_y = block.ny - block.ghost;
			}
			if (base_index_x != i || base_index_y != j)
			{
				base_node = base_index_x * (block.ny + 1) + base_index_y;
				bc_node = bc_index_x * (block.ny + 1) + bc_index_y;

				//cacluate mirror point
				node2d[inode][0] = 2 * node2d[bc_node][0] - node2d[base_node][0];
				node2d[inode][1] = 2 * node2d[bc_node][1] - node2d[base_node][1];
			}
		}
	}
}
