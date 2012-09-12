/**
 * $Id: object2shape_test.cpp 120 2012-02-28 05:41:05Z asarje $
 */

#include "object2shape.hpp"
#include "object2hdf5.h"

#include <iomanip>

int main(int narg, char** args) {
	if(narg != 3) {
		std::cout << "Please give object filename and output filename." << std::endl;
		return 0;
	} // if

	int rank, num_procs;
	MPI_Init(&narg, &args);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	o2s_converter my_convert(args[1], args[2], MPI_COMM_WORLD, false);

/*	double *shape_def = NULL;
	unsigned int num_triangles = 0;
	h5_shape_reader(args[2], &shape_def, &num_triangles, MPI_COMM_WORLD);

	std::cout << "NUM TRIANGLES = " << num_triangles << std::endl;
	for(int i = 0; i < num_triangles * 7; i += 7) {
		std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6)
			<< shape_def[i] << "\t"
			<< shape_def[i + 1] << "\t"
			<< shape_def[i + 2] << "\t"
			<< shape_def[i + 3] << "\t"
			<< shape_def[i + 4] << "\t"
			<< shape_def[i + 5] << "\t"
			<< shape_def[i + 6] << std::endl;
	} // for
	if(shape_def != NULL) free(shape_def);
*/

	MPI_Finalize();
	return 0;
} // main()
