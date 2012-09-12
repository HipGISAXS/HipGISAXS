/**
 * $Id: shape2hdf5_test.cpp 120 2012-02-28 05:41:05Z asarje $
 */

#include "shape2hdf5.hpp"
#include "object2hdf5.h"

#include <iomanip>

int main(int narg, char** args) {
	if(narg != 3) {
		std::cout << "Please give shape filename and output hdf5 filename." << std::endl;
		return 0;
	} // if

	int rank, num_procs;
	MPI_Init(&narg, &args);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	shape2hdf5_converter my_convert(args[1], args[2], MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
} // main()
