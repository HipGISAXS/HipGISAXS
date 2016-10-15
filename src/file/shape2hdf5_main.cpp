/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: shape2hdf5_main.cpp
 *  Created: Sep 12, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <file/shape2hdf5.hpp>
#include <file/object2hdf5.h>

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
