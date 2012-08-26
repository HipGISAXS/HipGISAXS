/**
 * $Id: object2hdf5.h 38 2012-08-09 23:01:20Z asarje $
 */

#ifndef _OBJECT2HDF5_H_
#define _OBJECT2HDF5_H_

//#include <mpi.h>
#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

//void s2h_converter(double** shape_def, unsigned int num_triangles, char* hdf5_filename, MPI_Comm comm);
void h5_shape_reader(const char* hdf5_filename, double** shape_def, unsigned int* num_triangles/*, MPI_Comm comm*/);

#ifdef __cplusplus
}
#endif

#endif // _OBJECT2HDF5_H_
