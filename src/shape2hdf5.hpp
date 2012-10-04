/**
  *  $Id: shape2hdf5.hpp 38 2012-08-09 23:01:20Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  */

#ifndef _SHAPE2HDF5_H_
#define _SHAPE2HDF5_H_

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

extern "C" {
void s2h_converter(double** shape_def, unsigned int num_triangles, char* hdf5_filename, MPI_Comm comm);
}

class shape2hdf5_converter {
	public:
		shape2hdf5_converter(char* filename, char* outfilename, MPI_Comm comm);
		~shape2hdf5_converter() {
			if(filename_ != NULL) delete filename_;
			if(outfilename_ != NULL) delete outfilename_;
			if(shape_def_ != NULL) delete[] shape_def_;
		} // ~converter()

	private:
		int load_shape(char* filename, std::vector<double> &shape_def);
		void convert(char* outfilename, std::vector<double> shape_def);

		std::string *filename_;
		std::string *outfilename_;
		double* shape_def_;
		MPI_Comm comm_;
}; // class shape2hdf5_converter

#endif // _SHAPE2HDF5_H_
