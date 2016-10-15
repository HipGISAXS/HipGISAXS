/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: object2hdf5.c
 *  Created: Aug 25, 2012
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

#include <string.h>

#include <file/hdf5shape_reader.h>

/**
 * reads a hdf5 file and converts to shape def format
 */
void h5_shape_reader(const char* hdf5_filename, double** shape_def, unsigned int* num_triangles/*, MPI_Comm comm*/) {
	hid_t file_id, dataset, dataspace, datatype, class, order;

	file_id = H5Fopen(hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen(file_id, "shape_def", H5P_DEFAULT);
	dataspace = H5Dget_space(dataset);
	datatype = H5Dget_type(dataset);
	class = H5Tget_class(datatype);
	order = H5Tget_order(datatype);
	int size = H5Tget_size(datatype);

	// check for datatype size ...

	int num_dim = H5Sget_simple_extent_ndims(dataspace);
	hsize_t* dims = (hsize_t*) malloc(sizeof(hsize_t) * num_dim);
	hsize_t* max_dims = (hsize_t*) malloc(sizeof(hsize_t) * num_dim);
	herr_t status = H5Sget_simple_extent_dims(dataspace, dims, max_dims);

	//printf("DATATYPE = %d, CLASS = %d, SIZE = %d\n", datatype, class, size);
	//printf("Dimensions = %d: %d x %d \n", num_dim, dims[0], dims[1]);

	(*num_triangles) = dims[0];

	(*shape_def) = (double *) malloc(dims[0] * dims[1] * size);
	if((*shape_def) == NULL) {
		fprintf(stderr, "error: cannot allocate memory to read shape definition data\n");
		*num_triangles = 0;
	} else {
		memset(&(**shape_def), 0, dims[0] * dims[1] * size);
		status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, *shape_def);
		if(status < 0) {
			fprintf(stderr, "error: shape definition data reading failed\n");
			*num_triangles = 0;
		} // if
	} // if-else

	free(max_dims);
	free(dims);
	H5Sclose(dataspace);
	H5Tclose(datatype);
	H5Dclose(dataset);
	H5Fclose(file_id);
} // h5_shape_reader()
