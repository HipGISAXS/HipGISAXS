/**
 * $Id: ff.cu 40 2012-08-21 20:43:52Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 */

//#ifndef _FF_CU_
//#define _FF_CU_

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cuComplex.h>
#include <omp.h>

#include "typedefs.hpp"
#include "constants.hpp"
#include "parameters.hpp"
#include "reduction.cuh"

namespace hig {

	/**
	 * Templatized class for computing Form Factor in either single or double precision on a single node/GPU.
	 */
	template<typename float_t, typename complex_t>
	class GFormFactor {
		public:
			GFormFactor(int block_cuda):
				block_cuda_(block_cuda), block_cuda_t_(0), block_cuda_y_(0), block_cuda_z_(0) { }
			GFormFactor(int block_cuda_t, int block_cuda_y, int block_cuda_z): block_cuda_(0),
				block_cuda_t_(block_cuda_t), block_cuda_y_(block_cuda_y), block_cuda_z_(block_cuda_z) { }
			~GFormFactor() { }

			/* original */
			unsigned int compute_form_factor(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					complex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					complex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
	#ifdef FINDBLOCK
					, const int, const int, const int, const int
	#endif
					);

			/* with double buffering */
			unsigned int compute_form_factor_db(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					complex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					complex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
	#ifdef FINDBLOCK
					, const int, const int, const int, const int
	#endif
					);

			/* with double buffering and optimized memory (incomplete ... )*/
			unsigned int compute_form_factor_db_mem(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					complex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					complex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
	#ifdef FINDBLOCK
					, const int, const int, const int, const int
	#endif
					);

		private:
			// for kernel 1
			int block_cuda_;

			// for kernel 2 and up
			int block_cuda_t_;
			int block_cuda_y_;
			int block_cuda_z_;

//			unsigned int read_shape_surface_file(char* filename, std::vector<float_t>& shape_def);
			void compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
					unsigned long int estimated_device_mem_need, unsigned long int device_mem_avail,
					unsigned int& b_nqx, unsigned int& b_nqy, unsigned int& b_nqz,
					unsigned int& b_num_triangles
	#ifdef FINDBLOCK
					, const int, const int, const int, const int
	#endif
					);
			void move_to_main_ff(complex_t* fq_buffer,
					unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
					unsigned int nqx, unsigned int nqy, unsigned int nqz,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z,
					complex_t* ff);
	}; // class GFormFactor

	/**
	 * Some forward declarations.
	 */
	/* K1: default kernel, with t decompostion, no shared memory */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*);
	/* K2: kernel with t, y, z decomposition, no shared memory */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*);
	/* K3: kernel with t, y, z decomposition, static shared mem for input */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*);
	/* K4: kernel with t, y, z decomposition, dynamic shared mem for input, static for output */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2(const float_t*, const float_t*, const complex_t*,
							const float_t*, const short int*,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							complex_t*);
	/* K5: kernel with t, y, z decomposition, dynamic shared mem for input, static for output, memopt? ... */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_mem(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*);
	/* K6: kernel with K3 and blocked y, z (probably incomplete ...) */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared_subblock(float_t*, float_t*, complex_t*, float_t*,
										short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*);
	/* K7: kernel with K1 (?) and blocked y, z (incomplete ...) */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_2(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*);
	/* K8: K4 with reduction combined in it */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_red(const float_t*, const float_t*, const complex_t*,
							const float_t*, const short int*,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							complex_t*,	complex_t*);
	__device__ cuFloatComplex compute_fq(float s, cuFloatComplex qt_d, cuFloatComplex qn_d);
	__device__ cuDoubleComplex compute_fq(double s, cuDoubleComplex qt_d, cuDoubleComplex qn_d);
	__device__ void compute_z(cuFloatComplex temp_z, float nz, float z,
				cuFloatComplex& qz2, cuFloatComplex& qzn, cuFloatComplex& qzt);
	__device__ void compute_z(cuDoubleComplex temp_z, double nz, double z,
				cuDoubleComplex& qz2, cuDoubleComplex& qzn, cuDoubleComplex& qzt);
	__device__ void compute_x(float temp_x, float qy2, cuFloatComplex qz2, float nx, float qyn,
				cuFloatComplex qzn, float x, float qyt,	cuFloatComplex qzt,
				cuFloatComplex& qn_d, cuFloatComplex& qt_d);
	__device__ void compute_x(double temp_x, double qy2, cuDoubleComplex qz2, double nx, double qyn,
				cuDoubleComplex qzn, double x, double qyt, cuDoubleComplex qzt,
				cuDoubleComplex& qn_d, cuDoubleComplex& qt_d);

	/*template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new(float_t*, float_t*, float_t*, float_t*, short int*,
					unsigned int, unsigned int, unsigned int, unsigned int,
					unsigned int, unsigned int, unsigned int, unsigned int,
					unsigned int, unsigned int, unsigned int, unsigned int,
					complex_t*);
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_2(float_t*, float_t*, float_t*, float_t*,
					unsigned int, unsigned int, unsigned int, unsigned int,
					unsigned int, unsigned int, unsigned int, unsigned int,
					unsigned int, unsigned int, unsigned int, unsigned int,
					const int, const int, complex_t*);
	*/

	/**
	 * Wrapper for the GFormFactor class function.
	 * This is required because with templates, everything is required to be in the
	 * same file, and with cuda, the whole .cu file cannot be included in other files!
	 */
	template<typename float_t, typename complex_t>
	unsigned int compute_form_factor_wrap(int rank, std::vector<float_t> &shape_def,
							unsigned int num_triangles,
							std::vector<short int> &axes,
							float_t* &qx_h, int nqx,
							float_t* &qy_h, int nqy,
							complex_t* &qz_h, int nqz,
							complex_t* &ff,
							float_t& kernel_time, float_t& red_time, float_t& mem_time
	#ifdef FINDBLOCK
							, const int block_x, const int block_y
							, const int block_z, const int block_t
	#endif
	#ifdef KERNEL2
							, unsigned int cuda_t, unsigned int cuda_y, unsigned int cuda_z
	#else // default kernel
							, unsigned int cuda_t
	#endif
							) {
	#ifdef KERNEL2
		GFormFactor<float_t, complex_t> gff(cuda_t, cuda_y, cuda_z);
	#else // default kernel
		GFormFactor<float_t, complex_t> gff(cuda_t);
	#endif
		//return gff.compute_form_factor(rank, shape_def, axes, ff, qx_h, nqx, qy_h, nqy, qz_h, nqz,
		//			kernel_time, red_time, mem_time
		return gff.compute_form_factor_db(rank, shape_def, axes, ff, qx_h, nqx, qy_h, nqy, qz_h, nqz,
					kernel_time, red_time, mem_time
	#ifdef FINDBLOCK
					, block_x, block_y, block_z, block_t
	#endif
					);
	} // compute_form_factor_wrap()

	// Instantiations for float and double
	template unsigned int compute_form_factor_wrap<float, cuFloatComplex>(
							int, std::vector<float>&, unsigned int, std::vector<short int>&,
							float*&, int,
							float*&, int,
							cuFloatComplex*&, int,
							cuFloatComplex*&,
							float& kernel_time, float& red_time, float& mem_time
	#ifdef FINDBLOCK
							, const int, const int, const int, const int
	#endif
	#ifdef KERNEL2
							, unsigned int cuda_t, unsigned int cuda_y, unsigned int cuda_z
	#else // default kernel
							, unsigned int cuda_t
	#endif
							);
	template unsigned int compute_form_factor_wrap<double, cuDoubleComplex>(
							int, std::vector<double>&, unsigned int, std::vector<short int>&,
							double*&, int,
							double*&, int,
							cuDoubleComplex*&, int,
							cuDoubleComplex*&,
							double& kernel_time, double& red_time, double& mem_time
	#ifdef FINDBLOCK
							, const int, const int, const int, const int
	#endif
	#ifdef KERNEL2
							, unsigned int cuda_t, unsigned int cuda_y, unsigned int cuda_z
	#else // default kernel
							, unsigned int cuda_t
	#endif
							);
/*	template unsigned int compute_form_factor_wrap<double, std::complex<double> >(
							int, std::vector<double>&, unsigned int, std::vector<short int>&,
							double*&, int,
							double*&, int,
							std::complex<double>*&, int,
							std::complex<double>*&,
							double& kernel_time, double& red_time, double& mem_time
	#ifdef FINDBLOCK
							, const int, const int, const int, const int
	#endif
	#ifdef KERNEL2
							, unsigned int cuda_t, unsigned int cuda_y, unsigned int cuda_z
	#else // default kernel
							, unsigned int cuda_t
	#endif
							); */

	/**
	 * The main host function called from outside, as part of the API for a single node.
	 */
	template<typename float_t, typename complex_t>
	unsigned int GFormFactor<float_t, complex_t>::compute_form_factor(int rank,
						std::vector<float_t> &shape_def, std::vector<short int>& axes,
						complex_t* &ff,
						float_t* &qx_h, int nqx,
						float_t* &qy_h, int nqy,
						complex_t* &qz_h, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
	#ifdef FINDBLOCK
						, const int block_x, const int block_y, const int block_z, const int block_t
	#endif
						) {
		float kernel_time, reduce_time, total_kernel_time = 0.0, total_reduce_time = 0.0,
				temp_mem_time = 0.0, total_mem_time = 0.0;

		//cudaSetDevice(2);		////////////////////////////////////////////////////////////////////////////
		//if(rank == 1) cudaSetDevice(2);

		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
		//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
		//cudaDeviceSetCacheConfig(cudaFuncCachePreferEqual);

		//cudaFuncCache cache_conf;
		//cudaDeviceGetCacheConfig(&cache_conf);
		//if(rank == 0) std::cout << "Device cache config: " << cache_conf << std::endl;

		cudaEvent_t begin_event, end_event;
		cudaEvent_t total_begin_event, total_end_event;
		cudaEvent_t mem_begin_event, mem_end_event;
		cudaEventCreate(&total_begin_event); cudaEventCreate(&total_end_event);
		cudaEventCreate(&mem_begin_event); cudaEventCreate(&mem_end_event);

		cudaEventRecord(total_begin_event, 0);

		cudaEventRecord(mem_begin_event, 0);
		unsigned long int total_qpoints = nqx * nqy * nqz;
		unsigned long int host_mem_usage = ((unsigned long int) nqx + nqy + nqz) * sizeof(float_t);
		size_t device_mem_avail, device_mem_total, device_mem_used;

		// allocate memory for the final FF 3D matrix
		ff = new (std::nothrow) complex_t[total_qpoints]();	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(complex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(complex_t);

		// read the shape file
#ifndef KERNEL2
		unsigned int num_triangles = shape_def.size() / 7;
#else // KERNEL2
		unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
#endif // KERNEL2
		if(num_triangles < 1) return 0;

		cudaError_t err;
		float_t *qx_d, *qy_d;
		complex_t *qz_d;
		if(cudaMalloc((void **) &qx_d, nqx * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qx_d. "
						<< "Size = " << nqx * sizeof(float_t) << " B" << std::endl;
			delete[] ff;
			return 0;
		} // if
		if(cudaMalloc((void **) &qy_d, nqy * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qy_d. "
						<< "Size = " << nqy * sizeof(float_t) << " B" << std::endl;
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		if(cudaMalloc((void **) &qz_d, nqz * sizeof(complex_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qz_d. "
						<< "Size = " << nqz * sizeof(complex_t) << " B" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpy(qx_d, qx_h, nqx * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qz_d, qz_h, nqz * sizeof(complex_t), cudaMemcpyHostToDevice);

		float_t *shape_def_d;
		float_t *shape_def_h = &shape_def[0];
#ifndef KERNEL2
		err = cudaMalloc((void **) &shape_def_d, 7 * num_triangles * sizeof(float_t));
		if(err != cudaSuccess) {
			std::cerr << "Device memory allocation failed for shape_def_d. "
						<< "Size = " << 7 * num_triangles * sizeof(float_t) << " B" << std::endl;
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpy(shape_def_d, shape_def_h, 7 * num_triangles * sizeof(float_t), cudaMemcpyHostToDevice);
#else // KERNEL2
		err = cudaMalloc((void **) &shape_def_d, T_PROP_SIZE_ * num_triangles * sizeof(float_t));
		if(err != cudaSuccess) {
			std::cerr << "Device memory allocation failed for shape_def_d. "
						<< "Size = " << T_PROP_SIZE_ * num_triangles * sizeof(float_t) << " B" << std::endl;
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpy(shape_def_d, shape_def_h, T_PROP_SIZE_ * num_triangles * sizeof(float_t),
					cudaMemcpyHostToDevice);
#endif // KERNEL2

		short int *axes_h = &axes[0];
		short int *axes_d;
		err = cudaMalloc((void **) &axes_d, 3 * sizeof(short int));
		cudaMemcpy(axes_d, axes_h, 3 * sizeof(short int), cudaMemcpyHostToDevice);
		
		unsigned long int matrix_size = (unsigned long int) nqx * nqy * nqz * num_triangles;
		
		cudaMemGetInfo(&device_mem_avail, &device_mem_total);
		device_mem_used = device_mem_total - device_mem_avail;
		size_t estimated_device_mem_need = matrix_size * sizeof(complex_t) + device_mem_used;
		if(rank == 0) {
			std::cout << "++       Available device memory: "
						<< (float) device_mem_avail / 1024 / 1024
						<< " MB" << std::endl;
			std::cout << "++  Estimated device memory need: "
						<< (float) estimated_device_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		//cudaEventRecord(mem_end_event, 0);
		//cudaEventSynchronize(mem_end_event);
		//cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
		//total_mem_time += temp_mem_time;

		// do hyper blocking if too large to fit in device memory
		unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
		compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
							estimated_device_mem_need, device_mem_avail,
							b_nqx, b_nqy, b_nqz, b_num_triangles
	#ifdef FINDBLOCK
							, block_x, block_y, block_z, block_t
	#endif
							);

		//cudaEventRecord(mem_begin_event, 0);
		unsigned long int blocked_3d_matrix_size = (unsigned long int) b_nqx * b_nqy * b_nqz;
		unsigned long int blocked_matrix_size = (unsigned long int) blocked_3d_matrix_size * b_num_triangles;
		
		complex_t *fq_buffer = NULL, *ff_buffer = NULL, *fq_d = NULL, *ff_d = NULL;
	#ifdef GPUR
		size_t estimated_host_mem_need = host_mem_usage + blocked_3d_matrix_size * sizeof(complex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += blocked_3d_matrix_size * sizeof(complex_t);
		if(cudaMallocHost(&ff_buffer, blocked_3d_matrix_size * sizeof(complex_t)) != cudaSuccess) {
			std::cerr << "Memory allocation failed for ff_buffer. blocked_3d_matrix_size = "
						<< blocked_3d_matrix_size << std::endl
						<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
	#else
		size_t estimated_host_mem_need = host_mem_usage + blocked_matrix_size * sizeof(complex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += blocked_matrix_size * sizeof(complex_t);
		fq_buffer = new (std::nothrow) complex_t[blocked_matrix_size]();
		if(fq_buffer == NULL) {
			std::cerr << "Memory allocation failed for fq_buffer. blocked_matrix_size = "
						<< blocked_matrix_size << std::endl
						<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
	#endif // GPUR
		} else {
			if(cudaMalloc((void **) &fq_d, blocked_matrix_size * sizeof(complex_t)) != cudaSuccess) {
				std::cerr << "Device memory allocation failed for fq_d. "
							<< "Size = " << blocked_matrix_size * sizeof(complex_t) << " B" << std::endl;
			} else {
				if(cudaMalloc((void **) &ff_d, blocked_3d_matrix_size * sizeof(complex_t)) != cudaSuccess) {
					std::cerr << "Device memory allocation failed for ff_d. "
								<< "Size = " << blocked_3d_matrix_size * sizeof(complex_t) << " B"
								<< std::endl;
				} else {
					cudaMemGetInfo(&device_mem_avail, &device_mem_total);
					device_mem_used = device_mem_total - device_mem_avail;
					if(rank == 0) {
						std::cout << "++            Device memory used: "
									<< (float) device_mem_used / 1024 / 1024
									<< " MB" << std::endl;
						std::cout << "++              Host memory used: "
									<< (float) host_mem_usage / 1024 / 1024
									<< " MB" << std::endl << std::flush;
					} // if

					// compute the number of sub-blocks, along each of the 4 dimensions
					// formulate loops over each dimension, to go over each sub block
					unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
					unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
					unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
					unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);

					unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
					unsigned int curr_b_num_triangles = b_num_triangles;
					unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;

					if(rank == 0) {
						std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
									<< " x " << b_nqz << " x " << b_num_triangles << std::endl;
						std::cout << "++         Number of hyperblocks: " << num_blocks
									<< " [" << nb_x << ", " << nb_y << ", " << nb_z << ", " << nb_t << "]"
									<< std::endl;

						std::cout << "++        Kernel CUDA block size: ";
	#ifndef KERNEL2
						std::cout << block_cuda_ << std::endl;
	#else
						std::cout << block_cuda_t_ << " x " << block_cuda_y_ << " x "
									<< block_cuda_z_ << std::endl;
	#endif
						std::cout << "++     Reduction CUDA block size: "
									<< BLOCK_REDUCTION_X_ << " x " << BLOCK_REDUCTION_Y_ << " x "
									<< BLOCK_REDUCTION_Z_ << std::endl;
					} // if

					cudaEventRecord(mem_end_event, 0);
					cudaEventSynchronize(mem_end_event);
					cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
					total_mem_time += temp_mem_time;

	#ifdef VERBOSE
					unsigned int block_num = 0;
	#else
					if(rank == 0) std::cout << "-- Computing form factor ... " << std::flush;
	#endif
					cudaEventCreate(&begin_event);
					cudaEventCreate(&end_event);

					// compute for each block
										// TODO: parallelize all these loops - block distrubution ...
					curr_b_nqx = b_nqx;
					for(unsigned int ib_x = 0; ib_x < nb_x; ++ ib_x) {
						if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;	// in last iteration
						curr_b_nqy = b_nqy;
						for(unsigned int ib_y = 0; ib_y < nb_y; ++ ib_y) {
							if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;	// in last iteration
							curr_b_nqz = b_nqz;
							for(unsigned int ib_z = 0; ib_z < nb_z; ++ ib_z) {
								if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;	// in last iteration
								curr_b_num_triangles = b_num_triangles;
								for(unsigned int ib_t = 0; ib_t < nb_t; ++ ib_t) {
									if(ib_t == nb_t - 1)	// in last iteration
										curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

									cudaEventRecord(begin_event, 0);
	#ifdef VERBOSE
									if(rank == 0) {
										std::cout << "- Executing kernel for block " << ++ block_num
													<< "/" << num_blocks << " ... " << std::flush;
									} // if
	#endif
									// call the main kernel on the device

	#ifndef KERNEL2
									// Kernel 1: decompose along triangles
									unsigned int cuda_block_size = block_cuda_;
									unsigned int cuda_num_blocks =
										(unsigned int) ceil((float) curr_b_num_triangles / cuda_block_size);

									form_factor_kernel<float_t>
										<<< cuda_num_blocks, cuda_block_size >>> (
											qx_d, qy_d, qz_d, shape_def_d, axes_d,
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ib_x, ib_y, ib_z, ib_t,
											fq_d);
	#else // KERNEL2
									// Kernel 4
									unsigned int ff_t_blocks = (unsigned int)
														ceil((float) curr_b_num_triangles / block_cuda_t_);
									unsigned int ff_y_blocks = (unsigned int)
														ceil((float) curr_b_nqy / block_cuda_y_);
									unsigned int ff_z_blocks = (unsigned int)
														ceil((float) curr_b_nqz / block_cuda_z_);
									dim3 ff_grid_size(ff_t_blocks, ff_y_blocks, ff_z_blocks);
									dim3 ff_block_size(block_cuda_t_, block_cuda_y_, block_cuda_z_);

									// dynamic shared memory
									size_t dyn_shared_mem_size = sizeof(float_t) *
																	(T_PROP_SIZE_ * block_cuda_t_ +
																	curr_b_nqx + block_cuda_y_) +
																	sizeof(complex_t) * block_cuda_z_;
									size_t stat_shared_mem_size = (FQ_COPY_SIZE_ + BANK_OFF_) *
																	MAX_NUM_THREADS_;
									size_t tot_shared_mem_size = dyn_shared_mem_size + stat_shared_mem_size;
									if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
										std::cout << "++       Required shared memory: "
													<< (float) tot_shared_mem_size << " B"
													<< std::endl;
									} // if
									if(tot_shared_mem_size > 32768) {
										std::cerr << "error: too much shared memory requested!" << std::endl;
										exit(1);
									} // if

									form_factor_kernel_new_shared2<float_t>
										<<< ff_grid_size, ff_block_size, dyn_shared_mem_size >>> (
											qx_d, qy_d, qz_d, shape_def_d, axes_d,
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ib_x, ib_y, ib_z, ib_t,
											fq_d);
	#endif // KERNEL2
									cudaThreadSynchronize();
									cudaEventRecord(end_event, 0);
									cudaEventSynchronize(end_event);
									cudaEventElapsedTime(&kernel_time, begin_event, end_event);
									total_kernel_time += kernel_time;
	#ifdef VERBOSE
									if(rank == 0) {
										std::cout << "done in " << kernel_time << "ms."
													<< std::endl << std::flush;
									} // if
	#endif

									err = cudaGetLastError();
									if(err != cudaSuccess) {
										std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
													<< cudaGetErrorString(err) << std::endl;
									} else {
										// call the reduction kernel
	#ifdef VERBOSE
										if(rank == 0) {
											std::cout << "- Reducing block, " << std::flush;
										} // if
	#endif
	#ifndef GPUR
										cudaEventRecord(mem_begin_event, 0);
										cudaMemcpy(fq_buffer, fq_d, blocked_matrix_size * sizeof(complex_t),
													cudaMemcpyDeviceToHost);
										cudaEventRecord(mem_end_event, 0);
										cudaEventSynchronize(mem_end_event);
										cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
										total_mem_time += temp_mem_time;

										cudaEventRecord(begin_event, 0);

										reduction_kernel(curr_b_nqx, curr_b_nqy, curr_b_nqz,
												curr_b_num_triangles, blocked_matrix_size,
												b_nqx, b_nqy, b_nqz, num_triangles,
												nqx, nqy, nqz,
												ib_x, ib_y, ib_z, ib_t,
												fq_buffer, ff);

										cudaEventRecord(end_event, 0);
										cudaEventSynchronize(end_event);
										cudaEventElapsedTime(&reduce_time, begin_event, end_event);
	#else // GPUR
										cudaEventRecord(begin_event, 0);

										int block_r = BLOCK_REDUCTION_Z_;
										dim3 r_block_size(block_r, block_r, block_r);
										unsigned int rx_num_blocks = (unsigned int) ceil((float) curr_b_nqx /
																									block_r);
										unsigned int ry_num_blocks = (unsigned int) ceil((float) curr_b_nqy /
																									block_r);
										unsigned int rz_num_blocks = (unsigned int) ceil((float) curr_b_nqz /
																									block_r);
	#ifdef VERBOSE
										if(rank == 0) {
											std::cout << "[" <<  rx_num_blocks << " " << ry_num_blocks
														<< " " << rz_num_blocks << "] ... " << std::flush;
										} // if
	#endif
										dim3 r_grid_size(rx_num_blocks, ry_num_blocks, rz_num_blocks);

										reduction_kernel <<< r_grid_size, r_block_size >>> (fq_d,
												curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
												b_nqx, b_nqy, b_nqz, b_num_triangles,
												ff_d);

										cudaThreadSynchronize();
										cudaEventRecord(end_event, 0);
										cudaEventSynchronize(end_event);
										cudaEventElapsedTime(&reduce_time, begin_event, end_event);
										total_reduce_time += reduce_time;

										cudaEventRecord(mem_begin_event, 0);
										cudaMemcpy(ff_buffer, ff_d, b_nqx * b_nqy * b_nqz * sizeof(complex_t),
													cudaMemcpyDeviceToHost);

										// move ff_buffer to correct location in ff
										move_to_main_ff(ff_buffer,
												curr_b_nqx, curr_b_nqy, curr_b_nqz,
												b_nqx, b_nqy, b_nqz,
												nqx, nqy, nqz,
												ib_x, ib_y, ib_z, ff);
										cudaEventRecord(mem_end_event, 0);
										cudaEventSynchronize(mem_end_event);
										cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
										total_mem_time += temp_mem_time;

										//std::cout << "ff_buffer:" << std::endl;
										//for(int i = 0; i < curr_b_nqx * curr_b_nqy * curr_b_nqz; ++ i)
										//	std::cout << ff_buffer[i].x << "," << ff_buffer[i].y << "\t";
										//std::cout << std::endl;
	#endif // GPUR

	#ifdef VERBOSE
										if(rank == 0) {
											std::cout << "done in " << reduce_time << "ms." << std::endl;
										} // if
	#endif
									} // if-else
								} // for ib_t
							} // for ib_z
						} // for ib_y
					} // for ib_x
	#ifndef VERBOSE
					if(rank == 0) std::cout << "done." << std::endl;
	#endif
					cudaEventRecord(mem_begin_event, 0);
					cudaFree(ff_d);
				} // if-else
				cudaFree(fq_d);
			} // if-else
			if(fq_buffer != NULL) delete[] fq_buffer;
			if(ff_buffer != NULL) cudaFreeHost(ff_buffer);
		} // if-else

		//std::cout << "FF:" << std::endl;
		//for(int i = 0; i < nqx * nqy * nqz; ++ i) {
		//	std::cout << ff[i].x << "," << ff[i].y << "\t";
		//} // for
		//std::cout << std::endl;

		cudaFree(axes_d);
		cudaFree(shape_def_d);
		cudaFree(qz_d);
		cudaFree(qy_d);
		cudaFree(qx_d);
		cudaEventRecord(mem_end_event, 0);
		cudaEventSynchronize(mem_end_event);
		cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
		total_mem_time += temp_mem_time;

		float total_time;
		cudaEventRecord(total_end_event, 0);
		cudaEventSynchronize(total_end_event);
		cudaEventElapsedTime(&total_time, total_begin_event, total_end_event);

		pass_kernel_time = total_kernel_time;
		red_time = total_reduce_time;
		mem_time = total_mem_time;

		return num_triangles;
	} // GFormFactor::compute_form_factor()


	/**
	 * The main host function called from outside, This one uses double buffering
	 */
	template<typename float_t, typename complex_t>
	unsigned int GFormFactor<float_t, complex_t>::compute_form_factor_db(int rank,
						std::vector<float_t> &shape_def, std::vector<short int> &axes,
						complex_t* &ff,
						float_t* &qx_h, int nqx,
						float_t* &qy_h, int nqy,
						complex_t* &qz_h, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
	#ifdef FINDBLOCK
						, const int block_x, const int block_y, const int block_z, const int block_t
	#endif
						) {
		cudaSetDevice(2);		////////////////////////////////////////////////////////////////////////////
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	
		float kernel_time = 0.0, reduce_time = 0.0;
		float total_kernel_time = 0.0, total_reduce_time = 0.0;
		float temp_mem_time = 0.0, total_mem_time = 0.0;
		cudaEvent_t start, stop;
		cudaEvent_t total_begin_event, total_end_event;
		cudaEvent_t mem_begin_event, mem_end_event;
		cudaEventCreate(&total_begin_event); cudaEventCreate(&total_end_event);
		cudaEventCreate(&mem_begin_event); cudaEventCreate(&mem_end_event);
	
		cudaEventRecord(total_begin_event, 0);
	
		cudaEventRecord(mem_begin_event, 0);
		unsigned long int total_qpoints = nqx * nqy * nqz;
		unsigned long int host_mem_usage = ((unsigned long int) nqx + nqy + nqz) * sizeof(float_t);
		size_t device_mem_avail, device_mem_total, device_mem_used;
	
		// allocate memory for the final FF 3D matrix
		ff = new (std::nothrow) complex_t[total_qpoints]();	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(complex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(complex_t);
	
		// read the shape file
	#ifndef KERNEL2
		unsigned int num_triangles = shape_def.size() / 7;
	#else // KERNEL2
		unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
	#endif // KERNEL2
		if(num_triangles < 1) return 0;
	
		cudaError_t err;
		float_t *qx_d, *qy_d;
		complex_t *qz_d;
		if(cudaMalloc((void **) &qx_d, nqx * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qx_d. "
						<< "Size = " << nqx * sizeof(float_t) << " B" << std::endl;
			delete[] ff;
			return 0;
		} // if
		if(cudaMalloc((void **) &qy_d, nqy * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qy_d. "
						<< "Size = " << nqy * sizeof(float_t) << " B" << std::endl;
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		if(cudaMalloc((void **) &qz_d, nqz * sizeof(complex_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qz_d. "
						<< "Size = " << nqz * sizeof(float_t) << " B" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpyAsync(qx_d, qx_h, nqx * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qz_d, qz_h, nqz * sizeof(complex_t), cudaMemcpyHostToDevice);
	
		float_t *shape_def_d;
		float_t *shape_def_h = &shape_def[0];
	#ifndef KERNEL2
		err = cudaMalloc((void **) &shape_def_d, 7 * num_triangles * sizeof(float_t));
		if(err != cudaSuccess) {
			std::cerr << "Device memory allocation failed for shape_def_d. "
						<< "Size = " << 7 * num_triangles * sizeof(float_t) << " B"
						<< std::endl;
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpyAsync(shape_def_d, shape_def_h, 7 * num_triangles * sizeof(float_t),
						cudaMemcpyHostToDevice);
	#else // KERNEL2
		err = cudaMalloc((void **) &shape_def_d, T_PROP_SIZE_ * num_triangles * sizeof(float_t));
		if(err != cudaSuccess) {
			std::cerr << "Device memory allocation failed for shape_def_d. "
						<< "Size = " << T_PROP_SIZE_ * num_triangles * sizeof(float_t) << " B"
						<< std::endl;
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpyAsync(shape_def_d, shape_def_h, T_PROP_SIZE_ * num_triangles * sizeof(float_t),
						cudaMemcpyHostToDevice);
	#endif // KERNEL2

		short int *axes_h = &axes[0];
		short int *axes_d;
		err = cudaMalloc((void**) &axes_d, 3 * sizeof(short int));
		cudaMemcpy(axes_d, axes_h, 3 * sizeof(short int), cudaMemcpyHostToDevice);
		
		unsigned long int matrix_size = (unsigned long int) nqx * nqy * nqz * num_triangles;
		
		cudaMemGetInfo(&device_mem_avail, &device_mem_total);
		device_mem_used = device_mem_total - device_mem_avail;
		size_t estimated_device_mem_need = matrix_size * sizeof(complex_t) + device_mem_used;
		if(rank == 0) {
			std::cout << "++       Available device memory: "
						<< (float) device_mem_avail / 1024 / 1024
						<< " MB" << std::endl;
			std::cout << "++  Estimated device memory need: "
						<< (float) estimated_device_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
	
		// do hyperblocking
		unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
		// this computes a hyperblock size
		compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
								estimated_device_mem_need, device_mem_avail,
								b_nqx, b_nqy, b_nqz, b_num_triangles
	#ifdef FINDBLOCK
								, block_x, block_y, block_z, block_t
	#endif
								);
	
		unsigned long int blocked_3d_matrix_size = (unsigned long int) b_nqx * b_nqy * b_nqz;
		unsigned long int blocked_matrix_size =
							(unsigned long int) blocked_3d_matrix_size * b_num_triangles;
		
		complex_t *ff_buffer = NULL, *fq_d = NULL, *ff_d = NULL;
		complex_t *ff_double_buff[2], *fq_double_d[2], *ff_double_d[2];
	
		size_t estimated_host_mem_need = host_mem_usage + 2 * blocked_3d_matrix_size * sizeof(complex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: "
						<< (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += 2 * blocked_3d_matrix_size * sizeof(complex_t);
		if(cudaMallocHost((void **) &ff_buffer, 2 * blocked_3d_matrix_size * sizeof(complex_t))
				!= cudaSuccess) {
			std::cerr << "Memory allocation failed for ff_buffer. blocked_3d_matrix_size = "
						<< blocked_3d_matrix_size << std::endl
						<< "Host memory usage = "
						<< (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
		} else {
			ff_double_buff[0] = ff_buffer;
			ff_double_buff[1] = ff_buffer + blocked_3d_matrix_size;
	
			if(cudaMalloc((void **) &fq_d, 2 * blocked_matrix_size * sizeof(complex_t)) != cudaSuccess) {
				std::cerr << "Device memory allocation failed for fq_d. "
							<< "Size = " << 2 * blocked_matrix_size * sizeof(complex_t) << " B"
							<< std::endl;
			} else {
				fq_double_d[0] = fq_d;
				fq_double_d[1] = fq_d + blocked_matrix_size;
	
				if(cudaMalloc((void **) &ff_d, 2 * blocked_3d_matrix_size * sizeof(complex_t))
						!= cudaSuccess) {
					std::cerr << "Device memory allocation failed for ff_d. "
								<< "Size = " << 2 * blocked_3d_matrix_size * sizeof(complex_t) << " B"
								<< std::endl;
				} else {
					ff_double_d[0] = ff_d;
					ff_double_d[1] = ff_d + blocked_3d_matrix_size;
	
					cudaMemGetInfo(&device_mem_avail, &device_mem_total);
					device_mem_used = device_mem_total - device_mem_avail;
					if(rank == 0) {
						std::cout << "++     Actual device memory used: "
									<< (float) device_mem_used / 1024 / 1024
									<< " MB" << std::endl;
						std::cout << "++       Actual host memory used: "
									<< (float) host_mem_usage / 1024 / 1024
									<< " MB" << std::endl << std::flush;
					} // if
	
					// compute the number of sub-blocks, along each of the 4 dimensions
					// formulate loops over each dimension, to go over each sub block
					unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
					unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
					unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
					unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);
	
					unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
					unsigned int curr_b_num_triangles = b_num_triangles;
					unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;
	
					if(rank == 0) {
						std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
									<< " x " << b_nqz << " x " << b_num_triangles << std::endl;
						std::cout << "++         Number of hyperblocks: " << num_blocks
									<< " [" << nb_x << " x " << nb_y << " x "
									<< nb_z << " x " << nb_t << "]"
									<< std::endl;
	
						std::cout << "++     FF kernel CUDA block size: ";
	#ifndef KERNEL2
						std::cout << block_cuda_ << std::endl;
	#else // KERNEL2
						std::cout << block_cuda_t_ << " x " << block_cuda_y_ << " x "
									<< block_cuda_z_ << std::endl;
	#endif // KERNEL2
						std::cout << "++     Reduction CUDA block size: "
	#ifndef REDUCTION2
									<< std::min((const unsigned int)BLOCK_REDUCTION_X_, b_nqx) << " x "
	#else // REDUCTION2
									<< (BLOCK_REDUCTION_T_ >> 1) << " x "
	#endif // REDUCTION2
									<< BLOCK_REDUCTION_Y_ << " x "
									<< BLOCK_REDUCTION_Z_ << std::endl;
					} // if
	
					cudaEventRecord(mem_end_event, 0);
					cudaEventSynchronize(mem_end_event);
					cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
					total_mem_time += temp_mem_time;
	
					cudaEventCreate(&start); cudaEventCreate(&stop);
	
					if(rank == 0) std::cout << "-- Computing form factor ... " << std::flush;
	
					// compute for each block with double buffering
	
					cudaStream_t stream[2];
					cudaStreamCreate(&(stream[0]));
					cudaStreamCreate(&(stream[1]));
					int active = 0, passive;
					unsigned int prev_ib_x = 0, prev_ib_y = 0, prev_ib_z = 0, prev_ib_t = 0;
					unsigned int prev_cbnqx = curr_b_nqx, prev_cbnqy = curr_b_nqy, prev_cbnqz = curr_b_nqz;
	
					cudaDeviceSynchronize();
	
					curr_b_nqx = b_nqx;
	
					for(unsigned int ib_x = 0; ib_x < nb_x; ++ ib_x) {
						if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
						curr_b_nqy = b_nqy;
	
						for(unsigned int ib_y = 0; ib_y < nb_y; ++ ib_y) {
							if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
							curr_b_nqz = b_nqz;
	
							for(unsigned int ib_z = 0; ib_z < nb_z; ++ ib_z) {
								if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
								curr_b_num_triangles = b_num_triangles;
	
								for(unsigned int ib_t = 0; ib_t < nb_t; ++ ib_t) {
									if(ib_t == nb_t - 1)
										curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;
	
									passive = 1 - active;
	
									if(ib_x + ib_y + ib_z + ib_t != 0) {
										cudaStreamSynchronize(stream[passive]);
										cudaMemcpyAsync(ff_double_buff[passive], ff_double_d[passive],
														b_nqx * b_nqy * b_nqz * sizeof(complex_t),
														cudaMemcpyDeviceToHost, stream[passive]);
									} // if
	
									// call the main kernel on the device
									cudaEventRecord(start, 0);
	#ifndef KERNEL2
									// Kernel 1: decompose along triangles
									unsigned int cuda_block_size = block_cuda_;
									unsigned int cuda_num_blocks = (unsigned int)
													ceil((float) curr_b_num_triangles / cuda_block_size);
	
									size_t dyna_shmem_size = 0;
									size_t stat_shmem_size = 0;
									size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
									if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
										std::cout << std::endl
													<< "++              FF shared memory: "
													<< (float) total_shmem_size
													<< " B" << std::endl;
									} // if
									if((float) total_shmem_size > 49152) {
										std::cerr << "error: too much shared memory requested!" << std::endl;
										exit(1);
									} // if

									form_factor_kernel<float_t>
									<<< cuda_num_blocks, cuda_block_size, 0, stream[active] >>> (
										qx_d, qy_d, qz_d, shape_def_d, axes_d,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										fq_double_d[active]);
	#else // KERNEL2
									// Kernel 4
									unsigned int ff_t_blocks = (unsigned int)
													ceil((float) curr_b_num_triangles / block_cuda_t_);
									unsigned int ff_y_blocks = (unsigned int)
													ceil((float) curr_b_nqy / block_cuda_y_);
									unsigned int ff_z_blocks = (unsigned int)
													ceil((float) curr_b_nqz / block_cuda_z_);
									dim3 ff_grid_size(ff_t_blocks, ff_y_blocks, ff_z_blocks);
									dim3 ff_block_size(block_cuda_t_, block_cuda_y_, block_cuda_z_);
	
									size_t dyna_shmem_size = sizeof(float_t) *
																(T_PROP_SIZE_ * block_cuda_t_ +
																b_nqx + block_cuda_y_) +
																sizeof(complex_t) * block_cuda_z_;
									size_t stat_shmem_size = 0;
									size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
									if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
										std::cout << std::endl
													<< "++              FF shared memory: "
													<< (float) total_shmem_size
													<< " B" << std::endl;
									} // if
									if((float) total_shmem_size > 49152) {
										std::cerr << "error: too much shared memory requested!" << std::endl;
										exit(1);
									} // if

									err = cudaGetLastError();

									form_factor_kernel_new_shared2<float_t>
									<<< ff_grid_size, ff_block_size, dyna_shmem_size, stream[active]>>> (
											qx_d, qy_d, qz_d, shape_def_d, axes_d,
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ib_x, ib_y, ib_z, ib_t,
											fq_double_d[active]);
	#endif // KERNEL2
									cudaEventRecord(stop, 0);
	
									if(err != cudaSuccess) {
										std::cerr << "CUDA error before main ff kernel ["
													<< __FILE__ << ":" << __LINE__ << "]: "
													<< cudaGetErrorString(err) << std::endl;
									} else {
										float temp_kernel;
										cudaStreamSynchronize(stream[active]);
										cudaStreamWaitEvent(stream[active], stop, 0);
										cudaEventElapsedTime(&temp_kernel, start, stop);
										kernel_time += temp_kernel;
										
										// call the reduction kernel
#ifndef REDUCTION2
										// original reduction kernel
										// perform decomposition along x, y, z
										// each thread adds all t's
										int block_rx =
												std::min((const unsigned int)BLOCK_REDUCTION_X_, b_nqx);
										int block_ry = BLOCK_REDUCTION_Y_;
										int block_rz = BLOCK_REDUCTION_Z_;
										dim3 r_block_size(block_rx, block_ry, block_rz);

										unsigned int rx_num_blocks = (unsigned int)
														ceil((float) curr_b_nqx / block_rx);
										unsigned int ry_num_blocks = (unsigned int)
														ceil((float) curr_b_nqy / block_ry);
										unsigned int rz_num_blocks = (unsigned int)
														ceil((float) curr_b_nqz / block_rz);
										dim3 r_grid_size(rx_num_blocks, ry_num_blocks, rz_num_blocks);
	
										dyna_shmem_size = 0;
										stat_shmem_size = 0;
										total_shmem_size = dyna_shmem_size + stat_shmem_size;
										if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
											std::cout << "++       Reduction shared memory: "
														<< (float) total_shmem_size
														<< " B" << std::endl;
										} // if
										if((float) total_shmem_size > 49152) {
											std::cerr << "error: too much shared memory requested!"
														<< std::endl;
											exit(1);
										} // if

										cudaEventRecord(start, 0);	// temp ...

										err = cudaGetLastError();

										reduction_kernel
										<<< r_grid_size, r_block_size, 0, stream[active] >>> (
											fq_double_d[active],
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ff_double_d[active]);
#else // REDUCTION2
										// the parallel reduction kernel
										// perform decomposition along t, y, z
										// threads cooperate to do reduction in parallel
										int block_rt = b_num_triangles >> 1; // == BLOCK_REDUCTION_T_ / 2
										int block_ry = BLOCK_REDUCTION_Y_;
										int block_rz = BLOCK_REDUCTION_Z_;
										dim3 r_block_size(block_rt, block_ry, block_rz);
										unsigned int rt_num_blocks = (unsigned int)
														ceil((float) (b_num_triangles >> 1) / block_rt);
										unsigned int ry_num_blocks = (unsigned int)
														ceil((float) curr_b_nqy / block_ry);
										unsigned int rz_num_blocks = (unsigned int)
														ceil((float) curr_b_nqz / block_rz);
										dim3 r_grid_size(rt_num_blocks, ry_num_blocks, rz_num_blocks);
	
										dyna_shmem_size = 0;
										stat_shmem_size = sizeof(complex_t) * BLOCK_REDUCTION_T_ *
															BLOCK_REDUCTION_Y_ * BLOCK_REDUCTION_Z_;
										total_shmem_size = dyna_shmem_size + stat_shmem_size;
										if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
											std::cout << "++       Reduction shared memory: "
														<< (float) total_shmem_size
														<< " B" << std::endl;
										} // if
										if((float) total_shmem_size > 49152) {
											std::cerr << "error: too much shared memory requested!"
														<< std::endl;
											exit(1);
										} // if

										cudaEventRecord(start, 0);	// temp ...

										err = cudaGetLastError();

										reduction_kernel_parallel
										<<< r_grid_size, r_block_size, 0, stream[active] >>> (
											fq_double_d[active],
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ff_double_d[active]);
#endif // REDUCTION2
										if(err != cudaSuccess) {
											std::cerr << "CUDA error before reduction kernel ["
														<< __FILE__ << ":" << __LINE__ << "]: "
														<< cudaGetErrorString(err) << std::endl;
										} // if
										cudaEventRecord(stop, 0);	// temp ...
	
										if(ib_x + ib_y + ib_z + ib_t != 0) {
											// move ff_buffer to correct location in ff
											cudaStreamSynchronize(stream[passive]);
											move_to_main_ff(ff_double_buff[passive],
															prev_cbnqx, prev_cbnqy, prev_cbnqz,
															b_nqx, b_nqy, b_nqz,
															nqx, nqy, nqz,
															prev_ib_x, prev_ib_y, prev_ib_z, ff);
										} // if
	
										float temp_reduce;
										cudaStreamSynchronize(stream[active]);
										cudaStreamWaitEvent(stream[active], stop, 0);
										cudaEventElapsedTime(&temp_reduce, start, stop);
										reduce_time += temp_reduce;
	
										active = 1 - active;
										prev_ib_x = ib_x;
										prev_ib_y = ib_y;
										prev_ib_z = ib_z;
										prev_ib_t = ib_t;
										prev_cbnqx = curr_b_nqx;
										prev_cbnqy = curr_b_nqy;
										prev_cbnqz = curr_b_nqz;
									} // if-else
								} // for ib_t
							} // for ib_z
						} // for ib_y
					} // for ib_x
	
					// for the last part
					passive = 1 - active;
					cudaStreamSynchronize(stream[passive]);
					cudaMemcpyAsync(ff_double_buff[passive], ff_double_d[passive],
									b_nqx * b_nqy * b_nqz * sizeof(complex_t),
									cudaMemcpyDeviceToHost, stream[passive]);
					cudaStreamSynchronize(stream[passive]);
					move_to_main_ff(ff_double_buff[passive],
									curr_b_nqx, curr_b_nqy, curr_b_nqz,
									b_nqx, b_nqy, b_nqz,
									nqx, nqy, nqz,
									nb_x - 1, nb_y - 1, nb_z - 1, ff);
	
					cudaDeviceSynchronize();
					cudaStreamDestroy(stream[0]);
					cudaStreamDestroy(stream[1]);
	
					cudaFree(ff_d);
				} // if-else
				cudaFree(fq_d);
			} // if-else
			if(ff_buffer != NULL) cudaFreeHost(ff_buffer);
		} // if-else
	
		cudaFree(shape_def_d);
		cudaFree(qz_d);
		cudaFree(qy_d);
		cudaFree(qx_d);
	
		float total_time;
		cudaEventRecord(total_end_event, 0);
		cudaEventSynchronize(total_end_event);
		cudaEventElapsedTime(&total_time, total_begin_event, total_end_event);
	
		//if(rank == 0) std::cout << "TOTAL TIME: " << total_time / 1000 << std::endl;
	
		pass_kernel_time = kernel_time;
		red_time = reduce_time;
		mem_time = total_mem_time;
	
		return num_triangles;
	} // GFormFactor::compute_form_factor_db()

	/**
	 * Function to compute the decomposition block size
	 * Improve it later ...
	 */
	template<typename float_t, typename complex_t>
	void GFormFactor<float_t, complex_t>::compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
			unsigned long int estimated_device_mem_need, unsigned long int device_mem_avail,
			unsigned int& b_nqx, unsigned int& b_nqy, unsigned int& b_nqz, unsigned int& b_num_triangles
	#ifdef FINDBLOCK
			, const int block_x, const int block_y, const int block_z, const int block_t
	#endif
			) {
		b_nqx = (unsigned int) nqx; b_nqy = (unsigned int) nqy; b_nqz = (unsigned int) nqz;
		b_num_triangles = (unsigned int) num_triangles;

	#ifdef FINDBLOCK
		b_nqx = (b_nqx > block_x) ? block_x : b_nqx;
		b_nqy = (b_nqy > block_y) ? block_y : b_nqy;
		b_nqz = (b_nqz > block_z) ? block_z : b_nqz;
		b_num_triangles = (b_num_triangles > block_t) ? block_t : b_num_triangles;
	#else
		b_nqx = (b_nqx > BLOCK_X_) ? BLOCK_X_ : b_nqx;
		b_nqy = (b_nqy > BLOCK_Y_) ? BLOCK_Y_ : b_nqy;
		b_nqz = (b_nqz > BLOCK_Z_) ? BLOCK_Z_ : b_nqz;
		b_num_triangles = (b_num_triangles > BLOCK_T_) ? BLOCK_T_ : b_num_triangles;
	#endif

		unsigned long int estimate = (unsigned long int) b_nqx * b_nqy * b_nqz * b_num_triangles *
										sizeof(complex_t);
		unsigned int i = 0;
		while((estimate + DEVICE_RESERVE_MEM_) > device_mem_avail) {
			-- b_nqx; -- b_nqy; -- b_nqz;
			estimate = (unsigned long int) b_nqx * b_nqy * b_nqz * b_num_triangles * sizeof(complex_t);
			++ i;
		} // if
	} // GFormFactor::compute_hyperblock_size()


	/**
	 * Function to read the input shape file.
	 */
	/*template<typename float_t, typename complex_t>
	unsigned int GFormFactor<float_t, complex_t>::read_shape_surface_file(char* filename,
			std::vector<float_t> &shape_def) {
		std::ifstream f(filename);
		if(!f.is_open()) {
			std::cerr << "Cannot open file " << filename << std::endl;
			exit(1);
		} // if
		float_t s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;

		while(true) {
			f >> s;
			if(f.eof() || !f.good()) break;
			f >> nx; f >> ny; f >> nz;
			f >> cx; f >> cy; f >> cz;
			shape_def.push_back(s);
			shape_def.push_back(nx);
			shape_def.push_back(ny);
			shape_def.push_back(nz);
			shape_def.push_back(cx);
			shape_def.push_back(cy);
			shape_def.push_back(cz);
		} // while

		f.close();
		return shape_def.size() / 7;
	} // GFormFactor::read_shape_surface_file() */


	/**
	 * Function to move a computed block of ff to its right place in ff.
	 */
	template<typename float_t, typename complex_t>
	void GFormFactor<float_t, complex_t>::move_to_main_ff(complex_t* ff_buffer,
				unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
				unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
				unsigned int nqx, unsigned int nqy, unsigned int nqz,
				unsigned int ib_x, unsigned int ib_y, unsigned int ib_z,
				complex_t* ff) {
		unsigned long int temp1 = nqx * nqy;
		unsigned long int temp2 = curr_b_nqx * curr_b_nqy;
		unsigned long int base_i = nqx * nqy * ib_z * b_nqz + nqx * ib_y * b_nqy + ib_x * b_nqx;
		// make these copy operations contiguous ... (very low priority)
		for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
			unsigned long int start_i = base_i + temp1 * i_z;
			unsigned long int super_i = 0;
			for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
				super_i = start_i + nqx * i_y;
				for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
					unsigned long int final_i = super_i + i_x;
					unsigned long int block_i = temp2 * i_z + curr_b_nqx * i_y + i_x;
					ff[final_i].x = ff[final_i].x + ff_buffer[block_i].x;
					ff[final_i].y = ff[final_i].y + ff_buffer[block_i].y;
				} // for i_x
			} // for i_y
		} // for i_z
	} // move_to_main_ff()


	/**
	 * the main Form Factor kernel functions called from host.
	 */

	/* K1:
	 * original kernel, with 1D decompostion along t; no shared memory
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel(
					float_t* qx_d, float_t* qy_d, complex_t* qz_d,
					float_t* shape_def_d, short int* axes,
					unsigned int curr_nqx, unsigned int curr_nqy,
					unsigned int curr_nqz, unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy,
					unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y,
					unsigned int ib_z, unsigned int ib_t,
					complex_t* fq_d) {
		// each thread is responsible for a different triangle
		unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;

		if(i < curr_num_triangles) {
			unsigned int shape_off = (ib_t * b_num_triangles + i) * 7;
			float_t s = shape_def_d[shape_off];
			float_t nx = shape_def_d[shape_off + axes[0] + 1];
			float_t ny = shape_def_d[shape_off + axes[1] + 1];
			float_t nz = shape_def_d[shape_off + axes[2] + 1];
			float_t x = shape_def_d[shape_off + axes[0] + 4];
			float_t y = shape_def_d[shape_off + axes[1] + 4];
			float_t z = shape_def_d[shape_off + axes[2] + 4];

			unsigned long int xy_size = (unsigned long int) curr_nqx * curr_nqy;
			unsigned long int matrix_off = xy_size * curr_nqz * i;
			unsigned int start_z = b_nqz * ib_z;
			unsigned int start_y = b_nqy * ib_y;
			unsigned int start_x = b_nqx * ib_x;

			for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz; ++ i_z, ++ global_i_z) {
				unsigned long int off_start = matrix_off + xy_size * i_z;
				complex_t temp_z = qz_d[global_i_z];
				complex_t qz2, qzn, qzt;
				compute_z(temp_z, nz, z, qz2, qzn, qzt);

				for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy; ++ i_y, ++ global_i_y) {
					unsigned long int xy_off_start = (unsigned long int) curr_nqx * i_y;
					float_t temp_y = qy_d[global_i_y];
					float_t qy2 = temp_y * temp_y;
					float_t qyn = temp_y * ny;
					float_t qyt = temp_y * y;
					
					for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx; ++ i_x, ++ global_i_x) {
						unsigned long int off = off_start + xy_off_start + i_x;
						float_t temp_x = qx_d[global_i_x];
						complex_t qn_d, qt_d;
						compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);

						fq_d[off] = compute_fq(s, qt_d, qn_d);
					} // for z
				} // for y
			} // for x
		} // if
	} // GFormFactor::form_factor_kernel()

	/* K2:
	 * kernel with 3D decomposition along t, y, z; no shared memory
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new(
						float_t*, float_t*, complex_t*, float_t*, short int*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						complex_t*) {
		// ...
	} // form_factor_kernel_new()

	/* K3:
	 * kernel with 3D decomposition along t, y, z; static shared mem for input
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared(
						float_t*, float_t*, complex_t*, float_t*, short int*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						complex_t*) {
		// ...
	} // form_factor_kernel_new_shared()


	extern __shared__ float_t dynamic_shared[];

	/* K4: default kernel
	 * kernel with 3D decomposition along t, y, z; dynamic shared mem for input, none/static for output
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2(
						const float_t* qx, const float_t* qy, const complex_t* qz,
						const float_t* shape_def, const short int* axes,
						const unsigned int curr_nqx, const unsigned int curr_nqy,
						const unsigned int curr_nqz, const unsigned int curr_num_triangles,
						const unsigned int b_nqx, const unsigned int b_nqy,
						const unsigned int b_nqz, const unsigned int b_num_triangles,
						const unsigned int ib_x, const unsigned int ib_y,
						const unsigned int ib_z, const unsigned int ib_t,
						complex_t* fq) {
		unsigned int i_t = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_y = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int i_z = blockDim.z * blockIdx.z + threadIdx.z;
		unsigned int i_thread = blockDim.x * blockDim.y * threadIdx.z +
								blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y * blockDim.z;

		// sizes are:	shared_shape_def = T_PROP_SIZE_ * blockDim.x
		// 				shared_qx = curr_nqx	// the whole qx
		// 				shared_qy = blockDim.y
		// 				shared_qz = blockDim.z
		// make these read only ... ?
		// reversed to fix alignment:
		float_t *shared_shape_def = (float_t*) dynamic_shared;
		complex_t *shared_qz = (complex_t*) &shared_shape_def[T_PROP_SIZE_ * blockDim.x];
		float_t *shared_qy = (float_t*) &shared_qz[blockDim.z];
		float_t *shared_qx = (float_t*) &shared_qy[blockDim.y];

		unsigned int i_shared, base_offset, num_loads;

		// load triangles
		unsigned int shape_def_size = T_PROP_SIZE_ * blockDim.x;
		num_loads = __float2int_ru(__fdividef(__int2float_ru(shape_def_size), num_threads));
		base_offset = T_PROP_SIZE_ * (b_num_triangles * ib_t + blockDim.x * blockIdx.x);
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < shape_def_size) shared_shape_def[i_shared] = shape_def[base_offset + i_shared];
		} // for

		// load qx
		num_loads = __float2uint_ru(__fdividef(__uint2float_ru(curr_nqx), num_threads));
		base_offset = b_nqx * ib_x;             // all qx of this hyperblock need to be loaded
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < curr_nqx)
				shared_qx[i_shared] = qx[base_offset + i_shared];
		} // for

		// load qy
		unsigned int i_qy = b_nqy * ib_y + i_y;
		if(threadIdx.x == 0 && threadIdx.z == 0 && i_y < curr_nqy)
			shared_qy[threadIdx.y] = qy[i_qy];	// M: spread about access ...

		// load qz
		unsigned int i_qz = b_nqz * ib_z + i_z;
		if(threadIdx.x == 0 && threadIdx.y == 0 && i_z < curr_nqz)
			shared_qz[threadIdx.z] = qz[i_qz];	// M: spread about access ...

		__syncthreads();	// sync to make sure all data is loaded and available

		// shared output buffer (static)
		//__shared__ complex_t shared_fq[(FQ_COPY_SIZE_ + BANK_OFF_) * MAX_NUM_THREADS_];

		if(i_t < curr_num_triangles && i_y < curr_nqy && i_z < curr_nqz) {
			unsigned int shape_off = T_PROP_SIZE_ * threadIdx.x;
			// this may be improved by making all accesses contiguous by reorganizing shared mem data ...
			float_t s = shared_shape_def[shape_off];
			float_t nx = shared_shape_def[shape_off + 1];
			float_t ny = shared_shape_def[shape_off + 2];
			float_t nz = shared_shape_def[shape_off + 3];
			float_t x = shared_shape_def[shape_off + 4];
			float_t y = shared_shape_def[shape_off + 5];
			float_t z = shared_shape_def[shape_off + 6];

			float_t temp_y = shared_qy[threadIdx.y];
			float_t qy2 = temp_y * temp_y;
			float_t qyn = temp_y * ny;
			float_t qyt = temp_y * y;
			complex_t temp_z = shared_qz[threadIdx.z];
			complex_t qz2, qzn, qzt;
			compute_z(temp_z, nz, z, qz2, qzn, qzt);

			unsigned int fq_base = curr_nqx * curr_nqy * curr_nqz * i_t +
									curr_nqx * curr_nqy * i_z +	curr_nqx * i_y;
			for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
				float_t temp_x = shared_qx[i_x];
				complex_t qn_d, qt_d;
				compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
				complex_t fq_temp = compute_fq(s, qt_d, qn_d);
				unsigned int fq_index = fq_base + i_x;
				fq[fq_index] = fq_temp;
			} // for */

			// use shared memory for output
			/*int num_iter = __float2int_ru(__fdividef(__uint2float_ru(curr_nqx), FQ_COPY_SIZE_F_));
			unsigned int shared_base = (FQ_COPY_SIZE_ + BANK_OFF_) * 
										(blockDim.y * (blockDim.z * threadIdx.x + threadIdx.z) +
										threadIdx.y);
			bool thread_odd = i_thread & 1;
			int num_copy_iter = __float2int_ru(__fdividef(FQ_COPY_SIZE_F_, num_threads));
			int num_copy_rem = curr_nqx % FQ_COPY_SIZE_;
			//int num_copy_rem = __float2int_ru(remainderf(curr_nqx, FQ_COPY_SIZE_));
			//unsigned int bx = __umulhi(blockDim.x, blockIdx.x);
			//unsigned int by = __umulhi(blockDim.y, blockIdx.y);
			//unsigned int bz = __umulhi(blockDim.z, blockIdx.z);
			unsigned int bx = blockDim.x * blockIdx.x;
			unsigned int by = blockDim.y * blockIdx.y;
			unsigned int bz = blockDim.z * blockIdx.z;

			for(unsigned int i_x = 0, qx_i_x = 0; i_x < num_iter; ++ i_x) {
				for(int i_xx = 0; i_xx < FQ_COPY_SIZE_ && qx_i_x < curr_nqx; ++ i_xx) {
					float_t temp_x = shared_qx[qx_i_x ++];
					complex_t qn_d, qt_d;
					compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
					complex_t fq_temp = compute_fq(s, qt_d, qn_d);

					// odd-evenized to remove shared bank write conflicts
					if(thread_odd) shared_fq[shared_base + i_xx] = fq_temp;
					if(!thread_odd) shared_fq[shared_base + i_xx] = fq_temp;
				} // for

				__syncthreads();

				for(int i_tt = 0; i_tt < blockDim.x; ++ i_tt) {
					unsigned int temp1 = curr_nqz * (bx + i_tt);
					unsigned int stemp1 = blockDim.z * i_tt;
					for(int i_zz = 0; i_zz < blockDim.z; ++ i_zz) {
						unsigned int temp2 = curr_nqy * (temp1 + bz + i_zz);
						unsigned int stemp2 = blockDim.y * (stemp1 + i_zz);
						for(int i_yy = 0; i_yy < blockDim.y; ++ i_yy) {
							unsigned int temp3 = curr_nqx * (temp2 + by + i_yy);
							unsigned int stemp3 = (FQ_COPY_SIZE_ + BANK_OFF_) * (stemp2 + i_yy);
							for(int copy_i = 0; copy_i < num_copy_iter; ++ copy_i) {
								if(i_thread < FQ_COPY_SIZE_ && i_thread < curr_nqx &&
										FQ_COPY_SIZE_ * copy_i + num_copy_rem - 1 < curr_nqx) {
									unsigned int fq_index = temp3 + FQ_COPY_SIZE_ * copy_i + i_thread;
									unsigned int shared_index = stemp3 + num_threads * copy_i + i_thread;
									fq[fq_index] = shared_fq[shared_index];
								} // if
							} // for
						} // for y
					} // for z
				} // for t
			} // for x */
		} // if
	} // form_factor_kernel_new_shared2()

	/* K8:
	 * kernel with 3D decomposition along t, y, z; dynamic shared mem for input, none for output (K4)
	 * and includes reduction
	 * INCOMPLETE CANNOT DO IT YET !!!!
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_red(
						const float_t* qx, const float_t* qy, const complex_t* qz,
						const float_t* shape_def, const short int* axes,
						const unsigned int curr_nqx, const unsigned int curr_nqy,
						const unsigned int curr_nqz, const unsigned int curr_num_triangles,
						const unsigned int b_nqx, const unsigned int b_nqy,
						const unsigned int b_nqz, const unsigned int b_num_triangles,
						const unsigned int ib_x, const unsigned int ib_y,
						const unsigned int ib_z, const unsigned int ib_t,
						complex_t* fq, complex_t* ff) {
		unsigned int i_t = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_y = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int i_z = blockDim.z * blockIdx.z + threadIdx.z;
		unsigned int i_thread = blockDim.x * blockDim.y * threadIdx.z +
								blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y * blockDim.z;

		float_t *shared_shape_def = (float_t*) dynamic_shared;
		complex_t *shared_qz = (complex_t*) &shared_shape_def[T_PROP_SIZE_ * blockDim.x];
		float_t *shared_qy = (float_t*) &shared_qz[blockDim.z];
		float_t *shared_qx = (float_t*) &shared_qy[blockDim.y];

		unsigned int i_shared, base_offset, num_loads;

		// load triangles
		unsigned int shape_def_size = T_PROP_SIZE_ * blockDim.x;
		num_loads = __float2int_ru(__fdividef(__int2float_ru(shape_def_size), num_threads));
		base_offset = T_PROP_SIZE_ * (b_num_triangles * ib_t + blockDim.x * blockIdx.x);
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < shape_def_size) shared_shape_def[i_shared] = shape_def[base_offset + i_shared];
		} // for

		// load qx
		num_loads = __float2uint_ru(__fdividef(__uint2float_ru(curr_nqx), num_threads));
		base_offset = b_nqx * ib_x;             // all qx of this hyperblock need to be loaded
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < curr_nqx)
				shared_qx[i_shared] = qx[base_offset + i_shared];
		} // for

		// load qy
		unsigned int i_qy = b_nqy * ib_y + i_y;
		if(threadIdx.x == 0 && threadIdx.z == 0 && i_y < curr_nqy)
			shared_qy[threadIdx.y] = qy[i_qy];	// M: spread about access ...

		// load qz
		unsigned int i_qz = b_nqz * ib_z + i_z;
		if(threadIdx.x == 0 && threadIdx.y == 0 && i_z < curr_nqz)
			shared_qz[threadIdx.z] = qz[i_qz];	// M: spread about access ...

		__syncthreads();	// sync to make sure all data is loaded and available

		if(i_t < curr_num_triangles && i_y < curr_nqy && i_z < curr_nqz) {
			unsigned int shape_off = T_PROP_SIZE_ * threadIdx.x;
			// this may be improved by making all accesses contiguous by reorganizing shared mem data ...
			float_t s = shared_shape_def[shape_off];
			float_t nx = shared_shape_def[shape_off + 1];
			float_t ny = shared_shape_def[shape_off + 2];
			float_t nz = shared_shape_def[shape_off + 3];
			float_t x = shared_shape_def[shape_off + 4];
			float_t y = shared_shape_def[shape_off + 5];
			float_t z = shared_shape_def[shape_off + 6];

			float_t temp_y = shared_qy[threadIdx.y];
			float_t qy2 = temp_y * temp_y;
			float_t qyn = temp_y * ny;
			float_t qyt = temp_y * y;
			complex_t temp_z = shared_qz[threadIdx.z];
			complex_t qz2, qzn, qzt;
			compute_z(temp_z, nz, z, qz2, qzn, qzt);

			unsigned int ff_base = curr_nqx * curr_nqy * i_z +	curr_nqx * i_y;
			unsigned int fq_base = curr_nqx * curr_nqy * curr_nqz * i_t + ff_base;
			for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
				float_t temp_x = shared_qx[i_x];
				complex_t qn_d, qt_d;
				compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
				complex_t fq_temp = compute_fq(s, qt_d, qn_d);
				unsigned int fq_index = fq_base + i_x;
				fq[fq_index] = fq_temp;
				// assume that ff is already initialized to 0
				// this does not work out!!!
				//fq_temp = cuCmulf(fq_temp, make_cuFloatComplex(0.0, -1.0));
				//unsigned int ff_index = ff_base + i_x;
				//// following should be atomic
				//ff[ff_index] = atomicAdd(ff[ff_index], fq_temp);
			} // for
		} // if
	} // form_factor_kernel_new_shared2_red()

	/* K5:
	 * kernel with 3D decomposition along t, y, z; dynamic shared mem for input, static for output
	 * and some memopt? ...
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_mem(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }

	/* K6:
	 * kernel with K3 and blocked along y, z
	 * INCOMPLETE ...
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared_subblock(float_t*, float_t*, complex_t*, float_t*,
										short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }

	/* K7:
	 * kernel with K1 (?) and blocked along y, z
	 * INCOMPLETE ...
	 */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_2(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }

	// decompose along y, z and t dimensions
	/*template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new(float_t* qx_d, float_t* qy_d, float_t* qz_d,
					float_t* shape_def_d, short int* axes,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					complex_t* fq_d) {
		// each thread is responsible for a different triangle, and y and z positions
		unsigned int i_t = threadIdx.x + blockDim.x * blockIdx.x;
		unsigned int i_y = threadIdx.y + blockDim.y * blockIdx.y;
		unsigned int i_z = threadIdx.z + blockDim.z * blockIdx.z;

		if(i_t < curr_num_triangles && i_y < curr_nqy && i_z < curr_nqz) {
			unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
			float_t s = shape_def_d[shape_off];
			float_t nx = shape_def_d[shape_off + axes[0] + 1];
			float_t ny = shape_def_d[shape_off + axes[1] + 1];
			float_t nz = shape_def_d[shape_off + axes[2] + 1];
			float_t x = shape_def_d[shape_off + axes[0] + 4];
			float_t y = shape_def_d[shape_off + axes[1] + 4];
			float_t z = shape_def_d[shape_off + axes[2] + 4];

			unsigned int global_i_z = b_nqz * ib_z + i_z;
			unsigned int global_i_y = b_nqy * ib_y + i_y;
			unsigned int start_x = b_nqx * ib_x;

			float_t temp_z = qz_d[global_i_z];
			float_t qz2 = temp_z * temp_z;
			float_t qzn = temp_z * nz;
			float_t qzt = temp_z * z;

			float_t temp_y = qy_d[global_i_y];
			float_t qy2 = temp_y * temp_y;
			float_t qyn = temp_y * ny;
			float_t qyt = temp_y * y;
			
			unsigned long int off_base = curr_nqx * curr_nqy * curr_nqz * i_t +
											curr_nqx * curr_nqy * i_z + curr_nqx * i_y;
	#pragma unroll 4
			for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx; ++ i_x, ++ global_i_x) {
				unsigned long int off = off_base + i_x;

				float_t temp_x = qx_d[global_i_x];

				// original
				//float_t q2 = temp_x * temp_x + qy2 + qz2;
				//float_t qn_d = (temp_x * nx + qyn + qzn) / q2;
				//float_t qt_d = temp_x * x + qyt + qzt;

				// new
				float_t q2 = fma(temp_x, temp_x, qy2 + qz2);
				float_t qn_d = fma(temp_x, nx, qyn + qzn) / q2;
				float_t qt_d = fma(temp_x, x, qyt + qzt);

				fq_d[off] = compute_fq(s, qt_d, qn_d);
			} // for x
		} // if
	} // GFormFactor::form_factor_kernel_new()
	*/

	/////// this does not give correct result ...
	/*template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_2(float_t* qx_d, float_t* qy_d, float_t* qz_d, float_t* shape_def_d,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					const int a_y, const int a_z,
					complex_t* fq_d) {
		// each thread is responsible for a different triangle, and subsets of y and z dimensions
		unsigned int i_dim_t = threadIdx.x + blockDim.x * blockIdx.x;
		unsigned int i_dim_y = threadIdx.y + blockDim.y * blockIdx.y;
		unsigned int i_dim_z = threadIdx.z + blockDim.z * blockIdx.z;

		if(i_dim_t < curr_num_triangles) {
			unsigned int shape_off = (ib_t * b_num_triangles + i_dim_t) * 7;
			float_t s = shape_def_d[shape_off];
			float_t nx = shape_def_d[shape_off + 1];
			float_t ny = shape_def_d[shape_off + 2];
			float_t nz = shape_def_d[shape_off + 3];
			float_t x = shape_def_d[shape_off + 4];
			float_t y = shape_def_d[shape_off + 5];
			float_t z = shape_def_d[shape_off + 6];

			unsigned int start_z = b_nqz * ib_z + a_z * i_dim_z;
			unsigned int start_y = b_nqy * ib_y + a_y * i_dim_y;
			unsigned int start_x = b_nqx * ib_x;

			for(unsigned int i_z = 0, global_i_z = start_z; i_z < a_z; ++ i_z, ++ global_i_z) {
				float_t temp_z = qz_d[global_i_z];
				float_t qz2 = temp_z * temp_z;
				float_t qzn = temp_z * nz;
				float_t qzt = temp_z * z;

				for(unsigned int i_y = 0, global_i_y = start_y; i_y < a_y; ++ i_y, ++ global_i_y) {
					float_t temp_y = qy_d[global_i_y];
					float_t qy2 = temp_y * temp_y;
					float_t qyn = temp_y * ny;
					float_t qyt = temp_y * y;
					
					for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx; ++ i_x, ++ global_i_x) {
						unsigned long int off = curr_nqx * curr_nqy * curr_nqz * i_dim_t +
												curr_nqx * curr_nqy * (a_z * i_dim_z + i_z) +
												curr_nqx * (a_y * i_dim_y + i_y) + i_x;
						float_t temp_x = qx_d[global_i_x];
						float_t q2 = temp_x * temp_x + qy2 + qz2;
						float_t qn_d = (temp_x * nx + qyn + qzn) / q2;
						float_t qt_d = temp_x * x + qyt + qzt;

						fq_d[off] = compute_fq(s, qt_d, qn_d);
					} // for x
				} // for y
			} // for z
		} // if
	} // GFormFactor::form_factor_kernel_new_2()
	*/

	/**
	 * Helper device functions
	 */

	/**
	 * For single precision
	 */

	__device__ cuFloatComplex compute_fq(float s, cuFloatComplex qt_d, cuFloatComplex qn_d) {
		cuFloatComplex v1 = cuCmulf(qn_d, make_cuFloatComplex(cosf(qt_d.x), sinf(qt_d.x)));
		float v2 = s * expf(qt_d.y);
		return cuCmulf(v1, make_cuFloatComplex(v2, 0.0f));
	} // compute_fq()

	__device__ void compute_z(cuFloatComplex temp_z, float nz, float z,
				cuFloatComplex& qz2, cuFloatComplex& qzn, cuFloatComplex& qzt) {
		qz2 = cuCmulf(temp_z, temp_z);
		qzn = cuCmulf(temp_z, make_cuFloatComplex(nz, 0.0f));
		qzt = cuCmulf(temp_z, make_cuFloatComplex(z, 0.0f));
	} // compute_z()

	__device__ void compute_x(float temp_x, float qy2, cuFloatComplex qz2, float nx, float qyn,
				cuFloatComplex qzn, float x, float qyt,	cuFloatComplex qzt,
				cuFloatComplex& qn_d, cuFloatComplex& qt_d) {
		cuFloatComplex q2 = cuCaddf(make_cuFloatComplex(fmaf(temp_x, temp_x, qy2), 0.0f), qz2);
		qn_d = cuCdivf(cuCaddf(make_cuFloatComplex(temp_x * nx, 0.0f),
								cuCaddf(make_cuFloatComplex(qyn, 0.0f), qzn)), q2);
		qt_d = cuCaddf(make_cuFloatComplex(fmaf(temp_x, x, qyt), 0.0f), qzt);
	} // compute_x()

	/**
	 * For double precision
	 */

	__device__ cuDoubleComplex compute_fq(double s, cuDoubleComplex qt_d, cuDoubleComplex qn_d) {
		cuDoubleComplex v1 = cuCmul(qn_d, make_cuDoubleComplex(cos(qt_d.x), sin(qt_d.x)));
		double v2 = s * exp(qt_d.y);
		return cuCmul(v1, make_cuDoubleComplex(v2, 0.0));
	} // compute_fq()

	__device__ void compute_z(cuDoubleComplex temp_z, double nz, double z,
				cuDoubleComplex& qz2, cuDoubleComplex& qzn, cuDoubleComplex& qzt) {
		qz2 = cuCmul(temp_z, temp_z);
		qzn = cuCmul(temp_z, make_cuDoubleComplex(nz, 0.0));
		qzt = cuCmul(temp_z, make_cuDoubleComplex(z, 0.0));
	} // compute_z()

	__device__ void compute_x(double temp_x, double qy2, cuDoubleComplex qz2, double nx, double qyn,
				cuDoubleComplex qzn, double x, double qyt,	cuDoubleComplex qzt,
				cuDoubleComplex& qn_d, cuDoubleComplex& qt_d) {
		cuDoubleComplex q2 = cuCadd(make_cuDoubleComplex(fma(temp_x, temp_x, qy2), 0.0), qz2);
		qn_d = cuCdiv(cuCadd(make_cuDoubleComplex(temp_x * nx, 0.0),
								cuCadd(make_cuDoubleComplex(qyn, 0.0), qzn)), q2);
		qt_d = cuCadd(make_cuDoubleComplex(fma(temp_x, x, qyt), 0.0), qzt);
	} // compute_x()

} // namespace hig

//#endif // _FF_CU_
