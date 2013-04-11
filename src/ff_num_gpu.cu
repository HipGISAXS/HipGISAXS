/**
 * $Id: ff.cu 40 2012-08-21 20:43:52Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cuComplex.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "woo/timer/woo_cudatimers.hpp"

#include "typedefs.hpp"
#include "constants.hpp"
#include "parameters.hpp"
#include "reduction.cuh"

#include "ff_num_gpu.cuh"
#include "cu_complex_numeric.cuh"

namespace hig {

	/**
	 * Some forward declarations.
	 */

	/* K1: default kernel, with t decompostion, no shared memory */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*);
	/* K2: kernel with t, y, z decomposition, no shared memory */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*);
	/* K3: kernel with t, y, z decomposition, static shared mem for input */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*);
	/* K4: kernel with t, y, z decomposition, dynamic shared mem for input, static for output */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2(const float_t*, const float_t*, const cucomplex_t*,
							const float_t*, const short int*,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							cucomplex_t*);
	/* K5: kernel with t, y, z decomposition, dynamic shared mem for input, static for output, memopt? ... */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_mem(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*);
	/* K6: kernel with K3 and blocked y, z (probably incomplete ...) */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared_subblock(float_t*, float_t*, cucomplex_t*, float_t*,
										short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*);
	/* K7: kernel with K1 (?) and blocked y, z (incomplete ...) */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_2(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*);
	/* K8: K4 with reduction combined in it */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_red(const float_t*, const float_t*, const cucomplex_t*,
							const float_t*, const short int*,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							cucomplex_t*,	cucomplex_t*);
	/* K9: combines ff K4 with reduction */
	__global__ void form_factor_kernel_fused(const float_t*, const float_t*, const cucomplex_t*,
							const float_t*, const short int*,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							const unsigned int, const unsigned int, const unsigned int, const unsigned int,
							cucomplex_t*);
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


	/**
	 * Wrapper for the NumericFormFactorG class function.
	 * This is required because with templates, everything is required to be in the
	 * same file, and with cuda, the whole .cu file cannot be included in other files!
	 */
	//template<typename float_t, typename complex_t>
/*	unsigned int compute_form_factor_gpu(int rank, std::vector<float_t> &shape_def,
							unsigned int num_triangles,
							std::vector<short int> &axes,
							float_t* &qx_h, int nqx,
							float_t* &qy_h, int nqy,
							cucomplex_t* &qz_h, int nqz,
							cucomplex_t* &ff,
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
		//NumericFormFactorG<float_t, complex_t> gff(cuda_t, cuda_y, cuda_z);
		NumericFormFactorG gff(cuda_t, cuda_y, cuda_z);
	#else // default kernel
		//NumericFormFactorG<float_t, complex_t> gff(cuda_t);
		NumericFormFactorG gff(cuda_t);
	#endif
		//return gff.compute_form_factor(rank, shape_def, axes, ff, qx_h, nqx, qy_h, nqy, qz_h, nqz,
		//			kernel_time, red_time, mem_time
		return gff.compute_form_factor_db(rank, shape_def, axes, ff, qx_h, nqx, qy_h, nqy, qz_h, nqz,
					kernel_time, red_time, mem_time
	#ifdef FINDBLOCK
					, block_x, block_y, block_z, block_t
	#endif
					);
	} // compute_form_factor_gpu()

	// Instantiations for float and double
	template unsigned int compute_form_factor_gpu<float, cuFloatComplex>(
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
	template unsigned int compute_form_factor_gpu<double, cuDoubleComplex>(
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
	template unsigned int compute_form_factor_wrap<double, std::complex<double> >(
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
							);
*/
	/**
	 * The main host function called from outside, as part of the API for a single node.
	 */
	//template<typename float_t, typename complex_t>
	//unsigned int NumericFormFactorG<float_t, complex_t>::compute_form_factor(int rank,
	unsigned int NumericFormFactorG::compute_form_factor(int rank,
						std::vector<float_t> &shape_def, std::vector<short int>& axes,
						cucomplex_t* &ff,
						float_t* &qx_h, int nqx,
						float_t* &qy_h, int nqy,
						cucomplex_t* &qz_h, int nqz,
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
		ff = new (std::nothrow) cucomplex_t[total_qpoints]();	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(cucomplex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(cucomplex_t);

		#ifndef KERNEL2
			unsigned int num_triangles = shape_def.size() / 7;
		#else
			unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
		#endif // KERNEL2
		if(num_triangles < 1) return 0;

		cudaError_t err;
		float_t *qx_d, *qy_d;
		cucomplex_t *qz_d;
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
		if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qz_d. "
						<< "Size = " << nqz * sizeof(cucomplex_t) << " B" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpy(qx_d, qx_h, nqx * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

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
		size_t estimated_device_mem_need = matrix_size * sizeof(cucomplex_t) + device_mem_used;
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
		
		cucomplex_t *fq_buffer = NULL, *ff_buffer = NULL, *fq_d = NULL, *ff_d = NULL;
	#ifdef GPUR
		size_t estimated_host_mem_need = host_mem_usage + blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(cudaMallocHost(&ff_buffer, blocked_3d_matrix_size * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "Memory allocation failed for ff_buffer. blocked_3d_matrix_size = "
						<< blocked_3d_matrix_size << std::endl
						<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
	#else
		size_t estimated_host_mem_need = host_mem_usage + blocked_matrix_size * sizeof(cucomplex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += blocked_matrix_size * sizeof(cucomplex_t);
		fq_buffer = new (std::nothrow) cucomplex_t[blocked_matrix_size]();
		if(fq_buffer == NULL) {
			std::cerr << "Memory allocation failed for fq_buffer. blocked_matrix_size = "
						<< blocked_matrix_size << std::endl
						<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
	#endif // GPUR
		} else {
			if(cudaMalloc((void **) &fq_d, blocked_matrix_size * sizeof(cucomplex_t)) != cudaSuccess) {
				std::cerr << "Device memory allocation failed for fq_d. "
							<< "Size = " << blocked_matrix_size * sizeof(cucomplex_t) << " B" << std::endl;
			} else {
				if(cudaMalloc((void **) &ff_d, blocked_3d_matrix_size * sizeof(cucomplex_t)) != cudaSuccess) {
					std::cerr << "Device memory allocation failed for ff_d. "
								<< "Size = " << blocked_3d_matrix_size * sizeof(cucomplex_t) << " B"
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
					if(rank == 0) std::cout << "-- Computing form factor on GPU ... " << std::flush;
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

									//form_factor_kernel<float_t>
									form_factor_kernel
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
																	sizeof(cucomplex_t) * block_cuda_z_;
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

									//form_factor_kernel_new_shared2<float_t>
									form_factor_kernel_new_shared2
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
										cudaMemcpy(fq_buffer, fq_d, blocked_matrix_size * sizeof(cucomplex_t),
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
										cudaMemcpy(ff_buffer, ff_d, b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
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
	} // NumericFormFactorG::compute_form_factor()


	/**
	 * The main host function called from outside, This one uses double buffering
	 */
	//template<typename float_t, typename complex_t>
	//unsigned int NumericFormFactorG<float_t, complex_t>::compute_form_factor_db(int rank,
	unsigned int NumericFormFactorG::compute_form_factor_db(int rank,
						std::vector<float_t> &shape_def, std::vector<short int> &axes,
						cucomplex_t* &ff,
						float_t* &qx_h, int nqx,
						float_t* &qy_h, int nqy,
						cucomplex_t* &qz_h, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
	#ifdef FINDBLOCK
						, const int block_x, const int block_y, const int block_z, const int block_t
	#endif
						) {
		//cudaSetDevice(2);		////////////////////////////////////////////////////////////////////////////
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	
		float kernel_time = 0.0, reduce_time = 0.0;
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
		ff = new (std::nothrow) cucomplex_t[total_qpoints]();	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(cucomplex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(cucomplex_t);
	
		// read the shape file
	#ifndef KERNEL2
		unsigned int num_triangles = shape_def.size() / 7;
	#else // KERNEL2
		unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
	#endif // KERNEL2
		if(num_triangles < 1) return 0;
	
		cudaError_t err;
		float_t *qx_d, *qy_d;
		cucomplex_t *qz_d;
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
		if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qz_d. "
						<< "Size = " << nqz * sizeof(float_t) << " B" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpyAsync(qx_d, qx_h, nqx * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
	
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
		size_t estimated_device_mem_need = matrix_size * sizeof(cucomplex_t) + device_mem_used;
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
		
		cucomplex_t *ff_buffer = NULL, *fq_d = NULL, *ff_d = NULL;
		cucomplex_t *ff_double_buff[2], *fq_double_d[2], *ff_double_d[2];
	
		size_t estimated_host_mem_need = host_mem_usage + 2 * blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: "
						<< (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += 2 * blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(cudaMallocHost((void **) &ff_buffer, 2 * blocked_3d_matrix_size * sizeof(cucomplex_t))
				!= cudaSuccess) {
			std::cerr << "Memory allocation failed for ff_buffer. blocked_3d_matrix_size = "
						<< blocked_3d_matrix_size << std::endl
						<< "Host memory usage = "
						<< (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
		} else {
			ff_double_buff[0] = ff_buffer;
			ff_double_buff[1] = ff_buffer + blocked_3d_matrix_size;
	
			if(cudaMalloc((void **) &fq_d, 2 * blocked_matrix_size * sizeof(cucomplex_t)) != cudaSuccess) {
				std::cerr << "Device memory allocation failed for fq_d. "
							<< "Size = " << 2 * blocked_matrix_size * sizeof(cucomplex_t) << " B"
							<< std::endl;
			} else {
				fq_double_d[0] = fq_d;
				fq_double_d[1] = fq_d + blocked_matrix_size;
	
				if(cudaMalloc((void **) &ff_d, 2 * blocked_3d_matrix_size * sizeof(cucomplex_t))
						!= cudaSuccess) {
					std::cerr << "Device memory allocation failed for ff_d. "
								<< "Size = " << 2 * blocked_3d_matrix_size * sizeof(cucomplex_t) << " B"
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
									<< std::min((unsigned int)BLOCK_REDUCTION_X_, b_nqx) << " x "
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
	
					if(rank == 0) std::cout << "-- Computing form factor on GPU ... " << std::flush;
	
					// compute for each block with double buffering
	
					cudaStream_t stream[2];
					cudaStreamCreate(&(stream[0]));
					cudaStreamCreate(&(stream[1]));
					int active = 0, passive;
					unsigned int prev_ib_x = 0, prev_ib_y = 0, prev_ib_z = 0;
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
														b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
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
									if(total_shmem_size > 49152) {
										std::cerr << "error: too much shared memory requested!" << std::endl;
										exit(1);
									} // if

									//form_factor_kernel<float_t>
									form_factor_kernel
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
																sizeof(cucomplex_t) * block_cuda_z_;
									size_t stat_shmem_size = 0;
									size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
									if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
										std::cout << std::endl
													<< "++              FF shared memory: "
													<< (float) total_shmem_size
													<< " B" << std::endl;
									} // if
									if(total_shmem_size > 49152) {
										std::cerr << "error: too much shared memory requested!" << std::endl;
										exit(1);
									} // if

									err = cudaGetLastError();

									form_factor_kernel_new_shared2
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
										int block_rx = std::min((unsigned int)BLOCK_REDUCTION_X_, b_nqx);
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
										if(total_shmem_size > 49152) {
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
										stat_shmem_size = sizeof(cucomplex_t) * BLOCK_REDUCTION_T_ *
															BLOCK_REDUCTION_Y_ * BLOCK_REDUCTION_Z_;
										total_shmem_size = dyna_shmem_size + stat_shmem_size;
										if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
											std::cout << "++       Reduction shared memory: "
														<< (float) total_shmem_size
														<< " B" << std::endl;
										} // if
										if(total_shmem_size > 49152) {
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
									b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
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
	} // NumericFormFactorG::compute_form_factor_db()


	/**
	 * The main host function.
	 * This one fuses ff and reduction and uses double buffering
	 */
	unsigned int NumericFormFactorG::compute_form_factor_db_fused(int rank,
						std::vector<float_t> &shape_def, std::vector<short int> &axes,
						cucomplex_t* &ff,
						float_t* &qx_h, int nqx,
						float_t* &qy_h, int nqy,
						cucomplex_t* &qz_h, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
						#ifdef FINDBLOCK
							, const int block_x, const int block_y, const int block_z, const int block_t
						#endif
						) {
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
	
		float kernel_time = 0.0, reduce_time = 0.0;
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
		ff = new (std::nothrow) cucomplex_t[total_qpoints]();	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(cucomplex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(cucomplex_t);
	
		// read the shape file
		#ifndef KERNEL2
			unsigned int num_triangles = shape_def.size() / 7;
		#else // KERNEL2
			unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
		#endif // KERNEL2
		if(num_triangles < 1) return 0;
	
		cudaError_t err;
		float_t *qx_d, *qy_d;
		cucomplex_t *qz_d;
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
		if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qz_d. "
						<< "Size = " << nqz * sizeof(float_t) << " B" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpyAsync(qx_d, qx_h, nqx * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
	
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
		size_t estimated_device_mem_need = matrix_size * sizeof(cucomplex_t) + device_mem_used;
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
		///////////////////////////////////////////////////////
		//b_nqy = 512;
		//b_nqz = 512;
		//b_num_triangles = 256;
		///////////////////////////////////////////////////////
	
		unsigned long int blocked_3d_matrix_size = (unsigned long int) b_nqx * b_nqy * b_nqz;
		
		cucomplex_t *ff_buffer = NULL, *ff_d = NULL;
		cucomplex_t *ff_double_buff[2], *ff_double_d[2];
	
		size_t estimated_host_mem_need = host_mem_usage + 2 * blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: "
						<< (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		host_mem_usage += 2 * blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(cudaMallocHost((void **) &ff_buffer, 2 * blocked_3d_matrix_size * sizeof(cucomplex_t))
				!= cudaSuccess) {
			std::cerr << "Memory allocation failed for ff_buffer. blocked_3d_matrix_size = "
						<< blocked_3d_matrix_size << std::endl
						<< "Host memory usage = "
						<< (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
		} else {
			ff_double_buff[0] = ff_buffer;
			ff_double_buff[1] = ff_buffer + blocked_3d_matrix_size;

			if(cudaMalloc((void **) &ff_d, 2 * blocked_3d_matrix_size * sizeof(cucomplex_t))
					!= cudaSuccess) {
				std::cerr << "Device memory allocation failed for ff_d. "
							<< "Size = " << 2 * blocked_3d_matrix_size * sizeof(cucomplex_t) << " B"
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

				////////////////////////////////////////////////////////////////////////
				//block_cuda_y_ = 16;
				//block_cuda_z_ = 16;
				////////////////////////////////////////////////////////////////////////

				// decompose only along y and z (no t since reducing along t)

				if(rank == 0) {
					std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
								<< " x " << b_nqz << " x " << b_num_triangles << std::endl;
					std::cout << "++         Number of hyperblocks: " << num_blocks
								<< " [" << nb_x << " x " << nb_y << " x "
								<< nb_z << " x " << nb_t << "]"
								<< std::endl;

					std::cout << "++     FF kernel CUDA block size: ";
					std::cout << block_cuda_y_ << " x " << block_cuda_z_ << std::endl;
				} // if

				cudaEventRecord(mem_end_event, 0);
				cudaEventSynchronize(mem_end_event);
				cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
				total_mem_time += temp_mem_time;

				cudaEventCreate(&start); cudaEventCreate(&stop);

				if(rank == 0) std::cout << "-- Computing form factor on GPU ... " << std::flush;

				// compute for each block with double buffering

				cudaStream_t stream[2];
				cudaStreamCreate(&(stream[0]));
				cudaStreamCreate(&(stream[1]));
				int active = 0, passive;
				unsigned int prev_ib_x = 0, prev_ib_y = 0, prev_ib_z = 0;
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
													b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
													cudaMemcpyDeviceToHost, stream[passive]);
								} // if

								// call the main kernel on the device
								//unsigned int ff_t_blocks = (unsigned int)
								//				ceil((float) curr_b_num_triangles / block_cuda_t_);
								unsigned int ff_y_blocks = (unsigned int)
												ceil((float) curr_b_nqy / block_cuda_y_);
								unsigned int ff_z_blocks = (unsigned int)
												ceil((float) curr_b_nqz / block_cuda_z_);
								dim3 ff_grid_size(ff_y_blocks, ff_z_blocks);
								dim3 ff_block_size(block_cuda_y_, block_cuda_z_);
								size_t dyna_shmem_size = sizeof(float_t) *
															(T_PROP_SIZE_ * curr_b_num_triangles +
															b_nqx + block_cuda_y_) +
															sizeof(cucomplex_t) * block_cuda_z_;
								size_t stat_shmem_size = 0;
								size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
								if(rank == 0 && ib_x + ib_y + ib_z + ib_t == 0) {
									std::cout << std::endl
												<< "++              FF shared memory: "
												<< (float) total_shmem_size
												<< " B" << std::endl;
								} // if
								if(total_shmem_size > 49152) {
									std::cerr << "error: too much shared memory requested!" << std::endl;
									exit(1);
								} // if

								#pragma omp sections nowait
								{
								#pragma omp section
								{
									form_factor_kernel_fused
									<<< ff_grid_size, ff_block_size, dyna_shmem_size, stream[active] >>> (
											qx_d, qy_d, qz_d, shape_def_d, axes_d,
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ib_x, ib_y, ib_z, ib_t,
											ff_double_d[active]);
								}
								#pragma omp section
								{
									if(ib_x + ib_y + ib_z + ib_t != 0) {
										// move ff_buffer to correct location in ff
										cudaStreamSynchronize(stream[passive]);
										move_to_main_ff(ff_double_buff[passive],
														prev_cbnqx, prev_cbnqy, prev_cbnqz,
														b_nqx, b_nqy, b_nqz,
														nqx, nqy, nqz,
														prev_ib_x, prev_ib_y, prev_ib_z, ff);
									} // if
								}
								}

								active = 1 - active;
								prev_ib_x = ib_x;
								prev_ib_y = ib_y;
								prev_ib_z = ib_z;
								prev_cbnqx = curr_b_nqx;
								prev_cbnqy = curr_b_nqy;
								prev_cbnqz = curr_b_nqz;
							} // for ib_t
						} // for ib_z
					} // for ib_y
				} // for ib_x

				// for the last part
				passive = 1 - active;
				cudaStreamSynchronize(stream[passive]);
				cudaMemcpyAsync(ff_double_buff[passive], ff_double_d[passive],
								b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
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
	} // NumericFormFactorG::compute_form_factor_db_fused()


	/**
	 * The main host function.
	 * This one fuses ff and reduction and uses k buffering
	 */
	unsigned int NumericFormFactorG::compute_form_factor_kb_fused(int rank,
						std::vector<float_t> &shape_def, std::vector<short int> &axes,
						cucomplex_t* &ff,
						float_t* &qx_h, int nqx,
						float_t* &qy_h, int nqy,
						cucomplex_t* &qz_h, int nqz,
						int k,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
						#ifdef FINDBLOCK
							, const int block_x, const int block_y, const int block_z, const int block_t
						#endif
						) {
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);

		if(k < 3) {
			if(rank == 0) {
				std::cerr << "error: when using k-bufferring, k should be >= 3." << std::endl;
			} // if
			return 0;
		} // if
	
		float kernel_time = 0.0, reduce_time = 0.0;
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
		ff = new (std::nothrow) cucomplex_t[total_qpoints]();	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(cucomplex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(cucomplex_t);
	
		unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
		if(num_triangles < 1) return 0;
	
		cudaError_t err;
		float_t *qx_d, *qy_d;
		cucomplex_t *qz_d;
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
		if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "Device memory allocation failed for qz_d. "
						<< "Size = " << nqz * sizeof(float_t) << " B" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			delete[] ff;
			return 0;
		} // if
		cudaMemcpyAsync(qx_d, qx_h, nqx * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpyAsync(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
	
		float_t *shape_def_d;
		float_t *shape_def_h = &shape_def[0];
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

		short int *axes_h = &axes[0];
		short int *axes_d;
		err = cudaMalloc((void**) &axes_d, k * sizeof(short int));
		cudaMemcpyAsync(axes_d, axes_h, k * sizeof(short int), cudaMemcpyHostToDevice);
		
		unsigned long int matrix_size = (unsigned long int) nqx * nqy * nqz * num_triangles;
		
		cudaMemGetInfo(&device_mem_avail, &device_mem_total);
		device_mem_used = device_mem_total - device_mem_avail;
		size_t estimated_device_mem_need = matrix_size * sizeof(cucomplex_t) + device_mem_used;
		if(rank == 0) {
			std::cout << "++       Available device memory: "
						<< (float) device_mem_avail / 1024 / 1024
						<< " MB" << std::endl;
			//std::cout << "++  Estimated device memory need: "
			//			<< (float) estimated_device_mem_need / 1024 / 1024
			//			<< " MB" << std::endl;
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
		#ifdef FF_NUM_AUTOTUNE_HB
			// AUTOTUNING: find optimal CUDA block size
			// ... TODO: improve ...

			std::cout << "-- Autotuning hyperblock size ... " << std::endl;
			double min_time_hb = 1000000.0;
			unsigned int min_b_nqx = 1, min_b_nqy = 1, min_b_nqz = 1, min_b_num_triangles = 1;
			woo::CUDATimer at_kernel_timer, at_overhead_timer;
			at_overhead_timer.start();
			cucomplex_t* ff_temp;
			cudaMallocHost((void **) &ff_temp, nqx * nqy * nqz * sizeof(cucomplex_t));
			for(unsigned int b_nqx_i = 1; b_nqx_i <= nqx; ++ b_nqx_i) {
				for(unsigned int b_nqy_i = block_cuda_y_; b_nqy_i <= nqy; b_nqy_i += 10) {
					for(unsigned int b_nqz_i = block_cuda_z_; b_nqz_i <= nqz; b_nqz_i += 10) {
						for(unsigned int b_nt_i = 10; b_nt_i <= num_triangles; b_nt_i += 10) {
							at_kernel_timer.start();

							// compute the number of sub-blocks, along each of the 4 dimensions
							unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx_i);
							unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy_i);
							unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz_i);
							unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_nt_i);
							unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;

							unsigned int ff_y_blocks = (unsigned int) ceil((float) b_nqy_i / block_cuda_y_);
							unsigned int ff_z_blocks = (unsigned int) ceil((float) b_nqz_i / block_cuda_z_);
							dim3 ff_grid_size(ff_y_blocks, ff_z_blocks);
							dim3 ff_block_size(block_cuda_y_, block_cuda_z_);
							size_t dyna_shmem_size = sizeof(float_t) * (T_PROP_SIZE_ * b_nt_i +
														b_nqx_i + block_cuda_y_) +
														sizeof(cucomplex_t) * block_cuda_z_;
							size_t stat_shmem_size = 0;
							size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
							if(total_shmem_size > 49152) continue;

							unsigned int ff_size = b_nqx_i * b_nqy_i * b_nqz_i;

							form_factor_kernel_fused <<< ff_grid_size, ff_block_size, dyna_shmem_size >>> (
									qx_d, qy_d, qz_d, shape_def_d, axes_d,
									b_nqx_i, b_nqy_i, b_nqz_i, b_nt_i,
									b_nqx_i, b_nqy_i, b_nqz_i, b_nt_i,
									0, 0, 0, 0, ff_temp);

							at_kernel_timer.stop();
							double curr_time = at_kernel_timer.elapsed_msec();
							double tot_time = curr_time * num_blocks;
							std::cout << "## " << b_nqx_i << " x " << b_nqy_i << " x " << b_nqz_i
										<< " x " << b_nt_i << "\t" << num_blocks << " : "
										<< curr_time << "\t" << tot_time << std::endl;
							if(tot_time < min_time_hb) {
								min_time_hb = tot_time;
								min_b_nqx = b_nqx_i; min_b_nqy = b_nqy_i; min_b_nqz = b_nqz_i;
								min_b_num_triangles = b_nt_i;
							} // if

						} // for
					} // for
				} // for
			} // for
			cudaFree(ff_temp);
			at_overhead_timer.stop();

			b_nqx = min_b_nqx; b_nqy = min_b_nqy; b_nqz = min_b_nqz; b_num_triangles = min_b_num_triangles;
			if(rank == 0) {
				std::cout << "##    HBlock Autotuner overhead: " << at_overhead_timer.elapsed_msec()
							<< " ms." << std::endl;
			} // if
		#endif
	
		unsigned long int blocked_3d_matrix_size = (unsigned long int) b_nqx * b_nqy * b_nqz;
		
		cucomplex_t *ff_buffer = NULL, *ff_d = NULL;
		cucomplex_t *ff_double_buff[k], *ff_double_d[k];
	
		host_mem_usage += k * blocked_3d_matrix_size * sizeof(cucomplex_t);
		if(cudaMallocHost((void **) &ff_buffer, k * blocked_3d_matrix_size * sizeof(cucomplex_t))
				!= cudaSuccess) {
			std::cerr << "Memory allocation failed for ff_buffer. blocked_3d_matrix_size = "
						<< blocked_3d_matrix_size << std::endl
						<< "Host memory usage = "
						<< (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
		} else {
			for(int i_k = 0; i_k < k; ++ i_k) {
				ff_double_buff[i_k] = ff_buffer + i_k * blocked_3d_matrix_size;
			} // for

			if(cudaMalloc((void **) &ff_d, k * blocked_3d_matrix_size * sizeof(cucomplex_t))
					!= cudaSuccess) {
				std::cerr << "Device memory allocation failed for ff_d. "
							<< "Size = " << k * blocked_3d_matrix_size * sizeof(cucomplex_t) << " B"
							<< std::endl;
			} else {
				for(int i_k = 0; i_k < k; ++ i_k) {
					ff_double_d[i_k] = ff_d + i_k * blocked_3d_matrix_size;
				} // for

				cudaMemGetInfo(&device_mem_avail, &device_mem_total);
				device_mem_used = device_mem_total - device_mem_avail;
				if(rank == 0) {
					std::cout << "++           Device memory usage: "
								<< (float) device_mem_used / 1024 / 1024
								<< " MB" << std::endl;
					std::cout << "++             Host memory usage: "
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

					//std::cout << "++     FF kernel CUDA block size: ";
					//std::cout << block_cuda_y_ << " x " << block_cuda_z_ << std::endl;
				} // if

				cudaEventRecord(mem_end_event, 0);
				cudaEventSynchronize(mem_end_event);
				cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
				total_mem_time += temp_mem_time;


				#ifdef FF_NUM_AUTOTUNE
					// AUTOTUNING: find optimal CUDA block size
					// ... TODO: improve ...

					std::cout << "-- Autotuning CUDA block size ... " << std::flush;
					unsigned int CUDA_BLOCK_MAX_THREADS_ = 1024;
					double min_time = 1000000.0;
					unsigned int min_bc_y = 1, min_bc_z = 1;
					woo::CUDATimer at_kernel_timer, at_overhead_timer;
					at_overhead_timer.start();
					for(unsigned int bc_y = 1; bc_y <= std::min(b_nqy, CUDA_BLOCK_MAX_THREADS_); ++ bc_y) {
						for(unsigned int bc_z = 1; bc_z <= std::min(b_nqz, CUDA_BLOCK_MAX_THREADS_);
								++ bc_z) {
							// maximum number of threads per block = 1024
							if(bc_y * bc_z > CUDA_BLOCK_MAX_THREADS_) continue;
							if(bc_y * bc_z % 32 != 0) continue;
	
							at_kernel_timer.start();

							unsigned int ff_y_blocks = (unsigned int) ceil((float) b_nqy / bc_y);
							unsigned int ff_z_blocks = (unsigned int) ceil((float) b_nqz / bc_z);
							dim3 ff_grid_size(ff_y_blocks, ff_z_blocks);
							dim3 ff_block_size(bc_y, bc_z);
							size_t dyna_shmem_size = sizeof(float_t) * (T_PROP_SIZE_ * b_num_triangles +
															b_nqx + bc_y) + sizeof(cucomplex_t) * bc_z;
							size_t stat_shmem_size = 0;
							size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
							if(total_shmem_size > 49152) continue;

							form_factor_kernel_fused <<< ff_grid_size, ff_block_size, dyna_shmem_size >>> (
									qx_d, qy_d, qz_d, shape_def_d, axes_d,
									b_nqx, b_nqy, b_nqz, b_num_triangles,
									b_nqx, b_nqy, b_nqz, b_num_triangles,
									0, 0, 0, 0, ff_double_d[0]);

							at_kernel_timer.stop();
							double curr_time = at_kernel_timer.elapsed_msec();
							std::cout << "== " << bc_y << " x " << bc_z << " : " << curr_time << std::endl;
							if(curr_time < min_time) {
								min_time = curr_time;
								min_bc_y = bc_y; min_bc_z = bc_z;
							} // if
						} // for
					} // for
					at_overhead_timer.stop();
					std::cout << "done." << std::endl;

					block_cuda_y_ = min_bc_y; block_cuda_z_ = min_bc_z;
					if(rank == 0) {
						std::cout << "++           Autotuner overhead: " << at_overhead_timer.elapsed_msec()
									<< " ms." << std::endl;
					} // if
				#endif

				if(rank == 0) {
					std::cout << "++    Autotuned CUDA block size: " << block_cuda_y_ << " x "
								<< block_cuda_z_ << std::endl;
				} // if

				cudaEventCreate(&start); cudaEventCreate(&stop);

				if(rank == 0) std::cout << "-- Computing form factor on GPU ... " << std::flush;
				woo::CUDATimer fftimer;
				fftimer.start();

				// compute for each block with double buffering

				cudaStream_t stream[k];
				for(int i_k = 0; i_k < k; ++ i_k) {
					cudaStreamCreate(&(stream[i_k]));
				} // for
				int active = 0, passive;
				unsigned int prev_ib_x = 0, prev_ib_y = 0, prev_ib_z = 0;
				unsigned int prev_cbnqx = curr_b_nqx, prev_cbnqy = curr_b_nqy, prev_cbnqz = curr_b_nqz;
				unsigned int mid_ib_x[k - 2], mid_ib_y[k - 2], mid_ib_z[k - 2];
				unsigned int mid_cbnqx[k - 2], mid_cbnqy[k - 2], mid_cbnqz[k - 2];
				for(int i_k = 0; i_k < k - 2; ++ i_k) {
					mid_ib_x[i_k] = 0; mid_ib_y[i_k] = 0; mid_ib_z[i_k] = 0;
					mid_cbnqx[i_k] = curr_b_nqx; mid_cbnqy[i_k] = curr_b_nqy; mid_cbnqz[i_k] = curr_b_nqz;
				} // for

				cudaDeviceSynchronize();

				unsigned int hblock_iter = 0;
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

								passive = (active + 1) % k;

								if(hblock_iter > k - 2) {
									cudaStreamSynchronize(stream[passive]);
									cudaMemcpyAsync(ff_double_buff[passive], ff_double_d[passive],
													b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
													cudaMemcpyDeviceToHost, stream[passive]);
								} // if

								// call the ff kernel on the device

								unsigned int ff_y_blocks = (unsigned int)
												ceil((float) curr_b_nqy / block_cuda_y_);
								unsigned int ff_z_blocks = (unsigned int)
												ceil((float) curr_b_nqz / block_cuda_z_);
								dim3 ff_grid_size(ff_y_blocks, ff_z_blocks);
								dim3 ff_block_size(block_cuda_y_, block_cuda_z_);
								size_t dyna_shmem_size = sizeof(float_t) *
															(T_PROP_SIZE_ * curr_b_num_triangles +
															b_nqx + block_cuda_y_) +
															sizeof(cucomplex_t) * block_cuda_z_;
								size_t stat_shmem_size = 0;
								size_t total_shmem_size = dyna_shmem_size + stat_shmem_size;
								if(rank == 0 && hblock_iter == 0) {
									std::cout << std::endl
												<< "++              FF shared memory: "
												<< (float) total_shmem_size
												<< " B" << std::endl;
								} // if
								if(total_shmem_size > 49152) {
									std::cerr << "error: too much shared memory requested!" << std::endl;
									exit(1);
								} // if

								#pragma omp sections nowait
								{
									#pragma omp section
									form_factor_kernel_fused
									<<< ff_grid_size, ff_block_size, dyna_shmem_size, stream[active] >>> (
											qx_d, qy_d, qz_d, shape_def_d, axes_d,
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ib_x, ib_y, ib_z, ib_t,
											ff_double_d[active]);
									#pragma omp section
									if(hblock_iter > k - 2) {
										// move ff_buffer to correct location in ff
										cudaStreamSynchronize(stream[passive]);
										move_to_main_ff(ff_double_buff[passive],
														prev_cbnqx, prev_cbnqy, prev_cbnqz,
														b_nqx, b_nqy, b_nqz,
														nqx, nqy, nqz,
														prev_ib_x, prev_ib_y, prev_ib_z, ff);
									} // if
								} // pragma omp sections

								active = (active + 1) % k;

								prev_ib_x = mid_ib_x[0]; prev_ib_y = mid_ib_y[0]; prev_ib_z = mid_ib_z[0];
								prev_cbnqx = mid_cbnqx[0]; prev_cbnqy = mid_cbnqy[0];
								prev_cbnqz = mid_cbnqz[0];
								for(int i_k = 0; i_k < k - 3; ++ i_k) {
									mid_ib_x[i_k] = mid_ib_x[i_k + 1]; mid_ib_y[i_k] = mid_ib_y[i_k + 1];
									mid_ib_z[i_k] = mid_ib_z[i_k + 1];
									mid_cbnqx[i_k] = mid_cbnqx[i_k + 1]; mid_cbnqy[i_k] = mid_cbnqy[i_k + 1];
									mid_cbnqz[i_k] = mid_cbnqz[i_k + 1];
								} // for
								mid_ib_x[k - 3] = ib_x; mid_ib_y[k - 3] = ib_y; mid_ib_z[k - 3] = ib_z;
								mid_cbnqx[k - 3] = curr_b_nqx; mid_cbnqy[k - 3] = curr_b_nqy;
								mid_cbnqz[k - 3] = curr_b_nqz;

								++ hblock_iter;
							} // for ib_t
						} // for ib_z
					} // for ib_y
				} // for ib_x

				// for the last k - 1 parts

				for(int i_k = 0; i_k < k - 1; ++ i_k) {
					passive = (active + 1) % k;
					cudaStreamSynchronize(stream[passive]);
					cudaMemcpyAsync(ff_double_buff[passive], ff_double_d[passive],
									b_nqx * b_nqy * b_nqz * sizeof(cucomplex_t),
									cudaMemcpyDeviceToHost, stream[passive]);
					cudaStreamSynchronize(stream[passive]);
					move_to_main_ff(ff_double_buff[passive],
									prev_cbnqx, prev_cbnqy, prev_cbnqz,
									b_nqx, b_nqy, b_nqz,
									nqx, nqy, nqz,
									prev_ib_x, prev_ib_y, prev_ib_z, ff);
					active = (active + 1) % k;
					prev_ib_x = mid_ib_x[0];
					prev_ib_y = mid_ib_y[0];
					prev_ib_z = mid_ib_z[0];
					prev_cbnqx = mid_cbnqx[0];
					prev_cbnqy = mid_cbnqy[0];
					prev_cbnqz = mid_cbnqz[0];
					for(int j_k = 0; j_k < k - 3; ++ j_k) {
						mid_ib_x[j_k] = mid_ib_x[j_k + 1];
						mid_ib_y[j_k] = mid_ib_y[j_k + 1];
						mid_ib_z[j_k] = mid_ib_z[j_k + 1];
						mid_cbnqx[j_k] = mid_cbnqx[j_k + 1];
						mid_cbnqy[j_k] = mid_cbnqy[j_k + 1];
						mid_cbnqz[j_k] = mid_cbnqz[j_k + 1];
					} // for
					mid_ib_x[k - 3] = nb_x - (k - i_k - 1);
					mid_ib_y[k - 3] = nb_y - (k - i_k - 1);
					mid_ib_z[k - 3] = nb_z - (k - i_k - 1);
					mid_cbnqx[k - 3] = curr_b_nqx;
					mid_cbnqy[k - 3] = curr_b_nqy;
					mid_cbnqz[k - 3] = curr_b_nqz;
				} // for

				cudaDeviceSynchronize();
				for(int i_k = 0; i_k < k; ++ i_k) {
					cudaStreamDestroy(stream[i_k]);
				} // for
	
				cudaFree(ff_d);

				fftimer.stop();
				if(rank == 0) {
					std::cout << "**                FF kernel time: " << fftimer.elapsed_msec() << " ms."
								<< std::endl;
				} // if
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
	
		pass_kernel_time = kernel_time;
		red_time = reduce_time;
		mem_time = total_mem_time;
	
		return num_triangles;
	} // NumericFormFactorG::compute_form_factor_kb_fused()


	/**
	 * Function to compute the decomposition block size
	 * TODO Improve it later ...
	 */
	//template<typename float_t, typename complex_t>
	//void NumericFormFactorG<float_t, complex_t>::compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
	void NumericFormFactorG::compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
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
										sizeof(cucomplex_t);
		unsigned int i = 0;
		while((estimate + DEVICE_RESERVE_MEM_) > device_mem_avail) {
			-- b_nqx; -- b_nqy; -- b_nqz;
			estimate = (unsigned long int) b_nqx * b_nqy * b_nqz * b_num_triangles * sizeof(cucomplex_t);
			++ i;
		} // if
	} // NumericFormFactorG::compute_hyperblock_size()


	/**
	 * Function to read the input shape file.
	 */
	/*template<typename float_t, typename complex_t>
	unsigned int NumericFormFactorG<float_t, complex_t>::read_shape_surface_file(char* filename,
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
	} // NumericFormFactorG::read_shape_surface_file() */


	/**
	 * Function to move a computed block of ff to its right place in ff.
	 */
	//template<typename float_t, typename complex_t>
	//void NumericFormFactorG<float_t, complex_t>::move_to_main_ff(complex_t* ff_buffer,
	void NumericFormFactorG::move_to_main_ff(cucomplex_t* ff_buffer,
				unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
				unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
				unsigned int nqx, unsigned int nqy, unsigned int nqz,
				unsigned int ib_x, unsigned int ib_y, unsigned int ib_z,
				cucomplex_t* ff) {
		unsigned long int temp1 = nqx * nqy;
		unsigned long int temp2 = curr_b_nqx * curr_b_nqy;
		unsigned long int base_i = nqx * nqy * ib_z * b_nqz + nqx * ib_y * b_nqy + ib_x * b_nqx;
		// make these copy operations contiguous ... (very low priority)
		#pragma omp parallel for
		for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
			unsigned long int start_i = base_i + temp1 * i_z;
			unsigned long int super_i = 0;
			for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
				super_i = start_i + nqx * i_y;
				for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
					unsigned long int final_i = super_i + i_x;
					unsigned long int block_i = temp2 * i_z + curr_b_nqx * i_y + i_x;
					#pragma omp atomic
					ff[final_i].x += ff_buffer[block_i].x;
					#pragma omp atomic
					ff[final_i].y += ff_buffer[block_i].y;
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
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel(
					float_t* qx_d, float_t* qy_d, cucomplex_t* qz_d,
					float_t* shape_def_d, short int* axes,
					unsigned int curr_nqx, unsigned int curr_nqy,
					unsigned int curr_nqz, unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy,
					unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y,
					unsigned int ib_z, unsigned int ib_t,
					cucomplex_t* fq_d) {
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
				cucomplex_t temp_z = qz_d[global_i_z];
				cucomplex_t qz2, qzn, qzt;
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
						cucomplex_t qn_d, qt_d;
						compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);

						fq_d[off] = compute_fq(s, qt_d, qn_d);
					} // for z
				} // for y
			} // for x
		} // if
	} // NumericFormFactorG::form_factor_kernel()

	/* K2:
	 * kernel with 3D decomposition along t, y, z; no shared memory
	 */
	//template<typename float_t, typename cucomplex_t>
	__global__ void form_factor_kernel_new(
						float_t*, float_t*, cucomplex_t*, float_t*, short int*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cucomplex_t*) {
		// TODO ...
	} // form_factor_kernel_new()

	/* K3:
	 * kernel with 3D decomposition along t, y, z; static shared mem for input
	 */
	//template<typename float_t, typename cucomplex_t>
	__global__ void form_factor_kernel_new_shared(
						float_t*, float_t*, cucomplex_t*, float_t*, short int*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cucomplex_t*) {
		// TODO ...
	} // form_factor_kernel_new_shared()


	extern __shared__ float_t dynamic_shared[];

	/* K4: default kernel
	 * kernel with 3D decomposition along t, y, z; dynamic shared mem for input, none/static for output
	 */
	//template<typename float_t, typename cucomplex_t>
	__global__ void form_factor_kernel_new_shared2(
						const float_t* qx, const float_t* qy, const cucomplex_t* qz,
						const float_t* shape_def, const short int* axes,
						const unsigned int curr_nqx, const unsigned int curr_nqy,
						const unsigned int curr_nqz, const unsigned int curr_num_triangles,
						const unsigned int b_nqx, const unsigned int b_nqy,
						const unsigned int b_nqz, const unsigned int b_num_triangles,
						const unsigned int ib_x, const unsigned int ib_y,
						const unsigned int ib_z, const unsigned int ib_t,
						cucomplex_t* fq) {
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
		cucomplex_t *shared_qz = (cucomplex_t*) &shared_shape_def[T_PROP_SIZE_ * blockDim.x];
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
		//__shared__ cucomplex_t shared_fq[(FQ_COPY_SIZE_ + BANK_OFF_) * MAX_NUM_THREADS_];

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
			cucomplex_t temp_z = shared_qz[threadIdx.z];
			cucomplex_t qz2, qzn, qzt;
			compute_z(temp_z, nz, z, qz2, qzn, qzt);

			unsigned int fq_base = curr_nqx * curr_nqy * curr_nqz * i_t +
									curr_nqx * curr_nqy * i_z +	curr_nqx * i_y;
			for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
				float_t temp_x = shared_qx[i_x];
				cucomplex_t qn_d, qt_d;
				compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
				cucomplex_t fq_temp = compute_fq(s, qt_d, qn_d);
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
					cucomplex_t qn_d, qt_d;
					compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
					cucomplex_t fq_temp = compute_fq(s, qt_d, qn_d);

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


	__global__ void form_factor_kernel_fused(
						const float_t* qx, const float_t* qy, const cucomplex_t* qz,
						const float_t* shape_def, const short int* axes,
						const unsigned int curr_nqx, const unsigned int curr_nqy,
						const unsigned int curr_nqz, const unsigned int curr_num_triangles,
						const unsigned int b_nqx, const unsigned int b_nqy,
						const unsigned int b_nqz, const unsigned int b_num_triangles,
						const unsigned int ib_x, const unsigned int ib_y,
						const unsigned int ib_z, const unsigned int ib_t,
						cucomplex_t* ff) {
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int i_thread = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y;

		// sizes are:	shared_shape_def = T_PROP_SIZE_ * blockDim.x
		// 				shared_qx = curr_nqx	// the whole qx
		// 				shared_qy = blockDim.y
		// 				shared_qz = blockDim.z
		// reversed to fix alignment:
		float_t *shared_shape_def = (float_t*) dynamic_shared;
		cucomplex_t *shared_qz = (cucomplex_t*) &shared_shape_def[T_PROP_SIZE_ * curr_num_triangles];
		float_t *shared_qy = (float_t*) &shared_qz[blockDim.y];
		float_t *shared_qx = (float_t*) &shared_qy[blockDim.x];

		unsigned int i_shared, base_offset, num_loads;

		// load triangles
		unsigned int shape_def_size = T_PROP_SIZE_ * curr_num_triangles;
		num_loads = __float2int_ru(__fdividef(__int2float_ru(shape_def_size), num_threads));
		base_offset = T_PROP_SIZE_ * b_num_triangles * ib_t;
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < shape_def_size) shared_shape_def[i_shared] = shape_def[base_offset + i_shared];
		} // for

		// load qx
		num_loads = __float2uint_ru(__fdividef(__uint2float_ru(curr_nqx), num_threads));
		base_offset = b_nqx * ib_x;             // all qx of this hyperblock need to be loaded
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < curr_nqx) shared_qx[i_shared] = qx[base_offset + i_shared];
		} // for

		// load qy
		unsigned int i_qy = b_nqy * ib_y + i_y;
		if(threadIdx.y == 0 && i_y < curr_nqy)
			shared_qy[threadIdx.x] = qy[i_qy];	// M: spread about access ...

		// load qz
		unsigned int i_qz = b_nqz * ib_z + i_z;
		if(threadIdx.x == 0 && i_z < curr_nqz)
			shared_qz[threadIdx.y] = qz[i_qz];	// M: spread about access ...

		__syncthreads();	// sync to make sure all data is loaded and available

		cucomplex_t ff_tot = make_cuC((float_t) 0.0, (float_t) 0.0);
		if(i_y < curr_nqy && i_z < curr_nqz) {
			float_t temp_y = shared_qy[threadIdx.x];
			cucomplex_t temp_z = shared_qz[threadIdx.y];
			for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
				float_t temp_x = shared_qx[i_x];
				for(int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
					unsigned int shape_off = T_PROP_SIZE_ * i_t;
					float_t s = shared_shape_def[shape_off];
					float_t nx = shared_shape_def[shape_off + 1];
					float_t ny = shared_shape_def[shape_off + 2];
					float_t nz = shared_shape_def[shape_off + 3];
					float_t x = shared_shape_def[shape_off + 4];
					float_t y = shared_shape_def[shape_off + 5];
					float_t z = shared_shape_def[shape_off + 6];

					float_t qy2 = temp_y * temp_y;
					float_t qyn = temp_y * ny;
					float_t qyt = temp_y * y;
					cucomplex_t qz2, qzn, qzt;
					compute_z(temp_z, nz, z, qz2, qzn, qzt);

					cucomplex_t qn_d, qt_d;
					compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
					cucomplex_t fq_temp = compute_fq(s, qt_d, qn_d);
					ff_tot = ff_tot + fq_temp;
				} // for
				unsigned int i_ff = curr_nqx * curr_nqy * i_z + curr_nqx * i_y + i_x;
				ff[i_ff] = ff_tot;
			} // for
		} // if
	} // form_factor_kernel_fused()


	/* K8:
	 * kernel with 3D decomposition along t, y, z; dynamic shared mem for input, none for output (K4)
	 * and includes reduction
	 * INCOMPLETE CANNOT DO IT YET !!!!
	 */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_red(
						const float_t* qx, const float_t* qy, const cucomplex_t* qz,
						const float_t* shape_def, const short int* axes,
						const unsigned int curr_nqx, const unsigned int curr_nqy,
						const unsigned int curr_nqz, const unsigned int curr_num_triangles,
						const unsigned int b_nqx, const unsigned int b_nqy,
						const unsigned int b_nqz, const unsigned int b_num_triangles,
						const unsigned int ib_x, const unsigned int ib_y,
						const unsigned int ib_z, const unsigned int ib_t,
						cucomplex_t* fq, cucomplex_t* ff) {
		unsigned int i_t = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_y = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int i_z = blockDim.z * blockIdx.z + threadIdx.z;
		unsigned int i_thread = blockDim.x * blockDim.y * threadIdx.z +
								blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y * blockDim.z;

		float_t *shared_shape_def = (float_t*) dynamic_shared;
		cucomplex_t *shared_qz = (cucomplex_t*) &shared_shape_def[T_PROP_SIZE_ * blockDim.x];
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
			cucomplex_t temp_z = shared_qz[threadIdx.z];
			cucomplex_t qz2, qzn, qzt;
			compute_z(temp_z, nz, z, qz2, qzn, qzt);

			unsigned int ff_base = curr_nqx * curr_nqy * i_z +	curr_nqx * i_y;
			unsigned int fq_base = curr_nqx * curr_nqy * curr_nqz * i_t + ff_base;
			for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
				float_t temp_x = shared_qx[i_x];
				cucomplex_t qn_d, qt_d;
				compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
				cucomplex_t fq_temp = compute_fq(s, qt_d, qn_d);
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
	//template<typename float_t, typename cucomplex_t>
	__global__ void form_factor_kernel_new_shared2_mem(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*) { }

	/* K6:
	 * kernel with K3 and blocked along y, z
	 * INCOMPLETE ...
	 */
	//template<typename float_t, typename cucomplex_t>
	__global__ void form_factor_kernel_new_shared_subblock(float_t*, float_t*, cucomplex_t*, float_t*,
										short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*) { }

	/* K7:
	 * kernel with K1 (?) and blocked along y, z
	 * INCOMPLETE ...
	 */
	//template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_2(float_t*, float_t*, cucomplex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										cucomplex_t*) { }

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
	} // NumericFormFactorG::form_factor_kernel_new()
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
	} // NumericFormFactorG::form_factor_kernel_new_2()
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
