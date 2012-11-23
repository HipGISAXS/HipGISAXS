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

			unsigned int read_shape_surface_file(char* filename, std::vector<float_t>& shape_def);
			void compute_block_size(int nqx, int nqy, int nqz, int num_triangles,
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
		return gff.compute_form_factor(rank, shape_def, axes, ff, qx_h, nqx, qy_h, nqy, qz_h, nqz,
					kernel_time, red_time, mem_time
	#ifdef FINDBLOCK
					, block_x, block_y, block_z, block_t
	#endif
					);
	} // compute_form_factor_wrap()

	// Instantiations for float and double
	template unsigned int compute_form_factor_wrap<float, cuComplex>(
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
		unsigned int num_triangles = shape_def.size() / 7;
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

		short int *axes_h = &axes[0];
		short int *axes_d;
		err = cudaMalloc((void **) &axes_d, 3 * sizeof(short int));
		cudaMemcpy(axes_d, axes_h, 3 * sizeof(short int), cudaMemcpyHostToDevice);
		
		unsigned long int matrix_size = (unsigned long int) nqx * nqy * nqz * num_triangles;
		
		cudaMemGetInfo(&device_mem_avail, &device_mem_total);
		device_mem_used = device_mem_total - device_mem_avail;
		size_t estimated_device_mem_need = matrix_size * sizeof(complex_t) + device_mem_used;
		if(rank == 0) {
			std::cout << "++       Available device memory: " << (float) device_mem_avail / 1024 / 1024
						<< " MB" << std::endl;
			std::cout << "++  Estimated device memory need: " << (float) estimated_device_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if
		//cudaEventRecord(mem_end_event, 0);
		//cudaEventSynchronize(mem_end_event);
		//cudaEventElapsedTime(&temp_mem_time, mem_begin_event, mem_end_event);
		//total_mem_time += temp_mem_time;

		// do blocking if too large to fit in device memory
		unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
		compute_block_size(nqx, nqy, nqz, num_triangles,
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
								<< "Size = " << blocked_3d_matrix_size * sizeof(complex_t) << " B" << std::endl;
				} else {
					cudaMemGetInfo(&device_mem_avail, &device_mem_total);
					device_mem_used = device_mem_total - device_mem_avail;
					if(rank == 0) {
						std::cout << "++            Device memory used: " << (float) device_mem_used / 1024 / 1024
									<< " MB" << std::endl;
						std::cout << "++              Host memory used: " << (float) host_mem_usage / 1024 / 1024
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
									<< BLOCK_REDUCTION_ << " x " << BLOCK_REDUCTION_ << " x "
									<< BLOCK_REDUCTION_ << std::endl;
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
	#else
									// Kernel 4
									unsigned int ff_t_blocks = (unsigned int)
														ceil((float) curr_b_num_triangles / block_cuda_t_);
									unsigned int ff_y_blocks = (unsigned int)
														ceil((float) curr_b_nqy / block_cuda_y_);
									unsigned int ff_z_blocks = (unsigned int)
														ceil((float) curr_b_nqz / block_cuda_z_);
									dim3 ff_grid_size(ff_t_blocks, ff_y_blocks, ff_z_blocks);
									dim3 ff_block_size(block_cuda_t_, block_cuda_y_, block_cuda_z_);

									form_factor_kernel_new_shared2<float_t>
										<<< ff_grid_size, ff_block_size >>> (
											qx_d, qy_d, qz_d, shape_def_d, axes_d,
											curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
											b_nqx, b_nqy, b_nqz, b_num_triangles,
											ib_x, ib_y, ib_z, ib_t,
											fq_d);
	#endif
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

										int block_r = BLOCK_REDUCTION_;
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
	 * Function to compute the decomposition block size
	 * Improve it later ...
	 */
	template<typename float_t, typename complex_t>
	void GFormFactor<float_t, complex_t>::compute_block_size(int nqx, int nqy, int nqz, int num_triangles,
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
	} // GFormFactor::compute_block_size()


	/**
	 * Function to read the input shape file.
	 */
	template<typename float_t, typename complex_t>
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
	} // GFormFactor::read_shape_surface_file() 


	/**
	 * Function to move a computed block of ff to its right place in ff.
	 * overlap with GPU later ...
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
		// make these copy operations contiguous ...
	//	#pragma omp parallel
	//	{
	//		#pragma omp for
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
	//	} // omp parallel
	} // move_to_main_ff()


	/**
	 * the main Form Factor kernel functions called from host.
	 */

	/* K1: default kernel, with t decompostion, no shared memory */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel(float_t* qx_d, float_t* qy_d, complex_t* qz_d,
					float_t* shape_def_d, short int* axes,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
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
				//complex_t qz2 = temp_z * temp_z;
				//complex_t qzn = temp_z * nz;
				//complex_t qzt = temp_z * z;

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
						//complex_t q2 = temp_x * temp_x + qy2 + qz2;
						//complex_t qn_d = (temp_x * nx + qyn + qzn) / q2;
						//complex_t qt_d = temp_x * x + qyt + qzt;

						fq_d[off] = compute_fq(s, qt_d, qn_d);
					} // for z
				} // for y
			} // for x
		} // if
	} // GFormFactor::form_factor_kernel()

	/* K2: kernel with t, y, z decomposition, no shared memory */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }

	/* K3: kernel with t, y, z decomposition, static shared mem for input */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }


	extern __shared__ float_t dynamic_shared[];

	/* K4: kernel with t, y, z decomposition, dynamic shared mem for input, static for output */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2(const float_t* qx, const float_t* qy, const complex_t* qz,
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

		// shared buffers TODO ... change at the construction place too!!
		// sizes are:	shared_shape_def = T_PROP_SIZE_ * blockDim.x
		// 				shared_qx = curr_nqx	// the whole qx
		// 				shared_qy = blockDim.y
		// 				shared_qz = blockDim.z
		// make these read only ... ?
/*		float_t *shared_shape_def = (float_t*) dynamic_shared;
		float_t *shared_qx = (float_t*) &shared_shape_def[T_PROP_SIZE_ * blockDim.x];
		float_t *shared_qy = (float_t*) &shared_qx[curr_nqx];
		complex_t *shared_qz = (complex_t*) &shared_qy[blockDim.y];
		// shared output buffer (static)
		__shared__ complex_t shared_fq[(FQ_COPY_SIZE_ + BANK_OFF_) * MAX_NUM_THREADS_];

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
		if(threadIdx.x == 0 && threadIdx.z == 0)
			shared_qy[threadIdx.y] = qy[i_qy];	// M: spread about access ...

		// load qz
		unsigned int i_qz = b_nqz * ib_z + i_z;
		if(threadIdx.x == 0 && threadIdx.y == 0)
			shared_qz[threadIdx.z] = qz[i_qz];	// M: spread about access ...

		__syncthreads();

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
			float_t temp_z = shared_qz[threadIdx.z];
			float_t qy2 = temp_y * temp_y;
			float_t qyn = temp_y * ny;
			float_t qyt = temp_y * y;
			complex_t qz2 = temp_z * temp_z;
			complex_t qzn = temp_z * nz;
			complex_t qzt = temp_z * z;

			int num_iter = __float2int_ru(__fdividef(__uint2float_ru(curr_nqx), FQ_COPY_SIZE_F_));
			unsigned int shared_base = (FQ_COPY_SIZE_ + BANK_OFF_) * 
										(blockDim.y * (blockDim.z * threadIdx.x + threadIdx.z) +
										threadIdx.y);
			unsigned int block_base = curr_nqx * (curr_nqy * (curr_nqz * i_t + i_z) + i_y) + thread_idx;
			int thread_odd = thread_idx & 1;

			for(unsigned int i_x = 0, qx_i_x = 0; i_x < num_iter; ++ i_x) {
				unsigned int block_base_2 = block_base + i_x * FQ_COPY_SIZE_;
				for(int i_xx = 0; i_xx < FQ_COPY_SIZE_ && qx_i_x < curr_nqx; ++ i_xx) {
					float_t temp_x = shared_qx[qx_i_x ++];
					complex_t q2 = fma(temp_x, temp_x, qy2 + qz2);
					complex_t qt_d = fma(temp_x, x, qyt + qzt);
					complex_t qn_d = fma(temp_x, nx, qyn + qzn) / q2;
					complex_t fq_temp = compute_fq(s, qt_d, qn_d);

					// to remove shared conflicts
					if(thread_odd == 0) shared_fq[shared_base + i_xx] = fq_temp;
					if(thread_odd != 0) shared_fq[shared_base + i_xx] = fq_temp;
				} // for

				__syncthreads();

				if(i_threadx < FQ_COPY_SIZE_) {
			// this is likely wrong ... should it be ceil(fq_copy_size / num_threads) ... ???
					for(int i_ww = 0; i_ww < num_threads; ++ i_ww) {    // FIXIT: this is not entirely correct
																		// when num_threads < FQ_COPY_SIZE_
						fq[block_base_2 + i_ww * curr_nqx] =
							shared_fq[i_ww * (FQ_COPY_SIZE_ + BANK_OFF_) + i_thread];
					} // for
				} // if
			} // for x
		} // if
*/
	} // form_factor_kernel_new_shared2()

	/* K5: kernel with t, y, z decomposition, dynamic shared mem for input, static for output, memopt? ... */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared2_mem(float_t*, float_t*, complex_t*, float_t*, short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }

	/* K6: kernel with K3 and blocked y, z (probably incomplete ...) */
	template<typename float_t, typename complex_t>
	__global__ void form_factor_kernel_new_shared_subblock(float_t*, float_t*, complex_t*, float_t*,
										short int*,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										unsigned int, unsigned int, unsigned int, unsigned int,
										complex_t*) { }

	/* K7: kernel with K1 (?) and blocked y, z (incomplete ...) */
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
	 * For single precision.
	 */
	/*__device__ cuComplex compute_fq(float s, float qt_d, float qn_d) {
		float temp = s * qn_d;
		return make_cuFloatComplex(temp * cosf(qt_d), temp * sinf(qt_d));
	} // compute_fq()
	*/

	/**
	 * For double precision.
	 */
	/*__device__ cuDoubleComplex compute_fq(double s, double qt_d, double qn_d) {
		double temp = s * qn_d;
		return make_cuDoubleComplex(temp * cos(qt_d), temp * sin(qt_d));
	} // compute_fq()
	*/

	/**
	 * For single precision.
	 */
	__device__ cuFloatComplex compute_fq(float s, cuFloatComplex qt_d, cuFloatComplex qn_d) {
		cuFloatComplex v1 = cuCmulf(qn_d, make_cuFloatComplex(cosf(qt_d.x), sinf(qt_d.x)));
		float v2 = s * expf(qt_d.y);
		return cuCmulf(v1, make_cuFloatComplex(v2, 0.0f));
	} // compute_fq()

	/**
	 * For double precision.
	 */
	__device__ cuDoubleComplex compute_fq(double s, cuDoubleComplex qt_d, cuDoubleComplex qn_d) {
		cuDoubleComplex v1 = cuCmul(qn_d, make_cuDoubleComplex(cos(qt_d.x), sin(qt_d.x)));
		double v2 = s * exp(qt_d.y);
		return cuCmul(v1, make_cuDoubleComplex(v2, 0.0));
	} // compute_fq()

	//template <typename float_t, typename complex_t>
	__device__ void compute_z(cuFloatComplex temp_z, float nz, float z,
				cuFloatComplex& qz2, cuFloatComplex& qzn, cuFloatComplex& qzt) {
		qz2 = cuCmulf(temp_z, temp_z);
		qzn = cuCmulf(temp_z, make_cuFloatComplex(nz, 0.0f));
		qzt = cuCmulf(temp_z, make_cuFloatComplex(z, 0.0f));
	} // compute_z()

	//template <typename float_t, typename complex_t>
	__device__ void compute_z(cuDoubleComplex temp_z, double nz, double z,
				cuDoubleComplex& qz2, cuDoubleComplex& qzn, cuDoubleComplex& qzt) {
		qz2 = cuCmul(temp_z, temp_z);
		qzn = cuCmul(temp_z, make_cuDoubleComplex(nz, 0.0));
		qzt = cuCmul(temp_z, make_cuDoubleComplex(z, 0.0));
	} // compute_z()

	//template <typename float_t, typename complex_t>
	__device__ void compute_x(float temp_x, float qy2, cuFloatComplex qz2, float nx, float qyn,
				cuFloatComplex qzn, float x, float qyt,	cuFloatComplex qzt,
				cuFloatComplex& qn_d, cuFloatComplex& qt_d) {
		cuFloatComplex q2 = cuCaddf(make_cuFloatComplex(temp_x * temp_x + qy2, 0.0f), qz2);
		qn_d = cuCdivf(cuCaddf(make_cuFloatComplex(temp_x * nx, 0.0f),
								cuCaddf(make_cuFloatComplex(qyn, 0.0f), qzn)), q2);
		qt_d = cuCaddf(make_cuFloatComplex(temp_x * x + qyt, 0.0f), qzt);
	} // compute_x()

	//template <typename float_t, typename complex_t>
	__device__ void compute_x(double temp_x, double qy2, cuDoubleComplex qz2, double nx, double qyn,
				cuDoubleComplex qzn, double x, double qyt,	cuDoubleComplex qzt,
				cuDoubleComplex& qn_d, cuDoubleComplex& qt_d) {
		cuDoubleComplex q2 = cuCadd(make_cuDoubleComplex(temp_x * temp_x + qy2, 0.0), qz2);
		qn_d = cuCdiv(cuCadd(make_cuDoubleComplex(temp_x * nx, 0.0),
								cuCadd(make_cuDoubleComplex(qyn, 0.0), qzn)), q2);
		qt_d = cuCadd(make_cuDoubleComplex(temp_x * x + qyt, 0.0), qzt);
	} // compute_x()

} // namespace hig

//#endif // _FF_CU_
