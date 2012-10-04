/**
 * $Id: reduction.cu 38 2012-08-09 23:01:20Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 */

#include <iostream>
#include <omp.h>
#include <cuComplex.h>

#include "reduction.hpp"

#ifndef GPUR

/**
 * OpenMP: For single precision
 */
void reduction_kernel(unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
						unsigned int curr_b_num_triangles, unsigned long int blocked_matrix_size,
						unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
						unsigned int num_triangles,
						unsigned int nqx, unsigned int nqy, unsigned int nqz,
						unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
						cuComplex* fq_buffer, cuComplex* ff) {
	unsigned long int curr_b_xyz = (unsigned long int) curr_b_nqx * curr_b_nqy * curr_b_nqz;
	// reduction over all triangles
	#pragma omp parallel
	{
		if(omp_get_thread_num() == 0)
			std::cout << "[" << omp_get_num_threads() << " threads] ... " << std::flush;
		#pragma omp for
		for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
			for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
				for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
					unsigned long int temp_i = (unsigned long int) curr_b_nqx * curr_b_nqy * i_z +
												curr_b_nqx * i_y + i_x;
					cuComplex total = make_cuFloatComplex(0.0, 0.0);
					for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
						//if(curr_b_xyz * i_t + temp_i >= blocked_matrix_size)
						//	std::cout << "OMG OMG OMG!!!!! " << curr_b_xyz << " " << i_t << " "
						//				<< temp_i << std::endl << std::flush;
						total = cuCaddf(total, fq_buffer[curr_b_xyz * i_t + temp_i]);
					} // for i_t
					unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
												nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
					ff[super_i] = cuCaddf(ff[super_i], cuCmulf(total, make_cuFloatComplex(0.0, -1.0)));
				} // for i_z
			} // for i_y
		} // for i_x
	} // omp parallel
} // reduction_kernel()


/**
 * OpenMP: For double precision
 */
void reduction_kernel(unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
						unsigned int curr_b_num_triangles, unsigned long int blocked_matrix_size,
						unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
						unsigned int num_triangles,
						unsigned int nqx, unsigned int nqy, unsigned int nqz,
						unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
						cuDoubleComplex* fq_buffer, cuDoubleComplex* ff) {
	unsigned long int curr_b_xyz = (unsigned long int) curr_b_nqx * curr_b_nqy * curr_b_nqz;
	// reduction over all triangles
	#pragma omp parallel
	{
		if(omp_get_thread_num() == 0)
			std::cout << "[" << omp_get_num_threads() << " threads] ... " << std::flush;
		#pragma omp for
		for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
			for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
				for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
					unsigned long int temp_i = (unsigned long int) curr_b_nqx * curr_b_nqy * i_z +
												curr_b_nqx * i_y + i_x;
					cuDoubleComplex total = make_cuDoubleComplex(0.0, 0.0);
					for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
						//if(curr_b_xyz * i_t + temp_i >= blocked_matrix_size)
						//	std::cout << "OMG OMG OMG!!!!! " << curr_b_xyz << " " << i_t << " "
						//				<< temp_i << std::endl << std::flush;
						total = cuCadd(total, fq_buffer[curr_b_xyz * i_t + temp_i]);
					} // for i_t
					unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
												nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
					ff[super_i] = cuCadd(ff[super_i], cuCmul(total, make_cuDoubleComplex(0.0, -1.0)));
					//ff[super_i] = make_cuDoubleComplex(ib_y, i_y);
				} // for i_z
			} // for i_y
		} // for i_x
	} // omp parallel
} // reduction_kernel()

#else	// GPUR

/**
 * GPU: For single precision
 */
__global__ void reduction_kernel(cuComplex* fq_d,
									unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
									unsigned int curr_num_triangles,
									unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
									unsigned int b_num_triangles,
									cuComplex* ff_d) {
	// 3D block, where each thread is responsible for one point in the x, y, z grid
	unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;
	unsigned int z = threadIdx.z + blockDim.z * blockIdx.z;
	unsigned int i_ff = curr_nqx * curr_nqy * z + curr_nqx * y + x;
	unsigned int temp = curr_nqx * curr_nqy * curr_nqz;

	if(x < curr_nqx && y < curr_nqy && z < curr_nqz) {
		cuComplex total = make_cuFloatComplex(0.0, 0.0);
		for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
			unsigned int i_fq = temp * i_t + i_ff;
			total = cuCaddf(total, fq_d[i_fq]);
		} // for
		ff_d[i_ff] = cuCmulf(total, make_cuFloatComplex(0.0, -1.0));
	} // if
} // reduction_kernel()

/**
 * GPU: For double precision
 */
__global__ void reduction_kernel(cuDoubleComplex* fq_d,
									unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
									unsigned int curr_num_triangles,
									unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
									unsigned int b_num_triangles,
									cuDoubleComplex* ff_d) {
	// 3D block, where each thread is responsible for one point in the x, y, z grid
	unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
	unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;
	unsigned int z = threadIdx.z + blockDim.z * blockIdx.z;
	unsigned int i_ff = curr_nqx * curr_nqy * z + curr_nqx * y + x;
	unsigned int temp = curr_nqx * curr_nqy * curr_nqz;
	
	if(x < curr_nqx && y < curr_nqy && z < curr_nqz) {
		cuDoubleComplex total = make_cuDoubleComplex(0.0, 0.0);
		for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
			unsigned int i_fq = temp * i_t + i_ff;
			total = cuCadd(total, fq_d[i_fq]);
		} // for
		ff_d[i_ff] = cuCmul(total, make_cuDoubleComplex(0.0, -1.0));
	} // if
} // reduction_kernel()

#endif // GPUR
