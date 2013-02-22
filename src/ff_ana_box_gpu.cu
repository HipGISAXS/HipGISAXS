/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Fri 22 Feb 2013 09:45:21 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <complex>
#include <cuComplex.h>

#include "ff_ana_gpu.cuh"
#include "enums.hpp"
#include "cu_complex_numeric.cuh"
#include "cu_utilities.cuh"


namespace hig {

	/**
	 * box on gpu
	 */

	__global__ void form_factor_box_kernel(unsigned int, unsigned int, unsigned int,
									float_t*, float_t*, cucomplex_t*, float_t, float_t, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									float_t*, cucomplex_t*);


	bool AnalyticFormFactorG::compute_box(const float_t tau, const float_t eta,
									const std::vector<float_t>& x,
									const std::vector<float_t>& distr_x,
									const std::vector<float_t>& y,
									const std::vector<float_t>& distr_y,
									const std::vector<float_t>& z,
									const std::vector<float_t>& distr_z,
									const float_t* rot_h, const std::vector<float_t>& transvec,
									std::vector<complex_t>& ff) {
		unsigned int n_x = x.size(), n_distr_x = distr_x.size();
		unsigned int n_y = y.size(), n_distr_y = distr_y.size();
		unsigned int n_z = z.size(), n_distr_z = distr_z.size();
		const float_t *x_h = x.empty() ? NULL : &*x.begin();
		const float_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
		const float_t *y_h = y.empty() ? NULL : &*y.begin();
		const float_t *distr_y_h = distr_y.empty() ? NULL : &*distr_y.begin();
		const float_t *z_h = z.empty() ? NULL : &*z.begin();
		const float_t *distr_z_h = distr_z.empty() ? NULL : &*distr_z.begin();

		unsigned int grid_size = nqx_ * nqy_ * nqz_;

		// construct device buffers
		float_t *x_d, *distr_x_d;
		float_t *y_d, *distr_y_d;
		float_t *z_d, *distr_z_d;

		cudaMalloc((void**) &x_d, n_x * sizeof(float_t));
		cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(float_t));
		cudaMalloc((void**) &y_d, n_y * sizeof(float_t));
		cudaMalloc((void**) &distr_y_d, n_distr_y * sizeof(float_t));
		cudaMalloc((void**) &z_d, n_z * sizeof(float_t));
		cudaMalloc((void**) &distr_z_d, n_distr_z * sizeof(float_t));

		// copy data to device buffers
		cudaMemcpy(x_d, x_h, n_x * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, n_y * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, n_z * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_y_d, distr_y_h, n_distr_y * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_z_d, distr_z_h, n_distr_z * sizeof(float_t), cudaMemcpyHostToDevice);

		run_init(rot_h, transvec);

		unsigned int cuda_block_y = 16, cuda_block_z = 8;
		unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) nqy_ / cuda_block_y);
		unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) nqz_ / cuda_block_z);
		dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
		dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

		size_t shared_mem_size = (nqx_ + cuda_block_y) * sizeof(float_t) +
									cuda_block_z * sizeof(cucomplex_t);
		if(shared_mem_size > 49152) {
			std::cerr << "Too much shared memory requested!" << std::endl;
			return false;
		} // if

		// the kernel
		//form_factor_box_kernel <<< ff_grid_size, ff_block_size >>> (
		form_factor_box_kernel <<< ff_grid_size, ff_block_size, shared_mem_size >>> (
				nqx_, nqy_, nqz_, qx_, qy_, qz_, tau, eta, rot_,
				n_x, x_d, n_distr_x, distr_x_d,
				n_y, y_d, n_distr_y, distr_y_d,
				n_z, z_d, n_distr_z, distr_z_d,
				transvec_,
				ff_);

		cudaThreadSynchronize();
		cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess) {
			std::cerr << "error: box form factor kernel failed [" << __FILE__ << ":" << __LINE__ << "]: "
						<< cudaGetErrorString(err) << std::endl;
			return false;
		} else {
			//std::cout << "block size: " << cby << " x " << cbz << ". ";
			construct_output_ff(ff);
		} // if-else

		cudaFree(distr_z_d);
		cudaFree(z_d);
		cudaFree(distr_y_d);
		cudaFree(y_d);
		cudaFree(distr_x_d);
		cudaFree(x_d);

		return true;
	} // AnalyticFormFactorG::compute_box()


/*	__global__ void form_factor_box_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
									float_t *qx, float_t *qy, cucomplex_t *qz,
									float_t tau, float_t eta, float_t *rot,
									unsigned int n_x, float_t *x, unsigned int n_distr_x, float_t *distr_x,
									unsigned int n_y, float_t *y, unsigned int n_distr_y, float_t *distr_y,
									unsigned int n_z, float_t *z, unsigned int n_distr_z, float_t *distr_z,
									float_t *transvec, cucomplex_t *ff) {
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;
		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				cucomplex_t mqx, mqy, mqz;
				compute_meshpoints(qx[i_x], qy[i_y], qz[i_z], rot, mqx, mqy, mqz);
				cucomplex_t tempa = sin(eta) * mqx;
				cucomplex_t tempb = cos(eta) * mqy;
				cucomplex_t temp_qm = tan(tau) * (tempa + tempb);
				cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
				for(unsigned int p_z = 0; p_z < n_z; ++ p_z) {
					for(unsigned int p_y = 0; p_y < n_y; ++ p_y) {
						for(unsigned int p_x = 0; p_x < n_x; ++ p_x) {
							cucomplex_t temp4 = fq_inv(mqz + temp_qm, y[p_y]);
							cucomplex_t temp8 = temp4 * cuCsinc(mqy * z[p_z]) * cuCsinc(mqx * x[p_x]);
							float_t temp9 = 4.0 * distr_x[p_x] * distr_y[p_y] * distr_z[p_z] *
											z[p_z] * x[p_x];
							temp_ff = temp_ff + temp9 * temp8;
						} // for x
					} // for y
				} // for z
				cucomplex_t temp_e = cuCexpi(mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
				unsigned int curr_index = base_index + i_x;
				ff[curr_index] = temp_ff * temp_e;
			} // for x
		} // if
	} // form_factor_box_kernel()*/


	extern __shared__ float_t dynamic_shared[];

	__global__ void form_factor_box_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
									float_t *qx, float_t *qy, cucomplex_t *qz,
									float_t tau, float_t eta, float_t *rot,
									unsigned int n_x, float_t *x, unsigned int n_distr_x, float_t *distr_x,
									unsigned int n_y, float_t *y, unsigned int n_distr_y, float_t *distr_y,
									unsigned int n_z, float_t *z, unsigned int n_distr_z, float_t *distr_z,
									float_t *transvec, cucomplex_t *ff) {
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;

		// shared buffers:
		cucomplex_t* qz_s = (cucomplex_t*) dynamic_shared;
		float_t* qx_s = (float_t*) &qz_s[blockDim.y];
		float_t* qy_s = (float_t*) &qx_s[nqx];

		// load all qx
		unsigned int i_thread = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y;
		unsigned int num_loads = ceil((float_t) nqx / num_threads);
		for(int i = 0; i < num_loads; ++ i) {
			unsigned int index = i * num_threads + i_thread;
			if(index < nqx) qx_s[index] = qx[index];
			else ;	// nop
		} // for
		// load part of qy
		if(i_y < nqy && threadIdx.y == 0)	// first row of threads
			qy_s[threadIdx.x] = qy[i_y];
		// load part of qz
		if(i_z < nqz && threadIdx.x == 0)	// first column of threads
			qz_s[threadIdx.y] = qz[i_z];

		// TODO: also put x, y, z, distr_x, distr_y, distr_z in shared mem ...
		// TODO: also put transvec in shared mem ...

		// make sure everything is in place
		__syncthreads();

		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				cucomplex_t mqx, mqy, mqz;
				compute_meshpoints(qx_s[i_x], qy_s[threadIdx.x], qz_s[threadIdx.y], rot, mqx, mqy, mqz);
				cucomplex_t tempa = sin(eta) * mqx;
				cucomplex_t tempb = cos(eta) * mqy;
				cucomplex_t temp_qm = tan(tau) * (tempa + tempb);
				cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
				for(unsigned int p_z = 0; p_z < n_z; ++ p_z) {
					for(unsigned int p_y = 0; p_y < n_y; ++ p_y) {
						for(unsigned int p_x = 0; p_x < n_x; ++ p_x) {
							cucomplex_t temp4 = fq_inv(mqz + temp_qm, y[p_y]);
							cucomplex_t temp8 = temp4 * cuCsinc(mqy * z[p_z]) * cuCsinc(mqx * x[p_x]);
							float_t temp9 = 4.0 * distr_x[p_x] * distr_y[p_y] * distr_z[p_z] *
											z[p_z] * x[p_x];
							temp_ff = temp_ff + temp9 * temp8;
						} // for x
					} // for y
				} // for z
				cucomplex_t temp_e = cuCexpi(mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
				unsigned int curr_index = base_index + i_x;
				ff[curr_index] = temp_ff * temp_e;
			} // for x
		} // if
	} // form_factor_box_kernel()

} // namespace hig

