/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_rand_cylinder_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Sat 02 Mar 2013 11:40:19 AM PST
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
	 * random cylinders on gpu - for SAXS
	 */

	__global__ void form_factor_random_cylinders_kernel(unsigned int, unsigned int, unsigned int,
									float_t*, float_t*, cucomplex_t*, float_t, float_t, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									float_t*, cucomplex_t*);


	bool AnalyticFormFactorG::compute_random_cylinders(const float_t tau, const float_t eta,
									const std::vector<float_t>& h,
									const std::vector<float_t>& distr_h,
									const std::vector<float_t>& r,
									const std::vector<float_t>& distr_r,
									const float_t* rot, const std::vector<float_t>& transvec,
									std::vector<complex_t>& ff) {
		unsigned int n_h = h.size(), n_distr_h = distr_h.size();
		unsigned int n_r = r.size(), n_distr_r = distr_r.size();
		const float_t *h_h = h.empty() ? NULL : &*h.begin();
		const float_t *distr_h_h = distr_h.empty() ? NULL : &*distr_h.begin();
		const float_t *r_h = r.empty() ? NULL : &*r.begin();
		const float_t *distr_r_h = distr_r.empty() ? NULL : &*distr_r.begin();

		run_init(rot, transvec);

		// construct device buffers
		float_t *h_d, *distr_h_d;
		float_t *r_d, *distr_r_d;

		cudaMalloc((void**) &h_d, n_h * sizeof(float_t));
		cudaMalloc((void**) &distr_h_d, n_distr_h * sizeof(float_t));
		cudaMalloc((void**) &r_d, n_r * sizeof(float_t));
		cudaMalloc((void**) &distr_r_d, n_distr_r * sizeof(float_t));

		// copy data to device buffers
		cudaMemcpy(h_d, h_h, n_h * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(r_d, r_h, n_r * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_h_d, distr_h_h, n_distr_h * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_r_d, distr_r_h, n_distr_r * sizeof(float_t), cudaMemcpyHostToDevice);

		unsigned int cuda_block_y = 16, cuda_block_z = 8;
		unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) nqy_ / cuda_block_y);
		unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) nqz_ / cuda_block_z);
		dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
		dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

		// the kernel
		form_factor_random_cylinders_kernel <<< ff_grid_size, ff_block_size >>> (
					nqx_, nqy_, nqz_, qx_, qy_, qz_, tau, eta, rot_,
					n_h, h_d, n_distr_h, distr_h_d,
					n_r, r_d, n_distr_r, distr_r_d,
					transvec_,
					ff_);

		cudaThreadSynchronize();
		cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess) {
			std::cerr << "error: random cylinders form factor kernel failed ["
						<< __FILE__ << ":" << __LINE__ << "]: "
						<< cudaGetErrorString(err) << std::endl;
			return false;
		} else {
			construct_output_ff(ff);
		} // if-else

		cudaFree(distr_r_d);
		cudaFree(r_d);
		cudaFree(distr_h_d);
		cudaFree(h_d);

		return true;
	} // AnalyticFormFactorG::compute_random_cylinders()


	__global__ void form_factor_random_cylinders_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
									float_t *qx, float_t *qy, cucomplex_t *qz,
									float_t tau, float_t eta, float_t *rot,
									unsigned int n_h, float_t *h, unsigned int n_distr_h, float_t *distr_h,
									unsigned int n_r, float_t *r, unsigned int n_distr_r, float_t *distr_r,
									float_t *transvec, cucomplex_t *ff) {
		// why does this not depend on eta? ... qx, qy, transvec, rot ???
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;
		float_t dx = 0.001;
		unsigned int nx = (1.0 - dx) / dx + 1;
		cucomplex_t qz_val = qz[i_z];
		if(i_y < nqy && i_z < nqz) {
			cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
			for(unsigned int p_r = 0; p_r < n_r; ++ p_r) {
				for(unsigned int p_h = 0; p_h < n_h; ++ p_h) {
					float_t x_val = 0.0;
					cucomplex_t temp_ffx = make_cuC((float_t) 0.0, (float_t) 0.0);
					for(unsigned int p_x = 0; p_x < nx; ++ p_x, x_val += dx) {
						cucomplex_t temp1 = cuCsinc(qz_val * h[p_h] * x_val / (float_t) 2.0);
						cucomplex_t temp2 = qz_val * r[p_r] * sqrt((float_t) 1.0 - x_val * x_val);
						temp_ffx = temp_ffx + temp1 * cuCcbessj(temp2, 1) / temp2;
					} // for
					temp_ff = temp_ff + 4.0 * distr_r[p_r] * distr_h[p_h] * temp_ffx * temp_ffx;
				} // for h
			} // for r
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				unsigned int curr_index = base_index + i_x;
				ff[curr_index] = temp_ff;
			} // for x
		} // if
	} // form_factor_random_cylinders_kernel()

} // namespace hig

