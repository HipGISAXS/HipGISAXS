/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_rand_cylinder_gpu.cu
 *  Created: Oct 16, 2012
 *  Modified: Thu 26 Sep 2013 10:25:21 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <complex>
#include <cuComplex.h>

#include "ff_ana_gpu.cuh"
#include "../../common/enums.hpp"
#include "../../numerics/gpu/cu_complex_numeric.cuh"
#include "../../utils/gpu/cu_utilities.cuh"

namespace hig {

	/**
	 * random cylinders on gpu - for SAXS
	 */

	__global__ void form_factor_random_cylinders_kernel(unsigned int, unsigned int, unsigned int,
									float_t*, float_t*, cucomplex_t*, float_t, float_t, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									float_t*, cucomplex_t*);


	__global__ void form_factor_random_cylinders_kernel_1D(unsigned int, unsigned int, unsigned int,
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

		run_init(rot, transvec, ff);

		// construct device buffers
		float_t *h_d, *distr_h_d;
		float_t *r_d, *distr_r_d;

		cudaMalloc((void**) &h_d, n_h * sizeof(float_t));
		cudaMalloc((void**) &distr_h_d, n_distr_h * sizeof(float_t));
		cudaMalloc((void**) &r_d, n_r * sizeof(float_t));
		cudaMalloc((void**) &distr_r_d, n_distr_r * sizeof(float_t));

		// decompose into blocks (to avoid kernel timeout, and reduce mem usage)
		// assuming all of the q-grid arrays are on the gpu (its done in init)

		// copy data to device buffers
		cudaMemcpy(h_d, h_h, n_h * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(r_d, r_h, n_r * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_h_d, distr_h_h, n_distr_h * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_r_d, distr_r_h, n_distr_r * sizeof(float_t), cudaMemcpyHostToDevice);

		size_t avail_device_mem, total_device_mem;
		cudaMemGetInfo(&avail_device_mem, &total_device_mem);
		unsigned int used_device_mem = total_device_mem - avail_device_mem;
		std::cout << "++           Device memory usage: " << (float) used_device_mem / 1024 / 1024
					<< " MB" << std::endl;

		unsigned int cuda_block_y = 16, cuda_block_z = 8;
		std::cout << "++        CUDA thread block size: " << cuda_block_y << " x " << cuda_block_z
					<< std::endl;

		unsigned int curr_b_nqx = b_nqx_, curr_b_nqy = b_nqy_, curr_b_nqz = b_nqz_;
		curr_b_nqz = b_nqz_;
		for(int ib_z = 0; ib_z < nb_z_; ++ ib_z) {
			if(ib_z == nb_z_ - 1) curr_b_nqz = nqz_ - b_nqz_ * ib_z;
			curr_b_nqy = b_nqy_;
			for(int ib_y = 0; ib_y < nb_y_; ++ ib_y) {
				if(ib_y == nb_y_ - 1) curr_b_nqy = nqy_ - b_nqy_ * ib_y;
				curr_b_nqx = b_nqx_;
				for(int ib_x = 0; ib_x < nb_x_; ++ ib_x) {
					if(ib_x == nb_x_ - 1) curr_b_nqx = nqx_ - b_nqx_ * ib_x;

					unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) curr_b_nqy / cuda_block_y);
					unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) curr_b_nqz / cuda_block_z);
					dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
					dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

					// the kernel
					float_t* curr_qx = qx_ + b_nqx_ * ib_x;
					float_t* curr_qy = qy_ + b_nqy_ * ib_y;
					cucomplex_t* curr_qz = qz_ + b_nqz_ * ib_z;
					//form_factor_random_cylinders_kernel <<< ff_grid_size, ff_block_size >>> (
					form_factor_random_cylinders_kernel_1D <<< ceil((float_t) curr_b_nqz / 32) , 32 >>> (
								curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_qx, curr_qy, curr_qz,
								tau, eta, rot_,
								n_h, h_d, n_distr_h, distr_h_d,
								n_r, r_d, n_distr_r, distr_r_d,
								transvec_,
								ff_buff_d_);

					cudaError_t err = cudaGetLastError();
					cudaDeviceSynchronize();
					if(err != cudaSuccess) {
						std::cerr << "error: random cylinders form factor kernel failed ["
									<< __FILE__ << ":" << __LINE__ << "]: "
									<< cudaGetErrorString(err) << std::endl;
						return false;
					} else {
						err = cudaGetLastError();
						if(err != cudaSuccess) {
							std::cerr << "error: random cylinders device sync failed ["
										<< __FILE__ << ":" << __LINE__ << "]: "
										<< cudaGetErrorString(err) << std::endl;
							return false;
						} else {
							move_ff_buff_to_host_ff(ff, curr_b_nqx, curr_b_nqy, curr_b_nqz,
													b_nqx_, b_nqy_, b_nqz_, ib_x, ib_y, ib_z);
						} // if-else
					} // if-else
				} // for x
			} // for y
		} // for z

		//construct_output_ff(ff);

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
		unsigned int nx = (unsigned int) ((1.0 - dx) / dx + 1.0);
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
						if(!(temp2.x == 0 && temp2.y == 0))
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


	__global__ void form_factor_random_cylinders_kernel_1D(
									unsigned int nqx, unsigned int nqy, unsigned int nqz,
									float_t *qx, float_t *qy, cucomplex_t *qz,
									float_t tau, float_t eta, float_t *rot,
									unsigned int n_h, float_t *h, unsigned int n_distr_h, float_t *distr_h,
									unsigned int n_r, float_t *r, unsigned int n_distr_r, float_t *distr_r,
									float_t *transvec, cucomplex_t *ff) {
		unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;
		float_t dx = 0.001;
		unsigned int nx = (unsigned int) ((1.0 - dx) / dx + 1.0);
		cucomplex_t qz_val = qz[i_z];
		if(i_z < nqz) {
			cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
			for(unsigned int p_r = 0; p_r < n_r; ++ p_r) {
				for(unsigned int p_h = 0; p_h < n_h; ++ p_h) {
					float_t x_val = 0.0;
					cucomplex_t temp_ffx = make_cuC((float_t) 0.0, (float_t) 0.0);
					for(unsigned int p_x = 0; p_x < nx; ++ p_x, x_val += dx) {
						cucomplex_t temp1 = cuCsinc(qz_val * h[p_h] * x_val / (float_t) 2.0);
						cucomplex_t temp2 = qz_val * r[p_r] * sqrt((float_t) 1.0 - x_val * x_val);
						if(!(temp2.x == 0 && temp2.y == 0))
							temp_ffx = temp_ffx + temp1 * cuCcbessj(temp2, 1) / temp2;
						else
							temp_ffx = temp_ffx + temp1 * make_cuC((float_t) 0.5, (float_t) 0.0);
					} // for
					temp_ff = temp_ff + 4.0 * distr_r[p_r] * distr_h[p_h] * temp_ffx * temp_ffx;
				} // for h
			} // for r
			unsigned int base_index = nqx * nqy * i_z;
			for(unsigned int i_y = 0; i_y < nqy; ++ i_y) {
				for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
					unsigned int curr_index = base_index + nqx * i_y + i_x;
					ff[curr_index] = temp_ff;
				} // for x
			} // for y
		} // if
	} // form_factor_random_cylinders_kernel_1D()

} // namespace hig

