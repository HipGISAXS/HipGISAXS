/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_cylinder_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Wed 20 Feb 2013 01:04:05 PM PST
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
	 * cylinder on gpu
	 */

	__global__ void form_factor_cylinder_kernel(unsigned int, unsigned int, unsigned int,
									float_t*, float_t*, cucomplex_t*, float_t, float_t, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									unsigned int, float_t*, unsigned int, float_t*,
									float_t*, cucomplex_t*);


	bool AnalyticFormFactorG::compute_cylinder(const float_t tau, const float_t eta,
									const std::vector<float_t>& h,
									const std::vector<float_t>& distr_h,
									const std::vector<float_t>& r,
									const std::vector<float_t>& distr_r,
									//const float_t* qx_h, const float_t* qy_h, const cucomplex_t* qz_h,
									const float_t* rot_h, const std::vector<float_t>& transvec,
									std::vector<complex_t>& ff) {
		unsigned int n_h = h.size(), n_distr_h = distr_h.size();
		unsigned int n_r = r.size(), n_distr_r = distr_r.size();
		//unsigned int n_transvec = transvec.size();
		const float_t *h_h = h.empty() ? NULL : &*h.begin();
		const float_t *distr_h_h = distr_h.empty() ? NULL : &*distr_h.begin();
		const float_t *r_h = r.empty() ? NULL : &*r.begin();
		const float_t *distr_r_h = distr_r.empty() ? NULL : &*distr_r.begin();
		//const float_t *transvec_h = transvec.empty() ? NULL : &*transvec.begin();

		unsigned int grid_size = nqx_ * nqy_ * nqz_;

		run_init(rot_h, transvec);

		// construct device buffers
		//float_t *qx_d, *qy_d;
		//cucomplex_t *qz_d, *ff_d;
		float_t *h_d, *distr_h_d;
		float_t *r_d, *distr_r_d;
		//float_t *transvec_d, *rot_d;

		//cudaMalloc((void**) &qx_d, nqx_ * sizeof(float_t));
		//cudaMalloc((void**) &qy_d, nqy_ * sizeof(float_t));
		//cudaMalloc((void**) &qz_d, nqz_ * sizeof(cucomplex_t));
		//cudaMalloc((void**) &ff_d, grid_size * sizeof(cucomplex_t));
		cudaMalloc((void**) &h_d, n_h * sizeof(float_t));
		cudaMalloc((void**) &distr_h_d, n_distr_h * sizeof(float_t));
		cudaMalloc((void**) &r_d, n_r * sizeof(float_t));
		cudaMalloc((void**) &distr_r_d, n_distr_r * sizeof(float_t));
		//cudaMalloc((void **) &transvec_d, n_transvec * sizeof(float_t));
		//cudaMalloc((void **) &rot_d, 9 * sizeof(float_t));

		// copy data to device buffers
		//cudaMemcpy(qx_d, qx_h, nqx_ * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(qy_d, qy_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(qz_d, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
		cudaMemcpy(h_d, h_h, n_h * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(r_d, r_h, n_r * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_h_d, distr_h_h, n_distr_h * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_r_d, distr_r_h, n_distr_r * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(transvec_d, transvec_h, n_transvec * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(rot_d, rot_h, 9 * sizeof(float_t), cudaMemcpyHostToDevice);

		unsigned int cuda_block_y = 16, cuda_block_z = 8;
		unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) nqy_ / cuda_block_y);
		unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) nqz_ / cuda_block_z);
		dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
		dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

		/*size_t shared_mem_size = (nqx_ + cuda_block_y) * sizeof(float_t) +
									cuda_block_z * sizeof(cucomplex_t);
		if(shared_mem_size > 49152) {
			std::cerr << "Too much shared memory requested!" << std::endl;
			return false;
		} // if*/

		// the kernel
		form_factor_cylinder_kernel <<< ff_grid_size, ff_block_size >>> (
					nqx_, nqy_, nqz_, qx_, qy_, qz_, tau, eta, rot_,
					n_h, h_d, n_distr_h, distr_h_d,
					n_r, r_d, n_distr_r, distr_r_d,
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
			/*cucomplex_t* ff_h = new (std::nothrow) cucomplex_t[grid_size];
			// copy result to host
			cudaMemcpy(ff_h, ff_, grid_size * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
			ff.clear();
			ff.reserve(grid_size);
			for(unsigned int i = 0; i < grid_size; ++ i) {
				ff.push_back(complex_t(ff_h[i].x, ff_h[i].y));
			} // for
			delete[] ff_h;*/
			construct_output_ff(ff);
		} // if-else

		//cudaFree(rot_d);
		//cudaFree(transvec_d);
		cudaFree(distr_r_d);
		cudaFree(r_d);
		cudaFree(distr_h_d);
		cudaFree(h_d);
		//cudaFree(ff_d);
		//cudaFree(qz_d);
		//cudaFree(qy_d);
		//cudaFree(qx_d);

		return true;
	} // AnalyticFormFactorG::compute_box()


	__global__ void form_factor_cylinder_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
									float_t *qx, float_t *qy, cucomplex_t *qz,
									float_t tau, float_t eta, float_t *rot,
									unsigned int n_h, float_t *h, unsigned int n_distr_h, float_t *distr_h,
									unsigned int n_r, float_t *r, unsigned int n_distr_r, float_t *distr_r,
									float_t *transvec, cucomplex_t *ff) {
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;
		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				cucomplex_t mqx = make_cuC(qy[i_y] * rot[0] + qx[i_x] * rot[1] + qz[i_z].x * rot[2],
											qz[i_z].y * rot[2]);
				cucomplex_t mqy = make_cuC(qy[i_y] * rot[3] + qx[i_x] * rot[4] + qz[i_z].x * rot[5],
											qz[i_z].y * rot[5]);
				cucomplex_t mqz = make_cuC(qy[i_y] * rot[6] + qx[i_x] * rot[7] + qz[i_z].x * rot[8],
											qz[i_z].y * rot[8]);
				cucomplex_t qpar = cuCsqrt(mqx * mqx + mqy * mqy);
				cucomplex_t temp1 = sin(eta) * mqx;
				cucomplex_t temp2 = cos(eta) * mqy;
				cucomplex_t temp_qm = tan(tau) * (temp1 + temp2);
				cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
				for(unsigned int p_r = 0; p_r < n_r; ++ p_r) {
					for(unsigned int p_h = 0; p_h < n_h; ++ p_h) {
						cucomplex_t temp3 = fq_inv(mqz + temp_qm, h[p_h]);
						cucomplex_t temp4 = qpar * r[p_r];
						cucomplex_t temp5 = cuCcbessj(temp4, 1);
						cucomplex_t temp6 = (temp5 / temp4) * temp3;
						temp_ff = temp_ff + distr_r[p_r] * distr_h[p_h] * 2 * PI_ * r[p_r] * r[p_r] * temp6;
					} // for h
				} // for r
				cucomplex_t temp_e = cuCexpi(mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
				unsigned int curr_index = base_index + i_x;
				ff[curr_index] = temp_ff * temp_e;
			} // for x
		} // if
	} // form_factor_cylinder_kernel()


/*	__device__ cucomplex_t fq_inv(cucomplex_t value, float_t y) {
		cucomplex_t temp1 = value * y / (float_t) 2.0;
		cucomplex_t temp = 2.0 * cuCexpi(temp1) * cuCsin(temp1) / value;
		if(cuCabsf(temp) < 1e-14) temp = make_cuC(y, (float_t) 0.0);
		return temp;
	} // fq_inv()
*/

} // namespace hig

