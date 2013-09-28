/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_sphere_gpu.cu
 *  Created: Oct 16, 2012
 *  Modified: Thu 26 Sep 2013 10:25:56 AM PDT
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

	// make gpu kernels members of the class ...

	// forward declarations of gpu kernels
	__global__ void form_factor_sphere_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
					float_t* qx, float_t* qy, cucomplex_t* mesh_qz, float_t* rot,
					unsigned int n_r, float_t* r, unsigned int n_distr_r, float_t* distr_r,
					float_t* transvec, cucomplex_t* ff);
	__device__ void ff_sphere_kernel_compute_tff(float r, float distr_r,
								cuFloatComplex q, cuFloatComplex mqz, cuFloatComplex& f);
	__device__ void ff_sphere_kernel_compute_tff(double r, double distr_r,
								cuDoubleComplex q, cuDoubleComplex mqz, cuDoubleComplex& f);
	__device__ cuFloatComplex ff_sphere_kernel_compute_ff(cuFloatComplex temp_f,
								cuFloatComplex mqx, cuFloatComplex mqy, cuFloatComplex mqz,
								float tx, float ty, float tz);
	__device__ cuDoubleComplex ff_sphere_kernel_compute_ff(cuDoubleComplex temp_f,
								cuDoubleComplex mqx, cuDoubleComplex mqy, cuDoubleComplex mqz,
								double tx, double ty, double tz);


	/**
	 * sphere host function
	 */
	bool AnalyticFormFactorG::compute_sphere(const std::vector<float_t>& r,
											const std::vector<float_t>& distr_r,
											const float_t* rot_h,
											const std::vector<float_t>& transvec,
											std::vector<complex_t>& ff) {
		unsigned int n_r = r.size();
		unsigned int n_distr_r = distr_r.size();
		const float_t *r_h = r.empty() ? NULL : &*r.begin();
		const float_t *distr_r_h = distr_r.empty() ? NULL : &*distr_r.begin();

		// construct device buffers
		float_t *r_d, *distr_r_d;
		if(cudaMalloc((void **) &r_d, n_r * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for r_d" << std::endl;
			return false;
		} // if
		if(cudaMalloc((void **) &distr_r_d, n_distr_r * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for n_distr_r_d" << std::endl;
			cudaFree(r_d);
			return false;
		} // if

		// copy data to device buffers
		cudaMemcpy(r_d, r_h, n_r * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_r_d, distr_r_h, n_distr_r * sizeof(float_t), cudaMemcpyHostToDevice);

		run_init(rot_h, transvec);

		//for(int cby = 2; cby < 129; cby += 2) {
		//for(int cbz = 2; cbz < 129; cbz += 2) {
		// decompose computations and construct and call the kernel
		// decomposing along y and z
		// note that (cuda x y z != dim x y z)
		unsigned int cuda_block_y = 16, cuda_block_z = 6;
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
		form_factor_sphere_kernel <<< ff_grid_size, ff_block_size, shared_mem_size >>> (
				nqx_, nqy_, nqz_, qx_, qy_, qz_, rot_,
				n_r, r_d, n_distr_r, distr_r_d, transvec_, //n_transvec, transvec_d,
				ff_);

		cudaThreadSynchronize();
		cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess) {
			std::cerr << "error: form factor kernel failed [" << __FILE__ << ":" << __LINE__ << "]: "
						<< cudaGetErrorString(err) << std::endl;
		} else {
			//std::cout << "block size: " << cby << " x " << cbz << ". ";
			//std::cout << "Analytical Sphere Kernel completed in " << kernel_time << " ms." << std::endl;
			construct_output_ff(ff);
		} // if-else

		//std::cout << "GPU memory time: " << mem_time << " ms." << std::endl;
		//} // for cbz
		//} // for cby

		cudaFree(distr_r_d);
		cudaFree(r_d);

		return true;
	} // AnalyticFormFactorG::compute_sphere()


	// sphere gpu kernel
/*	__global__ void form_factor_sphere_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
					float_t* qx, float_t* qy, cucomplex_t* qz, float_t* rot,
					unsigned int n_r, float_t* r, unsigned int n_distr_r, float_t* distr_r,
					unsigned int n_transvec, float_t* transvec, cucomplex_t* ff) {
		// decomposition is along y and z
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;
		// compute
		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				unsigned int index = base_index + i_x;
				// computing mesh values on the fly instead of storing them
				cucomplex_t temp_mqx = make_cuC(qy[i_y] * rot[0] + qx[i_x] * rot[1] + qz[i_z].x * rot[2],
												qz[i_z].y * rot[2]);
				cucomplex_t temp_mqy = make_cuC(qy[i_y] * rot[3] + qx[i_x] * rot[4] + qz[i_z].x * rot[5],
												qz[i_z].y * rot[5]);
				cucomplex_t temp_mqz = make_cuC(qy[i_y] * rot[6] + qx[i_x] * rot[7] + qz[i_z].x * rot[8],
												qz[i_z].y * rot[8]);
				cucomplex_t q = cuCnorm3(temp_mqx, temp_mqy, temp_mqz);
				cucomplex_t temp_f = make_cuC((float_t)0.0, (float_t)0.0);
				for(unsigned int i_r = 0; i_r < n_r; ++ i_r) {
					float_t temp_r = r[i_r];
					ff_sphere_kernel_compute_tff(temp_r, distr_r[i_r], q, temp_mqz, temp_f);

				} // for i_r
				ff[index] = ff_sphere_kernel_compute_ff(temp_f,	temp_mqx, temp_mqy, temp_mqz,
													transvec[0], transvec[1], transvec[2]);
			} // for x
		} // if
	} // form_factor_sphere_kernel() 
*/
	extern __shared__ float_t dynamic_shared[];

	// sphere gpu kernel
	__global__ void form_factor_sphere_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
					float_t* qx, float_t* qy, cucomplex_t* qz, float_t* rot,
					unsigned int n_r, float_t* r, unsigned int n_distr_r, float_t* distr_r,
					/*unsigned int n_transvec,*/ float_t* transvec, cucomplex_t* ff) {
		// decomposition is along y and z
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;

		// shared buffers
		//float_t* qx_s = (float_t*) dynamic_shared;
		//float_t* qy_s = (float_t*) &qx_s[nqx];
		//cucomplex_t* qz_s = (cucomplex_t*) &qy_s[blockDim.x];
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

		__syncthreads();

		// compute
		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				// computing mesh values on the fly instead of storing them
				cucomplex_t mqx, mqy, mqz;
				compute_meshpoints(qx_s[i_x], qy_s[threadIdx.x], qz_s[threadIdx.y], rot, mqx, mqy, mqz);
				cucomplex_t q = cuCnorm3(mqx, mqy, mqz);
				cucomplex_t temp_f = make_cuC((float_t)0.0, (float_t)0.0);
				for(unsigned int i_r = 0; i_r < n_r; ++ i_r) {
					float_t temp4 = distr_r[i_r] * 4 * PI_ * pow(r[i_r], 3);
					if(cuCiszero(q)) {
						temp_f = temp_f + temp4 / (float_t) 3.0;
					} else {
						cucomplex_t temp1 = q * r[i_r];
						cucomplex_t temp2 = cuCsin(temp1) - temp1 * cuCcos(temp1);
						cucomplex_t temp3 = temp1 * temp1 * temp1;
						temp_f = temp_f + temp4 * (temp2 / temp3) * cuCexpi(mqz * r[i_r]);
					} // if-else
				} // for i_r
				ff[base_index + i_x] = ff_sphere_kernel_compute_ff(temp_f,	mqx, mqy, mqz,
													transvec[0], transvec[1], transvec[2]);
			} // for x
		} // if
	} // form_factor_sphere_kernel()


	__device__ void ff_sphere_kernel_compute_tff(float r, float distr_r,
										cuFloatComplex q, cuFloatComplex mqz, cuFloatComplex &f) {
		cuFloatComplex temp1 = make_cuFloatComplex(q.x * r, q.y * r);
		cuFloatComplex temp2 = cuCsubf(cuCsin(temp1), cuCmulf(temp1, cuCcos(temp1)));
		cuFloatComplex temp3 = cuCmulf(temp1, cuCmulf(temp1, temp1));
		float temp4 = distr_r * 4 * PI_ * r * r * r;
		cuFloatComplex temp5 = cuCmulf(cuCdivf(temp2, temp3),
									cuCexp(make_cuFloatComplex(-r * mqz.y, r * mqz.x)));
		cuFloatComplex tempff = make_cuFloatComplex(temp4 * temp5.x, temp4 * temp5.y);
		f = cuCaddf(f, tempff);
	} // ff_sphere_kernel_compute_tff()


	__device__ void ff_sphere_kernel_compute_tff(double r, double distr_r,
										cuDoubleComplex q, cuDoubleComplex mqz, cuDoubleComplex &f) {
		cuDoubleComplex temp1 = cuCmul(q, make_cuDoubleComplex(r, 0.0));
		cuDoubleComplex temp2 = cuCsub(cuCsin(temp1), cuCmul(temp1, cuCcos(temp1)));
		cuDoubleComplex temp3 = cuCmul(temp1, cuCmul(temp1, temp1));
		double temp4 = distr_r * 4 * PI_ * r * r * r;
		cuDoubleComplex temp5 = cuCmul(cuCdiv(temp2, temp3),
									cuCexp(make_cuDoubleComplex(-r * mqz.y, r * mqz.x)));
		cuDoubleComplex tempff = make_cuDoubleComplex(temp4 * temp5.x, temp4 * temp5.y);
		f = cuCadd(f, tempff);
	} // ff_sphere_kernel_compute_tff()


	__device__ cuFloatComplex ff_sphere_kernel_compute_ff(cuFloatComplex temp_f,
								cuFloatComplex mqx, cuFloatComplex mqy, cuFloatComplex mqz,
								float tx, float ty, float tz) {
			float rl = tx * mqx.x + ty * mqy.x + tz * mqz.x;
			float im = tx * mqx.y + ty * mqy.y + tz * mqz.y;
			cuFloatComplex temp1 = cuCexp(make_cuFloatComplex(-im, rl));
			return cuCmulf(temp_f, temp1);
	} // ff_sphere_kernel_compute_ff()


	__device__ cuDoubleComplex ff_sphere_kernel_compute_ff(cuDoubleComplex temp_f,
								cuDoubleComplex mqx, cuDoubleComplex mqy, cuDoubleComplex mqz,
								double tx, double ty, double tz) {
			double rl = tx * mqx.x + ty * mqy.x + tz * mqz.x;
			double im = tx * mqx.y + ty * mqy.y + tz * mqz.y;
			cuDoubleComplex temp1 = cuCexp(make_cuDoubleComplex(-im, rl));
			return cuCmul(temp_f, temp1);
	} // ff_sphere_kernel_compute_ff()


} // namespace hig

