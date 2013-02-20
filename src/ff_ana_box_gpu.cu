/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Wed 20 Feb 2013 01:06:30 PM PST
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
									/*unsigned int,*/ float_t*, cucomplex_t*);


	bool AnalyticFormFactorG::compute_box(const float_t tau, const float_t eta,
										const std::vector<float_t>& x,
										const std::vector<float_t>& distr_x,
										const std::vector<float_t>& y,
										const std::vector<float_t>& distr_y,
										const std::vector<float_t>& z,
										const std::vector<float_t>& distr_z,
										/*const float_t* qx_h, const float_t* qy_h,
										const cucomplex_t* qz_h,*/
										const float_t* rot_h, const std::vector<float_t>& transvec,
										std::vector<complex_t>& ff) {
		unsigned int n_x = x.size(), n_distr_x = distr_x.size();
		unsigned int n_y = y.size(), n_distr_y = distr_y.size();
		unsigned int n_z = z.size(), n_distr_z = distr_z.size();
		//unsigned int n_transvec = transvec.size();
		const float_t *x_h = x.empty() ? NULL : &*x.begin();
		const float_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
		const float_t *y_h = y.empty() ? NULL : &*y.begin();
		const float_t *distr_y_h = distr_y.empty() ? NULL : &*distr_y.begin();
		const float_t *z_h = z.empty() ? NULL : &*z.begin();
		const float_t *distr_z_h = distr_z.empty() ? NULL : &*distr_z.begin();
		//const float_t *transvec_h = transvec.empty() ? NULL : &*transvec.begin();

		unsigned int grid_size = nqx_ * nqy_ * nqz_;

		// construct device buffers
		//float_t *qx_d, *qy_d;
		//cucomplex_t *qz_d, *ff_d;
		float_t *x_d, *distr_x_d;
		float_t *y_d, *distr_y_d;
		float_t *z_d, *distr_z_d;
		//float_t *transvec_d, *rot_d;

		//cudaMalloc((void**) &qx_d, nqx_ * sizeof(float_t));
		//cudaMalloc((void**) &qy_d, nqy_ * sizeof(float_t));
		//cudaMalloc((void**) &qz_d, nqz_ * sizeof(cucomplex_t));
		//cudaMalloc((void**) &ff_d, grid_size * sizeof(cucomplex_t));
		cudaMalloc((void**) &x_d, n_x * sizeof(float_t));
		cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(float_t));
		cudaMalloc((void**) &y_d, n_y * sizeof(float_t));
		cudaMalloc((void**) &distr_y_d, n_distr_y * sizeof(float_t));
		cudaMalloc((void**) &z_d, n_z * sizeof(float_t));
		cudaMalloc((void**) &distr_z_d, n_distr_z * sizeof(float_t));
		//cudaMalloc((void **) &transvec_d, n_transvec * sizeof(float_t));
		//cudaMalloc((void **) &rot_d, 9 * sizeof(float_t));

		// copy data to device buffers
		//cudaMemcpy(qx_d, qx_h, nqx_ * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(qy_d, qy_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(qz_d, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
		cudaMemcpy(x_d, x_h, n_x * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(y_d, y_h, n_y * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(z_d, z_h, n_z * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_y_d, distr_y_h, n_distr_y * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_z_d, distr_z_h, n_distr_z * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(transvec_d, transvec_h, n_transvec * sizeof(float_t), cudaMemcpyHostToDevice);
		//cudaMemcpy(rot_d, rot_h, 9 * sizeof(float_t), cudaMemcpyHostToDevice);

		run_init(rot_h, transvec);

		unsigned int cuda_block_y = 16, cuda_block_z = 8;
		unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) nqy_ / cuda_block_y);
		unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) nqz_ / cuda_block_z);
		dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
		dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

		// the kernel
		form_factor_box_kernel <<< ff_grid_size, ff_block_size >>> (
				nqx_, nqy_, nqz_, qx_, qy_, qz_, tau, eta, rot_,
				n_x, x_d, n_distr_x, distr_x_d,
				n_y, y_d, n_distr_y, distr_y_d,
				n_z, z_d, n_distr_z, distr_z_d,
				/*n_transvec,*/ transvec_,
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
		cudaFree(distr_z_d);
		cudaFree(z_d);
		cudaFree(distr_y_d);
		cudaFree(y_d);
		cudaFree(distr_x_d);
		cudaFree(x_d);
		//cudaFree(ff_d);
		//cudaFree(qz_d);
		//cudaFree(qy_d);
		//cudaFree(qx_d);

		return true;
	} // AnalyticFormFactorG::compute_box()


	__global__ void form_factor_box_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
									float_t *qx, float_t *qy, cucomplex_t *qz,
									float_t tau, float_t eta, float_t *rot,
									unsigned int n_x, float_t *x, unsigned int n_distr_x, float_t *distr_x,
									unsigned int n_y, float_t *y, unsigned int n_distr_y, float_t *distr_y,
									unsigned int n_z, float_t *z, unsigned int n_distr_z, float_t *distr_z,
									/*unsigned int n_transvec,*/ float_t *transvec, cucomplex_t *ff) {
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
				cucomplex_t tempa = sin(eta) * mqx;
				cucomplex_t tempb = cos(eta) * mqy;
				cucomplex_t temp_qm = tan(tau) * (tempa + tempb);
				cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
				for(unsigned int p_z = 0; p_z < n_z; ++ p_z) {
					for(unsigned int p_y = 0; p_y < n_y; ++ p_y) {
						for(unsigned int p_x = 0; p_x < n_x; ++ p_x) {
							cucomplex_t temp4 = fq_inv(mqz + temp_qm, y[p_y]);
							cucomplex_t temp8 = temp4 * sinc(mqy * z[p_z]) * sinc(mqx * x[p_x]);
							float_t temp9 = 4.0 * distr_x[p_x] * distr_y[p_y] * distr_z[p_z] *
											z[p_z] * x[p_x];
							temp_ff = temp_ff + temp9 * temp8;
						} // for x
					} // for y
				} // for z
				cucomplex_t temp_e = cuCexp(mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
				unsigned int curr_index = base_index + i_x;
				ff[curr_index] = temp_ff * temp_e;
			} // for x
		} // if
	} // form_factor_box_kernel()


} // namespace hig

