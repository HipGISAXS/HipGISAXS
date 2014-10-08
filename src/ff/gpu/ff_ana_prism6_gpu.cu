/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_prism6_gpu.cu
 *  Created: Oct 16, 2012
 *  Modified: Wed 08 Oct 2014 12:17:48 PM PDT
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

#include <ff/gpu/ff_ana_gpu.cuh>
#include <common/enums.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>
#include <utils/gpu/cu_utilities.cuh>


namespace hig {

  /**
   * prism6 on gpu
   */

  __global__ void form_factor_prism6_kernel(unsigned int, unsigned int, unsigned int,
                  float_t*, float_t*, cucomplex_t*, float_t, float_t, float_t*,
                  unsigned int, float_t*, unsigned int, float_t*,
                  unsigned int, float_t*, unsigned int, float_t*,
                  float_t*, cucomplex_t*);


  bool AnalyticFormFactorG::compute_prism6(const float_t tau, const float_t eta,
                  const std::vector<float_t>& l,
                  const std::vector<float_t>& distr_l,
                  const std::vector<float_t>& h,
                  const std::vector<float_t>& distr_h,
                  const float_t* rot_h, const std::vector<float_t>& transvec,
                  std::vector<complex_t>& ff) {
    unsigned int n_l = l.size(), n_distr_l = distr_l.size();
    unsigned int n_h = h.size(), n_distr_h = distr_h.size();
    const float_t *l_h = l.empty() ? NULL : &*l.begin();
    const float_t *distr_l_h = distr_l.empty() ? NULL : &*distr_l.begin();
    const float_t *h_h = h.empty() ? NULL : &*h.begin();
    const float_t *distr_h_h = distr_h.empty() ? NULL : &*distr_h.begin();

    // construct device buffers
    float_t *l_d, *distr_l_d;
    float_t *h_d, *distr_h_d;

    cudaMalloc((void**) &l_d, n_l * sizeof(float_t));
    cudaMalloc((void**) &distr_l_d, n_distr_l * sizeof(float_t));
    cudaMalloc((void**) &h_d, n_h * sizeof(float_t));
    cudaMalloc((void**) &distr_h_d, n_distr_h * sizeof(float_t));

    // copy data to device buffers
    cudaMemcpy(l_d, l_h, n_l * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(h_d, h_h, n_h * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_l_d, distr_l_h, n_distr_l * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_h_d, distr_h_h, n_distr_h * sizeof(float_t), cudaMemcpyHostToDevice);

    run_init(rot_h, transvec);

    unsigned int cuda_block_y = 16, cuda_block_z = 8;
    unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) nqy_ / cuda_block_y);
    unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) nqz_ / cuda_block_z);
    dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
    dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

    // the kernel
    form_factor_prism6_kernel <<< ff_grid_size, ff_block_size >>> (
        nqx_, nqy_, nqz_, qx_, qy_, qz_, tau, eta, rot_,
        n_l, l_d, n_distr_l, distr_l_d,
        n_h, h_d, n_distr_h, distr_h_d,
        transvec_,
        ff_);

    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "error: prism6 form factor kernel failed [" << __FILE__ << ":" << __LINE__ << "]: "
            << cudaGetErrorString(err) << std::endl;
      return false;
    } else {
      //std::cout << "block size: " << cby << " x " << cbz << ". ";
      construct_output_ff(ff);
    } // if-else

    cudaFree(distr_h_d);
    cudaFree(h_d);
    cudaFree(distr_l_d);
    cudaFree(l_d);

    return true;
  } // AnalyticFormFactorG::compute_box()


  __global__ void form_factor_prism6_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                  float_t *qx, float_t *qy, cucomplex_t *qz,
                  float_t tau, float_t eta, float_t *rot,
                  unsigned int n_l, float_t *l, unsigned int n_distr_l, float_t *distr_l,
                  unsigned int n_h, float_t *h, unsigned int n_distr_h, float_t *distr_h,
                  float_t *transvec, cucomplex_t *ff) {
    unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
    unsigned int base_index = nqx * nqy * i_z + nqx * i_y;
    float_t sqrt3 = sqrt(3.0);
    if(i_y < nqy && i_z < nqz) {
      for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
        cucomplex_t mqx, mqy, mqz;
        compute_meshpoints(qx[i_x], qy[i_y], qz[i_z], rot, mqx, mqy, mqz);
        cucomplex_t qm = tan(tau) * (sin(eta) * mqx + cos(eta) * mqy);
        cucomplex_t temp1 = (4.0 * sqrt3) / (3.0 * mqy * mqy - mqx * mqx);
        cucomplex_t temp_ff = make_cuC((float_t) 0.0, (float_t) 0.0);
        for(unsigned int p_h = 0; p_h < n_h; ++ p_h) {
          for(unsigned int p_l = 0; p_l < n_l; ++ p_l) {
            cucomplex_t rmqx = l[p_l] * mqx / sqrt3;
            cucomplex_t rmqy = l[p_l] * mqy;
            cucomplex_t temp2 = rmqy * rmqy * cuCsinc(rmqx) * cuCsinc(rmqy);
            cucomplex_t temp3 = cuCcos(2.0 * rmqx);
            cucomplex_t temp4 = cuCcos(rmqy) * cuCcos(rmqx);
            cucomplex_t temp5 = temp1 * (temp2 + temp3 - temp4);
            cucomplex_t temp6 = fq_inv(mqz + qm, h[p_h]);
            temp_ff = temp_ff + temp5 * temp6;
          } // for x
        } // for y
        cucomplex_t temp_e = cuCexpi(mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
        unsigned int curr_index = base_index + i_x;
        ff[curr_index] = temp_ff * temp_e;
      } // for x
    } // if
  } // form_factor_box_kernel()


} // namespace hig

