/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_prism6_gpu.cu
 *  Created: Oct 16, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>
#include <complex>
#include <cuComplex.h>
#include <stdio.h>

#include <ff/gpu/ff_ana_gpu.cuh>
#include <common/enums.hpp>
#include <common/constants.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>
#include <utils/gpu/cu_utilities.cuh>


namespace hig {

  extern __constant__ real_t tau_d;
  extern __constant__ real_t eta_d;
  extern __constant__ real_t transvec_d[3];

  /** Form Factor of Prism6:
   *  ff = It is complicated ...
   */
  __device__  __inline__ cucomplex_t FormFactorPrism6(cucomplex_t qx, cucomplex_t qy, cucomplex_t qz, 
          real_t length, real_t height){

    cucomplex_t tmp = 3. * qy * qy - qx * qx;
    if (cuC_abs(tmp) < 1.0E-20)
        return make_cuC(REAL_ZERO_, REAL_ZERO_);

    // define complex units
    const cucomplex_t P_J = make_cuC(REAL_ONE_, REAL_ZERO_);
    const cucomplex_t N_J = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
    const real_t sqrt3 = 1.732050808;

    real_t L = 0.5 * length;
    real_t H = 0.5 * height;

    cucomplex_t t1 = 4. * sqrt3 / tmp;
    cucomplex_t t2 = qy * qy * L * L * cuCsinc(qx * L / sqrt3) * cuCsinc(qy * L);
    cucomplex_t t3 = cuCcos(2 * qx * L / sqrt3);
    cucomplex_t t4 = cuCcos(qy * L) * cuCcos(qx * L / sqrt3);
    cucomplex_t t5 = cuCsinc(qz * H) * cuCexp(P_J * qz * H);
    return (t1 * (t2 + t3 - t4) * t5);
  }
 
  __global__ void ff_prism6_kernel (unsigned int nqy, unsigned int nqz, 
          real_t * qx, real_t * qy, cucomplex_t * qz, cucomplex_t * ff,
          RotMatrix_t rot,
          int nx, real_t * x, real_t * distr_x,
          int ny, real_t * y, real_t * distr_y) {

    int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if (i_z < nqz){
      int i_y = i_z % nqy;
      cucomplex_t c_neg_unit = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC(REAL_ZERO_, REAL_ZERO_);
      for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
          real_t wght = distr_x[i] * distr_y[j];
          temp_ff = temp_ff + FormFactorPrism6(mqx, mqy, mqz, x[i], y[j]) * wght;
        }
      }
      cucomplex_t temp1 = transvec_d[0] * mqx + transvec_d[1] * mqy + transvec_d[2] * mqz;
      ff[i_z] =  temp_ff * cuCexpi(temp1);
    }
  } // ff_prism6_kernel()

   
  bool AnalyticFormFactorG::compute_prism6(const real_t tau, const real_t eta,
                  const std::vector<real_t>& x,
                  const std::vector<real_t>& distr_x,
                  const std::vector<real_t>& y,
                  const std::vector<real_t>& distr_y,
                  const RotMatrix_t & rot, const std::vector<real_t>& transvec,
                  std::vector<complex_t>& ff) {
    unsigned int n_x = x.size(), n_distr_x = distr_x.size();
    unsigned int n_y = y.size(), n_distr_y = distr_y.size();

    const real_t *x_h = x.empty() ? NULL : &*x.begin();
    const real_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
    const real_t *y_h = y.empty() ? NULL : &*y.begin();
    const real_t *distr_y_h = distr_y.empty() ? NULL : &*distr_y.begin();
    real_t transvec_h[3] = {transvec[0], transvec[1], transvec[2]};

    // construct device buffers
    real_t *x_d, *distr_x_d;
    real_t *y_d, *distr_y_d;

    cudaMalloc((void**) &x_d, n_x * sizeof(real_t));
    cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(real_t));
    cudaMalloc((void**) &y_d, n_y * sizeof(real_t));
    cudaMalloc((void**) &distr_y_d, n_distr_y * sizeof(real_t));

    // copy data to device buffers
    cudaMemcpy(x_d, x_h, n_x * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y_h, n_y * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_y_d, distr_y_h, n_distr_y * sizeof(real_t), cudaMemcpyHostToDevice);

    //run_init(rot_h, transvec);
    cudaMemcpyToSymbol(tau_d, &tau, sizeof(real_t), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(eta_d, &eta, sizeof(real_t), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol(transvec_d, transvec_h, 3*sizeof(real_t), 0, cudaMemcpyHostToDevice); 

    int num_threads = 256;
    int num_blocks =  nqz_ / num_threads + 1;
    dim3 ff_grid_size(num_blocks, 1, 1);
    dim3 ff_block_size(num_threads, 1, 1);
    //std::cerr << "Q-Grid size = " << nqz_ << std::endl;

    // the kernel
    ff_prism6_kernel <<<num_blocks, num_threads >>> (nqy_, nqz_, 
            qx_, qy_, qz_, ff_, rot,
            n_x, x_d, distr_x_d, 
            n_y, y_d, distr_y_d);
    
    construct_output_ff(ff);

    cudaFree(distr_y_d);
    cudaFree(y_d);
    cudaFree(distr_x_d);
    cudaFree(x_d);
   
    return true;
  } // AnalyticFormFactorG::compute_prism6()

} // namespace hig
