/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_gpu.cu
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

  /** Form Factor of Pyramid:
   *  ff = It is complicated ...
   */
  __device__  __inline__ cucomplex_t FormFactorPyramid(cucomplex_t qx, cucomplex_t qy, cucomplex_t qz, 
          real_t length, real_t width, real_t height, real_t angle){

    real_t a = angle * PI_ / 180.;
    real_t tan_a = cuRtan(a);
    if ((2*height/length >= tan_a) || (2*height/width >= tan_a)) 
      return make_cuC(REAL_ZERO_, REAL_ZERO_);

    cucomplex_t tmp = (qx * qy);
    if (cuC_abs(tmp) < 1.0E-20)
      return make_cuC(REAL_ZERO_, REAL_ZERO_);

    // define complex units
    const cucomplex_t P_J = make_cuC(REAL_ZERO_, REAL_ONE_);
    const cucomplex_t N_J = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);

    // compute q[1-4]
    cucomplex_t q1 = 0.5 * ((qx - qy)/tan_a + qz);
    cucomplex_t q2 = 0.5 * ((qx - qy)/tan_a - qz);
    cucomplex_t q3 = 0.5 * ((qx + qy)/tan_a + qz);
    cucomplex_t q4 = 0.5 * ((qx + qy)/tan_a - qz);

    // compute K[1-4]
    cucomplex_t K1 = cuCsinc(q1 * height) * cuCexp(N_J * q1 * height) + 
        cuCsinc(q2 * height) * cuCexp(P_J * q2 * height);
    cucomplex_t K2 = cuCsinc(q1 * height) * cuCexp(N_J * q1 * height) * N_J + 
        cuCsinc(q2 * height) * cuCexp(P_J * q2 * height) * P_J;
    cucomplex_t K3 = cuCsinc(q3 * height) * cuCexp(N_J * q3 * height) + 
        cuCsinc(q4 * height) * cuCexp(P_J * q4 * height);
    cucomplex_t K4 = cuCsinc(q3 * height) * cuCexp(N_J * q3 * height) * N_J + 
        cuCsinc(q4 * height) * cuCexp(P_J * q4 * height) * P_J;

    cucomplex_t t1 = K1 * cuCcos(0.5 * (qx * length - qy * width));
    cucomplex_t t2 = K2 * cuCsin(0.5 * (qx * length - qy * width));
    cucomplex_t t3 = K3 * cuCcos(0.5 * (qx * length + qy * width));
    cucomplex_t t4 = K4 * cuCsin(0.5 * (qx * length + qy * width));
    
    return ((height/tmp) * (t1 + t2 - t3 - t4));
  }
 
  __global__ void ff_pyramid_kernel (unsigned int nqy, unsigned int nqz, 
          real_t * qx, real_t * qy, cucomplex_t * qz, cucomplex_t * ff,
          RotMatrix_t rot,
          int nx, real_t * x, real_t * distr_x,
          int ny, real_t * y, real_t * distr_y,
          int nz, real_t * z, real_t * distr_z,
          int na, real_t * a, real_t * distr_a) {
    int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if (i_z < nqz){
      int i_y = i_z % nqy;
      cucomplex_t c_neg_unit = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC(REAL_ZERO_, REAL_ZERO_);
      for (int i = 0; i < nx; i++){
        for (int j = 0; j < ny; j++){
          for (int k = 0; k < nz; k++){
            for (int l = 0; l < na; l++){
              real_t wght = distr_x[i] * distr_y[j] * distr_z[k] * distr_a[l];
              temp_ff = temp_ff + FormFactorPyramid(mqx, mqy, mqz, x[i], y[j], z[k], a[l]) * wght;
            }
          }
        }
      }
      cucomplex_t temp1 = transvec_d[0] * mqx + transvec_d[1] * mqy + transvec_d[2] * mqz;
      ff[i_z] =  temp_ff * cuCexpi(temp1);
    }
  } // ff_pyramid_kernel()

   
  bool AnalyticFormFactorG::compute_pyramid(const real_t tau, const real_t eta,
                  const std::vector<real_t>& x,
                  const std::vector<real_t>& distr_x,
                  const std::vector<real_t>& y,
                  const std::vector<real_t>& distr_y,
                  const std::vector<real_t>& z,
                  const std::vector<real_t>& distr_z,
                  const std::vector<real_t>& a,
                  const std::vector<real_t>& distr_a,
                  const RotMatrix_t & rot, const std::vector<real_t>& transvec,
                  std::vector<complex_t>& ff) {
    unsigned int n_x = x.size(), n_distr_x = distr_x.size();
    unsigned int n_y = y.size(), n_distr_y = distr_y.size();
    unsigned int n_z = z.size(), n_distr_z = distr_z.size();
    unsigned int n_a = a.size(), n_distr_a = distr_a.size();
    const real_t *x_h = x.empty() ? NULL : &*x.begin();
    const real_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
    const real_t *y_h = y.empty() ? NULL : &*y.begin();
    const real_t *distr_y_h = distr_y.empty() ? NULL : &*distr_y.begin();
    const real_t *z_h = z.empty() ? NULL : &*z.begin();
    const real_t *distr_z_h = distr_z.empty() ? NULL : &*distr_z.begin();
    const real_t *a_h = a.empty() ? NULL : &*a.begin();
    const real_t *distr_a_h = distr_a.empty() ? NULL : &*distr_a.begin();
    real_t transvec_h[3] = {transvec[0], transvec[1], transvec[2]};

    // construct device buffers
    real_t *x_d, *distr_x_d;
    real_t *y_d, *distr_y_d;
    real_t *z_d, *distr_z_d;
    real_t *a_d, *distr_a_d;

    cudaMalloc((void**) &x_d, n_x * sizeof(real_t));
    cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(real_t));
    cudaMalloc((void**) &y_d, n_y * sizeof(real_t));
    cudaMalloc((void**) &distr_y_d, n_distr_y * sizeof(real_t));
    cudaMalloc((void**) &z_d, n_z * sizeof(real_t));
    cudaMalloc((void**) &distr_z_d, n_distr_z * sizeof(real_t));
    cudaMalloc((void**) &a_d, n_a * sizeof(real_t));
    cudaMalloc((void**) &distr_a_d, n_distr_a * sizeof(real_t));

    // copy data to device buffers
    cudaMemcpy(x_d, x_h, n_x * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y_h, n_y * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(z_d, z_h, n_z * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(a_d, a_h, n_a * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_y_d, distr_y_h, n_distr_y * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_z_d, distr_z_h, n_distr_z * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_a_d, distr_a_h, n_distr_a * sizeof(real_t), cudaMemcpyHostToDevice);

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
    ff_pyramid_kernel <<<num_blocks, num_threads >>> (nqy_, nqz_, 
            qx_, qy_, qz_, ff_, rot,
            n_x, x_d, distr_x_d, 
            n_y, y_d, distr_y_d,
            n_z, z_d, distr_z_d,
            n_a, a_d, distr_a_d);
    
    construct_output_ff(ff);

    cudaFree(distr_a_d);
    cudaFree(a_d);
    cudaFree(distr_z_d);
    cudaFree(z_d);
    cudaFree(distr_y_d);
    cudaFree(y_d);
    cudaFree(distr_x_d);
    cudaFree(x_d);
   
    return true;
  } // AnalyticFormFactorG::compute_pyramid()

} // namespace hig
