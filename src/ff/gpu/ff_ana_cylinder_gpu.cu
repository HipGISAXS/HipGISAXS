/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_cylinder_gpu.cu
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

  /** Form Factor of Cylinder:
   *  R : (real) Length of the cylinder
   *  H : (real) Width of the cylinder
   *  q : (complex) q-vector
   *  ff = 2 * pi * R^2 * H * J1(qpar * R)/(qpar * R) * sinc(qz * H) exp(j * qz * H/2) 
   */
  __device__  __inline__ cucomplex_t FormFactorCylinder (
          cucomplex_t qx, cucomplex_t qy, cucomplex_t qz, 
          real_t radius, real_t height){
    cucomplex_t qpar = cuCsqrt(qx * qx + qy * qy);

    real_t     temp1 = 2 * PI_ * radius * radius * height;
    cucomplex_t temp2; temp2.x = 0.0; temp2.y = 0.0;
    if((qpar.x*qpar.x+qpar.y*qpar.y)>CUTINY_) temp2 = cuCcbessj(qpar * radius, 1) / (qpar * radius);
    cucomplex_t temp3 = cuCsinc (0.5 * qz * height);
    cucomplex_t temp4 = cuCexpi(0.5 * qz * height);
    return (temp1 * temp2 * temp3 * temp4);
  }
 
  __global__ void ff_cylinder_kernel (unsigned int nqy, unsigned int nqz, 
          real_t * qx, real_t * qy, cucomplex_t * qz, cucomplex_t * ff,
          RotMatrix_t rot,
          int nr, real_t * r, real_t * distr_r,
          int nh, real_t * h, real_t * distr_h) {
    int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if (i_z < nqz){
      int i_y = i_z % nqy;
      cucomplex_t c_neg_unit = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC(REAL_ZERO_, REAL_ZERO_);
      for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nh; j++) {
          temp_ff = temp_ff + FormFactorCylinder(mqx, mqy, mqz, r[i], h[j]);
        }
      }
      cucomplex_t temp1 = transvec_d[0] * mqx + transvec_d[1] * mqy + transvec_d[2] * mqz;
      ff[i_z] =  temp_ff * cuCexpi(temp1);
    }
  } // ff_cylinder_kernel()

   
  bool AnalyticFormFactorG::compute_cylinder(const real_t tau, const real_t eta,
                  const std::vector<real_t>& r,
                  const std::vector<real_t>& distr_r,
                  const std::vector<real_t>& h,
                  const std::vector<real_t>& distr_h,
                  const RotMatrix_t & rot, const std::vector<real_t>& transvec,
                  std::vector<complex_t>& ff) {
    unsigned int n_r = r.size(), n_distr_r = distr_r.size();
    unsigned int n_h = h.size(), n_distr_h = distr_h.size();

    const real_t *r_h = r.empty() ? NULL : &*r.begin();
    const real_t *distr_r_h = distr_r.empty() ? NULL : &*distr_r.begin();
    const real_t *h_h = h.empty() ? NULL : &*h.begin();
    const real_t *distr_h_h = distr_h.empty() ? NULL : &*distr_h.begin();
    real_t transvec_h[3] = {transvec[0], transvec[1], transvec[2]};

    // construct device buffers
    real_t *r_d, *distr_r_d;
    real_t *h_d, *distr_h_d;

    cudaMalloc((void**) &r_d, n_r * sizeof(real_t));
    cudaMalloc((void**) &distr_r_d, n_distr_r * sizeof(real_t));
    cudaMalloc((void**) &h_d, n_h * sizeof(real_t));
    cudaMalloc((void**) &distr_h_d, n_distr_h * sizeof(real_t));

    // copy data to device buffers
    cudaMemcpy(r_d, r_h, n_r * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(h_d, h_h, n_h * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_r_d, distr_r_h, n_distr_r * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_h_d, distr_h_h, n_distr_h * sizeof(real_t), cudaMemcpyHostToDevice);

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
    ff_cylinder_kernel <<<num_blocks, num_threads >>> (nqy_, nqz_, 
            qx_, qy_, qz_, ff_, rot,
            n_r, r_d, distr_r_d, 
            n_h, h_d, distr_h_d);
    
    construct_output_ff(ff);

    cudaFree(distr_h_d);
    cudaFree(h_d);
    cudaFree(distr_r_d);
    cudaFree(r_d);
   
    return true;
  } // AnalyticFormFactorG::compute_cylinder()

} // namespace hig
