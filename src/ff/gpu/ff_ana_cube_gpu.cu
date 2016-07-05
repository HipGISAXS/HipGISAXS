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

  /** Form Factor of Cube:
   *  L : (real) Length of the cube 
   *  q : (complex) q-vector
   *  ff = L^3  exp(j * qz * L/2) * sinc(qx * L/2) * sinc (qy * L/2) * sinc(qz * L/2)
   */
  __device__  __inline__ cucomplex_t FormFactorCube(cucomplex_t qx, cucomplex_t qy, cucomplex_t qz, 
          real_t L){
    cucomplex_t temp1 = cuCsinc (0.5 * qx * L);
    cucomplex_t temp2 = cuCsinc (0.5 * qy * L);
    cucomplex_t temp3 = cuCsinc (0.5 * qz * L);
    cucomplex_t temp4 = L * L * L * cuCexpi(0.5 * qz * L);
    return (temp1 * temp2 * temp3 * temp4);
  }
 
  __global__ void ff_cube_kernel (unsigned int nqy, unsigned int nqz, 
          real_t * qx, real_t * qy, cucomplex_t * qz, cucomplex_t * ff,
          RotMatrix_t rot,
          int nx, real_t * x, real_t * distr_x) {
    int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if (i_z < nqz){
      int i_y = i_z % nqy;
      cucomplex_t c_neg_unit = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC(REAL_ZERO_, REAL_ZERO_);
      for (int i = 0; i < nx; i++){
        temp_ff = temp_ff + distr_x[i] * FormFactorCube(mqx, mqy, mqz, x[i]); 
      }
      cucomplex_t temp1 = transvec_d[0] * mqx + transvec_d[1] * mqy + transvec_d[2] * mqz;
      ff[i_z] =  temp_ff * cuCexpi(temp1);
    }
  } // ff_cube_kernel()

   
  bool AnalyticFormFactorG::compute_cube(const real_t tau, const real_t eta,
                  const std::vector<real_t>& x,
                  const std::vector<real_t>& distr_x,
                  const RotMatrix_t & rot, const std::vector<real_t>& transvec,
                  std::vector<complex_t>& ff) {
    unsigned int n_x = x.size(), n_distr_x = distr_x.size();
    const real_t *x_h = x.empty() ? NULL : &*x.begin();
    const real_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
    real_t transvec_h[3] = {transvec[0], transvec[1], transvec[2]};

    // construct device buffers
    real_t *x_d, *distr_x_d;

    cudaMalloc((void**) &x_d, n_x * sizeof(real_t));
    cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(real_t));

    // copy data to device buffers
    cudaMemcpy(x_d, x_h, n_x * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(real_t), cudaMemcpyHostToDevice);

    //run_init(rot_h, transvec);
    cudaMemcpyToSymbol(transvec_d, transvec_h, 3*sizeof(real_t), 0, cudaMemcpyHostToDevice); 

    int num_threads = 256;
    int num_blocks =  nqz_ / num_threads + 1;
    dim3 ff_grid_size(num_blocks, 1, 1);
    dim3 ff_block_size(num_threads, 1, 1);

    // the kernel
    ff_cube_kernel <<<num_blocks, num_threads >>> (nqy_, nqz_, 
            qx_, qy_, qz_, ff_, rot, n_x, x_d, distr_x_d);
    
    construct_output_ff(ff);

    cudaFree(distr_x_d);
    cudaFree(x_d);
   
    return true;
  } // AnalyticFormFactorG::compute_cube()

} // namespace hig
