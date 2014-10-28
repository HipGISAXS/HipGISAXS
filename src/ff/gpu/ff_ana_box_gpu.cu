/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_gpu.cu
 *  Created: Oct 16, 2012
 *  Modified: Wed 08 Oct 2014 12:17:47 PM PDT
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

  __constant__ float_t d_transvec[3];
  __constant__ float_t d_rot[9];
  __constant__ float_t d_tau;
  __constant__ float_t d_eta;

  /** Form Factor of Box:
   *  L : (real) Length of the box 
   *  W : (real) Width of the box
   *  H : (real) Height of the box
   *  q : (complex) q-vector
   *  ff = L * W * H exp(j * qz * H/2) * sinc(qx * L/2) * sinc (qy * W/2) * sinc(qz * H/2)
   */
  __device__  __inline__ cucomplex_t FormFactorBox (float_t L, float_t W, float_t H, 
      cucomplex_t qx, cucomplex_t qy, cucomplex_t qz)
  {
    cucomplex_t temp1, temp2, temp3, temp4;
    temp1 = cuCsinc (0.5 * qx * L);
    temp2 = cuCsinc (0.5 * qy * W);
    temp3 = cuCsinc (0.5 * qz * H);
    temp4 = L * W * H * cuCexpi(0.5 * qz * H);
    return (temp1 * temp2 * temp3 * temp4);
  }
 

  __global__ void form_factor_box_kernel (unsigned int nqy, unsigned int nqz,
                  float_t *qx, float_t *qy, cucomplex_t *qz,
                  unsigned int n_x, float_t *x, unsigned int n_distr_x, float_t *distr_x,
                  unsigned int n_y, float_t *y, unsigned int n_distr_y, float_t *distr_y,
                  unsigned int n_z, float_t *z, unsigned int n_distr_z, float_t *distr_z,
                  cucomplex_t *ff) {
    unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int idx = threadIdx.x;
    if ( i_z > nqz ) {
      unsigned int i_y = i_z % nqy;

      /*
      // shared buffers:
      __shared__ cucomplex_t qz_s[blockDim.x];
      __shared__ float_t qx_s[blockDim.x];
      __shared__ float_t qy_s[blockDim.x];

      // load all q-vectors
      qx_s[threadIdx.x] = qx[i_y];
      qy_s[threadIdx.x] = qy[i_y];
      qz_s[threadIdx.x] = qz[i_z];
      */

      // TODO: also put x, y, z, distr_x, distr_y, distr_z in shared mem ...
      // TODO: also put transvec in shared mem ...

      // make sure everything is in place
      __syncthreads();

      cucomplex_t mqx, mqy, mqz;
      //compute_meshpoints(qx_s[idx], qy_s[idx], qz_s[idx], d_rot, mqx, mqy, mqz);
      compute_meshpoints(qx[i_y], qy[i_y], qz[i_z], d_rot, mqx, mqy, mqz);
      cucomplex_t temp_ff = make_cuC ((float_t) 0, (float_t) 0);
      for(unsigned int p_z = 0; p_z < n_z; ++ p_z) {
          for(unsigned int p_y = 0; p_y < n_y; ++ p_y) {
            for(unsigned int p_x = 0; p_x < n_x; ++ p_x) {
              float_t w = distr_x[p_x] * distr_y[p_y] * distr_z[p_z];
              cucomplex_t tff = FormFactorBox (x[p_x], y[p_y], z[p_z], mqx, mqy, mqz);
              temp_ff = temp_ff +  tff * w;
            } // for x
          } // for y
        } // for z
        ff[i_z] = temp_ff * cuCexpi(mqx * d_transvec[0] + mqy * d_transvec[1] + mqz * d_transvec[2]);
    } // if
  } // form_factor_box_kernel()

   
  bool AnalyticFormFactorG::compute_box(const float_t tau, const float_t eta,
                  const std::vector<float_t>& x,
                  const std::vector<float_t>& distr_x,
                  const std::vector<float_t>& y,
                  const std::vector<float_t>& distr_y,
                  const std::vector<float_t>& z,
                  const std::vector<float_t>& distr_z,
                  const float_t* rot_h, const std::vector<float_t>& transvec,
                  std::vector<complex_t>& ff) {
    unsigned int n_x = x.size(), n_distr_x = distr_x.size();
    unsigned int n_y = y.size(), n_distr_y = distr_y.size();
    unsigned int n_z = z.size(), n_distr_z = distr_z.size();
    const float_t *x_h = x.empty() ? NULL : &*x.begin();
    const float_t *distr_x_h = distr_x.empty() ? NULL : &*distr_x.begin();
    const float_t *y_h = y.empty() ? NULL : &*y.begin();
    const float_t *distr_y_h = distr_y.empty() ? NULL : &*distr_y.begin();
    const float_t *z_h = z.empty() ? NULL : &*z.begin();
    const float_t *distr_z_h = distr_z.empty() ? NULL : &*distr_z.begin();

    // construct device buffers
    float_t *x_d, *distr_x_d;
    float_t *y_d, *distr_y_d;
    float_t *z_d, *distr_z_d;

    cudaMalloc((void**) &x_d, n_x * sizeof(float_t));
    cudaMalloc((void**) &distr_x_d, n_distr_x * sizeof(float_t));
    cudaMalloc((void**) &y_d, n_y * sizeof(float_t));
    cudaMalloc((void**) &distr_y_d, n_distr_y * sizeof(float_t));
    cudaMalloc((void**) &z_d, n_z * sizeof(float_t));
    cudaMalloc((void**) &distr_z_d, n_distr_z * sizeof(float_t));

    // copy data to device buffers
    cudaMemcpy(x_d, x_h, n_x * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y_h, n_y * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(z_d, z_h, n_z * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_x_d, distr_x_h, n_distr_x * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_y_d, distr_y_h, n_distr_y * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(distr_z_d, distr_z_h, n_distr_z * sizeof(float_t), cudaMemcpyHostToDevice);

    run_init(rot_h, transvec);
    cudaMemcpyToSymbol (d_transvec, transvec_, 3 * sizeof(float_t), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol (d_rot, rot_, 9 * sizeof(float_t), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol (d_tau, &tau, sizeof(float_t), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol (d_eta, &eta, sizeof(float_t), 0, cudaMemcpyHostToDevice);

    unsigned int cuda_block = 256;
    unsigned int cuda_num_blocks = (unsigned int) ceil((float_t) nqz_ / cuda_block);
    dim3 ff_grid_size(cuda_num_blocks, 1, 1);
    dim3 ff_block_size(cuda_block, 1, 1);

    size_t shared_mem_size = 4 * cuda_block * sizeof(float_t);
    if(shared_mem_size > 49152) { //??? should this limit come from hardware?
      std::cerr << "Too much shared memory requested!" << std::endl;
      return false;
    } // if

    // the kernel
    form_factor_box_kernel <<< ff_grid_size, ff_block_size >>> (
        nqy_, nqz_, qx_, qy_, qz_,
        n_x, x_d, n_distr_x, distr_x_d,
        n_y, y_d, n_distr_y, distr_y_d,
        n_z, z_d, n_distr_z, distr_z_d,
        ff_);

    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "error: box form factor kernel failed [" << __FILE__ << ":" << __LINE__ << "]: "
            << cudaGetErrorString(err) << std::endl;
      return false;
    } else {
      std::cout << "block size: " << 1 << " x " << 1 << ". ";
      construct_output_ff(ff);
    } // if-else

    cudaFree(distr_z_d);
    cudaFree(z_d);
    cudaFree(distr_y_d);
    cudaFree(y_d);
    cudaFree(distr_x_d);
    cudaFree(x_d);

    return true;
  } // AnalyticFormFactorG::compute_box()

} // namespace hig

