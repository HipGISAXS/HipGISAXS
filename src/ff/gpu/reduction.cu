/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: reduction.cu
 *  Created: Aug 25, 2012
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
#include <omp.h>
#include <cuComplex.h>

#include <ff/gpu/reduction.cuh>
#include <common/parameters.hpp>

namespace hig {

  #ifdef GPUR

  /**
   * GPU: For single precision
   */

  /* Original reduction kernel; with no shared memory
   */
  __global__ void reduction_kernel(cuFloatComplex* fq_d,
            unsigned int curr_nqx, unsigned int curr_nqy,
            unsigned int curr_nqz, unsigned int curr_num_triangles,
            unsigned int b_nqx, unsigned int b_nqy,
            unsigned int b_nqz, unsigned int b_num_triangles,
            cuFloatComplex* ff_d) {
    // 3D block, where each thread is responsible for one point in the x, y, z grid
    unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
    unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;
    unsigned int z = threadIdx.z + blockDim.z * blockIdx.z;
    unsigned int i_ff = curr_nqx * curr_nqy * z + curr_nqx * y + x;
    unsigned int temp = curr_nqx * curr_nqy * curr_nqz;

    if(x < curr_nqx && y < curr_nqy && z < curr_nqz) {
      cuFloatComplex total = make_cuFloatComplex(0.0, 0.0);
      for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
        unsigned int i_fq = temp * i_t + i_ff;
        total = cuCaddf(total, fq_d[i_fq]);
      } // for
      ff_d[i_ff] = cuCmulf(total, make_cuFloatComplex(0.0, -1.0));
    } // if
  } // reduction_kernel()

  #ifdef REDUCTION2
  /* GPU parallel reduction:
   * perform 3D decomposition along t, y, z, and reduce along t in parallel
   */
  __global__ void reduction_kernel_parallel(cuFloatComplex* fq_d,
            unsigned int curr_nqx, unsigned int curr_nqy,
            unsigned int curr_nqz, unsigned int curr_num_triangles,
            unsigned int b_nqx, unsigned int b_nqy,
            unsigned int b_nqz, unsigned int b_num_triangles,
            cuFloatComplex* ff_d) {
    // b_num_triangles should be a power of 2
    // there will be only one block in the t direction
    // => blockDim.x = curr_num_triangles / 2 = h_t / 2
    // => hyperblock's t dimension size should be power of 2
    unsigned int i_tt = threadIdx.x + blockDim.x * blockIdx.x; // == threadIdx.x
    unsigned int i_yy = threadIdx.y + blockDim.y * blockIdx.y;
    unsigned int i_zz = threadIdx.z + blockDim.z * blockIdx.z;
    unsigned int i_ff_base = curr_nqx * (curr_nqy * i_zz + i_yy);
    unsigned int temp = curr_nqx * curr_nqy * curr_nqz;
    unsigned int i_fq_base = temp * i_tt + i_ff_base;
    unsigned int stempyz = blockDim.y * blockDim.z;
    unsigned int stemp = curr_nqx * stempyz;
    //unsigned int i_sbase = stemp * threadIdx.x + curr_nqx * (blockDim.y * threadIdx.z + threadIdx.y);

    // blockDim.x == BLOCK_REDUCTION_T_ / 2 == b_num_triangles / 2
    // blockDim.y == BLOCK_REDUCTION_Y_
    // blockDim.z == BLOCK_REDUCTION_Z_
    __shared__ cuFloatComplex shared_fq[BLOCK_REDUCTION_T_ * BLOCK_REDUCTION_Y_ * BLOCK_REDUCTION_Z_];

    unsigned int red_steps = __ffs(b_num_triangles) - 1;
    unsigned int i_sfq = stempyz * threadIdx.x + blockDim.y * threadIdx.z + threadIdx.y;
    for(unsigned int i_xx = 0; i_xx < curr_nqx; ++ i_xx) {
      if(i_tt < curr_num_triangles && i_yy < curr_nqy && i_zz < curr_nqz) {
        /* // without shared memory
        unsigned int curr_num_threads = blockDim.x; // == b_num_triangles >> 1
        unsigned int i_fq = i_fq_base + i_xx;
        for(unsigned int step = 0; step < red_steps; ++ step) {
          unsigned int i_tt_pair = curr_num_threads + i_tt;
          unsigned int i_fq_pair = temp * curr_num_threads + i_fq;
          if(curr_num_threads < curr_num_triangles &&
              i_tt < curr_num_threads && i_tt_pair < curr_num_triangles) {
            fq_d[i_fq] = cuCaddf(fq_d[i_fq], fq_d[i_fq_pair]);
          } // if
          curr_num_threads = curr_num_threads >> 1;
        } // for
        // now just copy the first part of fq_d to ff_d
        unsigned int i_ff = i_ff_base + i_xx;
        ff_d[i_ff] = cuCmulf(fq_d[i_ff], make_cuFloatComplex(0.0, -1.0));*/

        // with shared memory
        unsigned int fq_index = i_fq_base + i_xx;
        shared_fq[i_sfq] = fq_d[fq_index];
        unsigned int fq_index_pair = fq_index + temp * blockDim.x;
        unsigned int i_sfq_pair = i_sfq + stemp * blockDim.x;
        if(blockDim.x + i_tt < curr_num_triangles) shared_fq[i_sfq_pair] = fq_d[fq_index_pair];
      } // if

      __syncthreads();

      unsigned int curr_num_threads = blockDim.x;
      for(unsigned int step = 0; step < red_steps; ++ step) {
        unsigned int i_tt_pair = curr_num_threads + threadIdx.x;
        unsigned int i_sfq_pair = stempyz * curr_num_threads + i_sfq;
        if(i_tt < curr_num_triangles && i_yy < curr_nqy && i_zz < curr_nqz) {
          if(curr_num_threads < curr_num_triangles &&
              threadIdx.x < curr_num_threads && i_tt_pair < curr_num_triangles) {
            shared_fq[i_sfq] = cuCaddf(shared_fq[i_sfq], shared_fq[i_sfq_pair]);
          } // if
        } // if

        __syncthreads();

        curr_num_threads = curr_num_threads >> 1;
      } // for
      if(i_tt < curr_num_triangles && i_yy < curr_nqy && i_zz < curr_nqz) {
        // now just copy the first part of fq_d to ff_d
        unsigned int i_ff = i_ff_base + i_xx;
        if(i_sfq < stempyz) ff_d[i_ff] = cuCmulf(shared_fq[i_sfq], make_cuFloatComplex(0.0, -1.0));
      } // if
    } // for
  } // reduction_kernel_parallel()
  #endif // REDUCTION2

  /**
   * GPU: For double precision
   */

  /* Original reduction kernel; with no shared memory
   */
  __global__ void reduction_kernel(cuDoubleComplex* fq_d,
            unsigned int curr_nqx, unsigned int curr_nqy,
            unsigned int curr_nqz, unsigned int curr_num_triangles,
            unsigned int b_nqx, unsigned int b_nqy,
            unsigned int b_nqz, unsigned int b_num_triangles,
            cuDoubleComplex* ff_d) {
    // 3D block, where each thread is responsible for one point in the x, y, z grid
    unsigned int x = threadIdx.x + blockDim.x * blockIdx.x;
    unsigned int y = threadIdx.y + blockDim.y * blockIdx.y;
    unsigned int z = threadIdx.z + blockDim.z * blockIdx.z;
    unsigned int i_ff = curr_nqx * curr_nqy * z + curr_nqx * y + x;
    unsigned int temp = curr_nqx * curr_nqy * curr_nqz;

    if(x < curr_nqx && y < curr_nqy && z < curr_nqz) {
      cuDoubleComplex total = make_cuDoubleComplex(0.0, 0.0);
      for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
        unsigned int i_fq = temp * i_t + i_ff;
        total = cuCadd(total, fq_d[i_fq]);
      } // for
      ff_d[i_ff] = cuCmul(total, make_cuDoubleComplex(0.0, -1.0));
    } // if
  } // reduction_kernel()

  #ifdef REDUCTION2
  /* GPU parallel reduction:
   * perform 3D decomposition along t, y, z, and reduce along t in parallel
   */
  __global__ void reduction_kernel_parallel(cuDoubleComplex* fq_d,
            unsigned int curr_nqx, unsigned int curr_nqy,
            unsigned int curr_nqz, unsigned int curr_num_triangles,
            unsigned int b_nqx, unsigned int b_nqy,
            unsigned int b_nqz, unsigned int b_num_triangles,
            cuDoubleComplex* ff_d) {
    // b_num_triangles should be a power of 2
    // there will be only one block in the t direction
    // => blockDim.x = curr_num_triangles / 2 = h_t / 2
    // => hyperblock's t dimension size should be power of 2
    unsigned int i_tt = threadIdx.x + blockDim.x * blockIdx.x; // == threadIdx.x
    unsigned int i_yy = threadIdx.y + blockDim.y * blockIdx.y;
    unsigned int i_zz = threadIdx.z + blockDim.z * blockIdx.z;
    unsigned int i_ff_base = curr_nqx * (curr_nqy * i_zz + i_yy);
    unsigned int temp = curr_nqx * curr_nqy * curr_nqz;
    unsigned int stempyz = blockDim.y * blockDim.z;
    unsigned int stemp = curr_nqx * stempyz;
    unsigned int i_fq_base = temp * i_tt + i_ff_base;
    //unsigned int i_sbase = stemp * threadIdx.x + curr_nqx * (blockDim.y * threadIdx.z + threadIdx.y);

    // blockDim.x == BLOCK_REDUCTION_X_ / 2 == b_num_triangles / 2
    // blockDim.y == BLOCK_REDUCTION_Y_
    // blockDim.z == BLOCK_REDUCTION_Z_
    __shared__ cuDoubleComplex shared_fq[BLOCK_REDUCTION_T_ * BLOCK_REDUCTION_Y_ * BLOCK_REDUCTION_Z_];

    unsigned int red_steps = __ffs(b_num_triangles) - 1;
    unsigned int i_sfq = stempyz * threadIdx.x + blockDim.y * threadIdx.z + threadIdx.y;
    if(i_tt < curr_num_triangles && i_yy < curr_nqy && i_zz < curr_nqz) {
      for(unsigned int i_xx = 0; i_xx < curr_nqx; ++ i_xx) {
        /* // without shared memory
        unsigned int curr_num_threads = blockDim.x; // == b_num_triangles >> 1
        unsigned int i_fq = i_fq_base + i_xx;
        for(unsigned int step = 0; step < red_steps; ++ step) {
          unsigned int i_tt_pair = curr_num_threads + i_tt;
          unsigned int i_fq_pair = temp * curr_num_threads + i_fq;
          if(curr_num_threads < curr_num_triangles &&
              i_tt < curr_num_threads && i_tt_pair < curr_num_triangles) {
            fq_d[i_fq] = cuCadd(fq_d[i_fq], fq_d[i_fq_pair]);
          } // if
          curr_num_threads = curr_num_threads >> 1;
        } // for
        // now just copy the first part of fq_d to ff_d
        unsigned int i_ff = i_ff_base + i_xx;
        ff_d[i_ff] = cuCmul(fq_d[i_ff], make_cuDoubleComplex(0.0, -1.0));*/

        // with shared memory
        unsigned int fq_index = i_fq_base + i_xx;
        shared_fq[i_sfq] = fq_d[fq_index];
        unsigned int fq_index_pair = fq_index + temp * blockDim.x;
        unsigned int i_sfq_pair = i_sfq + stemp * blockDim.x;
        if(blockDim.x + i_tt < curr_num_triangles) shared_fq[i_sfq_pair] = fq_d[fq_index_pair];
        __syncthreads();

        unsigned int curr_num_threads = blockDim.x;
        for(unsigned int step = 0; step < red_steps; ++ step) {
          unsigned int i_tt_pair = curr_num_threads + threadIdx.x;
          i_sfq_pair = stempyz * curr_num_threads + i_sfq;
          if(curr_num_threads < curr_num_triangles &&
              threadIdx.x < curr_num_threads && i_tt_pair < curr_num_triangles) {
            shared_fq[i_sfq] = cuCadd(shared_fq[i_sfq], shared_fq[i_sfq_pair]);
            __syncthreads();
          } // if
          curr_num_threads = curr_num_threads >> 1;
        } // for
        // now just copy the first part of fq_d to ff_d
        unsigned int i_ff = i_ff_base + i_xx;
        if(i_sfq < stempyz) ff_d[i_ff] = cuCmul(shared_fq[i_sfq], make_cuDoubleComplex(0.0, -1.0));
      } // for
    } // if
  } // reduction_kernel_parallel()
  #endif // REDUCTION2

  #else // GPUR

  /**
   * OpenMP: For single precision
   */
  void reduction_kernel(unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
              unsigned int curr_b_num_triangles, unsigned long int blocked_matrix_size,
              unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
              unsigned int num_triangles,
              unsigned int nqx, unsigned int nqy, unsigned int nqz,
              unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
              cuComplex* fq_buffer, cuComplex* ff) {
    unsigned long int curr_b_xyz = (unsigned long int) curr_b_nqx * curr_b_nqy * curr_b_nqz;
    // reduction over all triangles
    #pragma omp parallel
    {
      if(omp_get_thread_num() == 0)
        std::cout << "[" << omp_get_num_threads() << " threads] ... " << std::flush;
      #pragma omp for
      for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
        for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
          for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
            unsigned long int temp_i = (unsigned long int) curr_b_nqx * curr_b_nqy * i_z +
                          curr_b_nqx * i_y + i_x;
            cuComplex total = make_cuFloatComplex(0.0, 0.0);
            for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
              //if(curr_b_xyz * i_t + temp_i >= blocked_matrix_size)
              //  std::cout << "OMG OMG OMG!!!!! " << curr_b_xyz << " " << i_t << " "
              //        << temp_i << std::endl << std::flush;
              total = cuCaddf(total, fq_buffer[curr_b_xyz * i_t + temp_i]);
            } // for i_t
            unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
                          nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
            ff[super_i] = cuCaddf(ff[super_i], cuCmulf(total, make_cuFloatComplex(0.0, -1.0)));
          } // for i_z
        } // for i_y
      } // for i_x
    } // omp parallel
  } // reduction_kernel()


  /**
   * OpenMP: For double precision
   */
  void reduction_kernel(unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
              unsigned int curr_b_num_triangles, unsigned long int blocked_matrix_size,
              unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
              unsigned int num_triangles,
              unsigned int nqx, unsigned int nqy, unsigned int nqz,
              unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
              cuDoubleComplex* fq_buffer, cuDoubleComplex* ff) {
    unsigned long int curr_b_xyz = (unsigned long int) curr_b_nqx * curr_b_nqy * curr_b_nqz;
    // reduction over all triangles
    #pragma omp parallel
    {
      if(omp_get_thread_num() == 0)
        std::cout << "[" << omp_get_num_threads() << " threads] ... " << std::flush;
      #pragma omp for
      for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
        for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
          for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
            unsigned long int temp_i = (unsigned long int) curr_b_nqx * curr_b_nqy * i_z +
                          curr_b_nqx * i_y + i_x;
            cuDoubleComplex total = make_cuDoubleComplex(0.0, 0.0);
            for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
              //if(curr_b_xyz * i_t + temp_i >= blocked_matrix_size)
              //  std::cout << "OMG OMG OMG!!!!! " << curr_b_xyz << " " << i_t << " "
              //        << temp_i << std::endl << std::flush;
              total = cuCadd(total, fq_buffer[curr_b_xyz * i_t + temp_i]);
            } // for i_t
            unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
                          nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
            ff[super_i] = cuCadd(ff[super_i], cuCmul(total, make_cuDoubleComplex(0.0, -1.0)));
            //ff[super_i] = make_cuDoubleComplex(ib_y, i_y);
          } // for i_z
        } // for i_y
      } // for i_x
    } // omp parallel
  } // reduction_kernel()

  #endif // GPUR

} // namespace hig

