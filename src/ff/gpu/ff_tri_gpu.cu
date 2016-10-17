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
#include <complex>

#include <cuComplex.h>
#include <stdio.h>

#include <ff/gpu/ff_num_gpu.cuh>
#include <common/enums.hpp>
#include <common/constants.hpp>
#include <common/globals.hpp>
#include <common/parameters.hpp>
#include <numerics/numeric_utils.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>
#include <utils/gpu/cu_utilities.cuh>


namespace hig {

  __constant__ int num_loads;
  __constant__ int buf_size;

  __device__ __inline__ cucomplex_t cuC_dot(cucomplex_t * a, real_t * b){
      return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  }

  __device__ __inline__ cucomplex_t cuC_dot(cucomplex_t * a, vector3_t b){
      return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  }

  __device__ cucomplex_t FormFactorTriangle (triangle_t & tri, 
          cucomplex_t * mq){
    cucomplex_t cuCZERO = make_cuC(REAL_ZERO_, REAL_ZERO_);
    cucomplex_t jp = make_cuC(REAL_ZERO_, REAL_ONE_);
    cucomplex_t jn = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
    cucomplex_t ff = cuCZERO;

    // calculate q^2
    real_t q_sqr = cuCnorm3(mq[0], mq[1], mq[2]);

    // form vertices
    vector3_t vertex[3];
    vertex[0] = vector3_t(tri.v1[0], tri.v1[1], tri.v1[2]);
    vertex[1] = vector3_t(tri.v2[0], tri.v2[1], tri.v2[2]);
    vertex[2] = vector3_t(tri.v3[0], tri.v3[1], tri.v3[2]);

    // form edges
    vector3_t edge[3];
    edge[0] = vertex[1] - vertex[0];
    edge[1] = vertex[2] - vertex[1];
    edge[2] = vertex[0] - vertex[2];

    // calculate outward normal and area
    vector3_t n_t = cross(edge[0], edge[1]);
    real_t t_area = 0.5 * n_t.abs();
    n_t =  n_t / n_t.abs(); // normalize vector
 
    // dot (q, n_t)
    cucomplex_t q_dot_nt = cuC_dot(mq, n_t);

    // calculate projection
    real_t proj_tq = q_sqr - cuC_norm(q_dot_nt);

    // CASE 1
    if (abs(proj_tq) < CUTINY_){
      cucomplex_t q_dot_v = cuC_dot(mq, vertex[0]);

      // calculate Form-Factor
      ff = jp * q_dot_nt * t_area / q_sqr * cuCexp(jn * q_dot_v);
    } else {
      // iterate on each edge :
      for (int e = 0; e < 3; e++) {

        // edge normal
        vector3_t n_e = cross(edge[e], n_t);
        n_e = n_e / n_e.abs(); // normalize

        // dot(q, n_e)
        cucomplex_t q_dot_ne = cuC_dot(mq, n_e);

        // proj_ne
        real_t proj_eq = proj_tq - cuC_norm(q_dot_ne);

        // CASE 2
        if (abs(proj_eq) < CUTINY_){
          // q_dot_v
          cucomplex_t q_dot_v = cuC_dot(mq, vertex[e]);

          // calculate contribution of edge
          real_t f0 = edge[e].abs() / (q_sqr * proj_tq);
          cucomplex_t c0 = -1 * q_dot_nt * q_dot_ne;
          cucomplex_t c1 = cuCexp(jn * q_dot_v);
          ff = ff + (f0 * c0 * c1);
        } else {
          // CASE 3 (General case)
          real_t f0 = q_sqr * proj_tq * proj_eq;

          // dot(q, v_a) vertex a
          cucomplex_t q_dot_v = cuC_dot(mq, vertex[e]);

          // vertex normal a ... whatever that means :-)
          vector3_t n_v = edge[e] / edge[e].abs();

          // dot (q, n_v)
          cucomplex_t q_dot_nv = cuC_dot(mq, n_v);

          // calculate contribution of vertex a
          cucomplex_t c0 = jn * q_dot_nt * q_dot_ne * q_dot_nv;
          cucomplex_t c1 = cuCexp(jn * q_dot_v);
          ff = ff + c0 * c1 / f0;

          // dot(q, v_b) the other vertex in the edge
          int ep = (e+1) % 3; // it is cyclic
          q_dot_v = cuC_dot(mq, vertex[ep]);

          // dot(q, n_v)
          q_dot_nv = cuC_dot(mq, n_v * REAL_MINUS_ONE_);

          // calculate contribution of the other vertex
          c0 = jn * q_dot_nt * q_dot_ne * q_dot_nv;
          c1 = cuCexp(jn * q_dot_v);
          ff = ff + c0 * c1 / f0;
        }
      }
    }
    return ff;
  }

  __global__ void ff_poly_kernel (unsigned int nqy, unsigned int nqz,
                  real_t *qx, real_t *qy, cucomplex_t *qz,
                  int num_triangles, triangle_t * triangles,
                  RotMatrix_t rot, cucomplex_t * ff) {

    unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if ( i_z < nqz ) {
      unsigned int i_y = i_z % nqy;
      cucomplex_t ff_temp = make_cuC(REAL_ZERO_, REAL_ZERO_);
      cucomplex_t mq[3];
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mq[0], mq[1], mq[2]);

      for (int i=0; i < num_triangles; i++)
        ff_temp = ff_temp + FormFactorTriangle(triangles[i], mq);

/***************
      // declare shared memory
      extern __shared__ triangle_t shared_triangles[];
      int i_thread = threadIdx.x;
      int nthreads = blockDim.x;

      // load triangles into shared memory
      const int T_SIZE_ = sizeof(triangle_t);
      int curr_num_tri = buf_size / T_SIZE_;
      int repeats = curr_num_tri / nthreads;
      int curr_pos = 0;

      // start loading
      for (int i_load = 0; i_load < num_loads; i_load++) {
        int residual = num_triangles - curr_pos;
        if ( residual * T_SIZE_ > buf_size ) {
          for (int i = 0; i < repeats; i++) {
            int i_shared = i * nthreads + i_thread;
            shared_triangles[i_shared] = triangles[curr_pos + i_shared];
          }
          curr_pos += curr_num_tri;
        } else {
          repeats = residual / nthreads;
          for (int i = 0; i < repeats; i++) {
            int i_shared = i * nthreads + i_thread;
            shared_triangles[i_shared] = triangles[curr_pos + i_shared];
          }
          int nn = residual % nthreads;
          if (i_thread < nn) {
            int i_shared = repeats * nthreads + i_thread;
            shared_triangles[i_shared] = triangles[curr_pos + i_shared];
          }
          curr_num_tri = num_triangles - curr_pos;
        }

        // wait for everyone to join here
        __syncthreads();

        for (int i = 0; i < curr_num_tri; i++) {
          ff_temp =  ff_temp + FormFactorTriangle (shared_triangles[i], mq);
        }
      }
***************/

      ff[i_z] = ff_temp;
    } // if (i_z < nqz)
  } //ff_poly_kernel

  /**
   * The main host function called from outside, as part of the API for a single node.
   */
  unsigned int NumericFormFactorG::compute_exact_triangle(
            triangle_t * triangles, int num_triangles,
            cucomplex_t *& ff,
            int nqy, real_t * qx_h, real_t * qy_h, 
            int nqz, cucomplex_t * qz_h,
            RotMatrix_t & rot, real_t & kernel_time) { 
      
    cudaError_t err;
    real_t *qx_d, *qy_d;
    cucomplex_t *qz_d;
    cucomplex_t * ff_d;
    triangle_t * triangles_d;

    std::cout << "**        Computing FF on GPU ..." << std::endl;
    cudaEvent_t begin_event, end_event;
    cudaEventCreate(&begin_event); cudaEventCreate(&end_event);
    cudaEventRecord(begin_event, 0);

    ff = new (std::nothrow) cucomplex_t[nqz];
    if (ff == NULL){
      std::cerr << "Error: failed to allocate memory of size: " 
          << nqz * sizeof(cucomplex_t) << std::endl;
      return false;
    }
 
    // Allocate memory for device-side ff 
    if(cudaMalloc((void **) &ff_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
        std::cerr << "Device memory allocation failed for ff_d. " << std::endl;
        return 0;
    }
    // Allocate memory for qx
    if(cudaMalloc((void **) &qx_d, nqy * sizeof(real_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qx_d. "
            << "Size = " << nqy * sizeof(real_t) << " B" << std::endl;
      return 0;
    } 
    // Allocate memory for qy
    if(cudaMalloc((void **) &qy_d, nqy * sizeof(real_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qy_d. "
            << "Size = " << nqy * sizeof(real_t) << " B" << std::endl;
      return 0;
    } 
    // Allocate memory for qz
    if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qz_d. "
            << "Size = " << nqz * sizeof(cucomplex_t) << " B" << std::endl;
      return 0;
    } 
    // Allocate memory for triangles
    err = cudaMalloc((void **) &triangles_d, num_triangles * sizeof(triangle_t));
    if(err != cudaSuccess) {
      std::cerr << "Device memory allocation failed for triangles. "
            << "Size = " <<  num_triangles * sizeof(triangle_t) << " B" << std::endl;
      return 0;
    }

    // copy buffers to device memory
    cudaMemcpy(qx_d, qx_h, nqy * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_d, qy_h, nqy * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
    cudaMemcpy(triangles_d, triangles, num_triangles * sizeof(triangle_t), cudaMemcpyHostToDevice);

    // number of cuda threads and blocks
    int cuda_num_threads = 256;
    int cuda_num_blocks = nqz / cuda_num_threads + 1;
    
    // calculate shared memeory size
    int max_load = 49152; // Maximum shared memory
    int single_load = cuda_num_threads * sizeof(triangle_t);
    int shared_buf = (max_load / single_load) * single_load;
    int num_of_loads = num_triangles * sizeof(triangle_t) / shared_buf + 1;

    // copy the rotation matrix and other stuff to constant memory
    err = cudaMemcpyToSymbol(num_loads, &num_of_loads, sizeof(int), 0, cudaMemcpyHostToDevice);
    err = cudaMemcpyToSymbol(buf_size, &shared_buf, sizeof(int), 0, cudaMemcpyHostToDevice);

    // Kernel
    ff_poly_kernel <<< cuda_num_blocks, cuda_num_threads, shared_buf >>> (nqy, nqz, qx_d, qy_d, qz_d, 
            num_triangles, triangles_d, rot, ff_d);
    err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "Error: ff_poly_kernel failed miserably!" 
          << cudaGetErrorString(err) << std::endl;
      return 0;
    }

    err = cudaMemcpy (ff, ff_d, nqz * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
      std::cerr << "Error: failed to copy buffer from device" 
          << cudaGetErrorString(err) << std::endl;
      return 0;
    } 
 
    cudaFree(triangles_d);
    cudaFree(qz_d);
    cudaFree(qy_d);
    cudaFree(qx_d);
    cudaFree(ff_d);

    float total_time;
    cudaEventRecord(end_event, 0);
    cudaEventSynchronize(end_event);
    cudaEventElapsedTime(&total_time, begin_event, end_event);
    kernel_time = total_time;
    return num_triangles;

  } // NumericFormFactorG::compute_exact_triangle()


  __global__ void ff_tri_kernel (unsigned int nqy, unsigned int nqz,
                  real_t *qx, real_t *qy, cucomplex_t *qz,
                  int num_tri, real_t * shape_def, RotMatrix_t rot,
                  cucomplex_t * ff) {

    unsigned int i_thread = threadIdx.x;
    unsigned int nthreads = blockDim.x;
    unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i_z < nqz ) {
      unsigned int i_y = i_z % nqy;

      cucomplex_t ff_temp = make_cuC(REAL_ZERO_, REAL_ZERO_);

      // shared memory
      extern __shared__ real_t shared_shape_def[];

      // do the rotations
      cucomplex_t mqx, mqy, mqz;
      rot.rotate(qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      real_t q2 = cuCnorm3(mqx, mqy, mqz);

      // load triangles into shared memory
      int curr_pos = 0;
      int i_shared = 0;
      int curr_num_tri = buf_size / T_PROP_SIZE_;
      int repeats = buf_size / nthreads;
      for (int i_load = 0; i_load < num_loads; i_load++) {
        int residual = num_tri * T_PROP_SIZE_ - curr_pos;
        if ( residual > buf_size ) {
          for (int i = 0; i < repeats; i++) {
            i_shared = i * nthreads + i_thread;
            shared_shape_def[i_shared] = shape_def[curr_pos + i_shared];
          }
          curr_pos += buf_size;
        } else {
          repeats = residual / nthreads;
          for (int i = 0; i < repeats; i++) {
            int i_shared = i * nthreads + i_thread;
            shared_shape_def[i_shared] = shape_def[curr_pos + i_shared];
          }
          int nn = residual % nthreads;
          if (i_thread < nn) {
            i_shared = repeats * nthreads + i_thread;
            shared_shape_def[i_shared] = shape_def[curr_pos + i_shared];
          }
          curr_num_tri = (num_tri * T_PROP_SIZE_ - curr_pos) / T_PROP_SIZE_;
        }

        // wait for everyone to join here
        __syncthreads();

        for (int i = 0; i < curr_num_tri; i++) {
          unsigned int offset = i * T_PROP_SIZE_;
          real_t s =  shared_shape_def[offset];
          real_t nx = shared_shape_def[offset + 1];
          real_t ny = shared_shape_def[offset + 2];
          real_t nz = shared_shape_def[offset + 3];
          real_t x  = shared_shape_def[offset + 4];
          real_t y  = shared_shape_def[offset + 5];
          real_t z  = shared_shape_def[offset + 6];

          cucomplex_t qn = mqx * nx + mqy * ny + mqz * nz;
          cucomplex_t qt = mqx * x  + mqy * y  + mqz * z;
          cucomplex_t nj = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);
          ff_temp = ff_temp +  (nj * qn * s * cuCexpi(qt) / q2);
        } 
      } // for (i_load = 0; ...)
      ff[i_z] = ff_temp;
      
    } // if(i_z < nqz)
  } // formfactor_tri_kernel

  /**
   * The main host function called from outside, as part of the API for a single node.
   */
  unsigned int NumericFormFactorG::compute_approx_triangle(
            std::vector<real_t> & shape_def, 
            cucomplex_t * & ff,
            int nqy, real_t * qx_h, real_t * qy_h,
            int nqz, cucomplex_t * qz_h,
            RotMatrix_t & rot, real_t & kernel_time) { 
      

    int i; 
    cudaError_t err;
    real_t * qx_d, * qy_d;
    cucomplex_t *qz_d;
    cucomplex_t * ff_d;

    unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
    if(num_triangles < 1) return 0;

    std::cout << "**        Computing FF on GPU ..." << std::endl;
    cudaEvent_t begin_event, end_event;
    cudaEventCreate(&begin_event); cudaEventCreate(&end_event);
    cudaEventRecord(begin_event, 0);

    // calculate triangle load sizes
    unsigned int cuda_num_threads = 256;
    unsigned int cuda_num_blocks = nqz / cuda_num_threads + 1;
    // assuming shared mem = 48KB 
    int max_per_load = 49152 / sizeof(real_t);
    for (i = 0; i < num_triangles; i++)
        if ( i * cuda_num_threads * T_PROP_SIZE_ > max_per_load)
            break;
    int shared_buf = (i-1) * cuda_num_threads * T_PROP_SIZE_;
    int num_of_loads = num_triangles * T_PROP_SIZE_ / shared_buf + 1;

    // allocate memory for the final FF matrix
    ff = new (std::nothrow) cucomplex_t[nqz];  // allocate and initialize to 0
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << nqz * sizeof(cucomplex_t) << " b" << std::endl;
      return 0;
    } // if
    for (int i = 0; i < nqz; i++) { 
        ff[i].x = (real_t) 0.;
        ff[i].y = (real_t) 0.;
    }

    if(cudaMalloc((void **) &ff_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
        std::cerr << "Device memory allocation failed for ff_d. " << std::endl;
        return 0;
    }

    if(cudaMalloc((void **) &qx_d, nqy * sizeof(real_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qx_d. "
            << "Size = " << nqy * sizeof(real_t) << " B" << std::endl;
      return 0;
    } // if

    if(cudaMalloc((void **) &qy_d, nqy * sizeof(real_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qy_d. "
            << "Size = " << nqy * sizeof(real_t) << " B" << std::endl;
      return 0;
    } // if

    if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qz_d. "
            << "Size = " << nqz * sizeof(cucomplex_t) << " B" << std::endl;
      return 0;
    } // if

    cudaMemcpy(qx_d, qx_h, nqy * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_d, qy_h, nqy * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

    // copy triangles
    real_t * shape_def_d;
    real_t * shape_def_h = &shape_def[0];
    err = cudaMalloc((void **) &shape_def_d, T_PROP_SIZE_ * num_triangles * sizeof(real_t));
    if(err != cudaSuccess) {
      std::cerr << "Device memory allocation failed for triangles. "
            << "Size = " <<  num_triangles << std::endl;
      return 0;
    } // if

    err = cudaMemcpy(shape_def_d, shape_def_h, 
            T_PROP_SIZE_ * num_triangles * sizeof(real_t), cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
          << cudaGetErrorString(err) << std::endl;
      return 0;
    }


    /* copy symbols to constant memory */
     
    // copy rotation matrix, num_of_loads and shared_buffer_size to constant memeory
    err = cudaMemcpyToSymbol(num_loads, &num_of_loads, sizeof(int), 0, cudaMemcpyHostToDevice);
    err = cudaMemcpyToSymbol(buf_size, &shared_buf, sizeof(int), 0, cudaMemcpyHostToDevice);
    
    // Kernel
    ff_tri_kernel <<< cuda_num_blocks, cuda_num_threads, shared_buf * sizeof(real_t) >>> (nqy, nqz, qx_d, 
            qy_d, qz_d, num_triangles, shape_def_d, rot, ff_d);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error: kernel crashed with message: " << cudaGetErrorString(err) << std::endl;
        return 0;
    }

    err = cudaMemcpy (ff, ff_d, nqz * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error: Memcpy failed with message: " << cudaGetErrorString(err) << std::endl;
      exit(1);
    } 

    float total_time;
    cudaEventRecord(end_event, 0);
    cudaEventSynchronize(end_event);
    cudaEventElapsedTime(&total_time, begin_event, end_event);
    kernel_time = total_time;

    cudaFree(qz_d);
    cudaFree(qy_d);
    cudaFree(qx_d);
    cudaFree(ff_d);
    cudaFree(shape_def_d);
    return num_triangles;

  } // NumericFormFactorG::compute_approx_triangle()

} // namespace hig

