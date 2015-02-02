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
#include <stdio.h>

#include <ff/gpu/ff_num_gpu.cuh>
#include <common/enums.hpp>
#include <common/constants.hpp>
#include <common/parameters.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>
#include <utils/gpu/cu_utilities.cuh>


namespace hig {

  extern __constant__ float_t rot_d[9];
  __constant__ int num_loads;
  __constant__ int buf_size;
  extern __shared__ float_t dynamic_shared[];

  /* calculates outward normal to the surface, using cross-product.
   * Order of the edges is important (counter-clockwise)
   */
  __device__ __inline__ float_t surface_normal (triangle_t & tri, float_t * t_norm) {
      int i;
      float_t edge1[3], edge2[3];
      // form edge vectors
      for (i = 0; i < 3; i++ ) edge1[i] = tri.v2[i] - tri.v1[i];
      for (i = 0; i < 3; i++ ) edge2[i] = tri.v3[i] - tri.v2[i];

      // calculate cross-product
      cross_prod (edge1, edge2, t_norm);

      // normalize to get unit-vector
      float_t nn = norm2 (t_norm);
      for (i = 0; i < 3; i++) t_norm[i] /= nn;
      return 0.5 * nn;
  }

  /* calculates normal to edge which lies in plane of the polygon */
  __device__ __inline__ void edge_normal (float_t * vert1, float_t * vert2,
          float_t * t_nrom, float_t * e_norm) {
      int i;
      float_t edge[3];

      // form edge vector
      for (int i = 0; i < 3; i++) edge[i] = vert2[i] - vert1[i];

      //calulcate cross-product
      cross_prod (edge, t_nrom, e_norm);

      // normalize to get unit-vector
      float_t nn = norm2 (e_norm);
      for (i = 0; i < 3; i++) e_norm[i] /= nn;
  }

  /* calculates normal to the vertex (whatever that means) */
  __device__ __inline__ void vert_normal (float_t * vert1, float_t * vert2,
          float_t * v_norm) {
      int i;
      for (i = 0; i < 3; i++) v_norm[i] = vert1[i] - vert2[i];
      float_t nn = norm2 (v_norm);
      for (i = 0; i < 3; i++) v_norm[i] /= nn;
  }


  /* length of the edge */
  __device__ __inline__ float_t length (float_t * v1, float_t * v2) {
      float_t len = 0;
      for (int i = 0; i < 3; i++) len += (v1[i] - v2[i]) * (v1[i] - v2[i]);
      return sqrt(len);
  }

  __device__ cucomplex_t FormFactorTriangle (triangle_t & tri,
          cucomplex_t qx, cucomplex_t qy, cucomplex_t qz) {
      
      int i;
      cucomplex_t ff = make_cuC ((float_t) 0.0, (float_t) 0.0);
      cucomplex_t unitc = make_cuC ((float_t) 0., (float_t) 1.);

      float_t t_norm[3];
      float_t e_norm[3];
      float_t v1_norm[3];
      float_t v2_norm[3];

      // calculate outward normal
      float_t area_t = surface_normal (tri, t_norm);

      // calculate proj_tq norm
      float_t qnorm = cuCnorm3 (qx, qy, qz);
      float_t qnorm2 = qnorm * qnorm;
      cucomplex_t q_dot_nt = qx * t_norm[0] + qy * t_norm[1] + qz * t_norm[2];
      float_t temp = cuCabsolute (q_dot_nt);
      float_t proj_tq = qnorm2 - temp * temp;
      if (proj_tq < 1.0E-08) {
           ff = unitc * q_dot_nt * cuCexpi (-1.0 * (qx * tri.v1[0] + qy * tri.v1[1] + qz * tri.v1[2])) 
              * area_t / qnorm2;
           return ff;
      }

      float_t edge[4][3];
      edge[0][0] = tri.v1[0]; edge[0][1] = tri.v1[1]; edge[0][2] = tri.v1[2];
      edge[1][0] = tri.v2[0]; edge[1][1] = tri.v2[1]; edge[1][2] = tri.v2[2];
      edge[2][0] = tri.v3[0]; edge[2][1] = tri.v3[1]; edge[2][2] = tri.v3[2];
      edge[3][0] = tri.v1[0]; edge[2][1] = tri.v1[1]; edge[3][2] = tri.v1[2];

      // calculate form-factor
      float_t p1[3], p2[3];
      for (int e = 0; e < 3; e++) {
          for (i = 0; i < 3; i++) p1[i] = edge[e][i];
          for (i = 0; i < 3; i++) p2[i] = edge[e+1][i];
         
          // edge-normal 
          edge_normal (p1, p2, t_norm, e_norm);

          // calculate proj_eq
          cucomplex_t q_dot_ne = qx * e_norm[0] + qy * e_norm[1] + qz * e_norm[2];
          float_t temp = cuCabsolute (q_dot_ne);
          float_t proj_eq = proj_tq - temp * temp;
          if ( proj_eq < 1.0E-08 ) {
              float_t e_len  = length (p1, p2);
              ff = ff - cuCexpi (-1.0 * (qx * p1[0] + qy * p1[1] + qz * p1[2])) * 
                  q_dot_nt * q_dot_ne * e_len / qnorm2 / proj_tq;
              ff = make_cuC(0.0f, 1.0f);
          } else {
              vert_normal (p1, p2, v1_norm);
              vert_normal (p2, p1, v2_norm);
              cucomplex_t q_dot_v1 = qx * v1_norm[0] + qy * v1_norm[1] + qz * v1_norm[2];
              cucomplex_t q_dot_v2 = qx * v2_norm[0] + qy * v2_norm[1] + qz * v2_norm[2];
              cucomplex_t tmpff = unitc * q_dot_nt * q_dot_ne 
                  / qnorm2 / proj_tq / proj_eq;
              ff = ff - tmpff * (cuCexpi(-1.0 * (qx * p1[0] + qy * p1[1] + qz * p1[2])) +
                      cuCexpi (-1.0 * (qx * p2[0] + qy * p2[1] + qz * p2[2])));
          }
      }
      return ff;
  }

  __global__ void ff_poly_kernel (unsigned int nqy, unsigned int nqz,
                  float_t *qx, float_t *qy, cucomplex_t *qz,
                  int num_triangles, triangle_t * triangles,
                  cucomplex_t *ff) {
    unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    if ( i_z < nqz ) {
      unsigned int i_y = i_z % nqy;
      cucomplex_t ff_temp = make_cuC(ZERO, ZERO);
      cucomplex_t mqx, mqy, mqz;

      rotate_q (rot_d, qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      for (int i = 0; i < num_triangles; i++) {
        ff_temp =  ff_temp + FormFactorTriangle (triangles[i], mqx, mqy, mqz);
      }
      ff[i_z] = ff_temp;
    }
  } //ff_poly_kernel

  /**
   * The main host function called from outside, as part of the API for a single node.
   */
  unsigned int NumericFormFactorG::compute_exact_triangle(
            triangle_t * triangles, int num_triangles,
            cucomplex_t *& ff,
            int nqy, float_t * qx_h, float_t * qy_h, 
            int nqz, cucomplex_t * qz_h,
            float_t * rot, float_t & kernel_time) { 
      
    int i;
    cudaError_t err;
    float_t *qx_d, *qy_d;
    cucomplex_t *qz_d;
    cucomplex_t * ff_d;
    triangle_t * triangles_d;

    cudaEvent_t begin_event, end_event;
    cudaEventCreate(&begin_event); cudaEventCreate(&end_event);
    cudaEventRecord(begin_event, 0);

    // allocate memory for the final FF matrix
    ff = new (std::nothrow) cucomplex_t[nqz];  // allocate and initialize to 0
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << nqz * sizeof(cucomplex_t) << " b" << std::endl;
      return 0;
    } // if
    for (int i = 0; i < nqy; i++) { 
        ff[i].x = (float_t) 0.;
        ff[i].y = (float_t) 0.;
    }

    if(cudaMalloc((void **) &ff_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
        std::cerr << "Device memory allocation failed for ff_d. " << std::endl;
        delete [] ff;
        return 0;
    }

    if(cudaMalloc((void **) &qx_d, nqy * sizeof(float_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qx_d. "
            << "Size = " << nqy * sizeof(float_t) << " B" << std::endl;
      delete[] ff;
      cudaFree(ff_d);
      return 0;
    } // if

    if(cudaMalloc((void **) &qy_d, nqy * sizeof(float_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qy_d. "
            << "Size = " << nqy * sizeof(float_t) << " B" << std::endl;
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    } // if

    if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qz_d. "
            << "Size = " << nqz * sizeof(cucomplex_t) << " B" << std::endl;
      cudaFree(qy_d);
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    } // if

    cudaMemcpy(qx_d, qx_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

    // copy triangles
    err = cudaMalloc((void **) &triangles_d, num_triangles * sizeof(triangle_t));
    if(err != cudaSuccess) {
      std::cerr << "Device memory allocation failed for triangles. "
            << "Size = " <<  num_triangles * sizeof(triangle_t) << " B" << std::endl;
      cudaFree(qz_d);
      cudaFree(qy_d);
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    } // if
    err = cudaMemcpy(triangles_d, triangles, num_triangles * sizeof(triangle_t), cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
          << cudaGetErrorString(err) << std::endl;
      return 0;
    }

    // calculate shared memeory size
    int num_threads = 256;
    int num_blocks = nqz / num_threads + 1;
    int max_shared_mem = 49152 / sizeof(triangle_t);
    for (i = 0; i < num_triangles; i++) {
      if ( i * sizeof(triangle_t) * num_threads > max_shared_mem)
        break;
    }
    size_t shared_mem = (i-1) * sizeof(triangle_t) * num_threads;
    int num_of_loads = i;

    // copy the rotation matrix and other stuff to constant memory
    /*
    err = cudaMemcpyToSymbol (rot_d, rot, 9 * sizeof(float_t), 0, cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
          << cudaGetErrorString(err) << std::endl;
      return 0;
    }

    err = cudaMemcpyToSymbol(num_loads, &num_of_loads, sizeof(int), 0, cudaMemcpyHostToDevice);
    f(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
          << cudaGetErrorString(err) << std::endl;
      return 0;
    }
    */

    // Kernel
    ff_poly_kernel <<< num_blocks, num_threads >>> (nqy, nqz, qx_d, qy_d, qz_d, 
            num_triangles, triangles_d, ff_d);
    err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
          << cudaGetErrorString(err) << std::endl;
      return 0;
    }

    err = cudaMemcpy (ff, ff_d, nqz * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
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
                  float_t *qx, float_t *qy, cucomplex_t *qz,
                  int num_tri, float_t * shape_def,
                  cucomplex_t * ff) {

    unsigned int i_thread = threadIdx.x;
    unsigned int nthreads = blockDim.x;
    unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i_z < nqz ) {
      unsigned int i_y = i_z % nqy;

      cucomplex_t ff_temp = make_cuC(ZERO, ZERO);

      // shared memory
      extern __shared__ float_t shared_shape_def[];

      // do the rotations
      cucomplex_t mqx, mqy, mqz;
      rotate_q (rot_d, qx[i_y], qy[i_y], qz[i_z], mqx, mqy, mqz);
      float_t q = cuCnorm3(mqx, mqy, mqz);
      float_t q2 = q * q;

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
          unsigned int shape_off = i * T_PROP_SIZE_;
          float_t s = shared_shape_def[shape_off];
          float_t nx = shared_shape_def[shape_off + 1];
          float_t ny = shared_shape_def[shape_off + 2];
          float_t nz = shared_shape_def[shape_off + 3];
          float_t x = shared_shape_def[shape_off + 4];
          float_t y = shared_shape_def[shape_off + 5];
          float_t z = shared_shape_def[shape_off + 6];

          cucomplex_t qn = mqx * nx + mqy * ny + mqz * nz;
          cucomplex_t qt = mqx * x  + mqx * y  + mqz * z;
          cucomplex_t nj = make_cuC(ZERO, NEG_ONE);
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
            std::vector<float_t> & shape_def, 
            cucomplex_t * & ff,
            int nqy, float_t * qx_h, float_t * qy_h,
            int nqz, cucomplex_t * qz_h,
            float_t * rot_h, float_t & kernel_time) { 
      

    int i; 
    cudaError_t err;
    float_t * qx_d, * qy_d;
    cucomplex_t *qz_d;
    cucomplex_t * ff_d;

    unsigned int num_triangles = shape_def.size() / T_PROP_SIZE_;
    if(num_triangles < 1) return 0;

    cudaEvent_t begin_event, end_event;
    cudaEventCreate(&begin_event); cudaEventCreate(&end_event);
    cudaEventRecord(begin_event, 0);


    /***
    // Get all the compute 2+ GPUs on the node
    int num_dev = 0;
    std::vector<int> dev_id;
    cudaDeviceProp dev_prop;
    err = cudaGetDeviceCount (&num_dev);
    for (i = 0; i < num_dev; i++) {
      err = cudaGetDeviceProperties (&dev_prop, dev_id[i]);
      if ( dev_prop.major >= 2 )
          dev_id.push_back(i);
    }
    num_dev = dev_id.size();
    
    // distribute data amongst 'n' devices
    ***/

    // calculate triangle load sizes
    unsigned int cuda_num_threads = 256;
    unsigned int cuda_num_blocks = nqz / cuda_num_threads + 1;
    // assuming shared mem = 48KB 
    int max_per_load = 49152 / sizeof(float_t);
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
        ff[i].x = (float_t) 0.;
        ff[i].y = (float_t) 0.;
    }

    if(cudaMalloc((void **) &ff_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
        std::cerr << "Device memory allocation failed for ff_d. " << std::endl;
        delete [] ff;
        return 0;
    }

    if(cudaMalloc((void **) &qx_d, nqy * sizeof(float_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qx_d. "
            << "Size = " << nqy * sizeof(float_t) << " B" << std::endl;
      delete[] ff;
      cudaFree(ff_d);
      return 0;
    } // if

    if(cudaMalloc((void **) &qy_d, nqy * sizeof(float_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qy_d. "
            << "Size = " << nqy * sizeof(float_t) << " B" << std::endl;
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    } // if

    if(cudaMalloc((void **) &qz_d, nqz * sizeof(cucomplex_t)) != cudaSuccess) {
      std::cerr << "Device memory allocation failed for qz_d. "
            << "Size = " << nqz * sizeof(cucomplex_t) << " B" << std::endl;
      cudaFree(qy_d);
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    } // if

    cudaMemcpy(qx_d, qx_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_d, qy_h, nqy * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_d, qz_h, nqz * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

    // copy triangles
    float_t * shape_def_d;
    float_t * shape_def_h = &shape_def[0];
    err = cudaMalloc((void **) &shape_def_d, T_PROP_SIZE_ * num_triangles * sizeof(float_t));
    if(err != cudaSuccess) {
      std::cerr << "Device memory allocation failed for triangles. "
            << "Size = " <<  num_triangles << std::endl;
      cudaFree(qz_d);
      cudaFree(qy_d);
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    } // if

    err = cudaMemcpy(shape_def_d, shape_def_h, 
            T_PROP_SIZE_ * num_triangles * sizeof(float_t), cudaMemcpyHostToDevice);
    if(err != cudaSuccess) {
      std::cerr << "CUDA Error [" << __FILE__ << ":" << __LINE__ << "]: "
          << cudaGetErrorString(err) << std::endl;
      cudaFree(qz_d);
      cudaFree(qy_d);
      cudaFree(qx_d);
      cudaFree(ff_d);
      delete[] ff;
      return 0;
    }


    //for (i = 0; i < num_dev; i++) {
    //  err = cudaSetDevice (dev_id[i]);

      /* copy symbols to constant memory */
     
      // copy rotation matrix
      err = cudaMemcpyToSymbol(rot_d, rot_h, 9 * sizeof(float_t), 0, cudaMemcpyHostToDevice);

      // copy num_of_loads and shared_buffer_size to constant memeory
      err = cudaMemcpyToSymbol(num_loads, &num_of_loads, sizeof(int), 0, cudaMemcpyHostToDevice);
      err = cudaMemcpyToSymbol(buf_size, &shared_buf, sizeof(int), 0, cudaMemcpyHostToDevice);
    //}


    // Kernel
    ff_tri_kernel <<< cuda_num_blocks, cuda_num_threads, shared_buf * sizeof(float_t) >>> (nqy, nqz, qx_d, 
            qy_d, qz_d, num_triangles, shape_def_d, ff_d);
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
    return num_triangles;

  } // NumericFormFactorG::compute_approx_triangle()

} // namespace hig

