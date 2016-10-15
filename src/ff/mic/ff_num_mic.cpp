/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_mic.cpp
 *  Created: Apr 02, 2013
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
#include <cmath>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PROFILE_PAPI
#include <papi.h>
#endif

// for mkl vml functions
//#include <mkl_vml.h>

#include <woo/timer/woo_boostchronotimers.hpp>

#include <common/constants.hpp>
#include <common/parameters.hpp>
#include <common/mic/parameters_mic.hpp>
#include <common/mic/definitions_mic.hpp>
#include <numerics/mic/mic_complex_numeric.hpp>

// only for mic intrinsics
#include <numerics/mic/mic_avx_numerics.hpp>

#include <ff/mic/ff_num_mic.hpp>
  
namespace hig {

  __attribute__((target(mic)))
  int get_target_omp_num_threads() {
    int thread_count;
    #pragma omp parallel
    {
      #pragma omp single
      //thread_count = omp_get_num_threads();
      thread_count = omp_get_max_threads();
    }
    return thread_count;
  } // get_target_omp_num_threads()
  
  NumericFormFactorM::NumericFormFactorM() { }

  NumericFormFactorM::~NumericFormFactorM() { }

  // TODO: ...
  bool NumericFormFactorM::init() { return true; }

  /**
   * The main host function called from outside, as part of the API for a single node.
   * DO NOT USE THIS!!! IT DOES NOT HAVE ROTATION IN IT, AND HAS SEPARATED KERNELS!
   */
  unsigned int NumericFormFactorM::compute_form_factor(int rank,
            real_vec_t &shape_def_vec, complex_t* &ff,
            real_t* qx, int nqx, real_t* qy, int nqy, complex_t* qz, int nqz,
            real_t* rot,
            real_t& pass_kernel_time, real_t& red_time, real_t& mem_time
            #ifdef FINDBLOCK
              , const int block_x, const int block_y, const int block_z, const int block_t
            #endif
            ) {
    double kernel_time = 0.0, reduce_time = 0.0, total_kernel_time = 0.0, total_reduce_time = 0.0,
        temp_mem_time = 0.0, total_mem_time = 0.0;
    #ifdef _OPENMP
      if(rank == 0)
        std::cout << "++ Number of Host OpenMP threads: " << omp_get_max_threads() << std::endl;
    #endif

    int num_mic = 0;
    #ifdef __INTEL_OFFLOAD
      num_mic = _Offload_number_of_devices();
      if(rank == 0) {
        std::cout << "++         Number of Target MICs: " << num_mic << std::endl;
      } // if
      if(num_mic == 0) {
        std::cerr << "error: no Target MIC found!" << std::endl;
        return 0;
      } // if
    #else
      std::cerr << "warning: offloading to MIC not set! using host." << std::endl;
    #endif

    unsigned int num_triangles = shape_def_vec.size() / 7;
    if(num_triangles < 1) return 0;

    real_t* shape_def = &shape_def_vec[0];
  
    unsigned int total_qpoints = nqx * nqy * nqz;
  
    // allocate memory for the final FF 3D matrix
    ff = new (std::nothrow) complex_t[total_qpoints];  // allocate and initialize to 0
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << total_qpoints * sizeof(complex_t) << " b" << std::endl;
      return 0;
    } // if
    memset(ff, 0, total_qpoints * sizeof(complex_t));

    scomplex_t* qz_flat = new (std::nothrow) scomplex_t[nqz];
    if(qz_flat == NULL) {
      std::cerr << "Memory allocation failed for qz_flat. Size = "
            << nqz * sizeof(scomplex_t) << " b" << std::endl;
      return 0;
    } // if
    for(int i = 0; i < nqz; ++ i) {
      //qz_flat[i].x = qz[i].real(); qz_flat[i].y = qz[i].imag();
      qz_flat[i] = make_sC(qz[i].real(), qz[i].imag());
    } // for

    unsigned int host_mem_usage = ((unsigned int) nqx + nqy) * sizeof(real_t) +  // qx, qy
                    (unsigned int) nqz * sizeof(complex_t) +    // qz
                    (unsigned int) nqz * sizeof(scomplex_t) +    // qz_flat
                    (unsigned int) shape_def_vec.size() * sizeof(real_t) +  // shape_def
                    total_qpoints * sizeof(complex_t);        // ff

    // allocate memory buffers on the target and transfer data
    #pragma offload_transfer target(mic:0) \
                in(qx: length(nqx) MIC_ALLOC) \
                in(qy: length(nqy) MIC_ALLOC) \
                in(qz_flat: length(nqz) MIC_ALLOC) \
                in(shape_def: length(7 * num_triangles) MIC_ALLOC)
  
    unsigned int target_mem_usage = ((unsigned int) nqx + nqy) * sizeof(real_t) +
                    (unsigned int) nqz * sizeof(scomplex_t) +
                    (unsigned int) num_triangles * 7 * sizeof(real_t);

    unsigned int matrix_size = (unsigned int) nqx * nqy * nqz * num_triangles;

    // get some information about the target and display ...
    // TODO ...
    
    // do hyperblocking to use less memory
    unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
    compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
              b_nqx, b_nqy, b_nqz, b_num_triangles
              #ifdef FINDBLOCK
                , block_x, block_y, block_z, block_t
              #endif
              );

    if(rank == 0) {
      std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
            << " x " << b_nqz << " x " << b_num_triangles << std::endl;
    } // if
  
    unsigned int blocked_3d_matrix_size = (unsigned int) b_nqx * b_nqy * b_nqz;
    unsigned int blocked_matrix_size = (unsigned int) blocked_3d_matrix_size * b_num_triangles;

    // do target memory usage estimation ...
    // TODO ...
    
    size_t estimated_host_mem_need = host_mem_usage + blocked_matrix_size * sizeof(complex_t);
    if(rank == 0) {
      std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
            << " MB" << std::endl;
    } // if

    // allocate ff and fq buffers on host first
    scomplex_t *fq_buffer = new (std::nothrow) scomplex_t[blocked_matrix_size]();
    scomplex_t *ff_buffer = new (std::nothrow) scomplex_t[blocked_3d_matrix_size]();
    if(fq_buffer == NULL || ff_buffer == NULL) {
      std::cerr << "Memory allocation failed for f buffers. blocked_matrix_size = "
            << blocked_matrix_size << std::endl
            << "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
            << std::endl;
      delete[] ff;
      return 0;
    } // if

    host_mem_usage += (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(complex_t);

    if(rank == 0) {
      std::cout << "++              Host memory used: " << (float) host_mem_usage / 1024 / 1024
            << " MB" << std::endl << std::flush;
    } // if

    // allocate ff and fq buffers on target
    #pragma offload_transfer target(mic:0) \
              nocopy(fq_buffer: length(blocked_matrix_size) MIC_ALLOC) \
              nocopy(ff_buffer: length(blocked_3d_matrix_size) MIC_ALLOC)

    // display actual target memory used ...
    // TODO ...
  
    // compute the number of sub-blocks, along each of the 4 dimensions
    // formulate loops over each dimension, to go over each sub block
    unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
    unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
    unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
    unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);
  
    unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
    unsigned int curr_b_num_triangles = b_num_triangles;
    unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;
  
    if(rank == 0) {
      std::cout << "++         Number of hyperblocks: " << num_blocks
            << " [" << nb_x << " x " << nb_y << " x " << nb_z << " x " << nb_t << "]"
            << std::endl;
    } // if
  
    unsigned int block_num = 0;

    if(rank == 0) std::cout << "-- Computing form factor on MIC ... " << std::flush;

    // compute for each hyperblock
    curr_b_nqx = b_nqx;
    for(unsigned int ib_x = 0; ib_x < nb_x; ++ ib_x) {
      if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
      curr_b_nqy = b_nqy;
      for(unsigned int ib_y = 0; ib_y < nb_y; ++ ib_y) {
        if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
        curr_b_nqz = b_nqz;
        for(unsigned int ib_z = 0; ib_z < nb_z; ++ ib_z) {
          if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
          curr_b_num_triangles = b_num_triangles;
          for(unsigned int ib_t = 0; ib_t < nb_t; ++ ib_t) {
            if(ib_t == nb_t - 1)
              curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

            #ifdef VERBOSITY_3
              if(rank == 0) {
                std::cout << "-- Processing hyperblock " << ++ block_num << " of "
                    << num_blocks << ": Kernel... " << std::flush;
              } // if
            #endif

            // call the main kernel offloaded to MIC
            form_factor_kernel(qx, qy, qz_flat, shape_def,
                curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                b_nqx, b_nqy, b_nqz, b_num_triangles,
                ib_x, ib_y, ib_z, ib_t,
                fq_buffer);

            #ifdef VERBOSITY_3
              if(rank == 0) {
                std::cout << "done [" << temp_kernel_time << "s]. Reduction... "
                    << std::flush;
              } // if
            #endif

            // call the reduction kernel offloaded to MIC
            reduction_kernel(curr_b_nqx, curr_b_nqy, curr_b_nqz,
                curr_b_num_triangles, blocked_matrix_size,
                b_nqx, b_nqy, b_nqz, num_triangles,
                nqx, nqy, nqz,
                ib_x, ib_y, ib_z, ib_t,
                fq_buffer, ff_buffer);

            #ifdef VERBOSITY_3
              if(rank == 0) {
                std::cout << "done [" << temp_reduce_time << "s]." << std::endl;
              } // if
            #endif

            move_to_main_ff(ff_buffer, curr_b_nqx, curr_b_nqy, curr_b_nqz,
                    b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                    ib_x, ib_y, ib_z, ff);
          } // for ib_t
        } // for ib_z
      } // for ib_y
    } // for ib_x

    // free f buffers on target
    #pragma offload_transfer target(mic:0) \
                nocopy(fq_buffer: length(0) MIC_FREE) \
                nocopy(ff_buffer: length(0) MIC_FREE)
    // and host
    delete[] ff_buffer;
    delete[] fq_buffer;

    #pragma offload_transfer target(mic:0) \
                nocopy(shape_def: length(0) MIC_FREE) \
                nocopy(qx: length(0) MIC_FREE) \
                nocopy(qy: length(0) MIC_FREE) \
                nocopy(qz_flat: length(0) MIC_FREE)
  
    delete[] qz_flat;

    if(rank == 0) std::cout << "done." << std::endl;

    pass_kernel_time = total_kernel_time;
    red_time = total_reduce_time;
    mem_time = total_mem_time;
  
    return num_triangles;
  } // NumericFormFactorM::compute_form_factor()
    

  /**
   * the main Form Factor kernel function - for one hyperblock.
   */
  void NumericFormFactorM::form_factor_kernel(real_t* qx, real_t* qy, scomplex_t* qz_flat,
          real_t* shape_def,
          unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
          unsigned int curr_num_triangles,
          unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
          unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
          scomplex_t* fq_buffer) {
    if(fq_buffer == NULL || qx == NULL || qy == NULL || qz_flat == NULL) return;
  
    #pragma offload target(mic:0) \
            in(shape_def: length(0) MIC_REUSE) \
            in(qx: length(0) MIC_REUSE) \
            in(qy: length(0) MIC_REUSE) \
            in(qz_flat: length(0) MIC_REUSE) \
            out(fq_buffer: length(0) MIC_REUSE)
    {
      omp_set_num_threads(MIC_OMP_NUM_THREADS_);
      #pragma omp parallel
      {
        #pragma omp for nowait //schedule(auto)
        for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
          unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
          real_t s = shape_def[shape_off];
          real_t nx = shape_def[shape_off + 1];
          real_t ny = shape_def[shape_off + 2];
          real_t nz = shape_def[shape_off + 3];
          real_t x = shape_def[shape_off + 4];
          real_t y = shape_def[shape_off + 5];
          real_t z = shape_def[shape_off + 6];
  
          unsigned int xy_size = curr_nqx * curr_nqy;
          unsigned int matrix_off = xy_size * curr_nqz * i_t;
          unsigned int start_z = b_nqz * ib_z;
          unsigned int start_y = b_nqy * ib_y;
          unsigned int start_x = b_nqx * ib_x;
  
          for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
              ++ i_z, ++ global_i_z) {
            unsigned int off_start = matrix_off + xy_size * i_z;
            scomplex_t temp_z = qz_flat[global_i_z];
            scomplex_t qz2 = temp_z * temp_z;
            scomplex_t qzn = temp_z * nz;
            scomplex_t qzt = temp_z * z;
  
            for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
                ++ i_y, ++ global_i_y) {
              unsigned int xy_off_start = curr_nqx * i_y;
              real_t temp_y = qy[global_i_y];
              real_t qy2 = temp_y * temp_y;
              real_t qyn = temp_y * ny;
              real_t qyt = temp_y * y;
  
              for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
                  ++ i_x, ++ global_i_x) {
                unsigned int off = off_start + xy_off_start + i_x;
                real_t temp_x = qx[global_i_x];
                scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
                scomplex_t qt = temp_x * x + qyt + qzt;
                scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

                fq_buffer[off] = compute_fq(s, qt, qn);
              } // for z
            } // for y
          } // for x
        } // for t
      } // pragma omp parallel
    } // pragma offload
  } // NumericFormFactorM::form_factor_kernel()
  
  
  /**
   * Computational kernel function.
   */

  // single precision
  __attribute__((target(mic:0)))
  inline float2_t NumericFormFactorM::compute_fq(float s, float2_t qt, float2_t qn) {
    float2_t v1 = qn * make_sC(cosf(qt.x), sinf(qt.x));
    float v2 = s * exp(qt.y);
    return v1 * v2;
  } // NumericFormFactorM::compute_fq()
  

  // double precision
  __attribute__((target(mic:0)))
  inline double2_t NumericFormFactorM::compute_fq(double s, double2_t qt, double2_t qn) {
    double2_t v1 = qn * make_sC(cos(qt.x), sin(qt.x));
    double v2 = s * exp(qt.y);
    return v1 * v2;
  } // NumericFormFactorM::compute_fq()
  
  
  /**
   * Reduction kernel
   */
  void NumericFormFactorM::reduction_kernel(
      unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
      unsigned int curr_b_num_triangles, unsigned int blocked_matrix_size,
      unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
      unsigned int num_triangles, unsigned int nqx, unsigned int nqy, unsigned int nqz,
      unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
      scomplex_t* fq_buffer, scomplex_t* ff_buffer) {
    if(fq_buffer == NULL || ff_buffer == NULL) return;
    
    #pragma offload target(mic:0) \
            in(fq_buffer: length(0) MIC_REUSE) \
            out(ff_buffer: length(curr_b_nqx * curr_b_nqy * curr_b_nqz) MIC_REUSE)
    {
      unsigned int curr_b_nqxyz = curr_b_nqx * curr_b_nqy * curr_b_nqz;
      unsigned int curr_b_nqxy = curr_b_nqx * curr_b_nqy;
      scomplex_t temp_complex = make_sC((real_t) 0.0, (real_t) -1.0);

      omp_set_num_threads(MIC_OMP_NUM_THREADS_);
  
      // reduction over all triangles
      for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
        #pragma omp parallel
        {
          #pragma omp for nowait //schedule(auto)
          for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
            for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
              scomplex_t total = make_sC((real_t) 0.0, (real_t) 0.0);
              for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
                unsigned int i_fq = curr_b_nqxyz * i_t + curr_b_nqxy * i_z +
                          curr_b_nqx * i_y + i_x;
                total = total + fq_buffer[i_fq];
              } // for i_t
              unsigned int i_ff = curr_b_nqxy * i_z + curr_b_nqx * i_y + i_x;
              //ff_buffer[i_ff] = total * temp_complex;
              ff_buffer[i_ff] = make_sC(total.y, - total.x);
            } // for i_z
          } // for i_y
        } // pragma omp parallel
      } // for i_x
    } // pragma offload
  } // NumericFormFactorM::reduction_kernel()
  
  
  /**
   * The main host function called from outside, as part of the API for a single node.
   * Double buffered version.
   */
  unsigned int NumericFormFactorM::compute_form_factor_db(int rank,
            real_vec_t &shape_def_vec, complex_t* &ff,
            real_t* qx, int nqx, real_t* qy, int nqy, complex_t* qz, int nqz,
            real_t* rot,
            real_t& pass_kernel_time, real_t& red_time, real_t& mem_time
            #ifdef FINDBLOCK
              , const int block_x, const int block_y, const int block_z, const int block_t
            #endif
            ) {
    double kernel_time = 0.0, reduce_time = 0.0, total_kernel_time = 0.0, total_reduce_time = 0.0,
        temp_mem_time = 0.0, total_mem_time = 0.0;
    woo::BoostChronoTimer kerneltimer;

    #ifdef _OPENMP
      if(rank == 0)
        std::cout << "++        OpenMP threads on host: " << omp_get_max_threads() << std::endl;
    #endif

    int num_mic = 0;
    #ifdef __INTEL_OFFLOAD
      num_mic = _Offload_number_of_devices();
      if(rank == 0) {
        std::cout << "++         Number of Target MICs: " << num_mic << std::endl;
      } // if
      if(num_mic == 0) {
        std::cerr << "error: no Target MIC found!" << std::endl;
        return 0;
      } // if
    #else
      std::cerr << "warning: offloading to MIC not set! using host." << std::endl;
    #endif

    unsigned int num_triangles = shape_def_vec.size() / 7;
    if(num_triangles < 1) return 0;

    real_t* shape_def = &shape_def_vec[0];
  
    unsigned int total_qpoints = nqx * nqy * nqz;
  
    // allocate memory for the final FF 3D matrix
    ff = new (std::nothrow) complex_t[total_qpoints];  // allocate and initialize to 0
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << total_qpoints * sizeof(complex_t) << " b" << std::endl;
      return 0;
    } // if
    memset(ff, 0, total_qpoints * sizeof(complex_t));
    unsigned long int host_mem_usage = total_qpoints * sizeof(complex_t);

    //scomplex_t* qz_flat = new (std::nothrow) scomplex_t[nqz];
    scomplex_t* qz_flat = (scomplex_t*) _mm_malloc(nqz * sizeof(scomplex_t), 64);
    if(qz_flat == NULL) {
      std::cerr << "Memory allocation failed for qz_flat. Size = "
            << nqz * sizeof(scomplex_t) << " b" << std::endl;
      return 0;
    } // if
    for(int i = 0; i < nqz; ++ i) {
      //qz_flat[i].x = qz[i].real(); qz_flat[i].y = qz[i].imag();
      qz_flat[i] = make_sC(qz[i].real(), qz[i].imag());
    } // for

    host_mem_usage += ((unsigned int) nqx + nqy) * sizeof(real_t) +  // qx, qy
              (unsigned int) nqz * sizeof(complex_t) +        // qz
              (unsigned int) nqz * sizeof(scomplex_t) +       // qz_flat
              (unsigned int) shape_def_vec.size() * sizeof(real_t) + // shape_def
              total_qpoints * sizeof(complex_t);              // ff

    // allocate memory buffers on the target and transfer data asynchronously
    #pragma offload_transfer target(mic:0) \
                in(qx: length(nqx) MIC_ALLOC) \
                in(qy: length(nqy) MIC_ALLOC) \
                in(qz_flat: length(nqz) MIC_ALLOC) \
                in(shape_def: length(7 * num_triangles) MIC_ALLOC) \
                in(rot: length(9) MIC_ALLOC) \
                signal(shape_def)
    
    unsigned long int target_mem_usage = ((unsigned int) nqx + nqy) * sizeof(real_t) +
                      (unsigned int) nqz * sizeof(scomplex_t) +
                      (unsigned int) num_triangles * 7 * sizeof(real_t);

    unsigned int matrix_size = (unsigned int) nqx * nqy * nqz * num_triangles;

    // do hyperblocking to use less memory and computation - mem transfer overlap
    unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
    compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
                b_nqx, b_nqy, b_nqz, b_num_triangles
                #ifdef FINDBLOCK
                  , block_x, block_y, block_z, block_t
                #endif
                );
    // compute the number of sub-blocks, along each of the 4 dimensions
    // formulate loops over each dimension, to go over each sub block
    unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
    unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
    unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
    unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);
  
    unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
    unsigned int curr_b_num_triangles = b_num_triangles;
    unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;
    // for the progress indicator
    unsigned int progress_delta = 10; // display progress at 10% intervals
    unsigned int progress_step = (unsigned int) ceil((float) num_blocks / progress_delta);
  
    if(rank == 0) {
      std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
            << " x " << b_nqz << " x " << b_num_triangles << std::endl;
      std::cout << "++         Number of hyperblocks: " << num_blocks
            << " [" << nb_x << " x " << nb_y << " x " << nb_z << " x " << nb_t << "]"
            << std::endl;
    } // if
  
    unsigned int blocked_3d_matrix_size = (unsigned int) b_nqx * b_nqy * b_nqz;
    unsigned long int blocked_matrix_size = (unsigned int) blocked_3d_matrix_size * b_num_triangles;

    // allocate ff and fq buffers on host first (double buffering)
    #ifndef FF_NUM_MIC_SWAP
      scomplex_t *fq_buffer0, *fq_buffer1;
    #endif
    scomplex_t *ff_buffer0, *ff_buffer1;
    //fq_buffer0 = new (std::nothrow) scomplex_t[blocked_matrix_size]();
    //fq_buffer1 = new (std::nothrow) scomplex_t[blocked_matrix_size]();
    //ff_buffer0 = new (std::nothrow) scomplex_t[blocked_3d_matrix_size]();
    //ff_buffer1 = new (std::nothrow) scomplex_t[blocked_3d_matrix_size]();
    #ifdef MIC_PADDING
      // to make buffers multiple of 64B:
      unsigned int pad_fq = (unsigned int) ceil(blocked_matrix_size * sizeof(scomplex_t) / 64.0f) *
                          64 - blocked_matrix_size * sizeof(scomplex_t);
      unsigned int pad_ff = (unsigned int) ceil(blocked_3d_matrix_size * sizeof(scomplex_t) / 64.0f) *
                          64 - blocked_3d_matrix_size * sizeof(scomplex_t);
      unsigned int fq_bytes = blocked_3d_matrix_size * sizeof(scomplex_t) + pad_fq;
      unsigned int ff_bytes = blocked_matrix_size * sizeof(scomplex_t) + pad_ff;
      unsigned int fq_size = fq_bytes / sizeof(scomplex_t);
      unsigned int ff_size = ff_bytes / sizeof(scomplex_t);
      if(rank == 0) {
        std::cout << "++                Buffer padding: " << pad_ff << " B" << std::endl;
      } // if
      #ifndef FF_NUM_MIC_SWAP
        fq_buffer0 = (scomplex_t*) _mm_malloc(fq_bytes, 64);
        fq_buffer1 = (scomplex_t*) _mm_malloc(fq_bytes, 64);
      #endif
      ff_buffer0 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
      ff_buffer1 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
      // FIXME: mem usage ...
      host_mem_usage += 2 * (fq_bytes + ff_bytes);
      target_mem_usage += 2 * (fq_bytes + ff_bytes);
    #else
      #ifndef FF_NUM_MIC_SWAP
        fq_buffer0 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_matrix_size, 64);
        fq_buffer1 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_matrix_size, 64);
      #endif
      ff_buffer0 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
      ff_buffer1 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
      // FIXME: mem usage ...
      host_mem_usage += 2 * (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(scomplex_t);
      target_mem_usage += 2 * (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(scomplex_t);
    #endif
    if(ff_buffer0 == NULL || ff_buffer1 == NULL
      #ifndef FF_NUM_MIC_SWAP
        || fq_buffer0 == NULL || fq_buffer1 == NULL
      #endif
      ) {
      std::cerr << "Memory allocation failed for f buffers. blocked_matrix_size = "
            << blocked_matrix_size << std::endl
            << "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
            << std::endl;
      return 0;
    } // if

    if(rank == 0) {
      std::cout << "++             Host memory usage: " << (float) host_mem_usage / 1024 / 1024
            << " MB" << std::endl << std::flush;
      #ifdef __INTEL_OFFLOAD
        std::cout << "++           Target memory usage: " << (float) target_mem_usage / 1024 / 1024
              << " MB" << std::endl << std::flush;
      #endif
    } // if

    // allocate ff and fq buffers on target
    #ifdef MIC_PADDING
      #ifndef FF_NUM_MIC_SWAP
        #pragma offload_transfer target(mic:0) \
              nocopy(fq_buffer0: length(fq_size) MIC_ALLOC) \
              nocopy(fq_buffer1: length(fq_size) MIC_ALLOC) \
              nocopy(ff_buffer0: length(ff_size) MIC_ALLOC) \
              nocopy(ff_buffer1: length(ff_size) MIC_ALLOC)
      #else
        #pragma offload_transfer target(mic:0) \
              nocopy(ff_buffer0: length(ff_size) MIC_ALLOC) \
              nocopy(ff_buffer1: length(ff_size) MIC_ALLOC)
      #endif
    #else
      #ifndef FF_NUM_MIC_SWAP
        #pragma offload_transfer target(mic:0) \
              nocopy(fq_buffer0: length(blocked_matrix_size) MIC_ALLOC) \
              nocopy(fq_buffer1: length(blocked_matrix_size) MIC_ALLOC) \
              nocopy(ff_buffer0: length(blocked_3d_matrix_size) MIC_ALLOC) \
              nocopy(ff_buffer1: length(blocked_3d_matrix_size) MIC_ALLOC)
      #else
        #pragma offload_transfer target(mic:0) \
              nocopy(ff_buffer0: length(blocked_3d_matrix_size) MIC_ALLOC) \
              nocopy(ff_buffer1: length(blocked_3d_matrix_size) MIC_ALLOC)
      #endif
    #endif

    if(rank == 0) {
      #ifdef __INTEL_OFFLOAD
        int num = 0;
        #pragma offload target(mic:0)
        num = get_target_omp_num_threads();
        std::cout << "-- OpenMP threads on MIC/process: " << num << std::endl;
        std::cout << "-- Computing form factor on MIC (DB) ... " << std::flush;
      #else
        std::cout << "-- Computing form factor on CPU (fallback) ... " << std::flush;
      #endif
    } // if

    int curr_buffer_i = 0;  // buffer index
    unsigned int prev_b_nqx, prev_b_nqy, prev_b_nqz, prev_b_num_triangles;
    int prev_ib_x, prev_ib_y, prev_ib_z, prev_ib_t;

    unsigned int hblock_counter = 0;

    // wait for input transfer to complete
    #pragma offload_wait target(mic:0) wait(shape_def)

    kerneltimer.start();

    // compute for each hyperblock
    curr_b_nqx = b_nqx;
    for(int ib_x = 0; ib_x < nb_x; ++ ib_x) {
      if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
      curr_b_nqy = b_nqy;
      for(int ib_y = 0; ib_y < nb_y; ++ ib_y) {
        if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
        curr_b_nqz = b_nqz;
        for(int ib_z = 0; ib_z < nb_z; ++ ib_z) {
          if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
          curr_b_num_triangles = b_num_triangles;
          for(int ib_t = 0; ib_t < nb_t; ++ ib_t) {
            if(ib_t == nb_t - 1)
              curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

            unsigned int curr_b_nqxyz = curr_b_nqx * curr_b_nqy * curr_b_nqz;

            // double buffering
            if(curr_buffer_i == 0) {
              // call the main kernel offloaded to MIC
              #ifndef FF_MIC_OPT
                #pragma offload target(mic:0) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(fq_buffer0: length(0) MIC_REUSE)
                form_factor_kernel_db(qx, qy, qz_flat, shape_def,
                    curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                    b_nqx, b_nqy, b_nqz, b_num_triangles,
                    ib_x, ib_y, ib_z, ib_t,
                    fq_buffer0);

                // call the reduction kernel offloaded to MIC
                #ifdef MIC_PADDING
                  #pragma offload target(mic:0) signal(ff_buffer0) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(fq_buffer0: length(0) MIC_REUSE) \
                      out(ff_buffer0: length(ff_size) MIC_REUSE)
                #else
                  #pragma offload target(mic:0) signal(ff_buffer0) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(fq_buffer0: length(0) MIC_REUSE) \
                      out(ff_buffer0: length(curr_b_nqxyz) MIC_REUSE)
                #endif
                reduction_kernel_db(curr_b_nqx, curr_b_nqy, curr_b_nqz,
                    curr_b_num_triangles, blocked_matrix_size,
                    b_nqx, b_nqy, b_nqz, num_triangles,
                    nqx, nqy, nqz,
                    ib_x, ib_y, ib_z, ib_t,
                    fq_buffer0, ff_buffer0);

              #else  // FF_MIC_OPT

                #ifdef MIC_PADDING
                  #ifndef FF_NUM_MIC_SWAP
                    #pragma offload target(mic:0) signal(ff_buffer0) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(fq_buffer0: length(0) MIC_REUSE) \
                      out(ff_buffer0: length(ff_size) MIC_REUSE)
                  #else
                    #pragma offload target(mic:0) signal(ff_buffer0) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(ff_buffer0: length(ff_size) MIC_REUSE)
                  #endif
                #else
                  #ifndef FF_NUM_MIC_SWAP
                    #pragma offload target(mic:0) signal(ff_buffer0) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(fq_buffer0: length(0) MIC_REUSE) \
                      out(ff_buffer0: length(curr_b_nqxyz) MIC_REUSE)
                  #else
                    #pragma offload target(mic:0) signal(ff_buffer0) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(ff_buffer0: length(curr_b_nqxyz) MIC_REUSE)
                  #endif
                #endif
                #ifndef FF_NUM_MIC_SWAP
                  form_factor_kernel_opt(qx, qy, qz_flat, shape_def,
                      curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                      b_nqx, b_nqy, b_nqz, b_num_triangles,
                      ib_x, ib_y, ib_z, ib_t,
                      rot,
                      fq_buffer0, ff_buffer0);
                #else
                  form_factor_kernel_loopswap(qx, qy, qz_flat, shape_def,
                      curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                      b_nqx, b_nqy, b_nqz, b_num_triangles,
                      nqx, nqy, nqz, num_triangles,
                      ib_x, ib_y, ib_z, ib_t,
                      rot,
                      ff_buffer0);
                #endif
              #endif

              if(ib_x + ib_y + ib_z + ib_t != 0) {
                // wait for transfer of 1 to finish before moving
                #pragma offload_wait target(mic:0) wait(ff_buffer1)

                // move computed ff block from buffer to final ff
                move_to_main_ff(ff_buffer1, prev_b_nqx, prev_b_nqy, prev_b_nqz,
                        b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                        prev_ib_x, prev_ib_y, prev_ib_z, ff);
              } // if
            } else {
              // call the main kernel offloaded to MIC
              #ifndef FF_MIC_OPT
                #pragma offload target(mic:0) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(fq_buffer1: length(0) MIC_REUSE)
                form_factor_kernel_db(qx, qy, qz_flat, shape_def,
                    curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                    b_nqx, b_nqy, b_nqz, b_num_triangles,
                    ib_x, ib_y, ib_z, ib_t,
                    fq_buffer1);

                // call the reduction kernel offloaded to MIC
                #ifdef MIC_PADDING
                  #pragma offload target(mic:0) signal(ff_buffer1) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(fq_buffer1: length(0) MIC_REUSE) \
                      out(ff_buffer1: length(ff_size) MIC_REUSE)
                #else
                  #pragma offload target(mic:0) signal(ff_buffer1) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(fq_buffer1: length(0) MIC_REUSE) \
                      out(ff_buffer1: length(curr_b_nqxyz) MIC_REUSE)
                #endif
                reduction_kernel_db(curr_b_nqx, curr_b_nqy, curr_b_nqz,
                    curr_b_num_triangles, blocked_matrix_size,
                    b_nqx, b_nqy, b_nqz, num_triangles,
                    nqx, nqy, nqz,
                    ib_x, ib_y, ib_z, ib_t,
                    fq_buffer1, ff_buffer1);

              #else  // FF_MIC_OPT

                #ifdef MIC_PADDING
                  #ifndef FF_NUM_MIC_SWAP
                    #pragma offload target(mic:0) signal(ff_buffer1) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(fq_buffer1: length(0) MIC_REUSE) \
                      out(ff_buffer1: length(ff_size) MIC_REUSE)
                  #else
                    #pragma offload target(mic:0) signal(ff_buffer1) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(ff_buffer1: length(ff_size) MIC_REUSE)
                  #endif
                #else
                  #ifndef FF_NUM_MIC_SWAP
                    #pragma offload target(mic:0) signal(ff_buffer1) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(fq_buffer1: length(0) MIC_REUSE) \
                      out(ff_buffer1: length(curr_b_nqxyz) MIC_REUSE)
                  #else
                    #pragma offload target(mic:0) signal(ff_buffer1) \
                      in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                      in(b_nqx, b_nqy, b_nqz, num_triangles) \
                      in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                      in(shape_def: length(0) MIC_REUSE) \
                      in(rot: length(0) MIC_REUSE) \
                      in(qx: length(0) MIC_REUSE) \
                      in(qy: length(0) MIC_REUSE) \
                      in(qz_flat: length(0) MIC_REUSE) \
                      out(ff_buffer1: length(curr_b_nqxyz) MIC_REUSE)
                  #endif
                #endif
                #ifndef FF_NUM_MIC_SWAP
                  form_factor_kernel_opt(qx, qy, qz_flat, shape_def,
                      curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                      b_nqx, b_nqy, b_nqz, b_num_triangles,
                      ib_x, ib_y, ib_z, ib_t,
                      rot,
                      fq_buffer1, ff_buffer1);
                #else
                  form_factor_kernel_loopswap(qx, qy, qz_flat, shape_def,
                      curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                      b_nqx, b_nqy, b_nqz, b_num_triangles,
                      nqx, nqy, nqz, num_triangles,
                      ib_x, ib_y, ib_z, ib_t,
                      rot,
                      ff_buffer1);
                #endif
              #endif

              if(ib_x + ib_y + ib_z + ib_t != 0) {
                // wait for transfer of 0 to finish before moving
                #pragma offload_wait target(mic:0) wait(ff_buffer0)

                // move computed ff block from buffer to final ff
                move_to_main_ff(ff_buffer0, prev_b_nqx, prev_b_nqy, prev_b_nqz,
                        b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                        prev_ib_x, prev_ib_y, prev_ib_z, ff);
              } // if
            } // if-else

            // flip
            curr_buffer_i = 1 - curr_buffer_i;
            prev_b_nqx = curr_b_nqx; prev_b_nqy = curr_b_nqy; prev_b_nqz = curr_b_nqz;
            prev_b_num_triangles = curr_b_num_triangles;
            prev_ib_x = ib_x; prev_ib_y = ib_y; prev_ib_z = ib_z; prev_ib_t = ib_t;

            ++ hblock_counter;

            if(rank == 0) {
              float progress = ((float) hblock_counter / num_blocks) * 100;
              double intpart, fracpart;
              fracpart = modf(progress, &intpart);
              if(((unsigned int)intpart) % progress_delta == 0 && fracpart < 100.0 / num_blocks)
                std::cout << intpart << "% ";
            } // if
          } // for ib_t
        } // for ib_z
      } // for ib_y
    } // for ib_x

    // move the last part
    if(curr_buffer_i == 0) {
      // wait for transfer of 1 to finish before moving
      #pragma offload_wait target(mic:0) wait(ff_buffer1)

      // move computed ff block from buffer to final ff
      move_to_main_ff(ff_buffer1, prev_b_nqx, prev_b_nqy, prev_b_nqz,
              b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
              prev_ib_x, prev_ib_y, prev_ib_z, ff);
    } else {
      // wait for transfer of 0 to finish before moving
      #pragma offload_wait target(mic:0) wait(ff_buffer0)

      // move computed ff block from buffer to final ff
      move_to_main_ff(ff_buffer0, prev_b_nqx, prev_b_nqy, prev_b_nqz,
              b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
              prev_ib_x, prev_ib_y, prev_ib_z, ff);
    } // if-else

    kerneltimer.stop();

    // free f buffers on target
    #ifndef FF_NUM_MIC_SWAP
      #pragma offload_transfer target(mic:0) \
          nocopy(fq_buffer0: length(0) MIC_FREE) \
          nocopy(fq_buffer1: length(0) MIC_FREE) \
          nocopy(ff_buffer0: length(0) MIC_FREE) \
          nocopy(ff_buffer1: length(0) MIC_FREE)
    #else
      #pragma offload_transfer target(mic:0) \
          nocopy(ff_buffer0: length(0) MIC_FREE) \
          nocopy(ff_buffer1: length(0) MIC_FREE)
    #endif
    // and host
    //delete[] ff_buffer0;
    //delete[] ff_buffer1;
    //delete[] fq_buffer0;
    //delete[] fq_buffer1;
    _mm_free(ff_buffer0);
    _mm_free(ff_buffer1);
    #ifndef FF_NUM_MIC_SWAP
      _mm_free(fq_buffer0);
      _mm_free(fq_buffer1);
    #endif

    #pragma offload_transfer target(mic:0) \
        nocopy(shape_def: length(0) MIC_FREE) \
        nocopy(rot: length(0) MIC_FREE) \
        nocopy(qx: length(0) MIC_FREE) \
        nocopy(qy: length(0) MIC_FREE) \
        nocopy(qz_flat: length(0) MIC_FREE)
  
    //delete[] qz_flat;
    _mm_free(qz_flat);

    if(rank == 0) std::cout << "done." << std::endl;

    pass_kernel_time = kerneltimer.elapsed_msec();
    red_time = total_reduce_time;
    mem_time = total_mem_time;
  
    return num_triangles;
  } // NumericFormFactorM::compute_form_factor_db()
    

  /**
   * The main host function called from outside, as part of the API for a single node.
   * K-buffered version. This uses only the "loopswap" kernel.
   */
  unsigned int NumericFormFactorM::compute_form_factor_kb(int rank,
            real_t* shape_def, unsigned int num_triangles, complex_t* &ff,
            real_t* qx, int nqx, real_t* qy, int nqy, complex_t* qz, int nqz, int k,
            real_t* rot,
            real_t& kernel_time, real_t& red_time, real_t& mem_time
            #ifdef FINDBLOCK
              , const int block_x, const int block_y, const int block_z, const int block_t
            #endif
            ) {
    k = 3;     // this is only triple-buffering for now
          // (not sure how to implement k version due to target memory allocation style)

    double reduce_time = 0.0, total_reduce_time = 0.0,
        temp_mem_time = 0.0, total_mem_time = 0.0;
    woo::BoostChronoTimer kerneltimer;

    #ifdef _OPENMP
      if(rank == 0)
        std::cout << "++        OpenMP threads on host: " << omp_get_max_threads() << std::endl;
    #endif

    int num_mic = 0;
    #ifdef __INTEL_OFFLOAD
      num_mic = _Offload_number_of_devices();
      if(rank == 0) {
        std::cout << "++         Number of Target MICs: " << num_mic << std::endl;
      } // if
      if(num_mic == 0) {
        std::cerr << "error: no Target MIC found!" << std::endl;
        return 0;
      } // if
    #else
      std::cerr << "warning: offloading to MIC not set! using host." << std::endl;
    #endif

    if(num_triangles < 1) return 0;

    unsigned long int total_qpoints = nqx * nqy * nqz;
    unsigned long int host_mem_usage = 0;
  
    // allocate memory for the final FF 3D matrix
    ff = new (std::nothrow) complex_t[total_qpoints];  // allocate and initialize to 0
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << total_qpoints * sizeof(complex_t) << " b" << std::endl;
      return 0;
    } // if
    memset(ff, 0, total_qpoints * sizeof(complex_t));
    host_mem_usage += total_qpoints * sizeof(complex_t);

    scomplex_t* qz_flat = (scomplex_t*) _mm_malloc(nqz * sizeof(scomplex_t), 64);
    if(qz_flat == NULL) {
      std::cerr << "Memory allocation failed for qz_flat. Size = "
            << nqz * sizeof(scomplex_t) << " b" << std::endl;
      return 0;
    } // if
    for(int i = 0; i < nqz; ++ i) qz_flat[i] = make_sC(qz[i].real(), qz[i].imag());

    // FIXME:
    host_mem_usage += ((unsigned int) nqx + nqy) * sizeof(real_t) +  // qx, qy
              (unsigned int) nqz * sizeof(complex_t) +        // qz
              (unsigned int) nqz * sizeof(scomplex_t) +       // qz_flat
              (unsigned int) num_triangles * 8 * sizeof(real_t) + // shape_def
              total_qpoints * sizeof(complex_t);              // ff

    // size of shape_def:
    // each of 7 components is aligned to 64 bytes with padding at end
    unsigned int shape_padding = (16 - (num_triangles & 15)) & 15;
    unsigned int shape_size = (num_triangles + shape_padding) * 7;

    // allocate memory buffers on the target and transfer data asynchronously
    #pragma offload_transfer target(mic:0) \
                in(qx: length(nqx) MIC_ALLOC) \
                in(qy: length(nqy) MIC_ALLOC) \
                in(qz_flat: length(nqz) MIC_ALLOC) \
                in(rot: length(9) MIC_ALLOC) \
                in(shape_def: length(shape_size) MIC_ALLOC) \
                signal(shape_def)

    // FIXME:    
    unsigned long int target_mem_usage = ((unsigned int) nqx + nqy) * sizeof(real_t) +
                      (unsigned int) nqz * sizeof(scomplex_t) +
                      (unsigned int) num_triangles * 7 * sizeof(real_t);

    unsigned int matrix_size = (unsigned int) nqx * nqy * nqz * num_triangles;

    // do hyperblocking to use less memory and computation - mem transfer overlap
    unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
    #ifndef FF_NUM_MIC_AUTOTUNE_HB
      compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
                  b_nqx, b_nqy, b_nqz, b_num_triangles
                  #ifdef FINDBLOCK
                    , block_x, block_y, block_z, block_t
                  #endif
                  );
    #else
      // Autotuning to find optimal hyperblock size ... TODO: improve ...

      std::cout << "-- Autotuning hyperblock size ... " << std::endl;
      double min_time_hb = 100000000.0;
      unsigned int min_b_nqx = 1, min_b_nqy = 1, min_b_nqz = 1, min_b_num_triangles = 1;
      woo::BoostChronoTimer at_kernel_timer, at_overhead_timer;
      at_overhead_timer.start();
      scomplex_t* ff_temp;
      ff_temp = new (std::nothrow) scomplex_t[nqx * nqy * nqz];
      #pragma offload_transfer target(mic:0) in(ff_temp: length(nqx * nqy * nqz) MIC_ALLOC)

      for(unsigned int b_nqx_i = 1; b_nqx_i <= nqx; ++ b_nqx_i) {
        for(unsigned int b_nqy_i = 10; b_nqy_i <= nqy; b_nqy_i += 10) {
          for(unsigned int b_nqz_i = 10; b_nqz_i <= nqz; b_nqz_i += 10) {
            for(unsigned int b_nt_i = 100; b_nt_i <= num_triangles; b_nt_i += 100) {
              at_kernel_timer.start();

              // compute the number of sub-blocks, along each of the 4 dimensions
              unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx_i);
              unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy_i);
              unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz_i);
              unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_nt_i);
              unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;

              unsigned int ff_size = b_nqx_i * b_nqy_i * b_nqz_i;

              // call the main kernel offloaded to MIC
              #pragma offload target(mic:0) \
                in(b_nqx, b_nqy, b_nqz, num_triangles) \
                in(nqx, nqy, nqz) \
                in(shape_def: length(0) MIC_REUSE) \
                in(rot: length(0) MIC_REUSE) \
                in(qx: length(0) MIC_REUSE) \
                in(qy: length(0) MIC_REUSE) \
                in(qz_flat: length(0) MIC_REUSE) \
                nocopy(ff_temp: length(0) MIC_REUSE)
              form_factor_kernel_loopswap(qx, qy, qz_flat, shape_def,
                  b_nqx, b_nqy, b_nqz, b_num_triangles,
                  b_nqx, b_nqy, b_nqz, b_num_triangles,
                  nqx, nqy, nqz, num_triangles,
                  0, 0, 0, 0,
                  rot,
                  ff_temp);
  
              at_kernel_timer.stop();
              double curr_time = at_kernel_timer.elapsed_msec();
              double tot_time = curr_time * num_blocks;
              std::cout << "## " << b_nqx_i << " x " << b_nqy_i << " x " << b_nqz_i
                    << " x " << b_nt_i << "\t" << num_blocks << " : "
                    << curr_time << "\t" << tot_time << std::endl;
              if(tot_time < min_time_hb) {
                min_time_hb = tot_time;
                min_b_nqx = b_nqx_i; min_b_nqy = b_nqy_i; min_b_nqz = b_nqz_i;
                min_b_num_triangles = b_nt_i;
              } // if

            } // for
          } // for
        } // for
      } // for

      #pragma offload_transfer target(mic:0) nocopy(ff_temp: length(0) MIC_FREE)
      delete[] ff_temp;
      at_overhead_timer.stop();

      b_nqx = min_b_nqx; b_nqy = min_b_nqy; b_nqz = min_b_nqz; b_num_triangles = min_b_num_triangles;
      if(rank == 0) {
        std::cout << "##    HBlock Autotuner overhead: " << at_overhead_timer.elapsed_msec()
              << " ms." << std::endl;
      } // if
    #endif

    // compute the number of sub-blocks, along each of the 4 dimensions
    // formulate loops over each dimension, to go over each sub block
    unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
    unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
    unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
    unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);
  
    unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
    unsigned int curr_b_num_triangles = b_num_triangles;
    unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;
    // for the progress indicator
    unsigned int progress_delta = 10; // display progress at 10% intervals
    unsigned int progress_step = (unsigned int) ceil((float) num_blocks / progress_delta);
  
    if(rank == 0) {
      std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
            << " x " << b_nqz << " x " << b_num_triangles << std::endl;
      std::cout << "++         Number of hyperblocks: " << num_blocks
            << " [" << nb_x << " x " << nb_y << " x " << nb_z << " x " << nb_t << "]"
            << std::endl;
    } // if
  
    unsigned int blocked_3d_matrix_size = (unsigned int) b_nqx * b_nqy * b_nqz;
    unsigned long int blocked_matrix_size = (unsigned int) blocked_3d_matrix_size * b_num_triangles;

    // allocate ff buffers on host first (k buffering)
    scomplex_t *ff_buffer0, *ff_buffer1, *ff_buffer2;
    #ifdef MIC_PADDING  // padded ff and fq
      // to make buffers multiple of 64B:
      unsigned int pad_fq = (unsigned int) ceil(blocked_matrix_size * sizeof(scomplex_t) / 64.0f) *
                          64 - blocked_matrix_size * sizeof(scomplex_t);
      unsigned int pad_ff = (unsigned int) ceil(blocked_3d_matrix_size * sizeof(scomplex_t) / 64.0f) *
                          64 - blocked_3d_matrix_size * sizeof(scomplex_t);
      unsigned int fq_bytes = blocked_3d_matrix_size * sizeof(scomplex_t) + pad_fq;
      unsigned int ff_bytes = blocked_matrix_size * sizeof(scomplex_t) + pad_ff;
      unsigned int fq_size = fq_bytes / sizeof(scomplex_t);
      unsigned int ff_size = ff_bytes / sizeof(scomplex_t);
      if(rank == 0) {
        std::cout << "++                Buffer padding: " << pad_ff << " B" << std::endl;
      } // if
      ff_buffer0 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
      ff_buffer1 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
      ff_buffer2 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
      // FIXME: mem usage ...
      host_mem_usage += 2 * (fq_bytes + ff_bytes);
      target_mem_usage += 2 * (fq_bytes + ff_bytes);
    #else
      ff_buffer0 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
      ff_buffer1 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
      ff_buffer2 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
      // FIXME: mem usage ...
      host_mem_usage += 2 * (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(scomplex_t);
      target_mem_usage += 2 * (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(scomplex_t);
    #endif
    if(ff_buffer0 == NULL || ff_buffer1 == NULL || ff_buffer2 == NULL) {
      std::cerr << "Memory allocation failed for ff buffers. blocked_3d_matrix_size = "
            << blocked_3d_matrix_size << std::endl
            << "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
            << std::endl;
      return 0;
    } // if

    if(rank == 0) {
      std::cout << "++             Host memory usage: " << (float) host_mem_usage / 1024 / 1024
            << " MB" << std::endl << std::flush;
      #ifdef __INTEL_OFFLOAD
        std::cout << "++           Target memory usage: " << (float) target_mem_usage / 1024 / 1024
              << " MB" << std::endl << std::flush;
      #endif
    } // if

    #ifdef PROFILE_PAPI
      long long papi_total_cycles = 0, papi_total_inst = 0, papi_total_iis = 0, papi_total_stall = 0;
      double overall_ipc = 0.0;
    #endif

    // allocate ff and fq buffers on target
    #ifdef MIC_PADDING
      #pragma offload_transfer target(mic:0) \
            nocopy(ff_buffer0: length(ff_size) MIC_ALLOC) \
            nocopy(ff_buffer1: length(ff_size) MIC_ALLOC) \
            nocopy(ff_buffer2: length(ff_size) MIC_ALLOC)
    #else
      #pragma offload_transfer target(mic:0) \
            nocopy(ff_buffer0: length(blocked_3d_matrix_size) MIC_ALLOC) \
            nocopy(ff_buffer1: length(blocked_3d_matrix_size) MIC_ALLOC) \
            nocopy(ff_buffer2: length(blocked_3d_matrix_size) MIC_ALLOC)
    #endif

    if(rank == 0) {
      #ifdef __INTEL_OFFLOAD
        int num = 0;
        #pragma offload target(mic:0)
        num = get_target_omp_num_threads();
        std::cout << "-- OpenMP threads on MIC/process: " << num << std::endl;
        std::cout << "-- Computing form factor on MIC (KB) ... " << std::flush;
      #else
        std::cout << "-- Computing form factor on CPU (fallback) ... " << std::flush;
      #endif
    } // if

    int curr_buffer_i = 0;  // buffer index 0, 1 or 2
    int active = 0, passive;
    int prev_ib_x = 0, prev_ib_y = 0, prev_ib_z = 0, prev_ib_t = 0;
    unsigned int prev_cbnqx = curr_b_nqx, prev_cbnqy = curr_b_nqy, prev_cbnqz = curr_b_nqz;
    unsigned int mid_ib_x[k - 2], mid_ib_y[k - 2], mid_ib_z[k - 2];
    unsigned int mid_cbnqx[k - 2], mid_cbnqy[k - 2], mid_cbnqz[k - 2];
    for(int i_k = 0; i_k < k - 2; ++ i_k) {
      mid_ib_x[i_k] = 0; mid_ib_y[i_k] = 0; mid_ib_z[i_k] = 0;
      mid_cbnqx[i_k] = curr_b_nqx; mid_cbnqy[i_k] = curr_b_nqy; mid_cbnqz[i_k] = curr_b_nqz;
    } // for

    unsigned int hblock_counter = 0;

    // wait for input transfer to complete
    #pragma offload_wait target(mic:0) wait(shape_def)

    kerneltimer.start();

    // compute for each hyperblock
    curr_b_nqx = b_nqx;

    for(int ib_x = 0; ib_x < nb_x; ++ ib_x) {
      if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
      curr_b_nqy = b_nqy;

      for(int ib_y = 0; ib_y < nb_y; ++ ib_y) {
        if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
        curr_b_nqz = b_nqz;

        for(int ib_z = 0; ib_z < nb_z; ++ ib_z) {
          if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
          curr_b_num_triangles = b_num_triangles;

          for(int ib_t = 0; ib_t < nb_t; ++ ib_t) {
            if(ib_t == nb_t - 1)
              curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

            unsigned int curr_b_nqxyz = curr_b_nqx * curr_b_nqy * curr_b_nqz;

            #ifdef PROFILE_PAPI
              int papi_events[4] = { PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_TOT_IIS, PAPI_RES_STL };
              long long  papi_counter_values[4];
              PAPI_start_counters(papi_events, 4);
            #endif

            //  triple buffering
            passive = (active + 1) % k;

            switch(active) {
              case 0:

                // call the main kernel offloaded to MIC
                #ifdef MIC_PADDING
                  #pragma offload target(mic:0) signal(ff_buffer0) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(rot: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(ff_buffer0: length(ff_size) MIC_REUSE)
                #else
                  #pragma offload target(mic:0) signal(ff_buffer0) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(rot: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(ff_buffer0: length(curr_b_nqxyz) MIC_REUSE)
                #endif
                form_factor_kernel_loopswap(qx, qy, qz_flat, shape_def,
                    curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                    b_nqx, b_nqy, b_nqz, b_num_triangles,
                    nqx, nqy, nqz, num_triangles,
                    ib_x, ib_y, ib_z, ib_t,
                    rot,
                    ff_buffer0);
  
                if(hblock_counter > k - 2) {
                  // wait for transfer of 1 to finish before moving
                  #pragma offload_wait target(mic:0) wait(ff_buffer1)
  
                  // move computed ff block from buffer to final ff -- BLOCKING CALL!
                  move_to_main_ff(ff_buffer1, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                          b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                          prev_ib_x, prev_ib_y, prev_ib_z, ff);
                } // if

                break;

              case 1:

                // call the main kernel offloaded to MIC
                #ifdef MIC_PADDING
                  #pragma offload target(mic:0) signal(ff_buffer1) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(rot: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(ff_buffer1: length(ff_size) MIC_REUSE)
                #else
                  #pragma offload target(mic:0) signal(ff_buffer1) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(rot: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(ff_buffer1: length(curr_b_nqxyz) MIC_REUSE)
                #endif
                form_factor_kernel_loopswap(qx, qy, qz_flat, shape_def,
                    curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                    b_nqx, b_nqy, b_nqz, b_num_triangles,
                    nqx, nqy, nqz, num_triangles,
                    ib_x, ib_y, ib_z, ib_t,
                    rot,
                    ff_buffer1);

                if(hblock_counter > k - 2) {
                  // wait for transfer of 0 to finish before moving
                  #pragma offload_wait target(mic:0) wait(ff_buffer2)

                  // move computed ff block from buffer to final ff -- BLOCKING CALL!
                  move_to_main_ff(ff_buffer2, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                          b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                          prev_ib_x, prev_ib_y, prev_ib_z, ff);
                } // if

                break;

              case 2:

                // call the main kernel offloaded to MIC
                #ifdef MIC_PADDING
                  #pragma offload target(mic:0) signal(ff_buffer2) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(rot: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(ff_buffer2: length(ff_size) MIC_REUSE)
                #else
                  #pragma offload target(mic:0) signal(ff_buffer2) \
                    in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
                    in(b_nqx, b_nqy, b_nqz, num_triangles) \
                    in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
                    in(shape_def: length(0) MIC_REUSE) \
                    in(rot: length(0) MIC_REUSE) \
                    in(qx: length(0) MIC_REUSE) \
                    in(qy: length(0) MIC_REUSE) \
                    in(qz_flat: length(0) MIC_REUSE) \
                    out(ff_buffer2: length(curr_b_nqxyz) MIC_REUSE)
                #endif
                form_factor_kernel_loopswap(qx, qy, qz_flat, shape_def,
                    curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
                    b_nqx, b_nqy, b_nqz, b_num_triangles,
                    nqx, nqy, nqz, num_triangles,
                    ib_x, ib_y, ib_z, ib_t,
                    rot,
                    ff_buffer2);

                if(hblock_counter > k - 2) {
                  // wait for transfer of 0 to finish before moving
                  #pragma offload_wait target(mic:0) wait(ff_buffer0)

                  // move computed ff block from buffer to final ff
                  move_to_main_ff(ff_buffer0, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                          b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                          prev_ib_x, prev_ib_y, prev_ib_z, ff);
                } // if

                break;

              default:

                if(rank == 0) {
                  std::cerr << "error: something went terribly wrong in k buffering"
                        << std::endl;
                  return 0;
                } // if
            } // switch

            #ifdef PROFILE_PAPI
              PAPI_stop_counters(papi_counter_values, 4);
              papi_total_cycles += papi_counter_values[0];
              papi_total_inst += papi_counter_values[1];
              papi_total_iis += papi_counter_values[2];
              papi_total_stall += papi_counter_values[3];
            #endif

            active = (active + 1) % k;

            prev_ib_x = mid_ib_x[0]; prev_ib_y = mid_ib_y[0]; prev_ib_z = mid_ib_z[0];
            prev_cbnqx = mid_cbnqx[0]; prev_cbnqy = mid_cbnqy[0]; prev_cbnqz = mid_cbnqz[0];
            for(int i_k = 0; i_k < k - 3; ++ i_k) {
              mid_ib_x[i_k] = mid_ib_x[i_k + 1]; mid_ib_y[i_k] = mid_ib_y[i_k + 1];
              mid_ib_z[i_k] = mid_ib_z[i_k + 1];
              mid_cbnqx[i_k] = mid_cbnqx[i_k + 1]; mid_cbnqy[i_k] = mid_cbnqy[i_k + 1];
              mid_cbnqz[i_k] = mid_cbnqz[i_k + 1];
            } // for
            mid_ib_x[k - 3] = ib_x; mid_ib_y[k - 3] = ib_y; mid_ib_z[k - 3] = ib_z;
            mid_cbnqx[k - 3] = curr_b_nqx; mid_cbnqy[k - 3] = curr_b_nqy;
            mid_cbnqz[k - 3] = curr_b_nqz;

            ++ hblock_counter;

            if(rank == 0) {
              float progress = ((float) hblock_counter / num_blocks) * 100;
              double intpart, fracpart;
              fracpart = modf(progress, &intpart);
              if(((unsigned int)intpart) % progress_delta == 0 && fracpart < 100.0 / num_blocks)
                std::cout << intpart << "% " << std::flush;
            } // if
          } // for ib_t
        } // for ib_z
      } // for ib_y
    } // for ib_x

    // move the last k - 1 parts

    passive = (active + 1) % k;
    switch(active) {
      case 0:
        #pragma offload_wait target(mic:0) wait(ff_buffer1)
        move_to_main_ff(ff_buffer1, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                prev_ib_x, prev_ib_y, prev_ib_z, ff);
        break;
      case 1:
        #pragma offload_wait target(mic:0) wait(ff_buffer2)
        move_to_main_ff(ff_buffer2, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                prev_ib_x, prev_ib_y, prev_ib_z, ff);
        break;
      case 2:
        #pragma offload_wait target(mic:0) wait(ff_buffer0)
        move_to_main_ff(ff_buffer0, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                prev_ib_x, prev_ib_y, prev_ib_z, ff);
        break;
      default:
        if(rank == 0) {
          std::cerr << "error: something went terribly wrong in k buffering"
                << std::endl;
          return 0;
        } // if
    } // switch
    active = (active + 1) % k;
    prev_ib_x = mid_ib_x[0]; prev_ib_y = mid_ib_y[0]; prev_ib_z = mid_ib_z[0];
    prev_cbnqx = mid_cbnqx[0]; prev_cbnqy = mid_cbnqy[0]; prev_cbnqz = mid_cbnqz[0];
    for(int i_k = 0; i_k < k - 3; ++ i_k) {
      mid_ib_x[i_k] = mid_ib_x[i_k + 1]; mid_ib_y[i_k] = mid_ib_y[i_k + 1];
      mid_ib_z[i_k] = mid_ib_z[i_k + 1];
      mid_cbnqx[i_k] = mid_cbnqx[i_k + 1]; mid_cbnqy[i_k] = mid_cbnqy[i_k + 1];
      mid_cbnqz[i_k] = mid_cbnqz[i_k + 1];
    } // for
    mid_ib_x[k - 3] = nb_x - (k - 1); mid_ib_y[k - 3] = nb_y - (k - 1);
    mid_ib_z[k - 3] = nb_z - (k - 1);
    mid_cbnqx[k - 3] = curr_b_nqx; mid_cbnqy[k - 3] = curr_b_nqy;
    mid_cbnqz[k - 3] = curr_b_nqz;

    passive = (active + 1) % k;
    switch(active) {
      case 0:
        #pragma offload_wait target(mic:0) wait(ff_buffer1)
        move_to_main_ff(ff_buffer1, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                prev_ib_x, prev_ib_y, prev_ib_z, ff);
        break;
      case 1:
        #pragma offload_wait target(mic:0) wait(ff_buffer2)
        move_to_main_ff(ff_buffer2, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                prev_ib_x, prev_ib_y, prev_ib_z, ff);
        break;
      case 2:
        #pragma offload_wait target(mic:0) wait(ff_buffer0)
        move_to_main_ff(ff_buffer0, prev_cbnqx, prev_cbnqy, prev_cbnqz,
                b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
                prev_ib_x, prev_ib_y, prev_ib_z, ff);
        break;
      default:
        if(rank == 0) {
          std::cerr << "error: something went terribly wrong in k buffering"
                << std::endl;
          return 0;
        } // if
    } // switch

    kerneltimer.stop();

    // free f buffers on target
    #pragma offload_transfer target(mic:0) \
        nocopy(ff_buffer0: length(0) MIC_FREE) \
        nocopy(ff_buffer1: length(0) MIC_FREE) \
        nocopy(ff_buffer2: length(0) MIC_FREE)
    // and host
    _mm_free(ff_buffer0);
    _mm_free(ff_buffer1);
    _mm_free(ff_buffer2);

    #pragma offload_transfer target(mic:0) \
        nocopy(shape_def: length(0) MIC_FREE) \
        nocopy(rot: length(0) MIC_FREE) \
        nocopy(qx: length(0) MIC_FREE) \
        nocopy(qy: length(0) MIC_FREE) \
        nocopy(qz_flat: length(0) MIC_FREE)
  
    _mm_free(qz_flat);

    if(rank == 0) std::cout << "done." << std::endl;

    kernel_time = kerneltimer.elapsed_msec();

    #ifdef PROFILE_PAPI
    if(rank == 0) {
      std::cout << "++                  PAPI_TOT_CYC: " << papi_total_cycles << std::endl;
      std::cout << "++                  PAPI_TOT_INS: " << papi_total_inst << std::endl;
      std::cout << "++                  PAPI_TOT_IIS: " << papi_total_iis << std::endl;
      std::cout << "++                  PAPI_RES_STL: " << papi_total_stall << std::endl;
      std::cout << "++                           IPC: "
            << (double) papi_total_inst / papi_total_cycles << std::endl;
      //std::cout << "++             Instruction Reply: "
      //      << (double) (papi_total_iis - papi_total_inst) * 100 / papi_total_iis
      //      << " %" << std::endl;
      std::cout << "++                Stalled Cycles: "
            << (double) papi_total_stall * 100 / papi_total_cycles
            << " %" << std::endl;
    } // if
    #endif

    //if(rank == 0) {
    //  std::cout << "**                FF kernel time: " << kernel_time << " ms."
    //        << std::endl;
    //} // if

    red_time = total_reduce_time;
    mem_time = total_mem_time;
  
    return num_triangles;
  } // NumericFormFactorM::compute_form_factor_kb()
    

  /***
   * the following kernels are used in double buffering case
   * TODO ... clean this
   */

  __attribute__((target(mic:0)))
  void NumericFormFactorM::form_factor_kernel_db(real_t* qx, real_t* qy, scomplex_t* qz_flat,
          real_t* shape_def,
          unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
          unsigned int curr_num_triangles,
          unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
          unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
          scomplex_t* fq_buffer) {

    //#ifdef __INTEL_OFFLOAD
    //  omp_set_num_threads(MIC_OMP_NUM_THREADS_);
    //#endif
    #pragma omp parallel
    {
      #pragma omp for nowait //schedule(auto)
      for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
        unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
        real_t s = shape_def[shape_off];
        real_t nx = shape_def[shape_off + 1];
        real_t ny = shape_def[shape_off + 2];
        real_t nz = shape_def[shape_off + 3];
        real_t x = shape_def[shape_off + 4];
        real_t y = shape_def[shape_off + 5];
        real_t z = shape_def[shape_off + 6];

        unsigned int xy_size = curr_nqx * curr_nqy;
        unsigned int matrix_off = xy_size * curr_nqz * i_t;
        unsigned int start_z = b_nqz * ib_z;
        unsigned int start_y = b_nqy * ib_y;
        unsigned int start_x = b_nqx * ib_x;

        for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
            ++ i_z, ++ global_i_z) {
          unsigned int off_start = matrix_off + xy_size * i_z;
          scomplex_t temp_z = qz_flat[global_i_z];
          scomplex_t qz2 = temp_z * temp_z;
          scomplex_t qzn = temp_z * nz;
          scomplex_t qzt = temp_z * z;

          for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
              ++ i_y, ++ global_i_y) {
            unsigned int xy_off_start = curr_nqx * i_y;
            real_t temp_y = qy[global_i_y];
            real_t qy2 = temp_y * temp_y;
            real_t qyn = temp_y * ny;
            real_t qyt = temp_y * y;

            // ... not vectorized by compiler ... TODO
            for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
                ++ i_x, ++ global_i_x) {
              unsigned int off = off_start + xy_off_start + i_x;
              real_t temp_x = qx[global_i_x];
              scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
              scomplex_t qt = temp_x * x + qyt + qzt;
              scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

              fq_buffer[off] = compute_fq(s, qt, qn);
            } // for z
          } // for y
        } // for x
      } // for t
    } // pragma omp parallel
  } // NumericFormFactorM::form_factor_kernel_db()
  
  
  __attribute__((target(mic:0)))
  void NumericFormFactorM::reduction_kernel_db(
        unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
        unsigned int curr_b_num_triangles, unsigned int blocked_matrix_size,
        unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
        unsigned int num_triangles, unsigned int nqx, unsigned int nqy, unsigned int nqz,
        unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
        scomplex_t* fq_buffer, scomplex_t* ff_buffer) {

    unsigned int curr_b_nqxyz = curr_b_nqx * curr_b_nqy * curr_b_nqz;
    unsigned int curr_b_nqxy = curr_b_nqx * curr_b_nqy;
    scomplex_t temp_complex = make_sC((real_t) 0.0, (real_t) -1.0);

    //#ifdef __INTEL_OFFLOAD
    //  omp_set_num_threads(MIC_OMP_NUM_THREADS_);
    //#endif

    // reduction over all triangles
    for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
      #pragma omp parallel
      {
        #pragma omp for nowait //schedule(auto)
        for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
          for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
            scomplex_t total = make_sC((real_t) 0.0, (real_t) 0.0);
            // ... not vectorized by compiler ... TODO
            for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
              unsigned int i_fq = curr_b_nqxyz * i_t + curr_b_nqxy * i_z +
                        curr_b_nqx * i_y + i_x;
              total = total + fq_buffer[i_fq];
            } // for i_t
            unsigned int i_ff = curr_b_nqxy * i_z + curr_b_nqx * i_y + i_x;
            //ff_buffer[i_ff] = total * temp_complex;
            ff_buffer[i_ff] = make_sC(total.y, - total.x);
          } // for i_z
        } // for i_y
      } // pragma omp parallel
    } // for i_x
  } // NumericFormFactorM::reduction_kernel_db()
  
  
  /***
   * the following kernels are with optimizations (trials)
   * TODO ... clean this
   */

  #ifndef FF_NUM_MIC_SWAP

  __attribute__((target(mic:0)))
  void NumericFormFactorM::form_factor_kernel_opt(real_t* qx, real_t* qy, scomplex_t* qz_flat,
          real_t* shape_def,
          unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
          unsigned int curr_num_triangles,
          unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
          unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
          scomplex_t* fq_buffer, scomplex_t* ff_buffer) {

    unsigned int curr_nqxy = curr_nqx * curr_nqy;
    #pragma omp parallel
    {
      #pragma omp for nowait //schedule(auto)
      for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
        unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
        real_t s = shape_def[shape_off];
        real_t nx = shape_def[shape_off + 1];
        real_t ny = shape_def[shape_off + 2];
        real_t nz = shape_def[shape_off + 3];
        real_t x = shape_def[shape_off + 4];
        real_t y = shape_def[shape_off + 5];
        real_t z = shape_def[shape_off + 6];

        unsigned int xy_size = curr_nqxy;
        unsigned int matrix_off = xy_size * curr_nqz * i_t;
        unsigned int start_z = b_nqz * ib_z;
        unsigned int start_y = b_nqy * ib_y;
        unsigned int start_x = b_nqx * ib_x;

        for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
            ++ i_z, ++ global_i_z) {
          unsigned int off_start = matrix_off + xy_size * i_z;
          scomplex_t temp_z = qz_flat[global_i_z];
          scomplex_t qz2 = temp_z * temp_z;
          scomplex_t qzn = temp_z * nz;
          scomplex_t qzt = temp_z * z;

          for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
              ++ i_y, ++ global_i_y) {
            unsigned int xy_off_start = curr_nqx * i_y;
            real_t temp_y = qy[global_i_y];
            real_t qy2 = temp_y * temp_y;
            real_t qyn = temp_y * ny;
            real_t qyt = temp_y * y;

            // ... not vectorized by compiler ... TODO (-ipo vectorizes)
            for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
                ++ i_x, ++ global_i_x) {
              unsigned int off = off_start + xy_off_start + i_x;
              real_t temp_x = qx[global_i_x];
              scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
              scomplex_t qt = temp_x * x + qyt + qzt;
              scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

              fq_buffer[off] = compute_fq(s, qt, qn);
              /*unsigned int i_ff = curr_nqxy * i_z + curr_nqx * i_y + i_x;
              scomplex_t res = compute_fq(s, qt, qn);
              res = make_sC(res.y, - res.x);
              #pragma omp atomic
              ff_buffer[i_ff].x += res.x;
              #pragma omp atomic
              ff_buffer[i_ff].y += res.y;*/
            } // for z
          } // for y
        } // for x
      } // for t
    } // pragma omp parallel
    
    // doing reduction here itself
    // a separate reduction helps avoid use of atomics otherwise (which are SLOW!)
    unsigned int curr_nqxyz = curr_nqx * curr_nqy * curr_nqz;
    //scomplex_t temp_complex = make_sC((real_t) 0.0, (real_t) -1.0);
    #pragma omp parallel
    {
      #pragma omp for nowait collapse(3) //schedule(auto)
      for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
        for(unsigned int i_y = 0; i_y < curr_nqy; ++ i_y) {
          for(unsigned int i_z = 0; i_z < curr_nqz; ++ i_z) {
            scomplex_t total = make_sC((real_t) 0.0, (real_t) 0.0);
            // ... not vectorized by compiler ... TODO
            for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
              unsigned int i_fq = curr_nqxyz * i_t + curr_nqxy * i_z +
                        curr_nqx * i_y + i_x;
              total = total + fq_buffer[i_fq];
            } // for i_t
            unsigned int i_ff = curr_nqxy * i_z + curr_nqx * i_y + i_x;
            //ff_buffer[i_ff] = total * temp_complex;
            ff_buffer[i_ff] = make_sC(total.y, - total.x);
          } // for i_z
        } // for i_y
      } // for i_x
    } // pragma omp parallel
  } // NumericFormFactorM::form_factor_kernel_opt()

  #else  
  
  __attribute__((target(mic:0)))
  void NumericFormFactorM::form_factor_kernel_loopswap(
          real_t* qx, real_t* qy, scomplex_t* qz_flat,
          real_t* shape_def,
          unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
          unsigned int curr_num_triangles,
          unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
          unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
          unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
          real_t* rot,
          scomplex_t* ff_buffer) {

    if(nqx == 1) {  // call the optimized special case
      #ifdef __MIC__
        form_factor_kernel_loopswap_vec_nqx1(qx, qy, qz_flat, shape_def,
                      curr_nqx, curr_nqy, curr_nqz, curr_num_triangles,
                      b_nqx, b_nqy, b_nqz, b_num_triangles,
                      nqx, nqy, nqz, num_triangles,
                      ib_x, ib_y, ib_z, ib_t,
                      rot,
                      ff_buffer);
      #else
        form_factor_kernel_loopswap_nqx1(qx, qy, qz_flat, shape_def,
                      curr_nqx, curr_nqy, curr_nqz, curr_num_triangles,
                      b_nqx, b_nqy, b_nqz, b_num_triangles,
                      nqx, nqy, nqz, num_triangles,
                      ib_x, ib_y, ib_z, ib_t,
                      rot,
                      ff_buffer);
      #endif
      return;
    } // if

    const unsigned int curr_nqxy = curr_nqx * curr_nqy;
    const unsigned int curr_nqxyz = curr_nqxy * curr_nqz;

    const unsigned int xy_size = curr_nqxy;
    const unsigned int start_z = b_nqz * ib_z;
    const unsigned int start_y = b_nqy * ib_y;
    const unsigned int start_x = b_nqx * ib_x;

    #pragma omp parallel
    {
      #pragma omp for collapse(3) nowait
      for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
        for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
          for(int i_x = 0; i_x < curr_nqx; ++ i_x) {

            int global_i_z = start_z + i_z;
            int global_i_y = start_y + i_y;
            int global_i_x = start_x + i_x;
            scomplex_t total = make_sC((real_t) 0.0, (real_t) 0.0);

            // TODO: do blocking for cache ... ?
            for(int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
              unsigned int matrix_off = xy_size * curr_nqz * i_t;
              unsigned int off_start = matrix_off + xy_size * i_z;
              unsigned int xy_off_start = curr_nqx * i_y;
              unsigned int off = off_start + xy_off_start + i_x;

              unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
              real_t s = shape_def[shape_off];
              real_t nx = shape_def[shape_off + 1];
              real_t ny = shape_def[shape_off + 2];
              real_t nz = shape_def[shape_off + 3];
              real_t x = shape_def[shape_off + 4];
              real_t y = shape_def[shape_off + 5];
              real_t z = shape_def[shape_off + 6];

              scomplex_t temp_qz = qz_flat[global_i_z];
              real_t temp_qy = qy[global_i_y];
              real_t temp_qx = qx[global_i_x];

              scomplex_t temp_x = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
              scomplex_t temp_y = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
              scomplex_t temp_z = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;

              scomplex_t qz2 = temp_z * temp_z;
              scomplex_t qzn = temp_z * nz;
              scomplex_t qzt = temp_z * z;

              scomplex_t qy2 = temp_y * temp_y;
              scomplex_t qyn = temp_y * ny;
              scomplex_t qyt = temp_y * y;

              scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
              scomplex_t qt = temp_x * x + qyt + qzt;
              scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

              total = total + compute_fq(s, qt, qn);
            } // for t

            unsigned int i_ff = curr_nqxy * i_z + curr_nqx * i_y + i_x;
            ff_buffer[i_ff] = make_sC(total.y, - total.x);
          } // for x
        } // for y
      } // for z
    } // pragma omp parallel
  } // NumericFormFactorM::form_factor_kernel_loopswap()
  
  #endif // FF_NUM_MIC_SWAP


  // include optimized kernels for MIC
  #include "ff_num_mic_kernels.hpp"

  
  /**
   * Function to move a computed block of ff to its right place in ff.
   */
  void NumericFormFactorM::move_to_main_ff(scomplex_t* ff_buffer,
        unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
        unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
        unsigned int nqx, unsigned int nqy, unsigned int nqz,
        unsigned int ib_x, unsigned int ib_y, unsigned int ib_z,
        complex_t* ff) {
    unsigned int temp1 = nqx * nqy;
    unsigned int temp2 = curr_b_nqx * curr_b_nqy;
    unsigned int base_i = nqx * nqy * ib_z * b_nqz + nqx * ib_y * b_nqy + ib_x * b_nqx;
    // make these copy operations contiguous ... (very low priority)
    for(int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
      unsigned int start_i = base_i + temp1 * i_z;
      for(int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
        unsigned int super_i = start_i + nqx * i_y;
        for(int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
          unsigned int final_i = super_i + i_x;
          unsigned int block_i = temp2 * i_z + curr_b_nqx * i_y + i_x;
          ff[final_i] += complex_t(ff_buffer[block_i].x, ff_buffer[block_i].y);
        } // for i_x
      } // for i_y
    } // for i_z
  } // move_to_main_ff()
  
  
  /**
   * Function to compute the decomposition block size
   * TODO: Improve it later ...
   */
  void NumericFormFactorM::compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
      unsigned int& b_nqx, unsigned int& b_nqy, unsigned int& b_nqz, unsigned int& b_num_triangles
      #ifdef FINDBLOCK
        , const int block_x, const int block_y, const int block_z, const int block_t
      #endif
      ) {
    b_nqx = (unsigned int) nqx; b_nqy = (unsigned int) nqy; b_nqz = (unsigned int) nqz;
    b_num_triangles = (unsigned int) num_triangles;
  
    #ifdef FINDBLOCK
      b_nqx = (b_nqx > block_x) ? block_x : b_nqx;
      b_nqy = (b_nqy > block_y) ? block_y : b_nqy;
      b_nqz = (b_nqz > block_z) ? block_z : b_nqz;
      b_num_triangles = (b_num_triangles > block_t) ? block_t : b_num_triangles;
    #else
      b_nqx = (b_nqx > MIC_BLOCK_X_) ? MIC_BLOCK_X_ : b_nqx;
      b_nqy = (b_nqy > MIC_BLOCK_Y_) ? MIC_BLOCK_Y_ : b_nqy;
      b_nqz = (b_nqz > MIC_BLOCK_Z_) ? MIC_BLOCK_Z_ : b_nqz;
      b_num_triangles = (b_num_triangles > MIC_BLOCK_T_) ? MIC_BLOCK_T_ : b_num_triangles;
    #endif
  } // NumericFormFactorM::compute_hyperblock_size()


} // namespace hig
 : b_num_triangles;
    #endif
  } // NumericFormFactorM::compute_hyperblock_size()


} // namespace hig
g
