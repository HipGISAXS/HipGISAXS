/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_exact_kernels.cpp
 *  Created: Mar 22, 2016
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

#include <woo/timer/woo_boostchronotimers.hpp>
#include <numerics/numeric_utils.hpp>
#include <common/constants.hpp>
#include <common/parameters.hpp>
#include <common/cpu/parameters_cpu.hpp>
#include <ff/cpu/ff_num_cpu.hpp>

#ifdef INTEL_AVX
#include <numerics/cpu/avx_numerics.hpp>
#endif
  

namespace hig {

#ifdef INTEL_AVX
  
  complex_t NumericFormFactorC::form_factor_kernel_exact(real_t qx, real_t qy, complex_t qz,
                                                         RotMatrix_t& rot,
                                                         const triangle_t& tri) {
    std::vector<complex_t> mq = rot.rotate(qx, qy, qz);
    real_t q_sqr = 0.0;
    for(int i = 0; i < 3; ++ i) q_sqr += std::norm(mq[i]);              // calculate q^2
          
    // construct vertices
    std::vector<vector3_t> vertex;
    vertex.resize(3);
    vertex[0] = vector3_t(tri.v1[0], tri.v1[1], tri.v1[2]);
    vertex[1] = vector3_t(tri.v2[0], tri.v2[1], tri.v2[2]);
    vertex[2] = vector3_t(tri.v3[0], tri.v3[1], tri.v3[2]);

    // construct edges
    std::vector<vector3_t> edge;
    edge.resize(3);
    edge[0] = vertex[1] - vertex[0];
    edge[1] = vertex[2] - vertex[1];
    edge[2] = vertex[0] - vertex[2];

    // calculate projection of n_t on q
    vector3_t n_t = cross(edge[0], edge[1]);
    real_t t_area = 0.5 * n_t.abs();
    n_t = n_t / n_t.abs();

    complex_t q_dot_nt = CMPLX_ZERO_;
    for(int i = 0; i < 3; ++ i) q_dot_nt += mq[i] * n_t[i];             // dot(q, n_t)
    real_t proj_tq = q_sqr - std::norm(q_dot_nt);                       // proj_tq

    complex_t ff = CMPLX_ZERO_;
    if(std::abs(proj_tq) < TINY_) {               // case 1
      complex_t q_dot_v = CMPLX_ZERO_;
      for (int i=0; i<3; i++) q_dot_v += mq[i] * tri.v1[i];
      ff = CMPLX_ONE_ * q_dot_nt * t_area / q_sqr * std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
    } else {
      // iterate over edges to compute form-factor
      for(int e = 0; e < 3; ++ e) {
        vector3_t n_e = cross(edge[e], n_t);
        n_e = n_e / n_e.abs();                                          // edge normal
        complex_t q_dot_ne = CMPLX_ZERO_;
        for(int i = 0; i < 3; ++ i) q_dot_ne += mq[i] * n_e[i];         // dot (q, n_e)
        real_t proj_eq = proj_tq - std::norm(q_dot_ne);
        if(std::abs(proj_eq) < TINY_) {           // case 2
          complex_t q_dot_v = CMPLX_ZERO_;
          for(int i = 0; i < 3; ++ i) q_dot_v += mq[i] * vertex[e][i];
          real_t f0 = edge[e].abs() / (q_sqr * proj_tq);
          complex_t c0 = - q_dot_nt * q_dot_ne;
          complex_t c1 = std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
          ff += f0 * c0 * c1;
        } else {                                  // case 3 [general case]
          real_t f0 = q_sqr * proj_tq * proj_eq;                        // denominator
          complex_t q_dot_v = CMPLX_ZERO_;
          for(int i = 0; i < 3; ++ i) q_dot_v += mq[i] * vertex[e][i];  // dot(q, v_a) vertex a
          vector3_t n_v = edge[e] / edge[e].abs();                      // vertrex normal a
          complex_t q_dot_nv = CMPLX_ZERO_;
          for (int i=0; i<3; i++) q_dot_nv += mq[i] * n_v[i];           // dot(q, n_v)
          // calculate contribution of vertex a
          complex_t c0 = CMPLX_MINUS_ONE_ * q_dot_nt * q_dot_ne * q_dot_nv;
          complex_t c1 = std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
          ff +=  c0 * c1 / f0;
          q_dot_v = CMPLX_ZERO_;
          int ep = (e + 1) % 3;
          for(int i = 0; i < 3; ++ i) q_dot_v += mq[i] * vertex[ep][i]; // dot(q, v) vertex b
          q_dot_nv = CMPLX_ZERO_;
          for(int i = 0; i < 3; ++ i) q_dot_nv -= mq[i] * n_v[i];       // dot (q, n_v)
          // calculate contribution of vertex b
          c0 = CMPLX_MINUS_ONE_ * q_dot_nt * q_dot_ne * q_dot_nv;
          c1 = std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
          ff += c0 * c1 / f0;
        } // if-else
      } // for
    } // if-else

    return ff;
  } // NumericFormFactorC::form_factor_kernel_exact()


  /**
   * Exact integration: hand vectorized
   */
  unsigned int NumericFormFactorC::compute_exact_triangle_vec(const triangle_t* shape_def,
                                                              int num_triangles,
                                                              complex_t* &ff,
                                                              int nqy,
                                                              const real_t* qx, const real_t* qy,
                                                              int nqz, const complex_t* qz,
                                                              RotMatrix_t& rot,
                                                              real_t& compute_time) {
    if(num_triangles < 1) return 0;

    // TODO:
    // vectorization is on the level of triangles (triangle loop)
    // hence, qx, qy, qz dont need to be vectorized (unlike analytical ff)
    // steps:
    // convert shape_def to the vectorized data structure
    // loop over all q-points
    // loop over ceil(num_triangles/4) vector triangles
    // compute a vector temp_ff in each vector triangle iteration
    // sum vector temp_ff over all vector triangle iterations
    // reduce the vector temp_ff by summing all 4 entries into one scalar ff for a q-point

    unsigned long int total_qpoints = nqz;
    ff = new (std::nothrow) complex_t[total_qpoints];
    if(ff == NULL) {
      std::cerr << "error: memory allocation failed for ff of size "
                << total_qpoints * sizeof(complex_t) << " bytes" << std::endl;
      return 0;
    } // if
   
    woo::BoostChronoTimer timer;
    timer.start();

    #pragma omp parallel for schedule(runtime)
    for(int i_z = 0; i_z < nqz; ++ i_z) {
      int i_y = i_z % nqy; 
      complex_t ff_temp = CMPLX_ZERO_;
      for(int i_t = 0; i_t < num_triangles; ++ i_t) {
        ff_temp += form_factor_kernel_exact(qx[i_y], qy[i_y], qz[i_z], rot, shape_def[i_t]);
      } // for
      ff[i_z] = ff_temp;
    } // for

    timer.stop();
    compute_time = timer.elapsed_msec();

    return num_triangles;
  } // NumericFormFactorC::compute_exact_triangle()


#endif // INTEL_AVX

} // namespace hig
