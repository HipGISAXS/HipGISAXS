/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
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

//#ifdef INTEL_SB_AVX
//#include <numerics/cpu/avx_numerics.hpp>
//#elif defined __SSE3__
//#include <numerics/cpu/sse3_numerics.hpp>
//#endif
#include <numerics/numeric_utils.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

#include <common/constants.hpp>
#include <common/parameters.hpp>
#include <common/cpu/parameters_cpu.hpp>

#include <ff/cpu/ff_num_cpu.hpp>
  
namespace hig {
  
  complex_t FormFactorTriangle(real_t qx, real_t qy, complex_t qz,
          RotMatrix_t & rot, triangle_t & tri) {
    complex_t ff = CMPLX_ZERO_;
    complex_t unitc = CMPLX_ONE_;
    complex_t n_unitc = CMPLX_MINUS_ONE_;

    // do the rotation
    std::vector<complex_t> mq = rot.rotate (qx, qy, qz);

    // calculate q^2
    real_t q_sqr = 0.;
    for(int i=0; i<3; i++)  q_sqr += std::norm(mq[i]);
          
    // form vertices
    std::vector<vector3_t> vertex;
    vertex.resize(3);
    vertex[0] = vector3_t(tri.v1[0], tri.v1[1], tri.v1[2]);
    vertex[1] = vector3_t(tri.v2[0], tri.v2[1], tri.v2[2]);
    vertex[2] = vector3_t(tri.v3[0], tri.v3[1], tri.v3[2]);

    // form edges
    std::vector<vector3_t> edge;
    edge.resize(3);
    edge[0] = vertex[1] - vertex[0];
    edge[1] = vertex[2] - vertex[1];
    edge[2] = vertex[0] - vertex[2];


    // calculate projection of n_t on q
    vector3_t n_t = cross(edge[0], edge[1]);
    real_t t_area = 0.5 * n_t.abs();
    n_t = n_t / n_t.abs();

    // dot(q, n_t)
    complex_t q_dot_nt = CMPLX_ZERO_;
    for (int i=0; i<3; i++) q_dot_nt += mq[i] * n_t[i];

    // proj_tq
    real_t proj_tq = q_sqr - std::norm(q_dot_nt);

    // CASE 1
    if (std::abs(proj_tq) < TINY_){
        complex_t q_dot_v = CMPLX_ZERO_;
        for (int i=0; i<3; i++) q_dot_v += mq[i] * tri.v1[i];
        // calculate form-factor (Case 1)
        ff = unitc * q_dot_nt * t_area / q_sqr * std::exp(n_unitc * q_dot_v);
    } else {
        // iterate of each edge to compute form-factor
      for (int e = 0; e < 3; e++) {

        // edge-normal
        vector3_t n_e = cross(edge[e], n_t);
        n_e = n_e / n_e.abs();

        // dot (q, n_e)
        complex_t q_dot_ne = CMPLX_ZERO_;
        for (int i=0; i<3; i++) q_dot_ne += mq[i] * n_e[i];

        // proj_eq
        real_t proj_eq = proj_tq - std::norm(q_dot_ne);
        // CASE 2
        if (std::abs(proj_eq) < TINY_){
        // q_dot_v
          complex_t q_dot_v = CMPLX_ZERO_;
          for (int i=0; i<3; i++) q_dot_v += mq[i] * vertex[e][i];
          real_t f0 = edge[e].abs() / (q_sqr * proj_tq);
          complex_t c0 = - q_dot_nt * q_dot_ne;
          complex_t c1 = std::exp(n_unitc * q_dot_v);
          ff += f0 * c0 * c1;
        } else {
          // CASE 3 (General case)
          // denominator
          real_t   f0 = q_sqr * proj_tq * proj_eq;

          // dot(q, v_a) vertex a
          complex_t q_dot_v = CMPLX_ZERO_;
          for (int i=0; i<3; i++) q_dot_v += mq[i] * vertex[e][i];

          // vertrex-normal a
          vector3_t n_v = edge[e] / edge[e].abs();

          // dot(q, n_v)
          complex_t q_dot_nv = CMPLX_ZERO_;
          for (int i=0; i<3; i++) q_dot_nv += mq[i] * n_v[i];

          // calculate contribution of vertex a
          complex_t c0 = n_unitc * q_dot_nt * q_dot_ne * q_dot_nv;
          complex_t c1 = std::exp(n_unitc * q_dot_v);
          ff +=  c0 * c1 / f0;

          // dot(q, v) the other vertex in the edge
          q_dot_v = CMPLX_ZERO_;
          int ep = (e+1)%3;
          for (int i=0; i<3; i++) q_dot_v += mq[i] * vertex[ep][i];

          // dot (q, n_v)
          q_dot_nv = CMPLX_ZERO_;
          for (int i=0; i<3; i++) q_dot_nv -= mq[i] * n_v[i];

          // calculate contribution of the other vertex
          c0 = n_unitc * q_dot_nt * q_dot_ne * q_dot_nv;
          c1 = std::exp(n_unitc * q_dot_v);
          ff += c0 * c1 / f0;
        }
      }
    }
    return ff;
  }


  /**
   * Exact integration
   */
  unsigned int NumericFormFactorC::compute_exact_triangle(
          triangle_t * shape_def, int num_triangles,
          complex_t* &ff,
          int nqy, real_t * qx, real_t * qy, int nqz, complex_t * qz,
          RotMatrix_t & rot, real_t & compute_time) {

    if(num_triangles < 1) return 0;
    unsigned long int total_qpoints = nqz;
  
    // allocate memory for the final FF 3D matrix
    ff = new (std::nothrow) complex_t[total_qpoints];  // allocate and initialize to 0
    memset(ff, 0, total_qpoints * sizeof(complex_t));
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << total_qpoints * sizeof(complex_t) << " b" << std::endl;
      return 0;
    } // if
 
    #ifdef DEBUG
    int n = 10;
    for (int i_t=0; i_t < num_triangles; i_t++){
        complex_t foo = FormFactorTriangle(qx[n], qy[n], qz[n], rot, shape_def[i_t]);
        std::cout << foo << std::endl;
    }
    #endif // DEBUG
   
#pragma omp parallel for
    for (int i_z = 0; i_z < nqz; i_z++) {
        int i_y = i_z % nqy; 
        complex_t ff_temp = CMPLX_ZERO_;
        for (int i_t = 0; i_t < num_triangles; i_t++) {
            ff_temp += FormFactorTriangle (qx[i_y], qy[i_y], qz[i_z], rot, shape_def[i_t]);
        }
        ff[i_z] = ff_temp;
    }
    return num_triangles;
  }

  /**
   * Approximated integration
   */
  unsigned int NumericFormFactorC::compute_approx_triangle(
          real_vec_t &shape_def,
          complex_t *& ff,
          int nqy, real_t * qx, real_t * qy, 
          int nqz, complex_t * qz, RotMatrix_t & rot, real_t &comp_time){

    int num_triangles = shape_def.size() / CPU_T_PROP_SIZE_;
    if (num_triangles < 1) return 0;

    ff = new (std::nothrow) complex_t[nqz];
    if (ff == NULL){
      std::cerr << "Memory allocation failed for ff. Requested size: " << nqz << std::endl;
      return 0;
    }
    memset(ff, 0, nqz * sizeof(complex_t));

    woo::BoostChronoTimer timer;
    timer.start();

#pragma omp parallel for
    for (int i_z = 0; i_z < nqz; i_z++){
      int i_y = i_z % nqy;
      for (int i_t = 0; i_t < num_triangles; i_t++){
        int offset = i_t * CPU_T_PROP_SIZE_;
        real_t s  = shape_def[offset];
        real_t nx = shape_def[offset + 1];
        real_t ny = shape_def[offset + 2];
        real_t nz = shape_def[offset + 3];
        real_t x  = shape_def[offset + 4];
        real_t y  = shape_def[offset + 5];
        real_t z  = shape_def[offset + 6];

        // rotate q-vector
        std::vector<complex_t> mq = rot.rotate(qx[i_y], qy[i_y], qz[i_z]);

        real_t q2 = std::norm(mq[0]) + std::norm(mq[1]) + std::norm(mq[2]);
        complex_t qn = mq[0] * nx + mq[1] * ny + mq[2] * nz;
        complex_t qt = mq[0] * x  + mq[1] * y  + mq[2] * z;
        complex_t nj = CMPLX_MINUS_ONE_;
        complex_t np = CMPLX_ONE_;
        ff[i_z] += (nj * qn * s * std::exp(np * qt) / q2);
      }
    }
    timer.stop();
    comp_time = timer.elapsed_msec();
    return num_triangles;
  }
} // namespace hig
