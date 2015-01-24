/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
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
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PROFILE_PAPI
#include <papi.h>
#endif

#ifdef INTEL_SB_AVX
#include <numerics/cpu/avx_numerics.hpp>
#elif defined __SSE3__
#include <numerics/cpu/sse3_numerics.hpp>
#endif
#include <numerics/numeric_utils.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

#include <common/constants.hpp>
#include <common/parameters.hpp>
#include <common/cpu/parameters_cpu.hpp>

#include <ff/cpu/ff_num_cpu.hpp>
  
namespace hig {
  
  complex_t FormFactorTriangle(float_t qx, float_t qy, complex_t qz,
          float_t * rot, triangle_t & tri) {
    complex_t ff = C_ZERO;
    complex_t unitc = C_ONE;
    complex_t n_unitc = C_NEG_ONE;

    // do the rotation
    std::vector<complex_t> mq; mq.resize(3);
    mq[0] = rot[0] * qx + rot[1] * qy + rot[2] * qz;
    mq[1] = rot[3] * qx + rot[4] * qy + rot[5] * qz;
    mq[2] = rot[6] * qx + rot[7] * qy + rot[8] * qz;

    // calculate q^2
    float_t q_sqr = 0.;
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
    vector3_t n_t = cross_product(edge[0], edge[1]);
    float_t t_area = 0.5 * n_t.norm2();
    n_t = n_t / n_t.norm2();

    // dot(q, n_t)
    complex_t q_dot_nt = C_ZERO;
    for (int i=0; i<3; i++) q_dot_nt += mq[i] * n_t[i];

    // proj_tq
    float_t proj_tq = q_sqr - std::norm(q_dot_nt);

    // CASE 1
    if (std::abs(proj_tq) < 1.0E-06){
        complex_t q_dot_v = C_ZERO;
        for (int i=0; i<3; i++) q_dot_v += mq[i] * tri.v1[i];
        // calculate form-factor (Case 1)
        ff = unitc * q_dot_nt * t_area / q_sqr * std::exp(n_unitc * q_dot_v);
    } else {
        // iterate of each edge to compute form-factor
      for (int e = 0; e < 3; e++) {

        // edge-normal
        vector3_t n_e = cross_product(edge[e], n_t);
        n_e = n_e / n_e.norm2();

        // dot (q, n_e)
        complex_t q_dot_ne = C_ZERO;
        for (int i=0; i<3; i++) q_dot_ne += mq[i] * n_e[i];

        // proj_eq
        float_t proj_eq = proj_tq - std::norm(q_dot_ne);
        // CASE 2
        if (std::abs(proj_eq) < 1.0E-08) {
        // q_dot_v
          complex_t q_dot_v = C_ZERO;
          for (int i=0; i<3; i++) q_dot_v += mq[i] * vertex[e][i];
          float_t f0 = edge[e].norm2() / (q_sqr * proj_tq);
          complex_t c0 = - q_dot_nt * q_dot_ne;
          complex_t c1 = std::exp(n_unitc * q_dot_v);
          ff += f0 * c0 * c1;
        } else {
          // denominator
          float_t   f0 = q_sqr * proj_tq * proj_eq;

          // dot(q, v_a) vertex a
          complex_t q_dot_v = C_ZERO;
          for (int i=0; i<3; i++) q_dot_v += mq[i] * vertex[e][i];

          // vertrex-normal a
          vector3_t n_v = edge[e] / edge[e].norm2();

          // dot(q, n_v)
          complex_t q_dot_nv = C_ZERO;
          for (int i=0; i<3; i++) q_dot_nv += mq[i] * n_v[i];

          // calculate ff for vertex a
          complex_t c0 = n_unitc * q_dot_nt * q_dot_ne * q_dot_nv;
          complex_t c1 = std::exp(n_unitc * q_dot_v);
          ff +=  c0 * c1 / f0;

          // dot(q, v) vertex b
          q_dot_v = C_ZERO;
          int ep = (e+1)%3;
          for (int i=0; i<3; i++) q_dot_v += mq[i] * vertex[ep][i];

          // dot (q, n_v)
          q_dot_nv = C_ZERO;
          for (int i=0; i<3; i++) q_dot_nv += mq[i] * (NEG_ONE * n_v[i]);

          c0 = n_unitc * q_dot_nt * q_dot_ne * q_dot_nv;
          c1 = std::exp(n_unitc * q_dot_v);
          ff += c0 * c1 / f0;
        }
      }
    }
    return ff;
  }


  /**
   * The main host function called from outside, as part of the API for a single node.
   */
  unsigned int NumericFormFactorC::compute_ana_triangle (int rank,
            triangle_t * shape_def, int num_triangles,
            complex_t* &ff,
            float_t * qx, int nqx, float_t * qy, int nqy, complex_t * qz, int nqz,
            float_t * rot,
            float_t & compute_time) {

    double temp_mem_time = 0.0, total_mem_time = 0.0;
    #ifdef _OPENMP
      if(rank == 0)
        std::cout << "++      Number of OpenMP threads: " << omp_get_max_threads() << std::endl;
    #endif
  
    if(num_triangles < 1) return 0;
    unsigned long int total_qpoints = nqz;
    unsigned long int host_mem_usage = ((unsigned long int) nqx + nqy) * sizeof(float_t) +
                        nqz * sizeof(complex_t);
  
    // allocate memory for the final FF 3D matrix
    ff = new (std::nothrow) complex_t[total_qpoints];  // allocate and initialize to 0
    memset(ff, 0, total_qpoints * sizeof(complex_t));
    if(ff == NULL) {
      std::cerr << "Memory allocation failed for ff. Size = "
            << total_qpoints * sizeof(complex_t) << " b" << std::endl;
      return 0;
    } // if
    host_mem_usage += total_qpoints * sizeof(complex_t);
 
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
        complex_t ff_temp = C_ZERO;
        for (int i_t = 0; i_t < num_triangles; i_t++) {
            ff_temp += FormFactorTriangle (qx[i_y], qy[i_y], qz[i_z], rot, shape_def[i_t]);
        }
        ff[i_z] = ff_temp;
    }
    return num_triangles;
  }
  
} // namespace hig
