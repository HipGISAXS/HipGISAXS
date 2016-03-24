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
  
  avx_m256c_t NumericFormFactorC::form_factor_kernel_exact(real_t qx, real_t qy, complex_t qz,
                                                           RotMatrix_t& rot,
                                                           const avx_triangle_t& vt) {
    std::vector<complex_t> mq = rot.rotate(qx, qy, qz);
    real_t q_sqr = 0.0;
    for(int i = 0; i < 3; ++ i) q_sqr += std::norm(mq[i]);              // calculate q^2
          
    // construct vertices
    //std::vector<vector3_t> vertex;
    //vertex.resize(3);
    //vertex[0] = vector3_t(tri.v1[0], tri.v1[1], tri.v1[2]);
    //vertex[1] = vector3_t(tri.v2[0], tri.v2[1], tri.v2[2]);
    //vertex[2] = vector3_t(tri.v3[0], tri.v3[1], tri.v3[2]);
    // 3 avx vertices of avx triangle are vt.a, vt.b and vt.c
    avx_vertex_t vertex[3];
    vertex[0].x = vt.a.x; vertex[0].y = vt.a.y; vertex[0].z = vt.a.z;
    vertex[1].x = vt.b.x; vertex[1].y = vt.b.y; vertex[1].z = vt.b.z;
    vertex[2].x = vt.c.x; vertex[2].y = vt.c.y; vertex[2].z = vt.c.z;

    // construct edges
    //std::vector<vector3_t> edge;
    //edge.resize(3);
    //edge[0] = vertex[1] - vertex[0];
    //edge[1] = vertex[2] - vertex[1];
    //edge[2] = vertex[0] - vertex[2];
    // 3 avx edges are:
    avx_edge_t edge[3];
    edge[0].x = _mm256_sub_pd(vt.b.x, vt.a.x);
    edge[0].y = _mm256_sub_pd(vt.b.y, vt.a.y);
    edge[0].z = _mm256_sub_pd(vt.b.z, vt.a.z);
    edge[1].x = _mm256_sub_pd(vt.c.x, vt.b.x);
    edge[1].y = _mm256_sub_pd(vt.c.y, vt.b.y);
    edge[1].z = _mm256_sub_pd(vt.c.z, vt.b.z);
    edge[2].x = _mm256_sub_pd(vt.a.x, vt.c.x);
    edge[2].y = _mm256_sub_pd(vt.a.y, vt.c.y);
    edge[2].z = _mm256_sub_pd(vt.a.z, vt.c.z);

    // calculate projection of n_t on q
    //vector3_t n_t = cross(edge[0], edge[1]);
    //real_t t_area = 0.5 * n_t.abs();
    //n_t = n_t / n_t.abs();
    avx_edge_t n_t = avx_cross_epep(edge[0], edge[1]);
    avx_m256_t t_area = avx_mul_rrp(avx_abs_ep(n_t), 0.5);
    n_t = avx_div_eprp(n_t, avx_abs_ep(n_t));

    //complex_t q_dot_nt = CMPLX_ZERO_;
    //for(int i = 0; i < 3; ++ i) q_dot_nt += mq[i] * n_t[i];             // dot(q, n_t)
    //real_t proj_tq = q_sqr - std::norm(q_dot_nt);                       // proj_tq
    avx_m256c_t q_dot_nt = avx_dot_cep(mq, n_t);
    avx_m256_t proj_tq = avx_sub_rrp(q_sqr, avx_norm_cp(q_dot_nt));

    // TODO: handle mixed cases:
    avx_m256c_t ff = avx_setzero_cp();
    // ... avx_m256_t cmp = _mm256_cmp_pd(avx_abs_rp(proj_tq), AVX_EPSILON_, _CMP_GT_OS);
    int cond = 1;
    // ... int cond = _mm256_movemask_pd(cmp);
    // if(std::abs(proj_tq) < TINY_) {               // case 1
    if(cond == 0) {         // case 1 for when all are < epsilon
      // ...
      // complex_t q_dot_v = CMPLX_ZERO_;
      // for (int i=0; i<3; i++) q_dot_v += mq[i] * tri.v1[i];
      // ff = CMPLX_ONE_ * q_dot_nt * t_area / q_sqr * std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
    } else {
      // iterate over edges to compute form-factor
      for(auto ei = 0; ei < 3; ++ ei) {
        // vector3_t n_e = cross(edge[e], n_t);
        // n_e = n_e / n_e.abs();                                          // edge normal
        // complex_t q_dot_ne = CMPLX_ZERO_;
        // for(int i = 0; i < 3; ++ i) q_dot_ne += mq[i] * n_e[i];         // dot (q, n_e)
        //real_t proj_eq = proj_tq - std::norm(q_dot_ne);
        avx_edge_t n_e = avx_cross_epep(edge[ei], n_t);
        n_e = avx_div_eprp(n_e, avx_abs_ep(n_e));
        avx_m256c_t q_dot_ne = avx_dot_cep(mq, n_e);
        avx_m256_t proj_eq = avx_sub_rprp(proj_tq, avx_norm_cp(q_dot_ne));
        // ... avx_m256_t cmp = _mm256_cmp_pd(avx_abs_rp(proj_eq), AVX_EPSILON_, _CMP_GT_OS);
        int cond = 1;
        // ... int cond = _mm256_movemask_pd(cmp);
        // if(std::abs(proj_eq) < TINY_) {           // case 2
        if(cond == 0) {     // case 2 for when all are < epsilon
          // ...
          // complex_t q_dot_v = CMPLX_ZERO_;
          // for(int i = 0; i < 3; ++ i) q_dot_v += mq[i] * vertex[e][i];
          // real_t f0 = edge[e].abs() / (q_sqr * proj_tq);
          // complex_t c0 = - q_dot_nt * q_dot_ne;
          // complex_t c1 = std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
          // ff += f0 * c0 * c1;
        } else {            // case 3 the general case
          // real_t f0 = q_sqr * proj_tq * proj_eq;                        // denominator
          // complex_t q_dot_v = CMPLX_ZERO_;
          // for(int i = 0; i < 3; ++ i) q_dot_v += mq[i] * vertex[e][i];  // dot(q, v_a) vertex a
          // vector3_t n_v = edge[e] / edge[e].abs();                      // vertrex normal a
          // complex_t q_dot_nv = CMPLX_ZERO_;
          // for (int i=0; i<3; i++) q_dot_nv += mq[i] * n_v[i];           // dot(q, n_v)
          // calculate contribution of vertex a
          // complex_t c0 = CMPLX_MINUS_ONE_ * q_dot_nt * q_dot_ne * q_dot_nv;
          // complex_t c1 = std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
          // ff +=  c0 * c1 / f0;
          // q_dot_v = CMPLX_ZERO_;
          // int ep = (e + 1) % 3;
          // for(int i = 0; i < 3; ++ i) q_dot_v += mq[i] * vertex[ep][i]; // dot(q, v) vertex b
          // q_dot_nv = CMPLX_ZERO_;
          // for(int i = 0; i < 3; ++ i) q_dot_nv -= mq[i] * n_v[i];       // dot (q, n_v)
          // calculate contribution of vertex b
          // c0 = CMPLX_MINUS_ONE_ * q_dot_nt * q_dot_ne * q_dot_nv;
          // c1 = std::exp(CMPLX_MINUS_ONE_ * q_dot_v);
          // ff += c0 * c1 / f0;
          avx_m256_t f0 = avx_mul_rrp(q_sqr, avx_mul_rprp(proj_tq, proj_eq));
          avx_m256c_t q_dot_v = avx_dot_cep(mq, vertex[ei]);
          avx_edge_t n_v = avx_div_eprp(edge[ei], avx_abs_ep(edge[ei]));
          avx_m256c_t q_dot_nv = avx_dot_cep(mq, n_v);
          avx_m256c_t c0 = avx_mul_cpcp(avx_mul_rcp(AVX_CMPLX_MINUS_ONE_, q_dot_nt),
                                        avx_mul_cpcp(q_dot_ne, q_dot_nv));
          avx_m256c_t c1 = avx_exp_cp(avx_mul_rcp(AVX_CMPLX_MINUS_ONE_, q_dot_v));
          ff = avx_add_cpcp(ff, avx_div_cprp(avx_mul_cpcp(c0, c1), f0));
          q_dot_v = avx_dot_cep(mq, vertex[(ei + 1) % 3]);
          q_dot_nv = avx_dot_cep(mq, n_v);
          c0 = avx_mul_cpcp(avx_mul_rcp(AVX_CMPLX_MINUS_ONE_, q_dot_nt),
                            avx_mul_cpcp(q_dot_ne, q_dot_nv));
          c1 = avx_exp_cp(avx_mul_rcp(AVX_CMPLX_MINUS_ONE_, q_dot_v));
          ff = avx_add_cpcp(ff, avx_div_cprp(avx_mul_cpcp(c0, c1), f0));
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

    // NOTE: vectorization is on the level of triangles (triangle loop)
    //       hence, qx, qy, qz dont need to be vectorized (unlike analytical ff)
    // STEPS:
    //   convert shape_def to the vectorized data structure
    //   loop over all q-points
    //   loop over ceil(num_triangles/AVX_VEC_LEN_) vector triangles
    //   compute a vector temp_ff in each vector triangle iteration
    //   sum vector temp_ff over all vector triangle iterations
    //   reduce the vector temp_ff by summing all 4 entries into one scalar ff for a q-point

    // assuming double precision for now ... FIXME

    // allocate aligned memory for avx vector triangles
    int num_vtriangles = ceil(num_triangles / AVX_VEC_LEN_);
    int pad = AVX_VEC_LEN_ - num_triangles % AVX_VEC_LEN_;  // this is the padding in the last vtriangle
    avx_triangle_t* vtriangles = (avx_triangle_t*) _mm_malloc(num_vtriangles * sizeof(avx_triangle_t),
                                                              AVX_ALIGNMENT_);
    if(vtriangles == NULL) {
      std::cerr << "error: memory allocation failed for vtriangles. requested size: "
                << num_vtriangles * sizeof(avx_triangle_t) << " bytes" << std::endl;
      return 0;
    } // if

    // convert scalar triangles to vector triangles
    int endi = (pad > 0) ? num_vtriangles - 1 : num_vtriangles; // to handle padding separately
    for(auto i = 0; i < endi; ++ i) {
      __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_];
      int ti = i * AVX_VEC_LEN_;  // index for scalar triangles array
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v1[0];
      vtriangles[i].a.x = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v2[0];
      vtriangles[i].b.x = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v3[0];
      vtriangles[i].c.x = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v1[1];
      vtriangles[i].a.y = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v2[1];
      vtriangles[i].b.y = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v3[1];
      vtriangles[i].c.y = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v1[2];
      vtriangles[i].a.z = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v2[2];
      vtriangles[i].b.z = _mm256_load_pd(real);
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = shape_def[ti + j].v3[2];
      vtriangles[i].c.z = _mm256_load_pd(real);
    } // for
    if(pad > 0) { // handle the last vector set with padding
      __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_];
      for(auto j = 0; j < AVX_VEC_LEN_; ++ j) real[j] = 0.0;  // to make sure pad is 0 for all
      int i = num_triangles - 1;
      int ti = i * AVX_VEC_LEN_;
      int rem = AVX_VEC_LEN_ - pad;
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v1[0];
      vtriangles[i].a.x = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v2[0];
      vtriangles[i].b.x = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v3[0];
      vtriangles[i].c.x = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v1[1];
      vtriangles[i].a.y = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v2[1];
      vtriangles[i].b.y = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v3[1];
      vtriangles[i].c.y = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v1[2];
      vtriangles[i].a.z = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v2[2];
      vtriangles[i].b.z = _mm256_load_pd(real);
      for(auto j = 0; j < rem; ++ j) real[j] = shape_def[ti + j].v3[2];
      vtriangles[i].c.z = _mm256_load_pd(real);
    } // if

    unsigned long int num_qpoints = nqz;
    ff = new (std::nothrow) complex_t[num_qpoints];
    if(ff == NULL) {
      std::cerr << "error: memory allocation failed for ff of size "
                << num_qpoints * sizeof(complex_t) << " bytes" << std::endl;
      return 0;
    } // if
   
    woo::BoostChronoTimer timer;
    timer.start();

    // NOTE: each thread accesses all triangles. TODO ... improve this later

    #pragma omp parallel for schedule(runtime)
    for(int i_z = 0; i_z < nqz; ++ i_z) {
      int i_y = i_z % nqy; 
      avx_m256c_t ff_temp = AVX_ZERO_, temp;
      for(int i_t = 0; i_t < num_vtriangles; ++ i_t) {    // NOTE: padding is included
        temp = form_factor_kernel_exact(qx[i_y], qy[i_y], qz[i_z], rot, vtriangles[i_t]);
        ff_temp = avx_add_cpcp(ff_temp, temp);
      } // for
      avx_m256c_t sum = avx_hadd_cpcp(ff_temp, ff_temp);  // partial sums
      __m128d sum_high = _mm256_extractf128_pd(sum, 1);   // extract upper 128 bits
      ff[i_z] = _mm_cvtsd_f64(_mm_add_pd(_mm256_extractf128_pd(sum, 1), _mm256_extract128_pd(sum, 0)));
    } // for

    timer.stop();
    compute_time = timer.elapsed_msec();

    return num_triangles;
  } // NumericFormFactorC::compute_exact_triangle()


#endif // INTEL_AVX

} // namespace hig
