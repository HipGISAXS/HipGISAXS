/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: avx_numerics.hpp
 *  Created: Apr 19, 2013
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

#ifndef __AVX_NUMERICS_HPP__
#define __AVX_NUMERICS_HPP__

#ifdef INTEL_AVX

#include <cstring>

#include <common/typedefs.hpp>
#include <common/constants.hpp>
#include <numerics/numeric_utils.hpp>

namespace hig {

	// ////////////////////////////////////////////////////
	// intrinsic wrapper naming (only for floating-point)
	// _mm256_xxx_abc  => a = r|c, b = p|s, c = s|d
	// _mm256_xxx_abcd => a = r|c, b = r|c, c = p|s, d = s|d
	// r = real,             c = complex,
	// p = packed (vector),  s = scalar,
	// s = single-precision, d = double-precision
	// ////////////////////////////////////////////////////

#ifdef DOUBLEP    // TODO: improve and move it within each functions instead

  const avx_m256_t AVX_ZERO_ = _mm256_setzero_pd();
  const avx_m256_t AVX_ONE_ = _mm256_set1_pd(1.0);
  const avx_m256_t AVX_ALLONE_ = _mm256_set1_pd(0xFFFFFFFFFFFFFFFF);
  const avx_m256_t AVX_EPSILON_ = _mm256_set1_pd(REAL_EPSILON_);

  /**
   * complex load/store/initialize
   */

  static inline avx_m256c_t avx_load_cp(const complex_t* p) {
    // TODO: use _mm256_i64scatter_pd() instead to directly store into p ...
    avx_m256c_t v;
    __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_], imag[AVX_VEC_LEN_];
    for(int i = 0; i < AVX_VEC_LEN_; ++ i) {
      real[i] = p[i].real(); imag[i] = p[i].imag();
    } // for
    v.real = _mm256_load_pd(real);
    v.imag = _mm256_load_pd(imag);
    return v;
  } // avx_load_cp()

  static inline void avx_store_cp(complex_t* p, avx_m256c_t v) {
    // TODO: use _mm256_i64scatter_pd() instead to directly store into p ...
    __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_], imag[AVX_VEC_LEN_];
    _mm256_store_pd(real, v.real);
    _mm256_store_pd(imag, v.imag);
    for(int i = 0; i < AVX_VEC_LEN_; ++ i) p[i] = complex_t(real[i], imag[i]);
  } // avx_store_cp()

  static inline avx_m256c_t avx_setzero_cp(void) {
    avx_m256c_t v;
    v.real = _mm256_setzero_pd();
    v.imag = _mm256_setzero_pd();
    return v;
  } // avx_setzero_cp()

  /**
   * absolute value
   */

  static inline avx_m256_t avx_fabs_rp(avx_m256_t v) {
    avx_m256_t m = _mm256_set1_pd(-0.0);
    return _mm256_andnot_pd(m, v);
  } // avx_fabs_rp()

  /**
   * addition
   */

  static inline avx_m256_t avx_add_rrp(avx_m256_t a, avx_m256_t b) {
    return _mm256_add_pd(a, b);
  } // avx_add_rrp()

  static inline avx_m256c_t avx_add_crp(avx_m256c_t a, avx_m256_t b) {
    avx_m256c_t v;
    v.real = _mm256_add_pd(a.real, b);
    v.imag = a.imag;
    return v;
  } // avx_add_crp()

  static inline avx_m256c_t avx_add_ccp(avx_m256c_t a, avx_m256c_t b) {
    avx_m256c_t v;
    v.real = _mm256_add_pd(a.real, b.real);
    v.imag = _mm256_add_pd(a.imag, b.imag);
    return v;
  } // avx_add_ccp()

  /**
   * multiplication
   */

  static inline avx_m256c_t avx_mul_crp(avx_m256c_t a, real_t b) {
    avx_m256c_t v;
    //avx_m256_t temp = _mm256_broadcast_sd(&b);
    avx_m256_t temp = _mm256_set1_pd(b);
    v.real = _mm256_mul_pd(a.real, temp);
    v.imag = _mm256_mul_pd(a.imag, temp);
    return v;
  } // avx_mul_crp()

  static inline avx_m256c_t avx_mul_crp(complex_t a, avx_m256_t b) {
    avx_m256c_t v;
    avx_m256c_t temp;
    temp.real = _mm256_set1_pd(a.real());
    temp.imag = _mm256_set1_pd(a.imag());
    v.real = _mm256_mul_pd(temp.real, b);
    v.imag = _mm256_mul_pd(temp.imag, b);
    return v;
  } // avx_mul_crp()

  static inline avx_m256c_t avx_mul_ccp(complex_t a, avx_m256c_t b) {
    avx_m256c_t v;
    avx_m256c_t tempa;
    tempa.real = _mm256_set1_pd(a.real());
    tempa.imag = _mm256_set1_pd(a.imag());
    avx_m256_t temp1 = _mm256_mul_pd(tempa.real, b.real);
    avx_m256_t temp2 = _mm256_mul_pd(tempa.imag, b.imag);
    avx_m256_t temp3 = _mm256_mul_pd(tempa.real, b.imag);
    avx_m256_t temp4 = _mm256_mul_pd(tempa.imag, b.real);
    v.real = _mm256_sub_pd(temp1, temp2);
    v.imag = _mm256_add_pd(temp3, temp4);
    return v;
  } // avx_mul_ccp()

  static inline avx_m256c_t avx_mul_ccp(avx_m256c_t a, avx_m256c_t b) {
    avx_m256c_t v;
    avx_m256_t temp1 = _mm256_mul_pd(a.real, b.real);
    avx_m256_t temp2 = _mm256_mul_pd(a.imag, b.imag);
    avx_m256_t temp3 = _mm256_mul_pd(a.real, b.imag);
    avx_m256_t temp4 = _mm256_mul_pd(a.imag, b.real);
    v.real = _mm256_sub_pd(temp1, temp2);
    v.imag = _mm256_add_pd(temp3, temp4);
    return v;
  } // avx_mul_ccp()

  static inline avx_m256c_t avx_fma_rccp(real_t a, avx_m256c_t b, avx_m256c_t c) {
    // a * b + c
    avx_m256c_t v;
    //avx_m256_t temp = _mm256_broadcast_sd(&a);
    avx_m256_t temp = _mm256_set1_pd(a);
    //v.real = _mm256_fmadd_pd(temp, b.real, c.real);
    v.real = _mm256_add_pd(_mm256_mul_pd(temp, b.real), c.real);
    //v.imag = _mm256_fmadd_pd(temp, b.imag, c.imag);
    v.imag = _mm256_add_pd(_mm256_mul_pd(temp, b.imag), c.imag);
    return v;
  } // avx_fma_rccp()

  static inline avx_m256_t avx_mul_rrp(real_t a, avx_m256_t b) {
    return _mm256_mul_pd(b, _mm256_set1_pd(a));
  } // avx_mul_rrp()

  static inline avx_m256_t avx_mul_rrp(avx_m256_t a, real_t b) {
    return _mm256_mul_pd(a, _mm256_set1_pd(b));
  } // avx_mul_rrp()

  /**
   * division
   */

  static inline avx_m256c_t avx_div_ccp(avx_m256c_t a, avx_m256c_t b) {
    avx_m256c_t v;
    avx_m256_t temp1 = _mm256_mul_pd(a.real, b.real);
    avx_m256_t temp2 = _mm256_mul_pd(a.imag, b.imag);
    avx_m256_t temp3 = _mm256_mul_pd(a.real, b.imag);
    avx_m256_t temp4 = _mm256_mul_pd(a.imag, b.real);
    avx_m256_t temp5 = _mm256_hypot_pd(b.real, b.imag);
    v.real = _mm256_div_pd(_mm256_add_pd(temp1, temp2), temp5);
    v.imag = _mm256_div_pd(_mm256_sub_pd(temp4, temp3), temp5);
    return v;
  } // avx_div_ccp()

  static inline avx_m256_t avx_div_rrp(avx_m256_t a, real_t b) {
    return _mm256_div_pd(a, _mm256_set1_pd(b));
  } // avx_div_rrp()

  /**
   * square-root
   */

  static inline avx_m256c_t avx_sqrt_cp(avx_m256c_t v) {
    avx_m256c_t res;
    avx_m256_t x = v.real;
    avx_m256_t y = v.imag;
    // TODO: improve and implement the special cases ...
    //if (x == 0) {
    //  double t = sqrt(fabs(y) / 2);
    //  return make_cuC(t, y < 0 ? -t : t);
    //} else if (y == 0) {
    //  return x < 0 ? make_cuC(0, sqrt(fabs(x))) : make_cuC(sqrt(x), 0);
    //} else {
    avx_m256_t hypv = _mm256_hypot_pd(x, y);
    avx_m256_t absx = avx_fabs_rp(x);
    avx_m256_t absy = avx_fabs_rp(y);
    avx_m256_t t = _mm256_sqrt_pd(avx_mul_rrp(_mm256_add_pd(hypv, absx), 2.0));
    avx_m256_t u = avx_div_rrp(t, 2.0);
    avx_m256_t ybt = _mm256_div_pd(y, t);
    avx_m256_t absybt = _mm256_div_pd(absy, t);
    avx_m256_t maskx = _mm256_cmp_pd(x, _mm256_set1_pd(0.0), _CMP_GT_OS);
    avx_m256_t masky = _mm256_cmp_pd(y, _mm256_set1_pd(0.0), _CMP_LT_OS);

    avx_m256_t r1 = _mm256_and_pd(maskx, u);
    avx_m256_t r2 = _mm256_andnot_pd(maskx, absybt);
    res.real = _mm256_or_pd(r1, r2);

    avx_m256_t i1 = _mm256_and_pd(maskx, ybt);
    avx_m256_t i21 = _mm256_and_pd(masky, _mm256_mul_pd(u, _mm256_set1_pd(-1.0)));
    avx_m256_t i22 = _mm256_andnot_pd(masky, u);
    avx_m256_t i2 = _mm256_or_pd(i21, i22);
    i2 = _mm256_andnot_pd(maskx, i2);
    res.imag = _mm256_or_pd(i1, i2);

    return res;
    //} // if-else
  } // avx_sqrt_cp()

  /**
   * exponential
   */

  static inline avx_m256c_t avx_exp_cp(avx_m256c_t v) {
    avx_m256c_t res;
    avx_m256_t expa = _mm256_exp_pd(v.real);
    avx_m256_t cosb = _mm256_cos_pd(v.imag);
    avx_m256_t sinb = _mm256_sin_pd(v.imag);
    res.real = _mm256_mul_pd(expa, cosb);
    res.imag = _mm256_mul_pd(expa, sinb);
    return res;
  } // avx_exp_cp()

  /**
   * sinc
   */

  static inline avx_m256c_t avx_sinc_cp(avx_m256c_t v) {
    // TODO: properly vectorize ...
    // currently it converts to scalar, computes, and converts back to vector
    avx_m256c_t res;
    __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_], imag[AVX_VEC_LEN_];
    _mm256_store_pd(real, v.real);
    _mm256_store_pd(imag, v.imag);
    for(int i = 0; i < AVX_VEC_LEN_; ++ i) {
      complex_t temp = sinc(complex_t(real[i], imag[i]));
      real[i] = temp.real();
      imag[i] = temp.imag();
    } // for
    res.real = _mm256_load_pd(real);
    res.imag = _mm256_load_pd(imag);
    return res;
  } // avx_sinc_cp()

  /**
   * bessel function
   */

  static inline avx_m256c_t avx_cbessj_cp(avx_m256c_t v, int o) {
    // convert to scalar, compute, and convert back to vector
    avx_m256c_t res;
    __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_], imag[AVX_VEC_LEN_];
    _mm256_store_pd(real, v.real);
    _mm256_store_pd(imag, v.imag);
    for(int i = 0; i < AVX_VEC_LEN_; ++ i) {
      complex_t temp = cbesselj(complex_t(real[i], imag[i]), o);
      real[i] = temp.real();
      imag[i] = temp.imag();
    } // for
    res.real = _mm256_load_pd(real);
    res.imag = _mm256_load_pd(imag);
    return res;
  } // avx_cbessj_cp()

  static inline avx_m256c_t avx_cj1_cp(avx_m256c_t v) {
    const int MAXK = 1e5;
    #ifdef DOUBLEP
      const int digits = 15;    // double precision
    #else
      const int digits = 7;     // single precision
    #endif
    const real_t threshold = 1.73 * digits;

    const avx_m256_t r = v.real;
    avx_m256_t asympt_mask = _mm256_cmp_pd(r, _mm256_set1_pd(threshold), _CMP_GT_OS);
    // asympt_mask = mask for asymptotic method
    int imask = _mm256_movemask_pd(asympt_mask);
    // imask == 0 => all are LE threshold
    // imask == 0xF => all are GT threshold
    // imask, otherwise has mixed and use asymp_mask as a mask
    avx_m256c_t direct_res, asympt_res;

    if(imask < 0xF) {          // either imask == 0 (all <= threshold) or mixed
      // the direct method to compute j1
      avx_m256_t g = AVX_ONE_;  // (1.0 / gamma(nu + k + 1)), nu = 1, k = 0 => (1.0 / gamma(2)) = 1.0
      avx_m256_t r2 = _mm256_andnot_pd(asympt_mask,
                                       _mm256_mul_pd(_mm256_set1_pd(0.25), _mm256_mul_pd(r, r)));
      avx_m256_t r0 = g, temp = g;
      int k = 1, m = 1;
      for(; k < MAXK; ++ k) {
        m *= -1;
        temp = _mm256_mul_pd(_mm256_div_pd(r2, _mm256_set1_pd(k * (k + 1))), temp);
                                                      // (-1)^k * (z^2 / 4)^k / (k * (k + 1))
        r0 = _mm256_add_pd(r0, _mm256_mul_pd(_mm256_set1_pd(m), temp));
        //if(fabs(temp) < REAL_EPSILON_) break;
        avx_m256_t cmp = _mm256_cmp_pd(temp, AVX_EPSILON_, _CMP_GT_OS);
        int cond = _mm256_movemask_pd(cmp);
        if(cond == 0) break;
      } // for
      if(k == MAXK) {
        std::cerr << "error: bessel loop overflow occured in direct method" << std::endl;
        return v;
      } // if
      direct_res.real = _mm256_mul_pd(_mm256_set1_pd(0.5), _mm256_mul_pd(r0, r));  // (z/2)^nu
      direct_res.imag = AVX_ZERO_;
    } // if

    if(imask > 0) {          // either imask == 0xF (all > threshold) or mixed
      // asymptotic method to compute j1
      avx_m256_t ir = _mm256_and_pd(asympt_mask, _mm256_div_pd(AVX_ONE_, r));
      avx_m256_t irk = ir;
      avx_m256_t ak = AVX_ONE_;                 // a0 = 1
      avx_m256_t xk = ak;                       // x0 = a0 / 1
      avx_m256_t xk_sum = xk;
      ak = _mm256_set1_pd(0.375);               // a1 = 3 / 8
      avx_m256_t yk = _mm256_mul_pd(ak, irk);   // y0 = a1 / z
      avx_m256_t yk_sum = yk;
      int k = 1, m = 1;
      avx_m256_t t0, t1, t2, t3, t4, t5, t6, t7;
      for(; k < MAXK; ++ k) {
        m *= -1;
        t0 = _mm256_set1_pd(16. * k);
        irk = _mm256_mul_pd(irk, ir);
        t1 = _mm256_set1_pd(4. * k - 1.);
        t2 = _mm256_div_pd(_mm256_sub_pd(_mm256_set1_pd(4.0), _mm256_mul_pd(t1, t1)), t0);
        ak = _mm256_mul_pd(ak, t2);
        xk = _mm256_mul_pd(_mm256_set1_pd(m), _mm256_mul_pd(ak, irk));
        xk_sum = _mm256_add_pd(xk_sum, xk);
        irk = _mm256_mul_pd(irk, ir);
        t1 = _mm256_set1_pd(4. * k + 1.);
        t2 = _mm256_div_pd(_mm256_sub_pd(_mm256_set1_pd(4.0), _mm256_mul_pd(t1, t1)),
                           _mm256_add_pd(t0, _mm256_set1_pd(8.0)));
        ak = _mm256_mul_pd(ak, t2);
        yk = _mm256_mul_pd(_mm256_set1_pd(m), _mm256_mul_pd(ak, irk));
        yk_sum = _mm256_add_pd(yk_sum, yk);
        // if(fabs(xk) < REAL_EPSILON_ && fabs(yk) < REAL_EPSILON_) break;
        avx_m256_t xcmp = _mm256_cmp_pd(xk, AVX_EPSILON_, _CMP_GT_OS);
        avx_m256_t ycmp = _mm256_cmp_pd(yk, AVX_EPSILON_, _CMP_GT_OS);
        int xcond = _mm256_movemask_pd(xcmp);
        int ycond = _mm256_movemask_pd(ycmp);
        if(xcond == 0 && ycond == 0) break;
      } // for
      if(k == MAXK) {
        std::cerr << "error: Bessel loop overflow occured in asymptotic method" << std::endl;
        return v;
      } // if
      avx_m256_t w = _mm256_sub_pd(r, _mm256_set1_pd(0.75 * PI_));
      t0 = _mm256_mul_pd(r, _mm256_set1_pd(PI_));
      t1 = _mm256_div_pd(_mm256_set1_pd(2.0), t0);
      t0 = _mm256_sqrt_pd(t1);
      t3 = _mm256_sincos_pd(&t2, w);    // t3 = sin(w), t2 = cos(w)
      t4 = _mm256_mul_pd(t2, xk_sum);
      t5 = _mm256_mul_pd(t3, yk_sum);
      t6 = _mm256_sub_pd(t4, t5);
      asympt_res.real = _mm256_mul_pd(t0, t6);
      asympt_res.imag = AVX_ZERO_;
    } // if

    if(imask == 0x0) return direct_res;
    if(imask == 0xF) return asympt_res;

    // mixed case
    avx_m256c_t mixed_res;
    avx_m256_t ares = _mm256_and_pd(asympt_mask, asympt_res.real);
    avx_m256_t dres = _mm256_andnot_pd(asympt_mask, direct_res.real);
    mixed_res.real = _mm256_or_pd(ares, dres);
    mixed_res.imag = AVX_ZERO_;

    return mixed_res;
  } // avx_cj1_cp()


  static inline avx_m256c_t avx_cbesselj_cp(avx_m256c_t v, int o) {
    if(o == 1) return avx_cj1_cp(v);
    std::cerr << "error: Bessel functions of order != 1 are not supported" << std::endl;
    return v;
  } // avx_cbesselj_cp()


  /**
   * other avx based operations with structs/arrays of vectors
   * t: triangle, v: vertex, e: edge, p: point/position
   */

  static inline avx_edge_t avx_sub_vpvp(avx_vertex_t a, avx_vertex_t b) {
    avx_edge_t res;
    res.x = _mm256_sub_pd(a.x, b.x);
    res.y = _mm256_sub_pd(a.y, b.y);
    res.z = _mm256_sub_pd(a.z, b.z);
    return res;
  } // avx_sub_epep()

  static inline avx_edge_t avx_cross_epep(avx_edge_t a, avx_edge_t b) {
    avx_edge_t c;
    c.x = _mm256_sub_pd(_mm256_mul_pd(a.y, b.z), _mm256_mul_pd(a.z, b.y));
    c.y = _mm256_sub_pd(_mm256_mul_pd(a.z, b.x), _mm256_mul_pd(a.x, b.z));
    c.z = _mm256_sub_pd(_mm256_mul_pd(a.x, b.y), _mm256_mul_pd(a.y, b.x));
    return c;
  } // avx_cross_epep()

  static inline avx_m256_t avx_abs_ep(avx_edge_t e) {
    return _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(e.x, e.x),
                                        _mm256_add_pd(_mm256_mul_pd(e.y, e.y),
                                                      _mm256_mul_pd(e.z, e.z))));
  } // avx_abs_ep()

  static inline avx_edge_t avx_div_eprp(avx_edge_t e, avx_m256_t r) {
    avx_edge_t res;
    avx_m256_t rr = _mm256_div_pd(AVX_ONE_, r);
    res.x = _mm256_mul_pd(e.x, rr);
    res.y = _mm256_mul_pd(e.y, rr);
    res.z = _mm256_mul_pd(e.z, rr);
    return res;
  } // avx_div_eprp()

  static inline avx_m256c_t avx_dot_cep(const std::vector<complex_t>& a, avx_edge_t e) {
    return avx_add_ccp(avx_add_ccp(avx_mul_crp(a[0], e.x), avx_mul_crp(a[1], e.y)),
                                   avx_mul_crp(a[2], e.z));
  } // avx_dot_cep()

  static inline avx_m256_t avx_norm_cp(avx_m256c_t a) {
    return _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(a.real, a.real), _mm256_mul_pd(a.imag, a.imag)));
  } // avx_norm_cp()

  static inline avx_m256_t avx_sub_rrp(real_t a, avx_m256_t b) {
    return _mm256_sub_pd(_mm256_set1_pd(a), b);
  } // avx_sub_rrp()

  // AVX_CMPLX_MINUS_ONE_

  static inline avx_m256c_t avx_mul_rpcp(avx_m256_t a, avx_m256c_t b) {
    avx_m256c_t res;
    res.real = _mm256_mul_pd(a, b.real);
    res.imag = _mm256_mul_pd(a, b.imag);
    return res;
  } // avx_mul_rcp()

  static inline avx_m256c_t avx_div_cprp(avx_m256c_t a, avx_m256_t b) {
    avx_m256c_t res;
    res.real = _mm256_div_pd(a.real, b);
    res.imag = _mm256_div_pd(a.imag, b);
    return res;
  } // avx_div_cprp()

  static inline avx_m256c_t avx_hadd_cpcp(avx_m256c_t a, avx_m256c_t b) {
    avx_m256c_t res;
    res.real = _mm256_hadd_pd(a.real, b.real);
    res.imag = _mm256_hadd_pd(a.imag, b.imag);
    return res;
  } // avx_add_cpcp()

  static inline complex_t avx_hreduce_cp(avx_m256c_t v) {
    __m128d sum_high_real, sum_high_imag, sum_low_real, sum_low_imag;
    avx_m256c_t sum = avx_hadd_cpcp(v, v);                // partial sums
    sum_high_real = _mm256_extractf128_pd(sum.real, 1);   // extract upper 128 bits
    sum_high_imag = _mm256_extractf128_pd(sum.imag, 1);   // extract upper 128 bits
    sum_low_real = _mm256_extractf128_pd(sum.real, 0);    // extract lower 128 bits
    sum_low_imag = _mm256_extractf128_pd(sum.imag, 0);    // extract lower 128 bits
    return complex_t(_mm_cvtsd_f64(_mm_add_pd(sum_high_real, sum_low_real)),
                     _mm_cvtsd_f64(_mm_add_pd(sum_high_imag, sum_low_imag)));
  } // avx_hreduce_cp()

  // load triangles to avx triangles
  static inline avx_triangle_t avx_load_triangles(const triangle_t* shape_def, int pad = 0) {
    __AVX_ALIGNED__ real_t real[AVX_VEC_LEN_];
    memset(real, 0, AVX_VEC_LEN_ * sizeof(real_t)); // to make sure pad is 0 for all
    avx_triangle_t vt;
    int nt = AVX_VEC_LEN_ - pad;                    // load these many triangles from shape_def to vt
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v1[0]; vt.a.x = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v2[0]; vt.b.x = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v3[0]; vt.c.x = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v1[1]; vt.a.y = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v2[1]; vt.b.y = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v3[1]; vt.c.y = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v1[2]; vt.a.z = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v2[2]; vt.b.z = _mm256_load_pd(real);
    for(auto j = 0; j < nt; ++ j) real[j] = shape_def[j].v3[2]; vt.c.z = _mm256_load_pd(real);
    return vt;
  } // avx_load_triangles()

#else     // single precision

	/**
	 * load store
	 */

	static inline avx_m256_t avx_load_rps(real_t* p) {
		return _mm256_load_ps(p);
	} // avx_load_rps()

	static inline void avx_addstore_css(complex_t* p, avx_m256c_t v) {
		real_t real = _mm256_cvtss_f32(v.xvec);
		real_t imag = _mm256_cvtss_f32(v.yvec);
		(*p) += complex_t(real, imag);
	} // avx_store_css()


	/**
	 * set
	 */
	
	static inline avx_m256c_t avx_setzero_cps() {
		avx_m256c_t vec;
		vec.xvec = _mm256_setzero_ps();
		vec.yvec = _mm256_setzero_ps();
		return vec;
	} // avx_setzero_cps()

	static inline avx_m256_t avx_set1_rps(real_t a) {
		return _mm256_set1_ps(a);
	} // avx_load1_rps()

	static inline avx_m256c_t avx_set1_cps(complex_t a) {
		avx_m256c_t vec;
		vec.xvec = _mm256_set1_ps(a.real());
		vec.yvec = _mm256_set1_ps(a.imag());
		return vec;
	} // avx_set1_cps()


	/**
	 * addition
	 */

	static inline avx_m256_t avx_add_rrps(avx_m256_t a, avx_m256_t b) {
		return _mm256_add_ps(a, b);
	} // avx_add_rrps()

	static inline avx_m256c_t avx_add_rcps(avx_m256_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_add_ps(a, b.xvec);
		vec.yvec = b.yvec;
		return vec;
	} // avx_add_rcps()

	static inline avx_m256c_t avx_add_ccps(avx_m256c_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_add_ps(a.xvec, b.xvec);
		vec.yvec = _mm256_add_ps(a.yvec, b.yvec);
		return vec;
	} // avx_add_ccps()

	static inline avx_m256c_t avx_hadd_ccps(avx_m256c_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_hadd_ps(a.xvec, b.xvec);
		vec.yvec = _mm256_hadd_ps(a.yvec, b.yvec);
		return vec;
	} // avx_hadd_ccps()


	/**
	 * multiplication
	 */

	static inline avx_m256_t avx_mul_rrps(avx_m256_t a, avx_m256_t b) {
		return _mm256_mul_ps(a, b);
	} // avx_mul_rrps()

	static inline avx_m256c_t avx_mul_crps(avx_m256c_t a, avx_m256_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_mul_ps(a.xvec, b);
		vec.yvec = _mm256_mul_ps(a.yvec, b);
		return vec;
	} // avx_mul_cpps()

	static inline avx_m256c_t avx_mul_ccps(avx_m256c_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_sub_ps(_mm256_mul_ps(a.xvec, b.xvec), _mm256_mul_ps(a.yvec, b.yvec));
		vec.yvec = _mm256_add_ps(_mm256_mul_ps(a.xvec, b.yvec), _mm256_mul_ps(a.yvec, b.xvec));
		return vec;
	} // avx_mul_ccps()


	/**
	 * reciprocal
	 */

	static inline avx_m256c_t avx_rcp_cps(avx_m256c_t a) {
		avx_m256_t temp1 = _mm256_add_ps(_mm256_mul_ps(a.xvec, a.xvec), _mm256_mul_ps(a.yvec, a.yvec));
		avx_m256_t temp2 = _mm256_rcp_ps(temp1);
		avx_m256c_t vec;
		avx_m256_t neg_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
		vec.xvec = _mm256_mul_ps(a.xvec, temp2);
		vec.yvec = _mm256_xor_ps(_mm256_mul_ps(a.yvec, temp2), neg_mask);
		return vec;
	} // avx_rcp_cps()


	/**
	 * exponential
	 */

	static inline avx_m256_t avx_exp_rps(avx_m256_t a) {
		return avx_mathfun::newexp_ps(a);
	} // avx_exp_rps()

	static inline void avx_exp_rps_dual(avx_m256_t a1, avx_m256_t a2, avx_m256_t* exp1, avx_m256_t* exp2) {
		avx_mathfun::newexp_ps_dual(a1, a2, exp1, exp2);
	} // avx_exp_rps()


	/**
	 * trigonometry
	 */

	static inline avx_m256_t avx_sin_rps(avx_m256_t a) {
		//return amd_vrs4_sinf(a);
		return avx_mathfun::newsin_ps(a);
	} // avx_exp_rps()

	static inline avx_m256_t avx_cos_rps(avx_m256_t a) {
		//return amd_vrs4_cosf(a);
		return avx_mathfun::newcos_ps(a);
	} // avx_exp_rps()

	static inline void avx_sincos_rps(avx_m256_t a, avx_m256_t* s, avx_m256_t* c) {
		avx_mathfun::newsincos_ps(a, s, c);
		//*s = avx_mathfun::newsin_ps(a);
		//*c = avx_mathfun::newcos_ps(a);
	} // avx_exp_rps()

	static inline void avx_sincos_rps_dual(avx_m256_t a1, avx_m256_t a2, avx_m256_t* s1, avx_m256_t* s2,
											avx_m256_t* c1, avx_m256_t* c2) {
		avx_mathfun::newsincos_ps_dual(a1, a2, s1, s2, c1, c2);
	} // avx_exp_rps()

#endif // DOUBLEP

} // namespace hig

#endif // INTEL_AVX

#endif // __AVX_NUMERICS_HPP__
