/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: numeric_utils.cpp
 *  Created: Oct 08, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *          Dinesh Kumar <dkumar@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>
#include <cmath>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include <common/constants.hpp>
#include <numerics/numeric_utils.hpp>
#ifdef INTEL_AVX
#include <numerics/cpu/avx_numerics.hpp>
#endif // INTEL_AVX

namespace hig {

  // vector3 dot product
  real_t dot(vector3_t& a, vector3_t& b) {
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  } // dot()

  // vector3 cross product
  vector3_t cross(vector3_t & u, vector3_t & v) {
    vector3_t w;
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    return w;
  } // cross()

  // sinc function
  complex_t sinc(complex_t x) {
    if(std::abs(x) > 1.0E-14) return std::sin(x) / x;
    else return (1.- x * x / 6. + x * x * x * x / 120.);  // taylor series approx
  } // sinc()

  #ifdef FF_CPU_OPT
  void sinc_vec(const int VEC_LEN, complex_t* x, complex_t* res) {
    for(int i = 0; i < VEC_LEN; ++ i) res[i] = sinc(x[i]);
  } // sinc_vec()
  #endif // FF_CPU_OPT


  // besselj functions (of the first kind)

  // gamma function
  // double gamma(double x) {
  //   double coef[6];
  //   coef[0] = 76.18009173;
  //   coef[1] = -86.50532033;
  //   coef[2] = 24.01409822;
  //   coef[3] = -1.231739516;
  //   coef[4] = 0.120858003e-2;
  //   coef[5] = -0.536382e-5;
  //   double stp = 2.50662827465;
  //   double temp = x + 5.5;
  //   temp = (x + 0.5) * log(temp) - temp;
  //   double ser = 1.0;
  //   for(int j = 0; j < 6; ++ j) {
  //     x += 1.0;
  //     ser += coef[j] / x;
  //   } // for
  //   return exp(temp + log(stp * ser));
  // } // gamma()

  // /*--------------------------------------------------
  //                         inf.     (-z^2/4)^k
  //     Jnu(z) = (z/2)^nu x Sum  ------------------
  //                         k=0  k! x Gamma(nu+k+1)
  //   (nu must be >= 0). Here k=15.
  // ---------------------------------------------------*/
  // complex_t cbessj(complex_t zz, int order) {
  //   std::complex<long double> z = zz;
  //   std::complex<long double> temp1 = pow(z / (long double) 2.0, order);
  //   std::complex<long double> z2 = - z * z / (long double) 4.0;
  //   std::complex<long double> sum(0.0, 0.0);
  //   long double factorial_k = 1.0;
  //   std::complex<long double> pow_z2_k = 1.0;
  //   //for(int k = 0; k <= MAXK; ++ k, pow_z2_k *= z2) {
  //   for(int k = 0; k <= 200; ++ k, pow_z2_k *= z2) {
  //     if(k == 0) factorial_k = 1.0;  // base case
  //     else factorial_k *= k;        // compute k!
  //     std::complex<long double> temp2 =
  //           pow_z2_k / (factorial_k * (long double) gamma((double) order + k + 1.0));
  //     sum += temp2;
  //   } // for
  //   temp1 *= sum;
  //   return complex_t((real_t) temp1.real(), (real_t) temp1.imag());
  // } // cbessj()


  // using GNU math. assuming imaginary component is 0
  complex_t cbessj(complex_t zz, int order) {
    return complex_t(j1(zz.real()), 0.0);
  } // cbessj()

  #ifdef FF_CPU_OPT
  // using GNU math. assuming imaginary component is 0
  void cbessj_vec(const int VEC_LEN, complex_t* zz, int order, complex_t* res) {
    for(int i = 0; i < VEC_LEN; ++ i) res[i] = complex_t(j1(zz[i].real()), 0.0);
  } // cbessj()
  #endif // FF_CPU_OPT


  // using new local implementation below
  complex_t cbesselj(complex_t z, int order) {
    if(order == 1) return cj1(z);
    std::cerr << "error: Bessel functions for order != 1 are not supported" << std::endl;
    return complex_t(0., 0.);
  } // cbesselj()

  #ifdef FF_CPU_OPT
  #ifdef FF_CPU_OPT_AVX
  // using new local avx implementation (this is just a wrapper)
  avx_m256c_t cbesselj_vec(avx_m256c_t z, int order) {
    return avx_cbesselj_cp(z, order);
  } // cbesselj()
  #else
  // using new local implementation below
  complex_t cbesselj_vec(complex_t z, int order) {
    if(order == 1) return cj1(z);
    std::cerr << "error: Bessel functions for order != 1 are not supported" << std::endl;
    return complex_t(0., 0.);
  } // cbesselj()
  #endif // FF_CPU_OPT_AVX
  #endif // FF_CPU_OPT

  // Bessel function with MPFUN-fort's mpbesselj as reference.
  // MPFUN-fort is developed by David Bailey et al.

  // generic implementation: bessel fuction with order 1
  complex_t cj1(complex_t z) {
    int MAXK = 1e5; // max number of iterations in the following
    // z = r + i, r = real, i = imaginary
    real_t r = z.real();
    // single precision ~ 7 digits precision (24 * log(2) / log(10))
    // double precision ~ 15 digits precision (53 * log(2) / log(10))
    int digits = 7;
    #ifdef DOUBLEP
      digits = 15;
    #endif
    real_t threshold = 1.73 * digits;
    if(r > threshold) {
      // the asymptotic method to compute j1      
      // J1(z) = (2/(PIz))^(1/2)(cos(w)\SUM_0^inf(A(k)) - sin(w)\SUM_0^inf(B(k)))
      double ir = 1. / r;
      double irk = ir;
      double ak = 1.;                     // a0 = 1
      double xk = ak;                     // x0 = a0 / 1
      double xk_sum = xk;
      ak = 0.375;                         // a1 = 3 / 8
      double yk = ak * irk;               // y0 = a1 / z
      double yk_sum = yk;
      int k = 1, m = 1;
      for(; k < MAXK; ++ k) {
        m *= -1;
        double t0 = 16. * k;
        irk *= ir;
        double t1 = 4. * k - 1.;
        ak *= (4. - t1 * t1) / t0;        // a{2k}
        xk = m * ak * irk;
        xk_sum += xk;
        irk *= ir;
        double t2 = 4. * k + 1.;
        ak *= (4. - t2 * t2) / (t0 + 8.); // a{2k+1}
        yk = m * ak * irk;
        yk_sum += yk;
        if(fabs(xk) < REAL_EPSILON_ && fabs(yk) < REAL_EPSILON_) break;
      } // for
      if(k == MAXK) {
        std::cerr << "error: bessel loop overflow with fabs(" << xk << ") > " << REAL_EPSILON_ << " and/or fabs(" << yk << ") > " << REAL_EPSILON_ << std::endl;
        assert(k != MAXK);
      } // if
      double w = r - 0.75 * PI_;
      return complex_t(sqrt(2. / (r * PI_)) * (cos(w) * xk_sum - sin(w) * yk_sum), 0.); 
    } else {
      // the direct method to compute j1
      double g = 1. / gamma((unsigned int)(1 + 1));   // 1 / gamma(nu + k + 1),
                                                      // nu = 1, k = 0
      double r2 = 0.25 * r * r;        // (z^2 / 4)
      double r0 = g, temp = g;
      int k = 1;
      for(; k < MAXK; ++ k) {
        temp = -1. * (r2 / (k * (k + 1))) * temp;    // (-1)^k * (z^2 / 4)^k / (k * (k + 1))
        r0 += temp;
        if(fabs(temp) < REAL_EPSILON_) break;
      } // for
      if(k == MAXK) {
        std::cerr << "error: bessel loop overflow with fabs(" << temp << ") > " << REAL_EPSILON_ << std::endl;
        assert(k != MAXK);
      } // if
      return complex_t(r0 * 0.5 * r, 0.);    // (z/2)^nu
    } // if-else
  } // cj1()
  
  // special case of gamma function for x = 2
  inline real_t gamma(unsigned int x) {
    assert(x == 2);
    return 1.;        // gamma(2) = 1
  } // gamma()

} // namespace hig
