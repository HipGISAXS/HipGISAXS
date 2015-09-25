/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: numeric_utils.cpp
 *  Created: Oct 08, 2012
 *  Modified: Wed 08 Oct 2014 12:17:46 PM PDT
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

#include <numerics/numeric_utils.hpp>

#include <cmath>

namespace hig {

  double gamma(double x) {
    double coef[6];
    coef[0] = 76.18009173;
    coef[1] = -86.50532033;
    coef[2] = 24.01409822;
    coef[3] = -1.231739516;
    coef[4] = 0.120858003e-2;
    coef[5] = -0.536382e-5;
    double stp = 2.50662827465;
    double temp = x + 5.5;
    temp = (x + 0.5) * log(temp) - temp;
    double ser = 1.0;
    for(int j = 0; j < 6; ++ j) {
      x += 1.0;
      ser += coef[j] / x;
    } // for
    return exp(temp + log(stp * ser));
  } // gamma()

  /*--------------------------------------------------
                          inf.     (-z^2/4)^k
      Jnu(z) = (z/2)^nu x Sum  ------------------
                          k=0  k! x Gamma(nu+k+1)
    (nu must be >= 0). Here k=15.
  ---------------------------------------------------*/
/*  complex_t cbessj(complex_t zz, int order) {
    std::complex<long double> z = zz;
    std::complex<long double> temp1 = pow(z / (long double) 2.0, order);
    std::complex<long double> z2 = - z * z / (long double) 4.0;
    std::complex<long double> sum(0.0, 0.0);
    long double factorial_k = 1.0;
    std::complex<long double> pow_z2_k = 1.0;
    //for(int k = 0; k <= MAXK; ++ k, pow_z2_k *= z2) {
    for(int k = 0; k <= 200; ++ k, pow_z2_k *= z2) {
      if(k == 0) factorial_k = 1.0;  // base case
      else factorial_k *= k;        // compute k!
      std::complex<long double> temp2 =
            pow_z2_k / (factorial_k * (long double) gamma((double) order + k + 1.0));
      sum += temp2;
    } // for
    temp1 *= sum;
    return complex_t((real_t) temp1.real(), (real_t) temp1.imag());
  } // cbessj()
*/

  // temporary fix. assuming imaginary component is 0
  complex_t cbessj(complex_t zz, int order) {
    return complex_t(j1(zz.real()), 0.0);
  } // cbessj()


#ifdef FF_CPU_OPT
  // TODO ...
  void cbessj_vec(const int VEC_LEN, complex_t* zz, int order, complex_t* res) {
    for(int i = 0; i < VEC_LEN; ++ i) res[i] = complex_t(j1(zz[i].real()), 0.0);
  } // cbessj()
#endif // FF_CPU_OPT

  real_t dot(vector3_t& a, vector3_t& b) {
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  } // dot()

  vector3_t cross(vector3_t & u, vector3_t & v) {
    vector3_t w;
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    return w;
  } // cross()

  complex_t sinc(complex_t x) {
    if(std::abs(x) > 1.0E-14) return std::sin(x) / x;
    else return (1.- x*x/6. + x*x*x*x/120.);  // taylor series approx
  } // sinc()

} // namespace hig
