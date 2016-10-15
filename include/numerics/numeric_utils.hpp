/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: numeric_utils.hpp
 *  Created: Oct 08, 2012
 *  Modified: Wed 08 Oct 2014 12:13:09 PM PDT
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

#ifndef NUMERIC_UTILS__H
#define NUMERIC_UTILS__H

#include <complex>
#include <cmath>

#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <common/cudafy.hpp>

namespace hig {

  extern double gamma(double x);
  extern complex_t cbessj(complex_t zz, int order);

  /* sign of a real value */
  inline real_t sign(real_t val) {
    if (val < 0) return -1;
    if (val > 0) return 1;
    return 0;
  }

  INLINE
  vector3_t cross(vector3_t & u, vector3_t & v){
    vector3_t w;
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    return w;
  }

  INLINE
  real_t dot(vector3_t & a, vector3_t & b){
   return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
  }

  // sinc function
  inline complex_t sinc(complex_t x){
    if (std::abs(x) > 1.0E-14)
      return std::sin(x) / x;
    else
      // Taylor series approx
      return (1.- x*x/6. + x*x*x*x/120.);
  }

} // namespace hig


#endif // NUMERIC_UTILS__H
