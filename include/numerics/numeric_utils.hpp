/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: numeric_utils.hpp
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

#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <common/cudafy.hpp>

namespace hig {

  // vector products
  extern vector3_t cross(vector3_t & u, vector3_t & v);
  extern real_t dot(vector3_t & a, vector3_t & b);

  // special functions
  // extern double gamma(double x);
  extern real_t gamma(unsigned int x);
  extern complex_t cbessj(complex_t z, int order);
  extern complex_t cbesselj(complex_t z, int order);
  extern complex_t cj1(complex_t z);
  extern complex_t sinc(complex_t x);

#ifdef FF_CPU_OPT
  extern void cbessj_vec(const int, complex_t*, int, complex_t*);
  extern void sinc_vec(const int VEC_LEN, complex_t* x, complex_t* res);
#endif // FF_CPU_OPT

} // namespace hig
