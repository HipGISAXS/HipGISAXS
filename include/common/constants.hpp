/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: constants.hpp
 *  Created: Aug 25, 2012
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

#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

#include <limits>
#include <common/typedefs.hpp>

namespace hig {

  const real_t TINY_ = std::numeric_limits<real_t>::min();
  const real_t CUTINY_ = 1.0e-13;

  const unsigned int LIGHT_SPEED_ = 3e+8;   /* speed of light in m/s */
  const real_t PI_ = 3.141592653589793;     /* PI correct upto 15 decimal places */
  const real_t SQRT_2PI_ = 2.506628;

  const unsigned int MAX_DEPTH_ = 150;      /* maximum depth allowed */

  // real constants
  const real_t REAL_ZERO_ = (real_t) 0.;
  const real_t REAL_ONE_  = (real_t) 1.;
  const real_t REAL_MINUS_ONE_ = (real_t) -1.;

  // complex constants
  const complex_t CMPLX_ZERO_ = complex_t(REAL_ZERO_, REAL_ZERO_);
  const complex_t CMPLX_ONE_  = complex_t(REAL_ZERO_, REAL_ONE_);
  const complex_t CMPLX_MINUS_ONE_ = complex_t(REAL_ZERO_, REAL_MINUS_ONE_);

  // miscellaneous helpers
  // const int VERBOSE_LEVEL_ONE = 1;

  // fitting related constants
  const int MAX_ITER_REG_ = 10;

} // namespace

#endif // __CONSTANTS_HPP__
