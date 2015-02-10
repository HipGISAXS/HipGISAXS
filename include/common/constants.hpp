/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: constants.hpp
 *  Created: Aug 25, 2012
 *  Modified: Wed 08 Oct 2014 12:13:01 PM PDT
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

#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

#include <common/typedefs.hpp>

namespace hig {

#ifdef DOUBLEP
  const double TINY = 1.0E-18;
#else
  const float  TINY = 1.0E-18;
#endif

  const unsigned int LIGHT_SPEED_ = 3e+8;    /* speed of light in m/s */
  const real_t PI_ = 3.141592653589793;    /* PI correct upto 15 decimal places */
  const real_t SQRT_2PI_ = 2.506628;

  //const unsigned int MAX_DEPTH_ = 500;    /* maximum depth allowed */
  const unsigned int MAX_DEPTH_ = 150;    /* maximum depth allowed */

  // real constants
  const real_t ZERO = (real_t) 0.;
  const real_t ONE  = (real_t) 1.;
  const real_t NEG_ONE = (real_t) -1.;

  // complex constants
  const complex_t C_ZERO = complex_t(ZERO,ZERO);
  const complex_t C_ONE  = complex_t(ZERO, ONE);
  const complex_t C_NEG_ONE = complex_t(ZERO, NEG_ONE);

} // namespace

#endif // __CONSTANTS_HPP__
