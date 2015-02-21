/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: cuconstants.cuh
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

#ifndef __CUCONSTANTS_HPP__
#define __CUCONSTANTS_HPP__

#include <common/typedefs.hpp>
#include <common/constants.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>

namespace hig {

  const real_t CUTINY_ = 1e-20;

  // complex constants
  const cucomplex_t CUCMPLX_ZERO_ = make_cuC(REAL_ZERO_, REAL_ZERO_);
  const cucomplex_t CUCMPLX_ONE_  = make_cuC(REAL_ZERO_, REAL_ONE_);
  const cucomplex_t CUCMPLX_MINUS_ONE_ = make_cuC(REAL_ZERO_, REAL_MINUS_ONE_);

} // namespace

#endif // __CUCONSTANTS_HPP__
