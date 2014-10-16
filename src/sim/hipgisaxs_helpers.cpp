/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_helpers.cpp
 *  Created: Apr 01, 2014
 *  Modified: Wed 08 Oct 2014 12:17:46 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
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

#include <sim/hipgisaxs_main.hpp>
#include <numerics/convolutions.hpp>

namespace hig {

  bool HipGISAXS::gaussian_smearing(float_t*& data, float_t sigma) {
    return Convolutions::instance().convolution_gaussian_2d(data, nqy_, nqz_, sigma);
  } // HipGISAXS::gaussian_smearing()

} // namespace hig
