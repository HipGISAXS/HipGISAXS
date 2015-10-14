/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_helpers.cpp
 *  Created: Apr 01, 2014
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
#include <boost/math/special_functions/fpclassify.hpp>

#include <sim/hipgisaxs_main.hpp>
#include <numerics/convolutions.hpp>

namespace hig {

  bool HipGISAXS::gaussian_smearing(real_t*& data, real_t sigma) {
    return Convolutions::instance().convolution_gaussian_2d(data, ncol_, nrow_, sigma);
  } // HipGISAXS::gaussian_smearing()

  bool HipGISAXS::check_finite(real_t* arr, unsigned int size) {
    for(unsigned int i = 0; i < size; ++ i) {
      if(!(boost::math::isfinite)(arr[i])) {
        std::cerr << "** ARRAY CHECK ** array entry not finite: " << i << std::endl;
        return false;
      } // if
    } // for
    return true;
  } // HipGISAXS::check_nan_inf()

  bool HipGISAXS::check_finite(complex_t* arr, unsigned int size) {
    bool res = true;
    for(unsigned int i = 0; i < size; ++ i) {
      if(!(boost::math::isfinite)(arr[i].real())) {
        std::cerr << "** ARRAY CHECK ** array entry not finite: " << i << ".real()" << std::endl;
        res = false;
      } // if
      if(!(boost::math::isfinite)(arr[i].imag())) {
        std::cerr << "** ARRAY CHECK ** array entry not finite: " << i << ".imag()" << std::endl;
        res = false;
      } // if
    } // for
    return res;
  } // HipGISAXS::check_nan_inf()

} // namespace hig
