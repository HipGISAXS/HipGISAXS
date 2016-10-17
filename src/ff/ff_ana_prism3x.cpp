/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_prism3x.cpp
 *  Created: Jul 12, 2012
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

#include <boost/math/special_functions/fpclassify.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

#include <ff/ff_ana.hpp>
#include <common/enums.hpp>
#include <common/constants.hpp>
#include <model/shape.hpp>
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#include <numerics/numeric_utils.hpp>

namespace hig {

  /**
   * triangular grating in the x-direction
   */
  bool AnalyticFormFactor::compute_prism3x(shape_param_list_t& params, std::vector<complex_t>& ff,
                        real_t tau, real_t eta, vector3_t transvec) {
    std::vector <real_t> lx, distr_lx;
    std::vector <real_t> ly, distr_ly;
    std::vector <real_t> h, distr_h;
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      switch((*i).second.type()) {
        case param_xsize:
          param_distribution((*i).second, lx, distr_lx);
          break;
        case param_ysize:
          param_distribution((*i).second, ly, distr_ly);
          break;
        case param_height:
          param_distribution((*i).second, h, distr_h);
          break;
        case param_edge:
        case param_radius:
        case param_baseangle:
          std::cerr << "warning: ignoring unwanted values for shape type 'prism3x'" << std::endl;
          break;
        default:
          std::cerr << "warning: unknown parameters for shape type 'prism3x'. ignoring"
                << std::endl;
      } // switch
    } // for

    // on cpu
    std::cout << "-- Computing prism3x FF on CPU ..." << std::endl;
    // initialize ff
    ff.clear();  ff.resize (nqz_, CMPLX_ZERO_);

    real_t d = 0.85;  // FIXME: hardcoded? variable?
    real_t gamma = 0.0;  // FIXME: hardcoded? variable?
    complex_t i(0.0, 1.0);

    #pragma omp parallel for
    for(unsigned int j_z = 0; j_z < nqz_; ++ j_z) {
      unsigned int j_y = j_z % nqy_;
      real_t temp_qx = QGrid::instance().qx(j_y);
      real_t temp_qy = QGrid::instance().qy(j_y);
      complex_t temp_qz = QGrid::instance().qz_extended(j_z);
      std::vector<complex_t> mq = rot_.rotate(temp_qx, temp_qy, temp_qz);
      real_t sg = sin(gamma);
      real_t cg = cos(gamma);
      real_t qx_rot = temp_qx * cg + temp_qy * sg;
      real_t qy_rot = temp_qy * cg - temp_qx * sg;
      complex_t temp_ff(0.0, 0.0);
      for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {        // H
        for(unsigned int i_y = 0; i_y < ly.size(); ++ i_y) {    // L
          for(unsigned int i_x = 0; i_x < lx.size(); ++ i_x) {  // Lx
            real_t temp_lx = lx[i_x] * 2, temp_ly = ly[i_y] * 2;// multiply by 2 (why?)
            real_t a1 = h[i_h] / (d * temp_ly);
            real_t b1 = 0.0;
            real_t a2 = h[i_h] / ((d - 1) * temp_ly);
            real_t b2 = h[i_h] / (1 - d);
            real_t temp1 = qx_rot * temp_lx / 2;
            real_t fqx = temp_lx;
            if(boost::math::fpclassify(qx_rot) != FP_ZERO) {
              fqx *= sin(temp1) / temp1;
            } // if
            complex_t k1 = temp_qz * a1 + qy_rot;
            complex_t k2 = temp_qz * a2 + qy_rot;
            complex_t i1 = exp(i * temp_qz * b1) * integral_e(0, d * temp_ly, k1);
            complex_t i2 = exp(i * temp_qz * b2) * integral_e(d * temp_ly, temp_ly, k2);
            complex_t i3 = integral_e(0, temp_ly, qy_rot);
            complex_t iy;
            if(boost::math::fpclassify(temp_qz.real()) == FP_ZERO &&
                boost::math::fpclassify(temp_qz.imag()) == FP_ZERO) {
              if(boost::math::fpclassify(qy_rot) == FP_ZERO) {
                iy = h[i_h] * temp_ly / 2;
              } else {
                iy = integral_xe(0, d * temp_ly, a1, b1, qy_rot) +
                    integral_xe(d * temp_ly, temp_ly, a2, b2, qy_rot);
              } // if-else
            } else {
              iy = (- i / temp_qz) * (i1 + i2 + i3);
            } // if-else
            temp_ff += fqx * iy;
          } // for i_x
        } // for i_y
      } // for i_h
      complex_t temp7 = (mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2]);
      ff[j_z] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
    } // for z

    return true;
  } // AnalyticFormFactor::compute_prism3x()

} // namespace hig
