/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_trunc_cone.cpp
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
#include <common/constants.hpp>
#include <common/enums.hpp>
#include <model/shape.hpp>
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#include <numerics/numeric_utils.hpp>

namespace hig {

  /**
   * truncated cone
   */
  complex_t FormFactorTrucatedCone(complex_t qx, complex_t qy, complex_t qz, 
          real_t radius, real_t height, real_t angle){
      complex_t ff = CMPLX_ZERO_;
      if ( height / radius >= tan(angle) )
        return ff;
   
      // else compute form-factor
      complex_t q_par = sqrt(qx * qx + qy * qy);
     return ff;
  }

  bool AnalyticFormFactor::compute_truncated_cone(shape_param_list_t& params, real_t tau, real_t eta,
                          std::vector<complex_t>& ff, vector3_t transvec) {
    std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_cone"
                          << std::endl;
    return false;

    // this is not working ... something is not right in matlab code ...

    std::vector <real_t> h, distr_h;  // for h dimension: param_height
    std::vector <real_t> r, distr_r;  // for r dimension: param_radius
    std::vector <real_t> a, distr_a;  // for a angle: param_baseangle
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      if(!(*i).second.isvalid()) {
        std::cerr << "warning: ignoring invalid shape parameter" << std::endl;
        continue;
      } // if
      switch((*i).second.type()) {
        case param_edge:
        case param_xsize:
        case param_ysize:
          std::cerr << "warning: ignoring unwanted input parameters for 'truncated cone'" << std::endl;
          break;
        case param_baseangle:
          param_distribution((*i).second, a, distr_a);
          break;
        case param_height:
          param_distribution((*i).second, h, distr_h);
          break;
        case param_radius:
          param_distribution((*i).second, r, distr_r);
          break;
        default:
          std::cerr << "error: unknown parameter type given" << std::endl;
          return false;
      } // switch
    } // for

    if(h.size() < 1 || r.size() < 1 || a.size() < 1) {
      std::cerr << "error: missing parameters for truncated cone" << std::endl;
      return false;
    } // if

    // on cpu
    std::cout << "-- Computing truncated cone FF on CPU ..." << std::endl;

    ff.clear();
    ff.resize(nqz_, CMPLX_ZERO_);

    unsigned int nz = 40;  // FIXME: hard-coded ... what is this???
    for(unsigned z = 0; z < nqz_; ++ z) {
      unsigned y = z % nqy_;
      std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), QGrid::instance().qy(y),
                QGrid::instance().qz_extended(z));

      complex_t qpar = sqrt(mq[0] * mq[0] + mq[1] * mq[1]);
      complex_t temp_ff(0.0, 0.0);

      for(unsigned int i_a = 0; i_a < a.size(); ++ i_a) {
        real_t temp1 = tan(a[i_a]);
        for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
          real_t dz = h[i_h] / (real_t)(nz - 1);
          for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
            real_t z_val = 0.0;
            complex_t temp_ffz(0.0, 0.0);
            for(unsigned int i_z = 0; i_z < nz; ++ i_z, z_val += dz) {
              real_t rz = r[i_r] - z_val / temp1;
              complex_t temp2 = exp(complex_t(-(mq[2] * rz).imag(), (mq[2] * rz).real()));
              complex_t temp3 = cbessj(qpar * rz, 1) / (qpar * rz);
              temp_ffz += rz * rz * temp2 * temp3;

                if((boost::math::fpclassify(temp_ffz.real()) == FP_ZERO &&
                  boost::math::fpclassify(temp_ffz.imag()) == FP_ZERO) ||
                  !boost::math::isfinite(temp_ffz.real()) ||
                  !boost::math::isfinite(temp_ffz.imag())) {
                  std::cout << "--------------------here it is: " << a[i_a] << ", "
                        << r[i_r] << ", " << h[i_h] << ", " << y
                        << ", " << z << ", " << i_z << std::endl;
                  exit(1);
              } // if
            } // for
            temp_ff += 2 * PI_ * distr_r[i_r] * distr_h[i_h] * distr_a[i_a] * temp_ffz;
          } // for r
        } // for h
      } // for a
      complex_t temp1 = mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2];
      complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
      ff[z] = temp_ff * temp2;
    } // for z

  } // AnalyticFormFactor::compute_truncated_cone()

} // namespace hig
