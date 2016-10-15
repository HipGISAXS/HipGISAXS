/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_prism6.cpp
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
   * six faceted prism
   */
  bool AnalyticFormFactor::compute_prism6(shape_param_list_t& params, std::vector<complex_t>& ff,
                      real_t tau, real_t eta, vector3_t transvec) {
    std::vector<real_t> l, distr_l;
    std::vector<real_t> h, distr_h;
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      switch((*i).second.type()) {
        case param_edge:
          param_distribution((*i).second, l, distr_l);
          break;
        case param_height:
          param_distribution((*i).second, h, distr_h);
          break;
        case param_xsize:
        case param_ysize:
        case param_baseangle:
        case param_radius:
          std::cerr << "warning: ignoring unwanted parameters in prism6" << std::endl;
          break;
        default:
          std::cerr << "error: invalid parameter given for prism6" << std::endl;
          return false;
      } // switch
    } // for

#ifdef TIME_DETAIL_2
    woo::BoostChronoTimer maintimer;
    maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
    // on gpu
    #ifdef FF_VERBOSE
      std::cout << "-- Computing prism6 FF on GPU ..." << std::endl;
    #endif

    std::vector<real_t> transvec_v;
    transvec_v.push_back(transvec[0]);
    transvec_v.push_back(transvec[1]); 
    transvec_v.push_back(transvec[2]);
    gff_.compute_prism6(tau, eta, l, distr_l, h, distr_h, rot_, transvec_v, ff);
#else
    // on cpu
    std::cout << "-- Computing prism6 FF on CPU ..." << std::endl;

    ff.clear(); ff.resize(nqz_, CMPLX_ZERO_);
    real_t sqrt3 = sqrt(3.0);

    #pragma omp parallel for 
    for(unsigned int z = 0; z < nqz_; ++ z) {
      unsigned int y = z % nqy_;
      std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), 
              QGrid::instance().qy(y), QGrid::instance().qz_extended(z));
      complex_t qm = tan(tau) * (mq[0] * sin(eta) + mq[1] * cos(eta));
      complex_t temp1 = ((real_t) 4.0 * sqrt3) / (3.0 * mq[1] * mq[1] - mq[0] * mq[0]);
      complex_t temp_ff(0.0, 0.0);
      for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
        for(unsigned int i_l = 0; i_l < l.size(); ++ i_l) {
          complex_t rmqx = l[i_l] * mq[0] / sqrt3;
          complex_t rmqy = l[i_l] * mq[1];
          complex_t temp2 = rmqy * rmqy *  sinc(rmqx) * sinc(rmqy);
          complex_t temp3 = cos(2.0 * rmqx);
          complex_t temp4 = cos(rmqy) * cos(rmqx);
          complex_t temp5 = temp1 * (temp2 + temp3 - temp4);
          complex_t temp6 = fq_inv(mq[2] + qm, h[i_h]);
          temp_ff += temp5 * temp6;
        } // for l
      } // for h
      complex_t temp7 = (mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2]);
      ff[z] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
    } // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
    maintimer.stop();
    std::cout << "**        Prism6 FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2
    return true;
  } // AnalyticFormFactor::compute_prism6()

} // namespace hig
