/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_prism.cpp
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
   * prism - 3 face
   */
  bool AnalyticFormFactor::compute_prism(shape_param_list_t& params, std::vector<complex_t>& ff,
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
        case param_radius:
        case param_baseangle:
          std::cerr << "warning: ignoring unwanted parameters in prism" << std::endl;
          break;
        default:
          std::cerr << "error: invalid parameter given for prism" << std::endl;
          return false;
      } // switch
    } // for

    if(h.size() < 1 || l.size() < 1) {
      std::cerr << "error: need edge and height for prism shape" << std::endl;
      return false;
    } // if

    // why not doing for a range of r, h for this? ... ???

#ifdef TIME_DETAIL_2
    woo::BoostChronoTimer maintimer;
    maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
    // on gpu
    #ifdef FF_VERBOSE
      std::cout << "-- Computing prism3 FF on GPU ..." << std::endl;
    #endif

    std::vector<real_t> transvec_v;
    transvec_v.push_back(transvec[0]);
    transvec_v.push_back(transvec[1]);
    transvec_v.push_back(transvec[2]);
    gff_.compute_prism3(tau, eta, l, distr_l, h, distr_h, rot_, transvec_v, ff);
#else
    // on cpu
    std::cout << "-- Computing prism3 FF on CPU ..." << std::endl;

    real_t sqrt3 = sqrt(3.0);
    complex_t unitc(0, 1.0);
    ff.clear(); ff.resize(nqz_, CMPLX_ZERO_);

    #pragma omp parallel for 
    for(unsigned int z = 0; z < nqz_; ++ z) {
      unsigned int y = z % nqy_;
      std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), 
              QGrid::instance().qy(y), QGrid::instance().qz_extended(z));
      complex_t mqx = mq[0];
      complex_t mqy = mq[1];
      complex_t mqz = mq[2];

      complex_t temp1 = mqx * (mqx * mqx - 3.0 * mqy * mqy);
      complex_t temp2 = mqz + tan(tau) * (mqx * sin(eta) + mqy * cos(eta));
      complex_t temp_ff(0.0, 0.0);
      for(unsigned int i_l = 0; i_l < l.size(); ++ i_l) {
        for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
          complex_t temp4 = mqx * exp(unitc * mqy * l[i_l] * sqrt3);
          complex_t temp5 = mqx * cos(mqx * l[i_l]);
          complex_t temp6 = unitc * sqrt3 * mqy * sin(mqx * l[i_l]);
          complex_t temp7 = fq_inv(temp2, h[i_h]);
          complex_t temp8 = (temp4 - temp5 - temp6) * temp7;
          complex_t temp3 = mqy * l[i_l] / sqrt3;
          complex_t temp9 = 2.0 * sqrt3 * exp(complex_t(temp3.imag(), -temp3.real()));
          temp_ff += (temp9 / temp1) * temp8;
        } // for h
      } // for r
      complex_t temp10 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
      ff[z] = temp_ff * exp(complex_t(-temp10.imag(), temp10.real()));
    } // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
    maintimer.stop();
    std::cout << "**        Prism3 FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2

    return true;
  } // AnalyticFormFactor::compute_prism()

} // namespace hig
