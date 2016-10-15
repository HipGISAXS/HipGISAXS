/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_sphere.cpp
 *  Created: Jul 12, 2012
 *  Modified: Wed 22 Oct 2014 05:31:59 PM PDT
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
   * sphere
   */
  complex_t FormFactorSphere(complex_t qx, complex_t qy, complex_t qz, real_t radius){
      complex_t qval = std::sqrt(qx * qx + qy * qy + qz * qz);
      if (std::abs(qval) < TINY_){
        return CMPLX_ZERO_;
      }
      complex_t qR   = qval * radius;
      complex_t qR3   = qR * qR * qR;
      complex_t t1   = std::conj(qR3) / std::norm(qR3);
      complex_t c0   = (std::sin(qR) - qR * std::cos(qR)) * t1;
      complex_t c1   = std::exp(CMPLX_ONE_ * qz * radius);
      real_t   f0   = 4 * PI_ * radius * radius * radius;
      return (f0 * c0 * c1);            
  }

  bool AnalyticFormFactor::compute_sphere(shape_param_list_t& params, std::vector<complex_t> &ff,
                      vector3_t transvec) {
    std::vector<real_t> r, distr_r;
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      switch((*i).second.type()) {
        case param_edge:
        case param_xsize:
        case param_ysize:
        case param_height:
        case param_baseangle:
          std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
          break;
        case param_radius:
          param_distribution((*i).second, r, distr_r);
          break;
        default:
          std::cerr << "error: unknown or invalid parameter given for sphere" << std::endl;
          return false;
      } // switch
    } // for

    if(r.size() < 1) {
      std::cerr << "error: radius parameter required for sphere" << std::endl;
      return false;
    } // if

#ifdef TIME_DETAIL_2
    woo::BoostChronoTimer maintimer;
    maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
    // on gpu
    #ifdef FF_VERBOSE
      std::cerr << "-- Computing sphere FF on GPU ..." << std::endl;
    #endif

    std::vector<real_t> transvec_v;
    transvec_v.push_back(transvec[0]);
    transvec_v.push_back(transvec[1]);
    transvec_v.push_back(transvec[2]);

    gff_.compute_sphere(r, distr_r, rot_, transvec_v, ff);
#else
    // on cpu
    std::cerr << "-- Computing sphere FF on CPU ..." << std::endl;

    ff.clear(); ff.resize(nqz_, CMPLX_ZERO_);

    #pragma omp parallel for 
    for(unsigned int z = 0; z < nqz_; ++ z) {
      unsigned int y = z % nqy_;
      std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), QGrid::instance().qy(y),
              QGrid::instance().qz_extended(z));
      complex_t temp_ff = CMPLX_ZERO_;
      for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
          temp_ff += distr_r[i_r] * FormFactorSphere(mq[0], mq[1], mq[2], r[i_r]);
      } // for r
      complex_t temp1 = mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2];
      complex_t temp2 = std::exp(complex_t(-temp1.imag(), temp1.real()));
      ff[z] = temp_ff * temp2;
    } // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
    maintimer.stop();
    std::cerr << "**        Sphere FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2
    
    return true;
  } // AnalyticFormFactor::compute_sphere()
} // namespace hig
