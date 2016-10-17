/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_box.cpp
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
   * box
   */

  inline complex_t FormFactorBox(complex_t qx, complex_t qy, complex_t qz, 
          real_t length, real_t width, real_t height) {
      real_t vol = length * width * height;
      complex_t exp_val = std::exp(CMPLX_ONE_ * 0.5 * qz * height);
      complex_t sinc_val = sinc(0.5 * qx * length) * sinc(0.5 * qy * width) 
          * sinc(0.5 * qz * height);
      return (vol * exp_val * sinc_val);
  }

  bool AnalyticFormFactor::compute_box(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                    std::vector<complex_t>& ff,
                    ShapeName shape, shape_param_list_t& params,
                    real_t tau, real_t eta, vector3_t &transvec){
    std::vector <real_t> x, distr_x;  // for x dimension: param_xsize  param_edge
    std::vector <real_t> y, distr_y;  // for y dimension: param_ysize  param_edge
    std::vector <real_t> z, distr_z;  // for z dimension: param_height param_edge
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      if(!(*i).second.isvalid()) {
        std::cerr << "warning: invalid shape parameter found" << std::endl;
        continue;
      } // if
      switch((*i).second.type()) {
        case param_edge:
          param_distribution((*i).second, x, distr_x);  // x == RR, distr_x == RRD
          param_distribution((*i).second, y, distr_y);
          param_distribution((*i).second, z, distr_z);
          break;
        case param_xsize:
          param_distribution((*i).second, x, distr_x);
          break;
        case param_ysize:
          param_distribution((*i).second, y, distr_y);
          break;
        case param_height:
          param_distribution((*i).second, z, distr_z);
          break;
        case param_radius:
        case param_baseangle:
          std::cerr << "warning: ignoring unwanted values for shape type 'box'" << std::endl;
          break;
        default:
          std::cerr << "warning: unknown parameters for shape type 'box'. ignoring"
                << std::endl;
      } // switch
    } // for

    // check if x y z etc are set or not
    if(x.size() < 1 || y.size() < 1 || z.size() < 1) {
      std::cerr << "error: invalid or not enough box parameters given" << std::endl;
      return false;
    } // if

    #ifdef TIME_DETAIL_2
      woo::BoostChronoTimer maintimer;
      maintimer.start();
    #endif // TIME_DETAIL_2
    #ifdef FF_ANA_GPU
      // on gpu
      #ifdef FF_VERBOSE
        std::cerr << "-- Computing box FF on GPU ..." << std::endl;
      #endif

      std::vector<real_t> transvec_v;
      transvec_v.push_back(transvec[0]);
      transvec_v.push_back(transvec[1]);
      transvec_v.push_back(transvec[2]);
      gff_.compute_box(tau, eta, x, distr_x, y, distr_y, z, distr_z, rot_, transvec_v, ff);
    #else
      // on cpu
      std::cerr << "-- Computing box FF on CPU ..." << std::endl;
      // initialize ff
      ff.clear();  ff.resize(nqz, CMPLX_ZERO_);

      #pragma omp parallel for 
      for(unsigned int i = 0; i < nqz; ++ i) {
        unsigned int j = i % nqy;
        std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(j), 
                QGrid::instance().qy(j), QGrid::instance().qz_extended(i));
        complex_t temp_ff(0.0, 0.0);
        for(unsigned int i_z = 0; i_z < z.size(); ++ i_z) {
          for(unsigned int i_y = 0; i_y < y.size(); ++ i_y) {
            for(unsigned int i_x = 0; i_x < x.size(); ++ i_x) {
              real_t wght = distr_x[i_x] * distr_y[i_y] * distr_z[i_z];
              temp_ff += FormFactorBox(mq[0], mq[1], mq[2], x[i_x], y[i_y], z[i_z]) * wght;
            } // for i_x
          } // for i_y
        } // for i_z
        complex_t temp7 = (mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2]);
        /*if(!(boost::math::isfinite(temp7.real()) && boost::math::isfinite(temp7.imag()))) {
          std::cerr << "---------------- here it is ------ " << j_x << ", "
                << j_y << ", " << j_z << std::endl;
          exit(1);
       } // if*/
       ff[i] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
       /*
       if(!(boost::math::isfinite(ff[curr_index].real()) &&
          boost::math::isfinite(ff[curr_index].imag()))) {
        std::cerr << "******************* here it is ********* " << j_x << ", "
              << j_y << ", " << j_z << std::endl;
        exit(1);
       } //if 
            */
      } // for i
    #endif // FF_ANA_GPU
    #ifdef TIME_DETAIL_2
      maintimer.stop();
      std::cerr << "**           Box FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
    #endif // TIME_DETAIL_2

    return true;
  } // AnalyticFormFactor::compute_box()
} // namespace hig
