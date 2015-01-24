/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_pyramid.cpp
 *  Created: Jul 12, 2012
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
   * pyramid
   */

  complex_t FormFactorPyramid (complex_t qx, complex_t qy, complex_t qz,
      float_t length, float_t width, float_t height, float_t angle) {

    float_t a = angle * PI_ / 180.;
    float_t tan_a = std::tan(a);
    if ((2*height/length >= tan_a) || (2*height/width >= tan_a))
        return C_ZERO;

    complex_t tmp = qx * qy;
    if (std::abs(tmp) < 1.0E-20)
        return C_ZERO;

    const complex_t P_J = C_ONE;
    const complex_t N_J = C_NEG_ONE;

    // q1,q2,q3,q4
    complex_t q1 = 0.5 * ((qx - qy)/tan_a + qz);
    complex_t q2 = 0.5 * ((qx - qy)/tan_a - qz);
    complex_t q3 = 0.5 * ((qx + qy)/tan_a + qz);
    complex_t q4 = 0.5 * ((qx + qy)/tan_a - qz);

    // k1,k2,k3,k4
    complex_t k1 = sinc(q1 * height) * std::exp(N_J * q1 * height) + 
        sinc(q2 * height) * std::exp(P_J * q2 * height);
    complex_t k2 = sinc(q1 * height) * std::exp(N_J * q1 * height) * N_J + 
        sinc(q2 * height) * std::exp(P_J * q2 * height) * P_J;
    complex_t k3 = sinc(q3 * height) * std::exp(N_J * q3 * height) + 
        sinc(q4 * height) * std::exp(P_J * q4 * height);
    complex_t k4 = sinc(q3 * height) * std::exp(N_J * q3 * height) * N_J +
        sinc(q4 * height) * std::exp(P_J * q4 * height) * P_J;

    // sins and cosines 
    complex_t t1 = k1 * std::cos (0.5 * (qx * length - qy * width));
    complex_t t2 = k2 * std::sin (0.5 * (qx * length - qy * width));
    complex_t t3 = k3 * std::cos (0.5 * (qx * length + qy * width));
    complex_t t4 = k4 * std::sin (0.5 * (qx * length + qy * width));

    // formfactor
    return ((height / tmp) * (t1 + t2 - t3 - t4));
  }


  bool AnalyticFormFactor::compute_pyramid(shape_param_list_t& params,
                            std::vector<complex_t>& ff,
                            float_t tau, float_t eta,
                            vector3_t transvec) {
    std::vector<float_t> x, distr_x;
    std::vector<float_t> y, distr_y;
    std::vector<float_t> h, distr_h;
    std::vector<float_t> b, distr_b;
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      switch((*i).second.type()) {
        case param_xsize:
          param_distribution((*i).second, x, distr_x);
          break;
        case param_ysize:
          param_distribution((*i).second, y, distr_y);
          break;
        case param_height:
          param_distribution((*i).second, h, distr_h);
          break;
        case param_baseangle:
          param_distribution((*i).second, b, distr_b);
          break;
        case param_edge:
        case param_radius:
          std::cerr << "warning: ignoring unwanted parameters in pyramid" << std::endl;
          break;
        default:
          std::cerr << "error: unknown/invalid parameter given for pyramid" << std::endl;
          return false;
      } // switch
    } // for

    #ifdef TIME_DETAIL_2
      woo::BoostChronoTimer maintimer;
      maintimer.start();
    #endif

    #ifdef FF_ANA_GPU
      std::cout << "-- Computing pyramid FF on GPU ..." << std::endl;
      std::vector<float_t> transvec_v; 
      transvec_v.push_back(transvec[0]);
      transvec_v.push_back(transvec[1]);
      transvec_v.push_back(transvec[2]);
      gff_.compute_pyramid(tau, eta, x, distr_x, y, distr_y, h, distr_h, b, distr_b, rot_, transvec_v, ff);
    #else

      std::cout << "-- Computing pyramid FF on CPU ..." << std::endl;
      ff.clear(); ff.resize(nqz_, C_ZERO);

      #pragma omp parallel for 
      for(int i = 0; i < nqz_; i++) {
        int j = i % nqy_;
        complex_t mqx, mqy, mqz;
        compute_meshpoints(QGrid::instance().qx(j), QGrid::instance().qy(j), 
            QGrid::instance().qz_extended(i), rot_, mqx, mqy, mqz);
        
        complex_t temp_ff(0.0, 0.0);
        for(int i_x = 0; i_x < x.size(); ++ i_x) {
          for(int i_y = 0; i_y < y.size(); ++ i_y) {
            for(int i_h = 0; i_h < h.size(); ++ i_h) {
              for(int i_b = 0; i_b < b.size(); ++ i_b) {
                float_t bb = b[i_b] * PI_ / 180;
                float_t prob = distr_x[i_x] * distr_y[i_y] * distr_h[i_h] * distr_b[i_b];
                temp_ff += FormFactorPyramid(mqx, mqy, mqz, x[i_x], y[i_y], h[i_h], b[i_b]) * prob;
              } // for b
            } // for h
          } // for y
        } // for x
        complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
        ff[i] = temp_ff * exp(complex_t(0, 1) * temp1);
      } // for z
    #endif // FF_ANA_GPU

    #ifdef TIME_DETAIL_2
      maintimer.stop();
      std::cout << "** Trunc Pyramid FF compute time: " << maintimer.elapsed_msec() << " ms."
            << std::endl;
    #endif

    return true;
  } // AnalyticFormFactor::compute_pyramid()
} // namespace hig

