/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_trunc_pyramid.cpp
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
   * truncated pyramid
   */
  bool AnalyticFormFactor::compute_truncated_pyramid(shape_param_list_t& params,
                            std::vector<complex_t>& ff,
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
          std::cerr << "warning: ignoring unwanted parameters in truncated pyramid" << std::endl;
          break;
        default:
          std::cerr << "error: unknown/invalid parameter given for truncated pyramid" << std::endl;
          return false;
      } // switch
    } // for

    //float_t tr = 0.4;  // FIXME: hardcoded ... ???

    #ifdef TIME_DETAIL_2
      woo::BoostChronoTimer maintimer;
      maintimer.start();
    #endif
      std::cout << "-- Computing truncated pyramid FF on CPU ..." << std::endl;
      ff.clear(); ff.reserve(nqz_);
      for(int i = 0; i < nqz_; ++ i) ff.push_back(complex_t(0, 0));

      #pragma omp parallel for 
      for(int zi = 0; zi < nqz_; ++ zi) {
        int yi = zi % nqy_;
        complex_t mqx, mqy, mqz;
        compute_meshpoints(QGrid::instance().qx(yi), QGrid::instance().qy(yi), 
            QGrid::instance().qz_extended(zi), rot_, mqx, mqy, mqz);
        
        complex_t temp_ff(0.0, 0.0);
        for(int i_x = 0; i_x < x.size(); ++ i_x) {
          for(int i_y = 0; i_y < y.size(); ++ i_y) {
            for(int i_h = 0; i_h < h.size(); ++ i_h) {
              for(int i_b = 0; i_b < b.size(); ++ i_b) {
                float_t xx = x[i_x] ;// * (1 - tr);
                float_t hh = h[i_h];// * (1 - tr);
                float_t yy = y[i_y] ;  // * (1 - tr);
                float_t bb = b[i_b] * PI_ / 180;
                float_t prob = distr_x[i_x] * distr_y[i_y] * distr_h[i_h] * distr_b[i_b];
                temp_ff += prob * truncated_pyramid_core(mqx, mqy, mqz, xx, yy, hh, bb);
                //temp_ff += prob * ffTruncatedPyramid (mqx, mqy, mqz, xx, yy, hh, bb);
              } // for b
            } // for h
          } // for y
        } // for x

        complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
        ff[zi] = temp_ff * exp(complex_t(0, 1) * temp1);
      } // for z
    //#endif
    #ifdef TIME_DETAIL_2
      maintimer.stop();
      std::cout << "** Trunc Pyramid FF compute time: " << maintimer.elapsed_msec() << " ms."
            << std::endl;
    #endif

    return true;
  } // AnalyticFormFactor::compute_truncated_pyramid()



  complex_t AnalyticFormFactor::ffTruncatedPyramid (complex_t rqx, complex_t rqy, complex_t rqz,
      float_t x, float_t y, float_t h, float_t ang) {
    float_t tan_a = tan(ang);
    complex_t jp = complex_t (0 ,1);
    complex_t jm = complex_t (0, -1);

    // q1,q2,q3,q4
    complex_t q1 = 0.5 * (((rqx - rqy)/tan_a) + rqz);
    complex_t q2 = 0.5 * (((rqx - rqy)/tan_a) - rqz);
    complex_t q3 = 0.5 * (((rqx + rqy)/tan_a) + rqz);
    complex_t q4 = 0.5 * (((rqx + rqy)/tan_a) - rqz);

    // k1,k2,k3,k4
    complex_t k1 = exp (jm * q2 * h) * sinc (q2 * h) + exp(jp * q1 * h) * sinc(q1 * h);
    complex_t k2 = jp * exp (jm * q2 * h) * sinc (q2 * h) - jp * exp(jp * q1 * h) * sinc(q1 * h);
    complex_t k3 = exp (jm * q4 * h) * sinc (q4 * h) + exp(jp * q3 * h) * sinc(q3 * h);
    complex_t k4 = jp * exp (jm * q4 * h) * sinc (q4 * h) - jp * exp(jp * q3 * h) * sinc(q3 * h);

    // sins and cosines 
    complex_t t1 = k1 * cos (rqx * x * 0.5 - rqy * y * 0.5);
    complex_t t2 = k2 * sin (rqx * x * 0.5 - rqy * y * 0.5);
    complex_t t3 = k3 * cos (rqx * x * 0.5 + rqy * y * 0.5);
    complex_t t4 = k4 * sin (rqx * x * 0.5 + rqy * y * 0.5);

    // formfactor
    complex_t val = complex_t (0.0, 0.0);
    complex_t tmp = rqx * rqy;
    if (std::norm (tmp) > std::numeric_limits<float_t>::epsilon()) {
      val = (h / tmp) * (t1 + t2 + t3 + t4);
    }
    return val;
  }

  complex_t AnalyticFormFactor::truncated_pyramid_core(complex_t rqx, complex_t rqy, complex_t rqz, 
      float_t x, float_t y, float_t h, float_t bang) {
    float_t tanbang = tan(bang);
    complex_t qxy_s = (rqx + rqy) / tanbang;
    complex_t qxy_d = (rqx - rqy) / tanbang;
    complex_t q1 = (rqz + qxy_d) / (float_t) 2.0;
    complex_t q2 = (-rqz + qxy_d) / (float_t) 2.0;
    complex_t q3 = (rqz + qxy_s) / (float_t) 2.0;
    complex_t q4 = (-rqz + qxy_s) / (float_t) 2.0;

    /* error in sinc expression fixed 03/06/2014 -Slim   */
    complex_t sinc_eiq1 = sinc(q1 * h) * exp(  complex_t(0, 1) * q1 * h);
    complex_t sinc_eiq2 = sinc(q2 * h) * exp(- complex_t(0, 1) * q2 * h);
    complex_t sinc_eiq3 = sinc(q3 * h) * exp(  complex_t(0, 1) * q3 * h);
    complex_t sinc_eiq4 = sinc(q4 * h) * exp(- complex_t(0, 1) * q4 * h);

    complex_t k1 = sinc_eiq1 + sinc_eiq2;
    complex_t k2 = - complex_t(0, 1) * (sinc_eiq1 - sinc_eiq2);
    complex_t k3 = sinc_eiq3 + sinc_eiq4;
    complex_t k4 = - complex_t(0, 1) * (sinc_eiq3 - sinc_eiq4);

    complex_t qxqy = rqx * rqy;
    complex_t qxqy_conj = complex_t(qxqy.real() , - qxqy.imag() );
    float_t qxqy_amp2 = qxqy.real()*qxqy.real() + qxqy.imag() * qxqy.imag();
    complex_t qxqy_inv = qxqy_conj /qxqy_amp2;

    float mach_eps = std::numeric_limits<float>::epsilon() ; //move to constants.hpp    
    complex_t val = complex_t( 0.0 ,  0.0 );
    complex_t tmp =  rqx * rqy;

    if(std::abs (tmp) > mach_eps) {
      val = (h / (rqx * rqy)) * (
             cos(rqx * x - rqy * y) * k1 +
             sin(rqx * x - rqy * y) * k2 -
             cos(rqx * x + rqy * y) * k3 -
             sin(rqx * x + rqy * y) * k4
                 );
    }
     // std::cout << "FF = "  << val << std::endl;
    return val;
  } // AnalyticFormFactor::truncated_pyramid_core()

} // namespace hig

