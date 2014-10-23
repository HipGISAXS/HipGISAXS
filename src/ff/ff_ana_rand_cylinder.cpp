/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_rand_cylinder.cpp
 *  Created: Jul 12, 2012
 *  Modified: Wed 22 Oct 2014 05:31:38 PM PDT
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
#include <common/enums.hpp>
#include <model/shape.hpp>
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#include <numerics/numeric_utils.hpp>

namespace hig {

  /**
   * random cylinders along z - for SAXS
   */
  bool AnalyticFormFactor::compute_random_cylinders(shape_param_list_t& params,
      std::vector<complex_t>& ff,  float_t tau, float_t eta, vector3_t transvec) {
    std::vector<float_t> r, distr_r;
    std::vector<float_t> h, distr_h;
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      switch((*i).second.type()) {
        case param_edge:
        case param_xsize:
        case param_ysize:
        case param_baseangle:
          std::cerr << "warning: ignoring unwanted parameters in random cylinders"
                << std::endl;
          break;
        case param_height:
          param_distribution((*i).second, h, distr_h);
          break;
        case param_radius:
          param_distribution((*i).second, r, distr_r);
          break;
        default:
          std::cerr << "error: unknown or invalid parameter given for random cylinders"
                << std::endl;
          return false;
      } // switch
    } // for

    if(r.size() < 1 || h.size() < 1) {
      std::cerr << "error: both radius and height parameters are required for random cylinders"
            << std::endl;
      return false;
    } // if

    #ifdef TIME_DETAIL_2
      woo::BoostChronoTimer maintimer;
      maintimer.start();
    #endif // TIME_DETAIL_2

    #ifdef FF_ANA_GPU
      // on gpu
      #ifdef FF_VERBOSE
        std::cout << "-- Computing random cylinders FF on GPU ..." << std::endl;
      #endif

      std::vector<float_t> transvec_v;
      transvec_v.push_back(transvec[0]);
      transvec_v.push_back(transvec[1]);
      transvec_v.push_back(transvec[2]);
      gff_.compute_random_cylinders(tau, eta, h, distr_h, r, distr_r, rot_, transvec_v, ff);
    #else
      // on cpu
      std::cout << "-- Computing random cylinder FF on CPU ..." << std::endl;

      ff.clear();
      ff.reserve(nqx_ * nqy_ * nqy_);
      for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));
      complex_t unitc(0, 1);

      float_t dx = 0.001;    // FIXME: hard-coded ???
      unsigned int nx = (unsigned int) ((1.0 - dx) / dx + 1.0);

      #pragma omp parallel for
      for(unsigned int z = 0; z < nqz_; ++ z) {
        complex_t qz = QGrid::instance().qz_extended(z);
        complex_t temp_ff(0.0, 0.0);
        for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
          for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
            float_t x_val = 0.0;
            complex_t temp_ffx(0.0, 0.0);
            for(unsigned int i_x = 0; i_x < nx; ++ i_x, x_val += dx) {
              complex_t temp1 = sinc(qz * h[i_h] * x_val / (float_t) 2.0);
              complex_t temp2 = qz * r[i_r] * sqrt(1 - x_val * x_val);
              if(!(temp2.real() == 0 && temp2.imag() == 0))  // TODO improve fp comparison ...
                temp_ffx += temp1 * cbessj(temp2, 1) / temp2;
              /*if(!(boost::math::isfinite(temp_ffx.real()) &&
                  boost::math::isfinite(temp_ffx.imag()))) {
                std::cout << "&&&&&&&&&&&&& OHO OHOH OHOHO: "
                      << z << ", " << i_x <<  ", " << i_r << ", " << i_h
                      << ": " << r[i_r] << ", " << h[i_h] << ", (" << temp1.real()
                      << "," << temp1.imag() << "), (" << temp2.real() << ","
                      << temp2.imag() << ")" << ", (" << qz.real() << "," << qz.imag()
                      << ")" << std::endl;
              } // if*/
            } // for
            temp_ff += 4.0 * distr_r[i_r] * distr_h[i_h] * temp_ffx * temp_ffx;
          } // for h
        } // for r
        // copy to all x and y
        for(unsigned int y = 0; y < nqy_; ++ y) {
          for(unsigned int x = 0; x < nqx_; ++ x) {
            unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
            ff[index] = temp_ff;
          } // for x
        } // for y
      } // for z
    #endif // FF_ANA_GPU
    #ifdef TIME_DETAIL_2
      maintimer.stop();
      std::cout << "** Rand Cylinder FF compute time: " << maintimer.elapsed_msec() << " ms."
            << std::endl;
    #endif // TIME_DETAIL_2
    return true;

  } // AnalyticFormFactor::compute_random_cylinders()

} // namespace hig
