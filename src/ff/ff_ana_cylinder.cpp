/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_cylinder.cpp
 *  Created: Jul 12, 2012
 *  Modified: Wed 08 Oct 2014 06:20:29 PM PDT
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
   * cylinder
   */
  complex_t FormFactorCylinder(complex_t qpar, complex_t qz, real_t radius, real_t height){
    real_t vol = 2. * PI_ * radius * radius * height;
    complex_t sinc_val = sinc(0.5 * qz * height);
    complex_t expt_val = std::exp(C_ONE * 0.5 * qz * height);
    complex_t t1 = qpar * radius;
    complex_t bess_val = cbessj(t1, 1) * std::conj(t1) / std::norm(t1);
    return (vol * sinc_val * expt_val * bess_val);
  }  

  bool AnalyticFormFactor::compute_cylinder(shape_param_list_t& params, real_t tau, real_t eta,
                                            std::vector<complex_t>& ff, vector3_t transvec) {
    std::vector <real_t> h, distr_h;  // for h dimension: param_height
    std::vector <real_t> r, distr_r;  // for r dimension: param_radius
    for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
      if(!(*i).second.isvalid()) {
        std::cerr << "warning: ignoring invalid shape parameter" << std::endl;
        continue;
      } // if
      switch((*i).second.type()) {
        case param_edge:
        case param_xsize:
        case param_ysize:
        case param_baseangle:
          std::cerr << "warning: ignoring unwanted input parameters for 'cylinder'" << std::endl;
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

    if(h.size() < 1 || r.size() < 1) {
      std::cerr << "error: missing parameters for cylinder" << std::endl;
      return false;
    } // if

#ifdef TIME_DETAIL_2
    woo::BoostChronoTimer maintimer;
    maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
    // on gpu
    std::cout << "-- Computing cylinder FF on GPU ..." << std::endl;

    std::vector<real_t> transvec_v;
    transvec_v.push_back(transvec[0]);
    transvec_v.push_back(transvec[1]);
    transvec_v.push_back(transvec[2]);
    gff_.compute_cylinder(tau, eta, h, distr_h, r, distr_r, rot_, transvec_v, ff);
#else
    // on cpu
    std::cout << "-- Computing cylinder FF on CPU ..." << std::endl;

    ff.clear(); ff.resize(nqz_, complex_t(0.,0.));

    #pragma omp parallel for 
    for(unsigned z = 0; z < nqz_; ++ z) {
      unsigned y = z % nqy_; 
      complex_t mqx, mqy, mqz;
      compute_meshpoints(QGrid::instance().qx(y), QGrid::instance().qy(y), 
          QGrid::instance().qz_extended(z), rot_, mqx, mqy, mqz);
      complex_t qpar = sqrt(mqx * mqx + mqy * mqy);
      complex_t temp_ff(0.0, 0.0);
      for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
        for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
          temp_ff += FormFactorCylinder(qpar, mqz, r[i_r], h[i_h]);
        } // for h
      } // for r
      complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
      complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
      ff[z] = temp_ff * temp2;
    } // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
    maintimer.stop();
    std::cout << "**      Cylinder FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2

    return true;
  } // AnalyticFormFactor::compute_cylinder()

} // namespace hig
