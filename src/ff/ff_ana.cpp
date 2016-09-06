/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana.cpp
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

  // TODO: decompose into two init functions:
  //   one for overall (sets qgrid in gff_),
  //   other for each run/invocation (sets rotation matrices)
  bool AnalyticFormFactor::init(RotMatrix_t & rot, std::vector<complex_t> &ff) {
    nqx_ = QGrid::instance().nqx();
    nqy_ = QGrid::instance().nqy();
    nqz_ = QGrid::instance().nqz_extended();

    // first make sure there is no residue from any previous computations
    ff.clear();

    #ifdef FF_ANA_GPU
      gff_.init(nqy_, nqz_);
    #endif // FF_ANA_GPU

    // rotation matrices are new for each ff calculation
    rot_ = rot;

    return true;
  } // AnalyticFormFactor::init()


  void AnalyticFormFactor::clear() {
    nqx_ = nqy_ = nqz_ = 0;
    #ifdef FF_ANA_GPU
      gff_.clear();
    #endif // FF_ANA_GPU
  } // AnalyticFormFactor::clear()


  bool AnalyticFormFactor::compute(ShapeName shape, real_t tau, real_t eta, vector3_t transvec,
                                    std::vector<complex_t>& ff,
                                    shape_param_list_t& params, real_t single_layer_thickness,
                                    RotMatrix_t & rot
                                    #ifdef USE_MPI
                                      , woo::MultiNode& world_comm, std::string comm_key
                                    #endif
                                    ) {

    #ifdef FF_VERBOSE
      std::cerr << "-- Computing form factor analytically ... " << std::endl;
    #endif
//    #ifdef TIME_DETAIL_1
      woo::BoostChronoTimer compute_timer;
      compute_timer.start();
//    #endif // TIME_DETAIL_1

    switch(shape) {
      case shape_cube:            // cube 
        if(!compute_cube(nqx_, nqy_, nqz_, ff, params, tau, eta, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a cube"
                << std::endl;
          return false;
        } // if
        break;
      case shape_box:            // box
        if(!compute_box(nqx_, nqy_, nqz_, ff, shape, params, tau, eta, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a box"
                << std::endl;
          return false;
        } // if
        break;
      case shape_cylinder:        // standing cylinder
        if(!compute_cylinder(params, tau, eta, ff, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a cylinder"
                << std::endl;
          return false;
        } // if
        break;
      case shape_sphere:          // simple sphere
        if(!compute_sphere(params, ff, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a sphere"
                << std::endl;
          return false;
        } // if
        break;
      case shape_prism3:          // triangular prism (prism with 3 sides)
        if(!compute_prism(params, ff, tau, eta, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a prism3"
                << std::endl;
          return false;
        } // if
        break;
      case shape_prism6:          // hexagonal prism (prism with 6 sides)
        if(!compute_prism6(params, ff, tau, eta, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a prism6"
                << std::endl;
          return false;
        } // if
        break;
      case shape_sawtooth_down:      // downwards sawtooth
        if(!compute_sawtooth_down()) {
          std::cerr << "error: something went wrong while computing FF for a sawtooth down"
                << std::endl;
          return false;
        } // if
        break;
      case shape_sawtooth_up:        // upwards sawtooth
        if(!compute_sawtooth_up()) {
          std::cerr << "error: something went wrong while computing FF for a sawtooth up"
                << std::endl;
          return false;
        } // if
        break;
      case shape_prism3x:          // triangular grating in x direction
        if(!compute_prism3x(params, ff, tau, eta, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a prism3x"
                << std::endl;
          return false;
        } // if
        break;
      case shape_pyramid:          // pyramid
        if(!compute_pyramid(params, ff, tau, eta, transvec)){
          std::cerr << "error: something went wrong while computing FF for a pyramid"
                << std::endl;
          return false;
        } // if
        break;
      case shape_trunccone:        // truncated cone
        if(!compute_truncated_cone(params, tau, eta, ff, transvec)) {
          std::cerr << "error: something went wrong while computing FF for a truncated cone"
                << std::endl;
          return false;
        } // if
        break;
      default:
        std::cerr << "error: invalid shape. given shape is not supported" << std::endl;
        return false;
    } // switch

    #ifdef TIME_DETAIL_1
      compute_timer.stop();
      std::cerr << "**               FF compute time: " << compute_timer.elapsed_msec() << " ms."
            << std::endl;
    #endif // TIME_DETAIL_1

    return true;
  } // AnalyticFormFactor::compute()


  /**
   * matrix computation helpers
   */

  bool AnalyticFormFactor::mat_fq_inv_in(unsigned int x_size,
                      unsigned int y_size,
                      unsigned int z_size,
                      complex_vec_t& matrix, real_t y) {
    for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      *i = fq_inv(*i, y);
    } // for
    
    return true;
  } // AnalyticFormFactor::mat_fq_inv()


  bool AnalyticFormFactor::mat_fq_inv(unsigned int x_size, unsigned int y_size, unsigned int z_size,
                    const complex_vec_t& matrix, real_t y, complex_vec_t& result) {
    result.clear();
    for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back(fq_inv(*i, y));
    } // for
    return true;
  } // AnalyticFormFactor::mat_fq_inv()


  complex_t AnalyticFormFactor::fq_inv(complex_t value, real_t y) {
    complex_t unitc(0, 1.0);
    complex_t temp = 2.0 * exp(unitc * value * y / (real_t) 2.0) *
              sin(value * y / (real_t) 2.0) / value;
    if(std::abs(temp) <= 1e-14) temp = y;
    return temp;
  } // AnalyticFormFactor::fq_inv()


  bool AnalyticFormFactor::mat_sinc(unsigned int x_size, unsigned int y_size, unsigned int z_size,
                    const complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(std::vector<complex_t>::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back(sinc(*i));
    } // for
    return true;
  } // AnalyticFormFactor::mat_sinc()


  bool AnalyticFormFactor::mat_sinc_in(unsigned int x_size,
                      unsigned int y_size,
                      unsigned int z_size,
                      std::vector<complex_t>& matrix) {
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      *i = sinc(*i);
    } // for
    return true;
  } // AnalyticFormFactor::mat_sinc()


  /**
   * Other helpers
   */

  complex_t AnalyticFormFactor::sinc(complex_t value) {
    complex_t temp;
    if(fabs(value.real()) <= 1e-14 && fabs(value.imag()) <= 1e-14) temp = complex_t(1.0, 0.0);
    else temp = sin(value) / value;
    return temp;
  } // AnalyticFormFactor::sinc()


  /**
   * generate parameter distribution
   */
  bool AnalyticFormFactor::param_distribution(ShapeParam& param, std::vector<real_t>& dim,
                        std::vector<real_t>& dim_vals) {
    if(!param.isvalid()) {
      std::cerr << "error: invalid shape parameter encountered" << std::endl;
      return false;
    } // if

    if(param.nvalues() < 1) {
      std::cerr << "error: empty parameter found (nvalues = 0)" << std::endl;
      return false;
    } // if
    real_t pmax = param.max(), pmin = param.min();
    if(pmax < pmin) pmax = pmin;
    if(param.nvalues() > 1) {
      real_t step = fabs(pmax - pmin) / (param.nvalues() - 1);
      real_t curr = pmin;
      do {
        dim.push_back(curr);
        curr += step;
      } while(curr < param.max());  // assumes min < max ...
    } else {
      dim.push_back(pmin);
    } // if-else

    if(param.stat() == stat_none || param.stat() == stat_null) {  // just one value
      dim_vals.push_back(1.0);
    } else if(param.stat() == stat_uniform) {
      for(unsigned int i = 0; i < dim.size(); ++ i) {
        dim_vals.push_back(1.0);
      } // for
    } else if(param.stat() == stat_gaussian) {
      real_t mean = param.mean();
      if(!boost::math::isfinite(mean)) {
        mean = (dim[0] + dim[dim.size() - 1]) / 2;
      } // if
      for(unsigned int i = 0; i < dim.size(); ++ i) {
        dim_vals.push_back(exp(-1.0 * pow((dim[i] - mean), 2) / (2 * pow(param.deviation(), 2)))
                  / (sqrt(2 * PI_) * param.deviation()));
      } // for
    } else if(param.stat() == stat_random) {
      std::cerr << "uh-oh: random statistic has not been implemented yet" << std::endl;
      return false;
      // ...
    } else {
      std::cerr << "error: an invalid statistic value given for shape parameter" << std::endl;
      return false;
    } // if-else

    return true;
  } // AnalyticFormFactor::param_distribution()

} // namespace hig
