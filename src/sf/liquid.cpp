/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: liquid.cpp
 *
 *  Author: Dinesh Kumar <dkumar@lbl.gov>
 *          Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>
#include <sf/sf.hpp>
#include <model/qgrid.hpp> 
#include <utils/utilities.hpp>
#include <common/constants.hpp>

namespace hig {

  bool StructureFactor::paracrystal1d(real_t mean_dist, real_t stddev_dist, 
          real_t domain_size){
    
    real_t exp_v1 = std::exp(-2. * mean_dist / domain_size);
#pragma omp parallel for
    for (int i = 0; i < nz_; i++){
      int j = i % ny_;
      real_t qy = QGrid::instance().qy(j);
      real_t qysy = stddev_dist * qy;
      real_t exp_v2 = std::exp(-1.0 * qysy * qysy);
      real_t exp_v3 = std::exp(-0.5 * qysy * qysy);
      real_t tmp  = (1. - exp_v2 * exp_v1) / (1. - 2. * exp_v3 * exp_v1 * std::cos(qy * mean_dist)
              + exp_v2 * exp_v1);
      sf_[i] = complex_t(std::sqrt(tmp), 0.);
    }
    return true;
  } // paracrystal1d()

  bool StructureFactor::paracrystal2d_cubic(real_t mean_dist_x, real_t mean_dist_y, 
          real_t stddev_dist_x, real_t stddev_dist_y, real_t domain_size){
    real_t exp_v1x = std::exp(-2. * mean_dist_x / domain_size);
    real_t exp_v1y = std::exp(-2. * mean_dist_y / domain_size);

#pragma omp parallel for
    for (int i = 0; i < nz_; i++){
      int j = i % ny_;
      real_t qx = QGrid::instance().qx(j);
      real_t qy = QGrid::instance().qy(j);
      real_t qpar = std::sqrt(qx * qx + qy * qy);
      real_t cos_vx = std::cos(qpar * mean_dist_x);
      real_t cos_vy = std::cos(qpar * mean_dist_y);

      real_t qxsx = qpar * stddev_dist_x;
      real_t exp_v2x = std::exp(-1. * qxsx * qxsx);
      real_t exp_v3x = std::exp(-0.5 * qxsx * qxsx);

      real_t qysy = qpar * stddev_dist_y;
      real_t exp_v2y = std::exp(-1. * qysy * qysy);
      real_t exp_v3y = std::exp(-0.5 * qysy * qysy);

      real_t s1 = (1. - exp_v2x * exp_v1x) / 
          (1. -2. * exp_v3x * exp_v1x * cos_vx + exp_v2x * exp_v1x);
      real_t s2 = (1. - exp_v2y * exp_v1y) / 
          (1. -2. * exp_v3y * exp_v1y * cos_vy + exp_v2y * exp_v1y);
       real_t tmp = std::abs(s1 * s2);
      sf_[i] = complex_t(std::sqrt(tmp), 0.);
    }
    return true;
  } // paracrystal2d_cubic

  bool StructureFactor::paracrystal2d_hex(real_t mean_dist_x, real_t mean_dist_y,
          real_t stddev_dist_x, real_t stddev_dist_y, real_t domain_size){

    const real_t C1 = 0.5 * std::sqrt(3.);
    real_t exp_v1x = std::exp(-2. * mean_dist_x / domain_size);
    real_t exp_v1y = std::exp(-2. * mean_dist_y / domain_size);

#pragma omp parallel for
    for (int i = 0; i < nz_; i++){
      int j = i % ny_;
      real_t qx = QGrid::instance().qx(j);
      real_t qy = QGrid::instance().qy(j);
      real_t qpar = std::sqrt(qx*qx + qy*qy);
      real_t qpsx = qpar * stddev_dist_x;
      real_t qpsy = qpar * stddev_dist_y;

      real_t phi_x = std::cos(C1 * qx * mean_dist_x + 0.5 * qy * mean_dist_y); 
      real_t phi_y = std::cos(qy * mean_dist_y);
      real_t exp_v2x = std::exp(-1. * qpsx * qpsx);
      real_t exp_v2y = std::exp(-1. * qpsy * qpsy);
      real_t exp_v3x = std::exp(-0.5 * qpsx * qpsx);
      real_t exp_v3y = std::exp(-0.5 * qpsy * qpsy);
      real_t s1 = (1. - exp_v2x * exp_v1x)/
          (1. -2. * exp_v3x * exp_v1x * phi_x + exp_v2x * exp_v1x);
      real_t s2 = (1. -exp_v2y * exp_v1y)/
          (1. -2. * exp_v3y * exp_v1y * phi_y + exp_v2y * exp_v1y);
      real_t tmp = std::abs(s1 * s2);
      sf_[i] = complex_t(std::sqrt(tmp), 0.);
    }
    return true;
  }

  real_t StructureFactor::PYFcn(real_t alpha, real_t beta, real_t gamma, real_t x){
    real_t sinx = std::sin(x);
    real_t cosx = std::cos(x);
    real_t x2   = x * x;
    real_t x3   = x2 * x;
    real_t x4   = x3 * x;
    real_t x5   = x4 * x;
    real_t t1 = alpha * (sinx - x * cosx)/ x2;
    real_t t2 = beta * (2 * x * sinx + (2 - x2) * cosx - 2)/ x3;
    real_t t3 = gamma * (- x4 * cosx + 4 * ((3 * x2 - 6) * cosx
            + (x3 - 6 * x) * sinx + 6)) / x5 ;
    return (t1 + t2 + t3);
  }

  bool StructureFactor::percus_yevik(real_t dia, real_t volf, int numd){
    real_t tmp = std::pow(1 - volf, 4);
    real_t alpha = std::pow(1 + 2 * volf, 2) / tmp;
    real_t beta  = -6 * volf * std::pow(1 + 0.5 * volf, 2) / tmp;
    real_t gamma = 0.5 * volf * alpha;

    if (numd == 2){
#pragma omp parallel for
      for (int i = 0; i < nz_; i++){
        int j = i % ny_;
        real_t qx = QGrid::instance().qx(j);
        real_t qy = QGrid::instance().qy(j);
        real_t qval = std::sqrt(qx * qx + qy * qy);
        if ( qval < 1.0E-06 )
          sf_[i] = alpha / 3. + beta / 4. + gamma;
        else {
          real_t x = qval * dia;
          real_t tmp = 1./(1 + 24 * volf * PYFcn(alpha, beta, gamma, x) / x);
          sf_[i] = complex_t(std::sqrt(tmp), 0.);
        } // if-else
      } // for
    } else if (numd == 3) {
#pragma omp parallel for
      for (int i = 0; i < nz_; i++){
        int j = i % ny_;
        real_t qx = QGrid::instance().qx(j);
        real_t qy = QGrid::instance().qy(j);
        complex_t qz = QGrid::instance().qz_extended(i);
        real_t qval = std::sqrt(qx * qx + qy * qy + std::norm(qz));
        if (qval < 1.0E-06)
          sf_[i] = alpha / 3. + beta / 4. + gamma;
        else {
          real_t x = qval * dia;
          real_t tmp = 1./std::sqrt(1 + 24 * volf * PYFcn(alpha, beta, gamma, x) / x);
          sf_[i] = complex_t(std::sqrt(tmp), 0.);
        } // if-else
      } // for
    } // if(numd ==?)
    return true;
  } // percus_yevik()
} // namespace hig

