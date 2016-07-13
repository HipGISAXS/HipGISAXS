/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: qgrid.cpp
 *  Created: Jun 17, 2012
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

#include <fstream>
#include <algorithm>

#include <common/constants.hpp>
#include <model/qgrid.hpp>
#include <config/hig_input.hpp>
#include <utils/utilities.hpp>

namespace hig {

  /**
   * create Q-grid in reciprocal space
   */
  bool QGrid::create(const ComputeParams & params, real_t alpha_i, real_t k0, int mpi_rank) {
    vector2_t min_point = params.output_minpoint();
    vector2_t max_point = params.output_maxpoint();
    OutputRegionType type = params.output_region_type();
    std::vector<int> pixels = params.resolution();

    nrow_ = pixels[1];
    ncol_ = pixels[0];
    vector3_t qmax, qmin, step;
    qmin[0] = 0.0; qmax[0] = 0.0;

    if(type == region_pixels) {
      std::cerr << "error: output data in pixel-space is not implemented yet.\nPlease use qspace." << std::endl;
      return false;
    } else if(type == region_qspace) {

      /* calculate the angles */
      // calculate angles in vertical direction
      real_t dq = (max_point[1] - min_point[1])/(pixels[1]-1);
      alpha_.resize(pixels[1]);
      for (int i = 0; i < pixels[1]; i++){
        real_t qpt = min_point[1] + i * dq;
        alpha_[i] = std::asin(qpt/k0 - std::sin(alpha_i));
      }

      // calculate angles in horizontal direction
      dq = (max_point[0] - min_point[0])/(pixels[0]-1);
      real_t cos_ai = std::cos(alpha_i);
      real_t cos_af = std::cos(alpha_[0]);
      real_vec_t theta; theta.resize(pixels[0]);
      for (int i = 0; i < pixels[0]; i++){
        real_t qpt = min_point[0] + i * dq;
        real_t kf2    = std::pow(qpt / k0, 2);
        real_t tmp = (cos_af * cos_af + cos_ai * cos_ai - kf2)/(2 * cos_af * cos_ai);
        if (tmp > 1.0) tmp = 1.0;
        theta[i] = sgn(qpt) * std::acos(tmp);
      }

      /* calculate q-vectors on the Ewald sphere for every pixel
       * on the detector.
       */
      std::reverse(alpha_.begin(), alpha_.end());
      for (int i = 0; i < pixels[1]; i++) {
        for (int j = 0; j < pixels[0]; j++) {
          real_t tth = theta[j];
          real_t alf = alpha_[i];
          qx_.push_back (k0 * (std::cos(alf) * std::cos(tth) - std::cos(alpha_i)));
          qy_.push_back (k0 * (std::cos(alf) * std::sin(tth)));
          qz_.push_back (k0 * (std::sin(alf) + std::sin(alpha_i)));
        }
      }
    } else {
      std::cerr << "error: unknown output region type" << std::endl;
      return false;
    } // if-else

    if(mpi_rank == 0) {
      std::cerr << "**                  Q-grid range: ("
            << qx_[0] << ", " << qy_[0] << ", " << qz_.back()
            << ") x ("
            << qx_.back() << ", " << qy_.back() << ", " << qz_[0]
            << ")" << std::endl;
    } // if

    return true;
  } // QGrid::create()


  bool QGrid::create_qz_extended(real_t k0, real_t alpha_i, complex_t dnl_q) {

    // reset qz_extended_ 
    size_t imsize = nrow_ * ncol_;
    qz_extended_.clear();
    qz_extended_.reserve(4 * imsize);

    // incoming vectors
    real_t sin_ai = std::sin(alpha_i);
    real_t kzi_0 = -1 * k0 * sin_ai;
    complex_t kzi = -1 * k0 * std::sqrt(sin_ai * sin_ai - dnl_q);

    std::vector<complex_t> kzf;
    for(int i = 0; i < imsize; ++ i){
      real_t kzf_0 = qz_[i] + kzi_0;
      kzf.push_back(sgn(kzf_0) * std::sqrt(kzf_0 * kzf_0 - k0 * k0 * dnl_q));
    } // for
  
    // calculate 4 components
    for(int i = 0; i < imsize; ++ i) qz_extended_.push_back( kzf[i] - kzi);
    for(int i = 0; i < imsize; ++ i) qz_extended_.push_back(-kzf[i] - kzi);
    for(int i = 0; i < imsize; ++ i) qz_extended_.push_back( kzf[i] + kzi);
    for(int i = 0; i < imsize; ++ i) qz_extended_.push_back(-kzf[i] + kzi);

    return true;
  } // QGrid::create_qz_extended()


  /**
   * for fitting, update the qgrid to match reference data
   */
  bool QGrid::update(unsigned int nqy, unsigned int nqz,
                     real_t qminy, real_t qminz, real_t qmaxy, real_t qmaxz,
                     real_t freq, real_t alpha_i, real_t k0, int mpi_rank) {

    nrow_ = nqz; ncol_ = nqy;

    // calculate angles in vertical direction
    real_t dq = (qmaxz - qminz) / (nqz - 1);
    alpha_.clear();
    alpha_.resize(nqz);
    for(int i = 0; i < nqz; ++ i) {
      real_t qpt = qminz + i * dq;
      alpha_[i] = std::asin(qpt / k0 - std::sin(alpha_i));
    } // for

    // calculate angles in horizontal direction
    dq = (qmaxy - qminy) / (nqy - 1);
    real_t cos_ai = std::cos(alpha_i);
    real_t cos_af = std::cos(alpha_[0]);
    real_vec_t theta; theta.resize(nqy);
    for(int i = 0; i < nqy; ++ i) {
      real_t qpt = qminy + i * dq;
      real_t kf2 = std::pow(qpt / k0, 2);
      real_t tmp = (cos_af * cos_af + cos_ai * cos_ai - kf2) / (2 * cos_af * cos_ai);
      tmp = min(tmp, 1.0);
      theta[i] = sgn(qpt) * std::acos(tmp);
    } // for

    // calculate q-vectors on the Ewald sphere for every pixel on the detector
    qx_.clear(); qy_.clear(); qz_.clear();
    std::reverse(alpha_.begin(), alpha_.end());
    for (int i = 0; i < nqz; i++) {
      for (int j = 0; j < nqy; j++) {
        real_t tth = theta[j];
        real_t alf = alpha_[i];
        qx_.push_back(k0 * (std::cos(alf) * std::cos(tth) - std::cos(alpha_i)));
        qy_.push_back(k0 * (std::cos(alf) * std::sin(tth)));
        qz_.push_back(k0 * (std::sin(alf) + std::sin(alpha_i)));
      } // for
    } // for

    // sanity check
    if(qy_.size() != nqy || qz_.size() != nqz) {
      std::cerr << "error: mismatch in the needed qgrid size and the constructed one" << std::endl;
      return false;
    } // if

    if(mpi_rank == 0) {
      std::cerr << "**              New Q-grid range: ("
                << qx_[0] << ", " << qy_[0] << ", " << qz_.back() << ") x ("
                << qx_.back() << ", " << qy_.back() << ", " << qz_[0] << ")" << std::endl;
    } // if
    return true;
  } // QGrid::update()


  /**
   * converts given pixel into q-space point
   * TODO: This doesn't seem right, the pixel values and sample-detector distance
   * should come either from input file or pre-build models of common detectors.
   * Also need should check if coordinates need to be corrected for alpha_i
   */
  bool QGrid::pixel_to_kspace(vector2_t pixel, real_t k0, real_t alpha_i,
          real_t rho, real_t distance, vector2_t beam,
          vector3_t& qpoint) {
    //real_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / distance) - alpha_i;
    //real_t theta_f = atan(rho * (pixel[0] - beam[0]) / distance);
    real_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / 4350) - alpha_i;
    real_t theta_f = atan(rho * (pixel[0] - beam[0]) / 4350);

    qpoint[0] = k0 * (cos(theta_f) * cos(alpha_f) - cos(alpha_i));  // x
    qpoint[1] = k0 * cos(alpha_f) * sin(theta_f);          // y
    qpoint[2] = k0 * (sin(alpha_f) + sin(alpha_i));          // z

    return true;
  } // QGrid::pixel_to_kspace()


  void QGrid::save (const char* filename) {
    std::ofstream f(filename);
    for (int i = 0; i < qx_.size(); i++)
      f << qx_[i] << ", " << qy_[i] << ", " << qz_[i] << std::endl;
    f.close();
  }
  vector3_t QGrid::pixel_to_kspace(vector2_t pixel, real_t k0, real_t alpha_i,
          real_t rho, real_t distance, vector2_t beam) {
    // 'distance' is not used. currently it is hard-coded as 4350
    real_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / 4350) - alpha_i;
    real_t theta_f = atan(rho * (pixel[0] - beam[0]) / 4350);

    vector3_t qpoint;
    qpoint[0] = k0 * (cos(theta_f) * cos(alpha_f) - cos(alpha_i));  // x
    qpoint[1] = k0 * cos(alpha_f) * sin(theta_f);          // y
    qpoint[2] = k0 * (sin(alpha_f) + sin(alpha_i));          // z
    return qpoint;
  } // QGrid::pixel_to_kspace()

  bool QGrid::kspace_to_pixel() {
    // needed ... ?
    // do it later ...
    return true;
  } // QGrid::kspace_to_pixel()

} // namespace hig
