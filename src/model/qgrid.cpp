/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: qgrid.cpp
 *  Created: Jun 17, 2012
 *  Modified: Wed 08 Oct 2014 12:17:46 PM PDT
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

#include <fstream>
#include <common/constants.hpp>
#include <model/qgrid.hpp>
#include <config/hig_input.hpp>
#include <utils/utilities.hpp>

namespace hig {

  /**
   * create Q-grid in reciprocal space
   */
  bool QGrid::create(real_t freq, real_t alpha_i, real_t k0, int mpi_rank) {
                    // x and y are reversed in slim's code ? ...
    vector2_t total_pixels = HiGInput::instance().detector_total_pixels();
    vector2_t min_point = HiGInput::instance().param_output_minpoint();
    vector2_t max_point = HiGInput::instance().param_output_maxpoint();
    OutputRegionType type = HiGInput::instance().param_output_type();

    vector2_t beam = HiGInput::instance().detector_direct_beam();
    real_t pixel_size = HiGInput::instance().detector_pixel_size();
    real_t sd_distance = HiGInput::instance().detector_sd_distance();
    std::vector<int> pixels = HiGInput::instance().param_resolution();
    nrow_ = pixels[1];
    ncol_ = pixels[0];
    qmin_ = HiGInput::instance().param_output_minpoint();
    qmax_ = HiGInput::instance().param_output_maxpoint();

    /* one-pixel range in k-space
     * in order to resolve a pixel, dq must be <= qpixel and lower to resolve subpixels */
    vector3_t q0 = pixel_to_kspace(vector2_t(0, 0), k0, alpha_i, pixel_size, sd_distance, beam);
    vector3_t q1 = pixel_to_kspace(vector2_t(1, 1), k0, alpha_i, pixel_size, sd_distance, beam);

    vector3_t qmax, qmin, step;

    qmin[0] = 0.0; qmax[0] = 0.0;

    if(type == region_pixels) {
      vector3_t temp_qmin = pixel_to_kspace(min_point, k0, alpha_i, pixel_size, sd_distance, beam);
      vector3_t temp_qmax = pixel_to_kspace(max_point, k0, alpha_i, pixel_size, sd_distance, beam);
      qmax[1] = max(fabs(temp_qmax[1]), fabs(temp_qmin[1]));
      qmin[1] = min(fabs(temp_qmax[1]), fabs(temp_qmin[1]));
      //qmin[1] = - qmax[1];
      qmax[2] = max(temp_qmax[2], temp_qmin[2]);
      qmin[2] = min(temp_qmax[2], temp_qmin[2]);
      step[0] = fabs(qmax[0] - qmin[0]) / 2;      // why ... ?
      step[1] = fabs(q1[1] - q0[1]) / pixels[0];
      step[2] = fabs(q1[2] - q0[2]) / pixels[1];
      if(HiGInput::instance().experiment() == "saxs") {
        qmin[0] = qmin[1] = qmin[2] = 0;
        step[0] = 0; step[1] = 1;
        qmax[0] = qmax[1] = 0;
      } // if
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
      std::cout << "**                  Q-grid range: ("
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
    for (int i = 0; i < imsize; i++){
      real_t kzf_0 = qz_[i] + kzi_0;
      kzf.push_back(sgn(kzf_0)*std::sqrt(kzf_0 * kzf_0 - k0 * k0 * dnl_q));
    }
  
    // calculate 4 components
    for (int i = 0; i < imsize; i++) qz_extended_.push_back( kzf[i] - kzi);
    for (int i = 0; i < imsize; i++) qz_extended_.push_back(-kzf[i] - kzi);
    for (int i = 0; i < imsize; i++) qz_extended_.push_back( kzf[i] + kzi);
    for (int i = 0; i < imsize; i++) qz_extended_.push_back(-kzf[i] + kzi);

    return true;
  } // QGrid::create_qz_extended()


  /**
   * for fitting, update the qgrid to match reference data
   */
  bool QGrid::update(unsigned int nqy, unsigned int nqz,
            real_t qminy, real_t qminz, real_t qmaxy, real_t qmaxz,
            real_t freq, real_t alpha_i, real_t k0, int mpi_rank) {

    vector3_t qmax, qmin, step;
    qmin[0] = 0.0; qmax[0] = 0.0;  // x dimension has just 0.0
    qmin[1] = qminy; qmax[1] = qmaxy;
    qmin[2] = qminz; qmax[2] = qmaxz;

    step[1] = (qmax[1] - qmin[1]) / nqy;
    //step[2] = (qmax[2] - qmin[2]) / (nqz - 1);
    step[2] = (qmax[2] - qmin[2]) / nqz;

    qx_.clear(); qy_.clear(); qz_.clear();
    qx_.push_back(qmin[0]);
    //for(real_t val = qmin[1]; val <= qmax[1]; val += step[1]) qy_.push_back(val);
    //for(real_t val = qmin[2]; val <= qmax[2]; val += step[2]) qz_.push_back(val);
    real_t val = qmin[1]; for(int i = 0; i < nqy; ++ i, val += step[1]) qy_.push_back(val);
    val = qmin[2]; for(int i = 0; i < nqz; ++ i, val += step[2]) qz_.push_back(val);
    std::reverse(qz_.begin(), qz_.end());

    // sanity check
    if(qy_.size() != nqy || qz_.size() != nqz) {
      std::cerr << "error: mismatch in the needed qgrid size and the constructed one" << std::endl;
      return false;
    } // if

    if(mpi_rank == 0) {
      std::cout << "**              New Q-grid range: ("
            << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << ") x ("
            << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << ")" << std::endl;
      std::cout << "**               NQX x NQY x NQZ: " << qx_.size() << " x " << qy_.size()
            << " x " << qz_.size() << std::endl;
    } // if
    return true;
  } // QGrid::create()


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
