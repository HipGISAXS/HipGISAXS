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
  bool QGrid::create(float_t freq, float_t alpha_i, float_t k0, int mpi_rank) {
                    // x and y are reversed in slim's code ? ...
    vector2_t total_pixels = HiGInput::instance().detector_total_pixels();
    vector2_t min_point = HiGInput::instance().param_output_minpoint();
    vector2_t max_point = HiGInput::instance().param_output_maxpoint();
    OutputRegionType type = HiGInput::instance().param_output_type();

    vector2_t beam = HiGInput::instance().detector_direct_beam();
    float_t pixel_size = HiGInput::instance().detector_pixel_size();
    float_t sd_distance = HiGInput::instance().detector_sd_distance();
    std::vector<int> pixels = HiGInput::instance().param_resolution();
    nrow_ = pixels[1];
    ncol_ = pixels[0];

    /* one-pixel range in k-space
     * in order to resolve a pixel, dq must be <= qpixel and lower to resolve subpixels */
    vector3_t q0 = pixel_to_kspace(vector2_t(0, 0), k0, alpha_i, pixel_size, sd_distance, beam);
    vector3_t q1 = pixel_to_kspace(vector2_t(1, 1), k0, alpha_i, pixel_size, sd_distance, beam);

    vector3_t qmax, qmin, step;
    float_t theta[2], alpha[2];

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
      theta[0] = acos (1. - 0.5 * pow(min_point[0] / k0, 2));
      if ( min_point[0] < 0 ) theta[0] *= -1;
      theta[1] = acos (1. - 0.5 * pow(max_point[0] / k0, 2));
      if ( max_point[0] < 0 ) theta[1] *= -1;
      alpha[0] = acos (1. - 0.5 * pow(max_point[1] / k0, 2)) - alpha_i;
      alpha[1] = acos (1. - 0.5 * pow(min_point[1] / k0, 2)) - alpha_i;
      step[0] = (theta[1] - theta[0]) / pixels[0];
      step[1] = (alpha[1] - alpha[0]) / pixels[1];
    } else {
      std::cerr << "error: unknown output region type" << std::endl;
      return false;
    } // if-else


    /* calculate q-vectors on the Ewald sphere for every pixel
     * on the detector.
     */
    for (int i = 0; i < pixels[1]; i++) {
        for (int j = 0; j < pixels[0]; j++) {
            float_t tth = theta[0] + j * step[0];
            float_t alf = alpha[0] + i * step[1];
            qx_.push_back (k0 * (cos(alf) * cos(tth) - cos(alpha_i)));
            qy_.push_back (k0 * (cos(alf) * sin(tth)));
            qz_.push_back (k0 * (sin(alf) + sin(alpha_i)));
        }
    }

    if(mpi_rank == 0) {
      std::cout << "**                  Q-grid range: ("
            << qx_[0] << ", " << qy_[0] << ", " << qz_.back()
            << ") x ("
            << qx_.back() << ", " << qy_.back() << ", " << qz_[0]
            << ")" << std::endl;
    } // if

    return true;
  } // QGrid::create()


  bool QGrid::create_qz_extended(float_t k0, float_t kzi_0, complex_t kzi_q, complex_t dnl_q) {
    cqvec_t qz_temp0, qz_temp1, qz_temp2, qz_temp3;
    qz_extended_.clear(); qz_temp0.clear(); qz_temp1.clear(); qz_temp2.clear(); qz_temp3.clear();
    for(qvec_iter_t q = qz_.begin(); q != qz_.end(); ++ q) {
      float_t temp0 = (*q) + kzi_0;
      float_t temp1 = temp0 * temp0;
      complex_t temp2 = k0 * k0 * dnl_q;
      complex_t temp3 = sqrt(temp1 - temp2);
      complex_t temp4 = temp3 - kzi_q;
      qz_temp0.push_back(temp4);
      complex_t temp5 = - temp3 - kzi_q;
      qz_temp1.push_back(temp5);
      complex_t temp6 = temp3 + kzi_q;
      qz_temp2.push_back(temp6);
      complex_t temp7 = - temp3 + kzi_q;
      qz_temp3.push_back(temp7);
    } // for
    // the 4 vectors be concatenated instead of copying ...
    for(cqvec_iter_t q = qz_temp0.begin(); q != qz_temp0.end(); ++ q) qz_extended_.push_back(*q);
    for(cqvec_iter_t q = qz_temp1.begin(); q != qz_temp1.end(); ++ q) qz_extended_.push_back(*q);
    for(cqvec_iter_t q = qz_temp2.begin(); q != qz_temp2.end(); ++ q) qz_extended_.push_back(*q);
    for(cqvec_iter_t q = qz_temp3.begin(); q != qz_temp3.end(); ++ q) qz_extended_.push_back(*q);

    return true;
  } // QGrid::create_qz_extended()


  /**
   * for fitting, update the qgrid to match reference data
   */
  bool QGrid::update(unsigned int nqy, unsigned int nqz,
            float_t qminy, float_t qminz, float_t qmaxy, float_t qmaxz,
            float_t freq, float_t alpha_i, float_t k0, int mpi_rank) {

    vector3_t qmax, qmin, step;
    qmin[0] = 0.0; qmax[0] = 0.0;  // x dimension has just 0.0
    qmin[1] = qminy; qmax[1] = qmaxy;
    qmin[2] = qminz; qmax[2] = qmaxz;

    step[1] = (qmax[1] - qmin[1]) / nqy;
    //step[2] = (qmax[2] - qmin[2]) / (nqz - 1);
    step[2] = (qmax[2] - qmin[2]) / nqz;

    qx_.clear(); qy_.clear(); qz_.clear();
    qx_.push_back(qmin[0]);
    //for(float_t val = qmin[1]; val <= qmax[1]; val += step[1]) qy_.push_back(val);
    //for(float_t val = qmin[2]; val <= qmax[2]; val += step[2]) qz_.push_back(val);
    float_t val = qmin[1]; for(int i = 0; i < nqy; ++ i, val += step[1]) qy_.push_back(val);
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
  bool QGrid::pixel_to_kspace(vector2_t pixel, float_t k0, float_t alpha_i,
          float_t rho, float_t distance, vector2_t beam,
          vector3_t& qpoint) {
    //float_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / distance) - alpha_i;
    //float_t theta_f = atan(rho * (pixel[0] - beam[0]) / distance);
    float_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / 4350) - alpha_i;
    float_t theta_f = atan(rho * (pixel[0] - beam[0]) / 4350);

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
  vector3_t QGrid::pixel_to_kspace(vector2_t pixel, float_t k0, float_t alpha_i,
          float_t rho, float_t distance, vector2_t beam) {
    // 'distance' is not used. currently it is hard-coded as 4350
    float_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / 4350) - alpha_i;
    float_t theta_f = atan(rho * (pixel[0] - beam[0]) / 4350);

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
