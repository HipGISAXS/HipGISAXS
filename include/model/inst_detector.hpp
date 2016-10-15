/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: inst_detector.hpp
 *  Created: Jun 05, 2012
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

#ifndef __INST_DETECTOR_HPP__
#define __INST_DETECTOR_HPP__

#include <iostream>
#include <string>

#include <common/globals.hpp>

namespace hig {

  class DetectorParams {
    friend class HiGInput;

    private:
      std::string origin_;   // change to enum - tl tr bl br
      real_t pixel_size_;
      real_t sd_distance_;    /* sample to detector distance */
      vector2_t total_pixels_;  // make it vector3_t for 3d ? ...
      vector2_t direct_beam_;

    public:
      DetectorParams();
      ~DetectorParams();

      void init();
      void clear();

      void origin(std::string s) { origin_ = s; }

      void direct_beam(vector2_t v) { direct_beam_ = v; }
      void total_pixels(vector2_t v) { total_pixels_ = v; }
      void direct_beam(real_t v, real_t w) { direct_beam_[0] = v; direct_beam_[1] = w; }
      void total_pixels(real_t v, real_t w) { total_pixels_[0] = v; total_pixels_[1] = w; }

      void sd_distance(real_t d) { sd_distance_ = d; }
      void pixel_size(real_t d) { pixel_size_ = d; }

      bool update_param(const std::string&, real_t);

      void print() {
        std::cout << " origin_ = " << origin_ << std::endl
              << " pixel_size_ = " << pixel_size_ << std::endl
              << " sd_distance_ = " << sd_distance_ << std::endl
              << " total_pixels_ = [" << total_pixels_[0] << ", "
              << total_pixels_[1] << "]" << std::endl
              << " direct_beam_ = [" << direct_beam_[0] << ", "
              << direct_beam_[1] << "]" << std::endl
              << std::endl;
      } // print()

  }; // class DetectorParams

} // namespace hig

#endif // __INST_DETECTOR_HPP__
