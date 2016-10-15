/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: inst_detector.cpp
 *  Created: Jun 12, 2012
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

#include <model/inst_detector.hpp>
#include <common/default_values.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>

namespace hig {

  DetectorParams::DetectorParams() { }
  DetectorParams::~DetectorParams() { }

  void DetectorParams::init() {  // fill default values
    origin_ = DEFAULT_ORIGIN_;
    total_pixels_[0] = DEFAULT_TOTAL_PIXELS_Y_;
    total_pixels_[1] = DEFAULT_TOTAL_PIXELS_Z_;
    pixel_size_ = DEFAULT_PIXEL_SIZE_;
    sd_distance_ = DEFAULT_SDD_;
    direct_beam_[0] = DEFAULT_DIRECT_BEAM_Y_;
    direct_beam_[1] = DEFAULT_DIRECT_BEAM_Z_;  // should be center of the detector
  } // init()


  bool DetectorParams::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    switch(TokenMapper::instance().get_keyword_token(keyword)) {
      case instrument_detector_origin_token:
      case instrument_detector_totpix_token:
      case instrument_detector_dirbeam_token:
        std::cerr << "earning: immutable param in '" << str << "'. ignoring." << std::endl;
        break;

      case instrument_detector_pixsize_token:
        pixel_size(new_val);
        break;

      case instrument_detector_sdd_token:
        sd_distance(new_val);
        break;

      case error_token:
        std::cerr << "error: invalid keyword in '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword in '" << str << "'" << std::endl;
        return false;
    } // switch
    return true;
  } // DetectorParams::update_param()

} // namespace hig

