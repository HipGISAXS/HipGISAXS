/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: compute_params.cpp
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

#include <model/compute_params.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>


namespace hig {

  ComputeParams::ComputeParams() { }

  void ComputeParams::init() {  // sets all default values
    pathprefix_ = ".";
    runname_ = "run_" + timestamp();
    method_ = "dwba";
    output_region_.type_ = region_qspace;
    output_region_.minpoint_[0] = -1;
    output_region_.minpoint_[1] = -1;
    output_region_.maxpoint_[0] = -1;
    output_region_.maxpoint_[1] = -1;
    resolution_.push_back(1); resolution_.push_back(1);
    nslices_ = 0;
    correlation_ = structcorr_null;
    palette_ = "default";
  } // ComputeParams::init()


  bool ComputeParams::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    switch(TokenMapper::instance().get_keyword_token(keyword)) {
      case compute_path_token:
      case compute_runname_token:
      case compute_method_token:
      case compute_outregion_token:
      case compute_resolution_token:
      case compute_nslices_token:
        std::cerr << "earning: immutable param in '" << str << "'. ignoring." << std::endl;
        break;

      case error_token:
        std::cerr << "error: invalid keyword in '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword in '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // ComputeParams::update_param()


  std::string ComputeParams::timestamp() {  // make this an independent utility ...
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[16];

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, 16, "%Y%m%d_%H%M%S", timeinfo);

    return std::string(buffer);
  } // ComputeParams::timestamp()

} // namespace hig

