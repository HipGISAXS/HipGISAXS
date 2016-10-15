/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: inst_scattering.cpp
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

#include <model/inst_scattering.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>


namespace hig {

  ScatteringParams::ScatteringParams() { }
  ScatteringParams::~ScatteringParams() { }

  void ScatteringParams::init() {
    expt_ = "gisaxs";
    alpha_i_.min_ = alpha_i_.max_ = alpha_i_.step_ = 0.1;      // check
    inplane_rot_.min_ = inplane_rot_.max_ = inplane_rot_.step_ = 0;  // check
    tilt_.min_ = tilt_.max_ = tilt_.step_ = 0;            // check
    photon_.value_ = 10000; photon_.unit_ = "ev";
    polarization_ = "s";
    coherence_ = 300;
    spot_area_ = 0.01; //0.001;
    //smearing_[0] = smearing_[1] = smearing_[2] = 1;
    smearing_ = 1.0;
  } // ScatteringParams::init()


  bool ScatteringParams::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    std::string keyword2, rem_str2;
    switch(TokenMapper::instance().get_keyword_token(keyword)) {
      case instrument_scatter_expt_token:
      case instrument_scatter_polarize_token:
      case instrument_scatter_smearing_token:
        std::cerr << "warning: immutable param in '" << str << "'. ignoring." << std::endl;
        break;

      case instrument_scatter_alphai_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case min_token:
            alphai_min(new_val);
            break;

          case max_token:
            alphai_max(new_val);
            break;

          case step_token:
            alphai_step(new_val);
            break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << rem_str << "'" << std::endl;
            return false;

          default:
            std::cerr << "error: misplaced keyword in '" << rem_str << "'" << std::endl;
            return false;
        } // switch
        break;

      case instrument_scatter_inplanerot_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case min_token:
            inplane_rot_min(new_val);
            break;

          case max_token:
            inplane_rot_max(new_val);
            break;

          case step_token:
            inplane_rot_step(new_val);
            break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << rem_str << "'" << std::endl;
            return false;

          default:
            std::cerr << "error: misplaced keyword in '" << rem_str << "'" << std::endl;
            return false;
        } // switch
        break;

      case instrument_scatter_tilt_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case min_token:
            tilt_min(new_val);
            break;

          case max_token:
            tilt_max(new_val);
            break;

          case step_token:
            tilt_step(new_val);
            break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << rem_str << "'" << std::endl;
            return false;

          default:
            std::cerr << "error: misplaced keyword in '" << rem_str << "'" << std::endl;
            return false;
        } // switch
        break;

      case instrument_scatter_photon_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case instrument_scatter_photon_value_token:
            photon_value(new_val);
            break;

          case instrument_scatter_photon_unit_token:
            std::cerr << "warning: immutable param in '" << rem_str
                  << "'. ignoring." << std::endl;
            break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << rem_str << "'" << std::endl;
            return false;

          default:
            std::cerr << "error: misplaced keyword in '" << rem_str << "'" << std::endl;
            return false;
        } // switch
        break;

      case instrument_scatter_coherence_token:
        coherence(new_val);
        break;

      case instrument_scatter_spotarea_token:
        spot_area(new_val);
        break;

      case error_token:
        std::cerr << "error: invalid keyword in '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword in '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // ScatteringParams::update_param()


} // namespace hig

