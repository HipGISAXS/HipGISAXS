/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: layer.cpp
 *  Created: Jun 13, 2012
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

#include <iostream>
#include <iomanip>

#include <model/layer.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>


namespace hig {

  Layer::Layer() { }
  Layer::~Layer() { }

  void Layer::init() {
    clear();
  } // Layer::init()

  void Layer::clear() {
    key_.clear();
    order_ = 0;
    thickness_ = 0.0;
    refindex_.delta(0.0); refindex_.beta(0.0);
  } // Layer::clear()

  bool Layer::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
      std::string keyword2, rem_str2;
    switch(TokenMapper::instance().get_keyword_token(keyword)) {
      case key_token:
      case layer_order_token:
        std::cerr << "warning: immutable param in '" << str << "'. ignoring." << std::endl;
        break;

      case layer_thickness_token:
        thickness_ = new_val;
        break;

      case refindex_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case refindex_beta_token:
            refindex_.beta(new_val);
            break;

          case refindex_delta_token:
            refindex_.delta(new_val);
              break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << str << "'" << std::endl;
            return false;

          default:
            std::cerr << "error: misplaced keyword in '" << str << "'" << std::endl;
            return false;
        } // switch
        break;

      case error_token:
        std::cerr << "error: invalid keyword in '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword in '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // Layer::update_param()

  void Layer::print() {
    std::cout << " key_ = " << key_ << std::endl
          << " order_ = " << order_ << std::endl
          << " thickness_ = " << thickness_ << std::endl
          << " refindex_ = [" << refindex_.delta() << ", "
          << refindex_.beta() << "]" << std::endl << std::endl;
  } // Layer::print()

} // namespace hig

