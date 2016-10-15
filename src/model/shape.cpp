/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: shape.cpp
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

#include <iostream>

#include <model/shape.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>


namespace hig {

  /* ShapeParam */


  void ShapeParam::init() {
    type_ = param_error;  // this object is required from user
    stat_ = stat_none;
    type_name_.clear();
    max_ = min_ = p1_ = p2_ = 0.0;
    nvalues_ = 1;      // default nvalue
  } // ShapeParam::init()

  void ShapeParam::clear() {
    type_name_.clear();
    type_ = param_null;
    stat_ = stat_null;
    max_ = min_ = p1_ = p2_ = 0.0;
    nvalues_ = 1;
  } // ShapeParam::clear()

  void ShapeParam::print() {
    std::cout << "  type_name_ = " << type_name_ << std::endl
          << "  type_ = " << type_ << std::endl
          << "  stat_ = " << stat_ << std::endl
          << "  max_ = " << max_ << std::endl
          << "  min_ = " << min_ << std::endl
          << "  p1_ = " << p1_ << std::endl
          << "  p2_ = " << p2_ << std::endl
          << "  nvalues_ = " << nvalues_ << std::endl
          << "  isvalid_ = " << isvalid_ << std::endl
          << std::endl;
  } // ShapeParam::print()


  /* Shape */

  Shape::Shape() { init(); }
  Shape::~Shape() { }

  // not used
  Shape::Shape(const std::string& key, const ShapeName name, const vector3_t& origin,
          const real_t zrot, const real_t yrot, const real_t xrot, shape_param_list_t& param_list) :
          key_(key), name_(name), originvec_(origin),
          zrot_(zrot), yrot_(yrot), xrot_(xrot) {
    for(shape_param_iterator_t i = param_list.begin(); i != param_list.end(); i ++) {
      parse_param((*i).second);
      insert_param((*i));
    } // for
  } // Shape::Shape()

  // not used
  Shape::Shape(const std::string& key, const ShapeName name) : key_(key), name_(name) {
    originvec_[0] = originvec_[1] = originvec_[2] = 0.0;
    zrot_ = yrot_ = xrot_ = 0.;
  } // Shape::Shape()


  void Shape::init() {  // the user needs to provide the shape, no defaults
    key_.clear();
    params_.clear();
    name_ = shape_error;
    name_str_.clear();
    originvec_[0] = originvec_[1] = originvec_[2] = 0.0;
    zrot_ = yrot_ = xrot_ = 0.;
    
  } // Shape::init()


  void Shape::clear() {
    key_.clear();
    params_.clear();
    name_ = shape_null;
    name_str_.clear();
    originvec_[0] = originvec_[1] = originvec_[2] = 0.0;
    zrot_ = yrot_ = xrot_ = 0.;
  } // Shape::clear()


  bool Shape::parse_param(const ShapeParam& param) const {
    // do some micro error checking
    return true;
  } // Shape::parse_param()


  bool Shape::insert_param(const std::pair <std::string, ShapeParam>& param) {
    // TODO: check if it already exists ...
    params_[param.first] = param.second;
    return true;
  } // Shape::insert_param()


  bool Shape::insert_param(const std::string& type, const ShapeParam& param) {
    // TODO: check if it already exists ...
    params_[type] = param;
    return true;
  } // Shape::insert_param()


  void Shape::print() {
    std::cout << " key_ = " << key_ << std::endl
          << " name_ = " << name_ << std::endl
          << " name_str_ = " << name_str_ << std::endl
          << " originvec_ = [" << originvec_[0] << ", " << originvec_[1]
          << ", " << originvec_[2] << "]" << std::endl
          << " zrot_ = " << zrot_ << std::endl
          << " yrot_ = " << yrot_ << std::endl
          << " xrot_ = " << xrot_ << std::endl
          << " params_: " << params_.size() << std::endl;
    for(shape_param_iterator_t i = params_.begin(); i != params_.end(); i ++) {
      (*i).second.print();
    } // for
    std::cout << std::endl;
  } // Shape::print()


  /**
   * updates for fitting
   */


  bool ShapeParam::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    std::string keyword_name, key;
    if(!extract_keyword_name_and_key(keyword, keyword_name, key)) return false;
    switch(TokenMapper::instance().get_keyword_token(keyword_name)) {
      case min_token:
        min_ = new_val;
        break;

      case max_token:
        max_ = new_val;
        break;

      case shape_param_p1_token:
        p1_ = new_val;
        break;

      case shape_param_p2_token:
        p2_ = new_val;
        break;

      case shape_param_nvalues_token:
        nvalues_ = int(new_val);
        break;

      case stat_token:
        std::cerr << "warning: immutable param in '" << str << "'. ignoring." << std::endl;
        break;

      case error_token:
        std::cerr << "error: invalid keyword '" << keyword_name
              << "' in param '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword '" << keyword_name
              << "' in param '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // ShapeParam::update_param()


  bool Shape::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str, rem_str2;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    std::string keyword_name, key;
    if(!extract_keyword_name_and_key(keyword, keyword_name, key)) return false;
    std::string coord;
    switch(TokenMapper::instance().get_keyword_token(keyword_name)) {
      case key_token:
      case shape_name_token:
      case shape_originvec_token:
        //std::cerr << "warning: immutable param in '" << str << "'. ignoring." << std::endl;
        if(!extract_first_keyword(rem_str, coord, rem_str2)) return false;
        if(coord.compare("x") == 0 || coord.compare("0") == 0) {
          originvec_[0] = new_val;
        } else if(coord.compare("y") == 0 || coord.compare("1") == 0) {
          originvec_[1] = new_val;
        } else if(coord.compare("z") == 0 || coord.compare("2") == 0) {
          originvec_[2] = new_val;
        } else {
          std::cerr << "error: invalid keyword '" << coord << "' in param '"
                << str << "'" << std::endl;
          return false;
        } // if-else
        break;

      case shape_zrot_token:
        zrot_ = new_val;
        break;

      case shape_yrot_token:
        yrot_ = new_val;
        break;

      case shape_xrot_token:
        xrot_ = new_val;

      case shape_param_token:
        #ifdef __INTEL_COMPILER
          if(params_.count(key) == 0 ||
              !params_[key].update_param(rem_str, new_val)) return false;
        #else
          if(!params_.at(key).update_param(rem_str, new_val)) return false;
        #endif
        break;

      case error_token:
        std::cerr << "error: invalid keyword '" << keyword_name
              << "' in param '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword '" << keyword_name
              << "' in param '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // Shape::update_param()


} // namespace hig
