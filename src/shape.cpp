/***
  *  $Id: shape.cpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project:
  *
  *  File: shape.cpp
  *  Created: Jun 05, 2012
  *  Modified: Thu 30 Aug 2012 03:06:46 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>

#include "shape.hpp"

namespace hig {

	/* ShapeParam */

	void ShapeParam::init() {
		type_ = param_error;	// this object is required from user
		stat_ = stat_error;
		type_name_.clear();
		max_ = min_ = p1_ = p2_ = 0.0;
		nvalues_ = 0;
	} // ShapeParam::init()

	void ShapeParam::clear() {
		type_name_.clear();
		type_ = param_null;
		stat_ = stat_null;
		max_ = min_ = p1_ = p2_ = 0.0;
		nvalues_ = 0;
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
					const float_t ztilt, const float_t xyrot, shape_param_list_t& param_list) :
					key_(key), name_(name), originvec_(origin),
					ztilt_(ztilt), xyrotation_(xyrot) {
		for(shape_param_iterator_t i = param_list.begin(); i != param_list.end(); i ++) {
			parse_param((*i).second);
			insert_param((*i));
		} // for
	} // Shape::Shape()

	// not used
	Shape::Shape(const std::string& key, const ShapeName name) : key_(key), name_(name) {
		originvec_[0] = originvec_[1] = originvec_[2] = 0.0;
		ztilt_ = 0;
		xyrotation_ = 0;
	} // Shape::Shape()


	void Shape::init() {	// the user needs to provide the shape, no defaults
		key_.clear();
		params_.clear();
		name_ = shape_error;
		name_str_.clear();
		originvec_[0] = originvec_[1] = originvec_[2] = 0.0;
		ztilt_ = xyrotation_ = 0.0;
	} // Shape::init()


	void Shape::clear() {
		key_.clear();
		params_.clear();
		name_ = shape_null;
		name_str_.clear();
		originvec_[0] = originvec_[1] = originvec_[2] = 0.0;
		ztilt_ = xyrotation_ = 0.0;
	} // Shape::clear()


	bool Shape::parse_param(const ShapeParam& param) const {
		// do some micro error checking
		return true;
	} // Shape::parse_param()


	bool Shape::insert_param(const std::pair <std::string, ShapeParam>& param) {
		// check if it already exists ...
		params_[param.first] = param.second;
		return true;
	} // Shape::insert_param()


	bool Shape::insert_param(const std::string& type, const ShapeParam& param) {
		// check if it already exists ...
		params_[type] = param;
		return true;
	} // Shape::insert_param()


	void Shape::print() {
		std::cout << " key_ = " << key_ << std::endl
					<< " name_ = " << name_ << std::endl
					<< " name_str_ = " << name_str_ << std::endl
					<< " originvec_ = [" << originvec_[0] << ", " << originvec_[1]
					<< ", " << originvec_[2] << "]" << std::endl
					<< " ztilt_ = " << ztilt_ << std::endl
					<< " xyrotation_ = " << xyrotation_ << std::endl
					<< " params_: " << params_.size() << std::endl;
		for(shape_param_iterator_t i = params_.begin(); i != params_.end(); i ++) {
			(*i).second.print();
		} // for
		std::cout << std::endl;
	} // Shape::print()

} // namespace hig
