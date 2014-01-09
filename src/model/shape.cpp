/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: shape.cpp
 *  Created: Jun 05, 2012
 *  Modified: Wed 08 Jan 2014 12:56:49 PM PST
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

#include <iostream>

#include "shape.hpp"

namespace hig {

	/* ShapeParam */

	void ShapeParam::init() {
		type_ = param_error;	// this object is required from user
		stat_ = stat_none;
		type_name_.clear();
		max_ = min_ = p1_ = p2_ = 0.0;
		nvalues_ = 1;			// default nvalue
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
					<< " ztilt_ = " << ztilt_ << std::endl
					<< " xyrotation_ = " << xyrotation_ << std::endl
					<< " params_: " << params_.size() << std::endl;
		for(shape_param_iterator_t i = params_.begin(); i != params_.end(); i ++) {
			(*i).second.print();
		} // for
		std::cout << std::endl;
	} // Shape::print()


	/**
	 * updates for fitting
	 */


	bool ShapeParam::update_param(const std::string& str, float_t new_val) {
		std::string keyword, rem_str;
		if(!get_first_keyword(str, keyword, rem_str)) return false;
		std::string keyword_name, key;
		if(!get_keyword_name_and_key(keyword, keyword_name, key)) return false;
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


	bool Shape::update_param(const std::string& str, float_t new_val) {
		std::string keyword, rem_str;
		if(!get_first_keyword(str, keyword, rem_str)) return false;
		std::string keyword_name, key;
		if(!get_keyword_name_and_key(keyword, keyword_name, key)) return false;
		switch(TokenMapper::instance().get_keyword_token(keyword_name)) {
			case key_token:
			case shape_name_token:
			case shape_originvec_token:
				std::cerr << "warning: immutable param in '" << str << "'. ignoring." << std::endl;
				break;

			case shape_ztilt_token:
				ztilt_ = new_val;
				break;

			case shape_xyrot_token:
				xyrotation_ = new_val;
				break;

			case shape_param_token:
				if(!params_.at(key).update_param(rem_str, new_val)) return false;
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
