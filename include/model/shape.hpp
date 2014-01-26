/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: shape.hpp
 *  Created: Jun 05, 2012
 *  Modified: Sun 26 Jan 2014 09:59:45 AM PST
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

#ifndef _SHAPE_HPP_
#define _SHAPE_HPP_

#include <string>
#include <unordered_map>

#include <common/enums.hpp>
#include <common/globals.hpp>

namespace hig {

	// analytical shapes can possibly be pre-defined using class inheritance ...
	class ShapeParam {
		private:
			std::string type_name_;	// this is the key for a param
			ShapeParamType type_;	// change to this as the key ... ?
			StatisticType stat_;
			float_t max_;
			float_t min_;
			float_t p1_;				/* mean */
			float_t p2_;				/* deviation */
			int nvalues_;
			bool isvalid_;

		public:
			ShapeParam(): stat_(stat_none), isvalid_(false) { init(); };
			~ShapeParam() { }

			void init();
			void clear();

			/* getters */
			bool isvalid() const { return isvalid_; }
			const std::string& type_name() const { return type_name_; }
			ShapeParamType type() const { return type_; }
			float_t max() const { return max_; }
			float_t min() const { return min_; }
			StatisticType stat() const { return stat_; }
			float_t p1() const { return p1_; }
			float_t mean() const { return p1_; }
			float_t p2() const { return p2_; }
			float_t deviation() const { return p2_; }
			int nvalues() const { return nvalues_; }

			/* setters */
			void set(void) { isvalid_ = true; }
			void unset(void) { isvalid_ = false; }
			void type_name(std::string s) { type_name_ = s; }
			void type(ShapeParamType s) { type_ = s; }
			void stat(StatisticType s) { stat_ = s; }
			void max(float_t d) { max_ = d; }
			void min(float_t d) { min_ = d; }
			void p1(float_t d) { p1_ = d; }
			void mean(float_t d) { p1_ = d; }
			void p2(float_t d) { p2_ = d; }
			void deviation(float_t d) { p2_ = d; }
			void nvalues(float_t d) { nvalues_ = int(d); }

			/* modifiers (update) */
			bool update_param(const std::string&, float_t);

			/* copy constructor */
			ShapeParam& operator=(ShapeParam const& a) {
				type_name_ = a.type_name_;
				type_ = a.type_;
				max_ = a.max_;
				min_ = a.min_;
				stat_ = a.stat_;
				p1_ = a.p1_;
				p2_ = a.p2_;
				nvalues_ = a.nvalues_;
				isvalid_ = a.isvalid_;

				return *this;
			} // operator=()

			void print();

			friend class Shape;

	}; // class ShapeParam

	typedef std::unordered_map <std::string, ShapeParam> shape_param_list_t;
	typedef shape_param_list_t::iterator shape_param_iterator_t;

	class Shape {
		private:
			std::string key_;
			ShapeName name_;			/* name of a predefined shape, or custom file */
			std::string name_str_;		/* shape file name: used for custom shape only */
			vector3_t originvec_;
			float_t ztilt_;				/* shape tau */
			float_t xyrotation_;		/* shape eta */
			shape_param_list_t params_;	// a map of all params (key is the type)

			bool parse_param(const ShapeParam& param) const;

		public:
			Shape();
			Shape(const std::string& key, const ShapeName name);
			Shape(const std::string& key, const ShapeName name, const vector3_t& origin,
				const float_t ztilt, const float_t xyrot, shape_param_list_t& param_list);
			~Shape(void);

			void init();
			void clear();

			/* setters */

			void key(const std::string& s) { key_ = s; }
			void name(ShapeName s) { name_ = s; }
			void name_str(const std::string& s) { name_str_ = s; }

			void originvec(vector3_t vec) { originvec_ = vec; }
			void originvec(float_t a, float_t b, float_t c) {
				originvec_[0] = a; originvec_[1] = b; originvec_[2] = c; }

			void ztilt(float_t d) { ztilt_ = d; }
			void xyrotation(float_t d) { xyrotation_ = d; }

			Shape& operator=(const Shape& sh) {
				key_ = sh.key_;
				name_ = sh.name_;
				name_str_ = sh.name_str_;
				originvec_ = sh.originvec_;
				ztilt_ = sh.ztilt_;
				xyrotation_ = sh.xyrotation_;
				for(auto i = sh.params_.begin(); i != sh.params_.end(); i ++) {
					params_[(*i).first] = (*i).second;
				} // for
				return *this;
			} // operator=()

			bool insert_param(const std::pair<std::string, ShapeParam>&);
			bool insert_param(const std::string&, const ShapeParam&);

			/* getters */

			std::string key() const { return key_; }
			ShapeName name() const { return name_; }
			float_t ztilt() const { return ztilt_; }
			float_t xyrotation() const { return xyrotation_; }
			vector3_t originvec() const { return originvec_; }
			std::string filename() const { return name_str_; }
			shape_param_list_t& param_list() { return params_; }
			shape_param_iterator_t param_begin() { return params_.begin(); }
			shape_param_iterator_t param_end() { return params_.end(); }

			/* modifiers (updates) */
			bool update_param(const std::string&, float_t);

			void print();

	}; // class Shape

	typedef std::unordered_map <std::string, Shape> shape_list_t;	// ??????
	typedef shape_list_t::iterator shape_iterator_t;

} // namespace hig

#endif /* _SHAPE_HPP_ */
