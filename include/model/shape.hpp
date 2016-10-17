/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: shape.hpp
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

#ifndef __SHAPE_HPP__
#define __SHAPE_HPP__

#include <string>
#include <unordered_map>

#include <common/enums.hpp>
#include <common/globals.hpp>
#include <model/common.hpp>

namespace hig {

  // analytical shapes can possibly be pre-defined using class inheritance ...
  class ShapeParam {
    private:
      std::string type_name_;  // this is the key for a param
      ShapeParamType type_;  // change to this as the key ... ?
      StatisticType stat_;
      real_t max_;
      real_t min_;
      real_t p1_;        /* mean */
      real_t p2_;        /* deviation */
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
      real_t max() const { return max_; }
      real_t min() const { return min_; }
      StatisticType stat() const { return stat_; }
      real_t p1() const { return p1_; }
      real_t mean() const { return p1_; }
      real_t p2() const { return p2_; }
      real_t deviation() const { return p2_; }
      int nvalues() const { return nvalues_; }

      /* setters */
      void set(void) { isvalid_ = true; }
      void unset(void) { isvalid_ = false; }
      void type_name(std::string s) { type_name_ = s; }
      void type(ShapeParamType s) { type_ = s; }
      void stat(StatisticType s) { stat_ = s; }
      void max(real_t d) { max_ = d; }
      void min(real_t d) { min_ = d; }
      void p1(real_t d) { p1_ = d; }
      void mean(real_t d) { p1_ = d; }
      void p2(real_t d) { p2_ = d; }
      void deviation(real_t d) { p2_ = d; }
      void nvalues(real_t d) { nvalues_ = int(d); }

      /* modifiers (update) */
      bool update_param(const std::string&, real_t);

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
      ShapeName name_;      /* name of a predefined shape, or custom file */
      std::string name_str_;    /* shape file name: used for custom shape only */
      RefractiveIndex refindex_;
      vector3_t originvec_;
      real_t xrot_, yrot_, zrot_;
      shape_param_list_t params_;  // a map of all params (key is the type)

      bool parse_param(const ShapeParam& param) const;

    public:
      Shape();
      Shape(const std::string& key, const ShapeName name);
      Shape(const std::string& key, const ShapeName name, const vector3_t& origin,
        const real_t zrot, const real_t yrot, real_t xrot, shape_param_list_t& param_list);
      ~Shape(void);

      void init();
      void clear();

      /* setters */

      void key(const std::string& s) { key_ = s; }
      void name(ShapeName s) { name_ = s; }
      void name_str(const std::string& s) { name_str_ = s; }

      void originvec(vector3_t vec) { originvec_ = vec; }
      void originvec(real_t a, real_t b, real_t c) {
        originvec_[0] = a; originvec_[1] = b; originvec_[2] = c; }

      void refindex(RefractiveIndex & n) { refindex_ = n; }
      void refindex_beta(real_t b) { refindex_.beta(b); }
      void refindex_delta(real_t d) { refindex_.delta(d); }
      void zrot(real_t d) { zrot_ = d; }
      void yrot(real_t d) { yrot_ = d; }
      void xrot(real_t d) { xrot_ = d; }

      Shape& operator=(const Shape& sh) {
        key_ = sh.key_;
        name_ = sh.name_;
        name_str_ = sh.name_str_;
        originvec_ = sh.originvec_;
        refindex_ = sh.refindex_;
        zrot_ = sh.zrot_;
        yrot_ = sh.yrot_;
        xrot_ = sh.xrot_;
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
      real_t zrot() const { return zrot_; }
      real_t yrot() const { return yrot_; }
      real_t xrot() const { return xrot_; }
      vector3_t originvec() const { return originvec_; }
      RefractiveIndex refindex() const { return refindex_; }
      complex_t one_minus_n2() const { return refindex_.one_minus_n2(); }
      std::string filename() const { return name_str_; }
      shape_param_list_t& param_list() { return params_; }
      shape_param_iterator_t param_begin() { return params_.begin(); }
      shape_param_iterator_t param_end() { return params_.end(); }

      /* modifiers (updates) */
      bool update_param(const std::string&, real_t);

      void print();

  }; // class Shape

  typedef std::unordered_map <std::string, Shape> shape_list_t;
  typedef shape_list_t::iterator shape_iterator_t;

} // namespace hig

#endif // __SHAPE_HPP__
