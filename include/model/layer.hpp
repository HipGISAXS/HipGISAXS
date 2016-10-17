/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: layer.hpp
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

#ifndef __LAYER_HPP__
#define __LAYER_HPP__

#include <string>
#include <map>
#include <unordered_map>

#include <common/typedefs.hpp>
#include <model/common.hpp>

namespace hig {

  class Layer {
    private:
      std::string key_;      /* a unique key string */
      int order_;          /* layer order: -1(substrate) 0 1 2 3 ... */
      real_t thickness_;      /* layer thickness */
      real_t z_val_;        /* cumulative distance */
      RefractiveIndex refindex_;  /* layer's refractive index */

    public:
      Layer();
      ~Layer();

      void init();
      void clear();

      /* getters */
      std::string key() const { return key_; }
      RefractiveIndex& refindex() { return refindex_; }
      real_t thickness() const { return thickness_; }
      int order() const { return order_; }
      real_t z_val() const { return z_val_; }
      complex_t one_minus_n2() { return refindex_.one_minus_n2(); }

      /* setters */
      void key(std::string s) { key_ = s; }
      void refindex(RefractiveIndex d) { refindex_ = d; }
      void refindex_delta(real_t d) { refindex_.delta(d); }
      void refindex_beta(real_t d) { refindex_.beta(d); }
      void order(real_t d) { order_ = (int) d; }
      void thickness(real_t d) { thickness_ = d; }
      void z_val(real_t z) { z_val_ = z; }

      /* modifiers (updat) */
      bool update_param(const std::string&, real_t);

      void print();

  }; // class Layer

  typedef std::map <int, Layer> layer_list_t;    // need the layers sorted on order
  typedef layer_list_t::iterator layer_iterator_t;
  typedef layer_list_t::const_iterator layer_citerator_t;
  typedef std::unordered_map <std::string, int> layer_key_t;

  // Defilne multilayer
  class MultiLayer {
    private:
      std::vector<Layer> layers_;

    public:
      // constructors
      MultiLayer();
      
      // index operator
      Layer operator[] (int j) const {
        return layers_[j];
      }

      // hipgisaxs style init and clear
      bool init(const layer_list_t &);
      void clear();

      //gets
      int count(){ return layers_.size(); }
      
      // calculate Ts and Rs of an angle
      complex_vec_t parratt_recursion(real_t, real_t, int);

      // calculate transmission and reflection coefficents
      bool propagation_coeffs(complex_vec_t &, real_t, real_t, int);

      /** DEBUG **/
     void debug_multilayer();
  };
} // namespace hig

#endif // __LAYER_HPP__
