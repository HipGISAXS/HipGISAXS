/***
  *  $Id: layer.hpp 27 2012-07-15 05:37:29Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: layer.hpp
  *  Created: Jun 13, 2012
  *  Modified: Mon 01 Oct 2012 11:14:53 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _LAYER_HPP_
#define _LAYER_HPP_

#include <string>
#include <map>
#include <unordered_map>

#include "typedefs.hpp"
#include "common.hpp"

namespace hig {

	class Layer {
		private:
			std::string key_;			/* a unique key string */
			int order_;					/* layer order: -1(substrate) 0 1 2 3 ... */
			float_t thickness_;			/* layer thickness */
			float_t z_val_;				/* cumulative distance */
			RefractiveIndex refindex_;	/* layer's refractive index */

		public:
			Layer();
			~Layer();

			void init();
			void clear();

			/* getters */
			std::string key() const { return key_; }
			RefractiveIndex& refindex() { return refindex_; }
			float_t thickness() const { return thickness_; }
			int order() const { return order_; }
			float_t z_val() const { return z_val_; }

			/* setters */
			void key(std::string s) { key_ = s; }
			void refindex_delta(float_t d) { refindex_.delta(d); }
			void refindex_beta(float_t d) { refindex_.beta(d); }
			void order(float_t d) { order_ = (int) d; }
			void thickness(float_t d) { thickness_ = d; }
			void z_val(float_t z) { z_val_ = z; }

			void print();

	}; // class Layer

	typedef std::map <int, Layer> layer_list_t;		// need the layers sorted on order
	typedef layer_list_t::iterator layer_iterator_t;
	typedef layer_list_t::const_iterator layer_citerator_t;
	typedef std::unordered_map <std::string, int> layer_key_t;

} // namespace hig

#endif /* _LAYER_HPP_ */
