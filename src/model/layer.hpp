/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: layer.hpp
 *  Created: Jun 13, 2012
 *  Modified: Wed 08 Jan 2014 01:00:55 PM PST
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

#ifndef _LAYER_HPP_
#define _LAYER_HPP_

#include <string>
#include <map>
#include <unordered_map>

#include "../common/typedefs.hpp"
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

			/* modifiers (updat) */
			bool update_param(const std::string&, float_t);

			void print();

	}; // class Layer

	typedef std::map <int, Layer> layer_list_t;		// need the layers sorted on order
	typedef layer_list_t::iterator layer_iterator_t;
	typedef layer_list_t::const_iterator layer_citerator_t;
	typedef std::unordered_map <std::string, int> layer_key_t;

} // namespace hig

#endif /* _LAYER_HPP_ */
