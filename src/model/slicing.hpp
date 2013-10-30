/***
  *  Project:
  *
  *  File: slicing.hpp
  *  Created: May 13, 2013
  *  Modified: Wed 30 Oct 2013 11:03:57 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __SLICING_HPP__
#define __SLICING_HPP__

#include <vector>

#include "../common/typedefs.hpp"
#include "common.hpp"

namespace hig {

	class Slice {
		private:
			float_t z_top_;
			float_t z_bottom_;
			float_t thickness_;					// slice thickness in z direction
												// for now uniform slicing, TODO: adaptive slicing
			RefractiveIndex avg_refindex_;		// average ref index of the slice
			RefractiveIndex layer_refindex_;	// ref index of the layer the slice is in
												// TODO: for now assumption = slice is within single layer
			float_t slice_volume_;
			float_t sample_volume_;				// is this really needed ??? ...

		public:
			Slice();
			~Slice();

			void init();
			void clear();

			bool compute_slice_data();
			// ...

	}; // class Slice


	typedef std::vector <Slice> slice_list_t;
	typedef slice_list_t::iterator slice_iterator_t;
	typedef slice_list_t::const_iterator slice_citerator_t;


} // namespace hig

#endif // __SLICING_HPP__
