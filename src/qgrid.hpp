/***
  *  $Id: qgrid.hpp 33 2012-08-06 16:22:01Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: qgrid.hpp
  *  Created: Jun 17, 2012
  *  Modified: Mon 01 Oct 2012 11:16:28 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _QGRID_HPP_
#define _QGRID_HPP_

#include <vector>

#include "hig_input.hpp"

namespace hig {

	class QGrid {								// should these all be complex ... ???
												// use the new 'std::array<>' instead ... ?
		/* data types for qgrid data */
		typedef std::vector<float_t> qvec_t;
		typedef std::vector<complex_t> cqvec_t;
		/* iterators for qgrid data types */
		typedef qvec_t::const_iterator qvec_iter_t;
		typedef cqvec_t::const_iterator cqvec_iter_t;

		private:
			qvec_t qx_;
			qvec_t qy_;
			qvec_t qz_;
			cqvec_t qz_extended_;

			/* singleton */
			QGrid() { }
			QGrid(const QGrid&);
			QGrid& operator=(const QGrid&);

			bool pixel_to_kspace(vector2_t, float_t, float_t, float_t, float_t, vector2_t, vector3_t&);
			bool kspace_to_pixel();		// not implemented yet ...

		public:
			static QGrid& instance() {
				static QGrid qgrid;
				return qgrid;
			} // instance()

			/* create the Q-grid */
			bool create(float_t, float_t, float_t);
			bool create_qz_extended(float_t, float_t, complex_t, complex_t);
			bool create_test();

			/* sizes */
			int nqx() const { return qx_.size(); }
			int nqy() const { return qy_.size(); }
			int nqz() const { return qz_.size(); }
			int nqz_extended() const { return qz_extended_.size(); }

			/* value accessors */
			float_t qx(int i) const { return qx_[i]; }	// do some error checking also ...
			float_t qy(int i) const { return qy_[i]; }
			float_t qz(int i) const { return qz_[i]; }
			complex_t qz_extended(int i) const { return qz_extended_[i]; }
			float_t qx(unsigned int i) const { return qx_[i]; }	// do some error checking also ...
			float_t qy(unsigned int i) const { return qy_[i]; }
			float_t qz(unsigned int i) const { return qz_[i]; }
			complex_t qz_extended(unsigned int i) const { return qz_extended_[i]; }

			/* iterator helpers */
			qvec_iter_t qx_begin() const { return qx_.begin(); }
			qvec_iter_t qy_begin() const { return qy_.begin(); }
			qvec_iter_t qz_begin() const { return qz_.begin(); }
			cqvec_iter_t qz_extended_begin() const { return qz_extended_.begin(); }
			qvec_iter_t qx_end() const { return qx_.end(); }
			qvec_iter_t qy_end() const { return qy_.end(); }
			qvec_iter_t qz_end() const { return qz_.end(); }
			cqvec_iter_t qz_extended_end() const { return qz_extended_.end(); }

	}; // class QGrid

} // namespace hig

#endif /* QGRID_HPP_ */
