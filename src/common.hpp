/***
  *  $Id: common.hpp 33 2012-08-06 16:22:01Z asarje $
  *
  *  Project:
  *
  *  File: common.hpp
  *  Created: Jun 13, 2012
  *  Modified: Thu 02 Aug 2012 05:19:34 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _COMMON_HPP_
#define _COMMON_HPP_

namespace hig {

	class RefractiveIndex {
		private:
			float_t delta_;
			float_t beta_;

		public:
			RefractiveIndex(): delta_(0.0), beta_(0.0) { }
			RefractiveIndex(float_t delta, float_t beta): delta_(delta), beta_(beta) { }
			~RefractiveIndex() { }

			float_t delta() { return delta_; }
			float_t beta() { return beta_; }

			void delta(float_t a) { delta_ = a; }
			void beta(float_t a) { beta_ = a; }

			void init();
			void clear();
	}; // class RefractiveIndex

} // namespace hig

#endif /* _COMMON_HPP_ */
