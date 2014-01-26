/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: common.hpp
 *  Created: Jun 13, 2012
 *  Modified: Sun 26 Jan 2014 09:58:05 AM PST
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

#ifndef __MODEL_COMMON_HPP__
#define __MODEL_COMMON_HPP__

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

#endif /* __MODEL_COMMON_HPP__ */
