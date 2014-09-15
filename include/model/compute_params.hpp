/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: compute_params.hpp
 *  Created: Jun 05, 2012
 *  Modified: Sun 14 Sep 2014 09:06:11 AM PDT
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

#ifndef _COMPUTE_PARAMS_HPP_
#define _COMPUTE_PARAMS_HPP_

#include <iostream>
#include <string>

#include <common/globals.hpp>
#include <common/enums.hpp>

namespace hig {

	class ComputeParams {	// just one instance of this ... TODO ...
		friend class HiGInput;

		private:
			std::string pathprefix_;
			std::string runname_;
			std::string method_;	// TODO: ... change to enum - "dwba" ?
			StructCorrelationType correlation_;		/* grain/ensemble correlation type */
			struct OutputRegion {
				OutputRegionType type_;
				vector2_t minpoint_;
				vector2_t maxpoint_;
				OutputRegion() { }
			} output_region_;
			vector2_t resolution_;
			std::string palette_;
			unsigned int nslices_;

			std::string timestamp();

		public:
			ComputeParams();

			void init();
			void clear();

			/* setters */

			const std::string& pathprefix() const { return pathprefix_; }
			const std::string& runname() const { return runname_; }

			void pathprefix(std::string s) { pathprefix_ = s; }
			void runname(std::string s) { runname_ = s + "_" + timestamp(); }
			void method(std::string s) { method_ = s; }

			void output_region_type(OutputRegionType o) { output_region_.type_ = o; }
			void output_region_minpoint(vector2_t v) { output_region_.minpoint_ = v; }
			void output_region_maxpoint(vector2_t v) { output_region_.maxpoint_ = v; }
			void output_region_minpoint(float_t v, float_t w) {
				output_region_.minpoint_[0] = v; output_region_.minpoint_[1] = w; }
			void output_region_maxpoint(float_t v, float_t w) {
				output_region_.maxpoint_[0] = v; output_region_.maxpoint_[1] = w; }
			void resolution(float_t a, float_t b) { resolution_[0] = a; resolution_[1] = b; }

			void palette(std::string p) { palette_ = p; }
			void nslices(float_t d) { nslices_ = (unsigned int) d; }
			void structcorrelation(StructCorrelationType c) { correlation_ = c; }

			/* getters */
			//OutputRegion output_region() { return output_region_; }

			/* modifiers (updates) */
			bool update_param(const std::string&, float_t);

			void print() {
				std::cout << " pathprefix_ = " << pathprefix_ << std::endl
							<< " runname_ = " << runname_ << std::endl
							<< " method_ = " << method_ << std::endl
							<< " output_region_: type = " << output_region_.type_
							<< ", minpoint = [" << output_region_.minpoint_[0] << ", "
							<< output_region_.minpoint_[1] << "], maxpoint = ["
							<< output_region_.maxpoint_[0] << ", "
							<< output_region_.maxpoint_[1] << "]" << std::endl
							<< " resolution_ = [" << resolution_[0] << ", "
							<< resolution_[1] << "]" << std::endl
							<< " nslices_ = " << nslices_ << std::endl
							<< " palette_ = " << palette_ << std::endl
							<< std::endl;
			} // print()

	}; // class ComputeParams

} // namespace hig

#endif /* _COMPUTE_PARAMS_HPP_ */
