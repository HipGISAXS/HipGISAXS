/***
  *  $Id: compute_params.hpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS - High Performance GISAXS
  *
  *  File: compute_params.hpp
  *  Created: Jun 05, 2012
  *  Modified: Tue 21 Aug 2012 06:23:08 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _COMPUTE_PARAMS_HPP_
#define _COMPUTE_PARAMS_HPP_

#include <iostream>
#include <string>

#include "globals.hpp"
#include "enums.hpp"

namespace hig {

	class ComputeParams {	// just one instance of this
		friend class HiGInput;

		private:
			std::string pathprefix_;
			std::string runname_;
			std::string method_;	// change to enum - "dwba" ?
			struct OutputRegion {
				OutputRegionType type_;
				vector2_t minpoint_;
				vector2_t maxpoint_;
				OutputRegion() { }
			} output_region_;
			vector2_t resolution_;
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

			void nslices(float_t d) { nslices_ = (unsigned int) d; }

			/* getters */
			//OutputRegion output_region() { return output_region_; }

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
							<< std::endl;
			} // print()

	}; // class ComputeParams

} // namespace hig

#endif /* _COMPUTE_PARAMS_HPP_ */
