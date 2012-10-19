/***
  *  $Id: compute_params.cpp 27 2012-07-15 05:37:29Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: compute_params.cpp
  *  Created: Jun 05, 2012
  *  Modified: Tue 16 Oct 2012 02:26:45 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "compute_params.hpp"

namespace hig {

	ComputeParams::ComputeParams() { }

	void ComputeParams::init() {	// sets all default values
		pathprefix_ = ".";
		runname_ = "run_" + timestamp();
		method_ = "dwba";
		output_region_.type_ = region_qspace;
		output_region_.minpoint_[0] = -1;
		output_region_.minpoint_[1] = -1;
		output_region_.maxpoint_[0] = -1;
		output_region_.maxpoint_[1] = -1;
		resolution_[0] = 1;
		resolution_[1] = 1;
		nslices_ = 0;
	} // ComputeParams::init()


	std::string ComputeParams::timestamp() {	// make this an independent utility ...
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[16];

		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(buffer, 16, "%Y%m%d_%H%M%S", timeinfo);

		return std::string(buffer);
	} // ComputeParams::timestamp()

} // namespace hig
