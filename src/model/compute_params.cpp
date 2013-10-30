/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: compute_params.cpp
 *  Created: Jun 05, 2012
 *  Modified: Tue 16 Jul 2013 11:49:29 AM PDT
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
		correlation_ = structcorr_null;
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
