/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: inst_detector.cpp
 *  Created: Jun 12, 2012
 *  Modified: Tue 16 Jul 2013 11:51:07 AM PDT
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

#include "inst_detector.hpp"
#include "default_values.hpp"

namespace hig {

	DetectorParams::DetectorParams() { }
	DetectorParams::~DetectorParams() { }

	void DetectorParams::init() {	// fill default values
		origin_ = DEFAULT_ORIGIN_;
		total_pixels_[0] = DEFAULT_TOTAL_PIXELS_Y_;
		total_pixels_[1] = DEFAULT_TOTAL_PIXELS_Z_;
		pixel_size_ = DEFAULT_PIXEL_SIZE_;
		sd_distance_ = DEFAULT_SDD_;
		direct_beam_[0] = DEFAULT_DIRECT_BEAM_Y_;
		direct_beam_[1] = DEFAULT_DIRECT_BEAM_Z_;	// should be center of the detector
	} // init()

	/*void DetectorParams::init() {	// fill default values
		origin_ = "bl";
		total_pixels_[0] = total_pixels_[1] = 1000;
		pixel_size_ = 0.2;
		sd_distance_ = 2000;
		direct_beam_[0] = direct_beam_[1] = 0;	// should be center of detector
	} // init()*/

} // namespace hig
