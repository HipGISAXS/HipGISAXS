/***
  *  $Id: inst_detector.cpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: inst_detector.cpp
  *  Created: Jun 12, 2012
  *  Modified: Fri 12 Jul 2013 02:39:38 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
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
