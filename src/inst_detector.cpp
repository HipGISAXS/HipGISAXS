/***
  *  $Id: inst_detector.cpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project:
  *
  *  File: inst_detector.cpp
  *  Created: Jun 12, 2012
  *  Modified: Wed 22 Aug 2012 12:37:53 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "inst_detector.hpp"

namespace hig {

	DetectorParams::DetectorParams() { }
	DetectorParams::~DetectorParams() { }

	void DetectorParams::init() {	// fill default values
		origin_ = "bl";
		total_pixels_[0] = total_pixels_[1] = 1000;
		pixel_size_ = 0.2;
		sd_distance_ = 2000;
		direct_beam_[0] = direct_beam_[1] = 0;	// should be center of detector
	} // init()

} // namespace hig
