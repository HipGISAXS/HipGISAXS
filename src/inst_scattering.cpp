/***
  *  $Id: inst_scattering.cpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: inst_scattering.cpp
  *  Created: Jun 12, 2012
  *  Modified: Mon 01 Oct 2012 11:14:39 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "inst_scattering.hpp"

namespace hig {

	ScatteringParams::ScatteringParams() { }
	ScatteringParams::~ScatteringParams() { }

	void ScatteringParams::init() {
		expt_ = "gisaxs";
		alpha_i_.min_ = alpha_i_.max_ = alpha_i_.step_ = 0.1;			// check
		inplane_rot_.min_ = inplane_rot_.max_ = inplane_rot_.step_ = 0;	// check
		tilt_.min_ = tilt_.max_ = tilt_.step_ = 0;						// check
		photon_.value_ = 10000; photon_.unit_ = "ev";
		polarization_ = "s";
		coherence_ = 300;
		spot_area_ = 1;
		smearing_[0] = smearing_[1] = smearing_[2] = 1;
	} // ScatteringParams::init()


} // namespace hig
