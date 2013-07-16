/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: inst_scattering.cpp
 *  Created: Jun 12, 2012
 *  Modified: Tue 16 Jul 2013 11:51:09 AM PDT
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
		spot_area_ = 0.01; //0.001;
		smearing_[0] = smearing_[1] = smearing_[2] = 1;
	} // ScatteringParams::init()


} // namespace hig
