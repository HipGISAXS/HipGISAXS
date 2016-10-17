/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: parameters_mic.hpp
 *  Created: Apr 02, 2013
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __PARAMETERS_MIC_HPP__
#define __PARAMETERS_MIC_HPP__

namespace hig {

	/***
	 * tunable parameters for computations on MIC
	 */

	// hyperblock size
	const int MIC_BLOCK_X_ = 100;
	const int MIC_BLOCK_Y_ = 128;	//40;
	const int MIC_BLOCK_Z_ = 160;	//30;
	const int MIC_BLOCK_T_ = 1024;	// TODO: has to be multiple of 16 for now

	// max number of OMP threads on MIC
	__attribute__((target(mic)))
	const int MIC_OMP_NUM_THREADS_ = 240;

} // namespace hig

#endif // __PARAMETERS_MIC_HPP__
