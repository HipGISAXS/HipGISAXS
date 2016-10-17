/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: parameters_cpu.hpp
 *  Created: Feb 28, 2013
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

#ifndef __PARAMETERS_CPU_HPP__
#define __PARAMETERS_CPU_HPP__

namespace hig {

	/***
	 * parameters for computations on cpu
	 */

	const unsigned int CPU_BLOCK_X_ = 100;
	const unsigned int CPU_BLOCK_Y_ = 256; // 20;
	const unsigned int CPU_BLOCK_Z_ = 256; // 15;
	const unsigned int CPU_BLOCK_T_ = 256; // 2000;

	// padding shape definitions
	const unsigned int CPU_T_PROP_SIZE_ = 8;

} // namespace hig

#endif // __PARAMETERS_CPU_HPP_
