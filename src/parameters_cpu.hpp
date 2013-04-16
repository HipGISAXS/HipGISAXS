/***
  *  Project:
  *
  *  File: parameters_cpu.hpp
  *  Created: Feb 28, 2013
  *  Modified: Thu 28 Feb 2013 11:26:55 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __PARAMETERS_CPU_HPP__
#define __PARAMETERS_CPU_HPP__

namespace hig {

	/***
	 * parameters for computations on cpu
	 */

	const unsigned int CPU_BLOCK_X_ = 100;
	const unsigned int CPU_BLOCK_Y_ = 20;
	const unsigned int CPU_BLOCK_Z_ = 15;
	const unsigned int CPU_BLOCK_T_ = 2000;

	// padding shape definitions
	const unsigned int CPU_T_PROP_SIZE_ = 8;

} // namespace hig

#endif // __PARAMETERS_CPU_HPP_
