/***
  *  Project:
  *
  *  File: definitions_mic.hpp
  *  Created: Apr 02, 2013
  *  Modified: Tue 02 Apr 2013 01:05:28 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __DEFINITIONS_MIC_HPP__
#define __DEFINITIONS_MIC_HPP__

namespace hig {

	#define MIC_ALLOC	alloc_if(1) free_if(0)
	#define MIC_FREE	alloc_if(0) free_if(1)
	#define MIC_REUSE	alloc_if(0) free_if(0)

} // namespace hig

#endif // __DEFINITIONS_MIC_HPP__
