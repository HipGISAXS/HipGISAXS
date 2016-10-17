/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: definitions_mic.hpp
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

#ifndef __DEFINITIONS_MIC_HPP__
#define __DEFINITIONS_MIC_HPP__

namespace hig {

	#define MIC_ALLOC	alloc_if(1) free_if(0)
	#define MIC_FREE	alloc_if(0) free_if(1)
	#define MIC_REUSE	alloc_if(0) free_if(0)

} // namespace hig

#endif // __DEFINITIONS_MIC_HPP__
