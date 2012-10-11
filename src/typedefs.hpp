/***
  *  $Id: typedefs.hpp 30 2012-07-24 20:17:58Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Thu 04 Oct 2012 05:07:04 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _TYPEDEFS_HPP_
#define _TYPEDEFS_HPP_

#include <complex>
#include <cuComplex.h>
#include <vector>

namespace hig {

#ifdef DOUBLEP	// double precision
	typedef double					float_t;
	typedef std::complex<double>	complex_t;
	typedef cuDoubleComplex			cucomplex_t;
#else			// single precision
	typedef float					float_t;
	typedef std::complex<float>		complex_t;
	typedef cuFloatComplex			cucomplex_t;
#endif

	typedef std::vector<float_t> 	float_vec_t;
	typedef std::vector<complex_t>	complex_vec_t;

} // namespace


#endif /* _TYPEDEFS_HPP_ */
