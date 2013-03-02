/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Sat 02 Mar 2013 09:23:29 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _TYPEDEFS_HPP_
#define _TYPEDEFS_HPP_

#include <vector>
#include <complex>
#include <cuComplex.h>

namespace hig {

#ifdef DOUBLEP	// double precision
	typedef double						float_t;
	typedef cuDoubleComplex				cucomplex_t;
#else			// single precision
	typedef float						float_t;
	typedef cuFloatComplex				cucomplex_t;
#endif	// DOUBLEP

	// TODO: multiprecision? ...

	typedef std::complex<float_t>		complex_t;
	typedef std::vector<float_t> 		float_vec_t;
	typedef std::vector<complex_t>		complex_vec_t;
	typedef std::vector<cucomplex_t>	cucomplex_vec_t;

} // namespace


#endif /* _TYPEDEFS_HPP_ */
