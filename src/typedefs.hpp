/***
  *  $Id: typedefs.hpp 30 2012-07-24 20:17:58Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Mon 01 Oct 2012 11:19:44 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _TYPEDEFS_HPP_
#define _TYPEDEFS_HPP_

#include <complex>
#include <cuComplex.h>

namespace hig {

#ifdef DOUBLEP
	typedef double float_t;
	typedef std::complex <float_t> complex_t;
	typedef cuDoubleComplex cucomplex_t;
#else
	typedef float float_t;
	typedef std::complex<float_t> complex_t;
	typedef cuFloatComplex cucomplex_t;
#endif

} // namespace


#endif /* _TYPEDEFS_HPP_ */
