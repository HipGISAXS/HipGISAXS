/***
  *  $Id: typedefs.hpp 30 2012-07-24 20:17:58Z asarje $
  *
  *  Project:
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Tue 24 Jul 2012 12:35:18 PM PDT
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
