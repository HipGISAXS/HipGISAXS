/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Tue 02 Apr 2013 09:42:11 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _TYPEDEFS_HPP_
#define _TYPEDEFS_HPP_

#include <vector>
#include <complex>
#include <cuComplex.h>

namespace hig {

	#if defined USE_MIC
		typedef struct {	// serialized complex
			double x;
			double y; }		double2_t;

		typedef struct {	// serialized complex
			float x;
			float y; }		float2_t;
	#endif

	#ifdef DOUBLEP						// double precision
		typedef double					float_t;
		typedef cuDoubleComplex			cucomplex_t;
		#if defined USE_MIC
			typedef double2_t			scomplex_t;
		#endif
	#else								// single precision
		typedef float					float_t;
		typedef cuFloatComplex			cucomplex_t;
		#if defined USE_MIC
			typedef float2_t			scomplex_t;
		#endif
	#endif

	// TODO: handle multiprecision? ...


	typedef std::complex<float_t>		complex_t;
	typedef std::vector<float_t> 		float_vec_t;
	typedef std::vector<complex_t>		complex_vec_t;
	typedef std::vector<cucomplex_t>	cucomplex_vec_t;

} // namespace


#endif /* _TYPEDEFS_HPP_ */
