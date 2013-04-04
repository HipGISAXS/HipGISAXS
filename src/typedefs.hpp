/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Wed 03 Apr 2013 07:27:44 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _TYPEDEFS_HPP_
#define _TYPEDEFS_HPP_

#include <vector>
#include <complex>
#ifdef USE_GPU
	#include <cuComplex.h>
#endif

namespace hig {

	#if defined USE_MIC
		typedef struct {	// serialized complex
			double x;
			double y; }						double2_t;

		typedef struct {	// serialized complex
			float x;
			float y; }						float2_t;
	#endif

	#ifdef DOUBLEP						// double precision
		typedef double						float_t;
		#ifdef USE_GPU
			typedef cuDoubleComplex			cucomplex_t;
		#endif
		#if defined USE_MIC
			typedef double2_t				scomplex_t;
		#endif
	#else								// single precision
		typedef float						float_t;
		#ifdef USE_GPU
			typedef cuFloatComplex			cucomplex_t;
		#endif
		#if defined USE_MIC
			typedef float2_t				scomplex_t;
		#endif
	#endif

	typedef std::complex<float_t>			complex_t;
	typedef std::vector<float_t> 			float_vec_t;
	typedef std::vector<complex_t>			complex_vec_t;

	#ifdef USE_GPU
		typedef std::vector<cucomplex_t>	cucomplex_vec_t;
	#endif

	// TODO: handle multiprecision? ...


} // namespace


#endif /* _TYPEDEFS_HPP_ */
