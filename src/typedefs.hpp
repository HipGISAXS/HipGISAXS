/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: typedefs.hpp
  *  Created: Jul 08, 2012
  *  Modified: Thu 02 May 2013 11:51:18 PM PDT
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
#ifdef INTEL_SB_AVX
	#include <immintrin.h>
#elif defined __SSE3__
	#include <pmmintrin.h>
#elif defined USE_MIC
	#include <immintrin.h>
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


	// SSE/AVX vector types:

	#ifdef INTEL_SB_AVX
		typedef __m256						avx_m256_t;
		typedef struct {
			__m256 xvec;
			__m256 yvec;
		}									avx_m256c_t;
	#elif defined __SSE3__
		typedef __m128						sse_m128_t;
		typedef struct {
			__m128 xvec;
			__m128 yvec;
		}									sse_m128c_t;
	#endif

	#ifdef USE_MIC
		typedef __m512						mic_m512_t;
		typedef struct {
			__m512 xvec;
			__m512 yvec;
		}									mic_m512c_t;
	#endif

	// TODO: handle multiprecision? ...


} // namespace


#endif /* _TYPEDEFS_HPP_ */
