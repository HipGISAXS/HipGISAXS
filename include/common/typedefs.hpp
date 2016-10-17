/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: typedefs.hpp
 *  Created: Jul 08, 2012
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

#ifndef __HIG_TYPEDEFS_HPP__
#define __HIG_TYPEDEFS_HPP__

#include <vector>
#include <map>
#include <complex>
#include <string>

#ifdef USE_GPU
  #include <cuComplex.h>
//#elif defined USE_MIC
//  #include <immintrin.h>
//#elif defined INTEL_SB_AVX
//  #include <immintrin.h>
//#elif defined __SSE3__
//  #include <pmmintrin.h>
#endif

namespace hig {

  #if defined USE_MIC
    typedef struct {  // serialized complex
      double x;
      double y; }                           double2_t;

    typedef struct {  // serialized complex
      float x;
      float y; }                            float2_t;
  #endif

  #ifdef DOUBLEP      // double precision
    typedef double                          real_t;
    #ifdef USE_GPU
      typedef cuDoubleComplex               cucomplex_t;
    #endif
    #if defined USE_MIC
      typedef double2_t                     scomplex_t;
    #endif
  #else               // single precision
    typedef float                           real_t;
    #ifdef USE_GPU
      typedef cuFloatComplex                cucomplex_t;
    #endif
    #if defined USE_MIC
      typedef float2_t                      scomplex_t;
    #endif
  #endif

  typedef std::complex<real_t>             complex_t;
  typedef std::vector<real_t>              real_vec_t;
  typedef std::vector<complex_t>           complex_vec_t;
  typedef std::vector<unsigned int>        uint_vec_t;

  typedef std::pair <real_t, real_t>       real_pair_t;

  #ifdef USE_GPU
    typedef std::vector<cucomplex_t>        cucomplex_vec_t;
/*  #elif defined INTEL_SB_AVX
    // SSE/AVX vector types:
    typedef __m256                          avx_m256_t;
    typedef struct {
      __m256 xvec;
      __m256 yvec;
    }                                       avx_m256c_t;
  #elif defined __SSE3__
    typedef __m128                          sse_m128_t;
    typedef struct {
      __m128 xvec;
      __m128 yvec;
    }                                       sse_m128c_t; */
  #endif

/*  #ifdef USE_MIC
    typedef __m512                          mic_m512_t;
    typedef struct {
      __m512 xvec;
      __m512 yvec;
    }                                       mic_m512c_t;
  #endif */

  typedef std::string                       string_t;
  typedef std::map <std::string, real_t>    map_t;


  // Triangle with vertices counter-clockwise order
  typedef struct Triangle {
      real_t v1[3];
      real_t v2[3];
      real_t v3[3]; 
  }                                         triangle_t;
 

  // TODO: handle multiprecision? ...

} // namespace


#endif /* __HIG_TYPEDEFS_HPP__ */
