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
#elif defined INTEL_AVX
  #include <immintrin.h>
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

  typedef std::complex<real_t>              complex_t;
  typedef std::vector<real_t>               real_vec_t;
  typedef std::vector<complex_t>            complex_vec_t;
  typedef std::vector<unsigned int>         uint_vec_t;
  typedef std::pair <real_t, real_t>        real_pair_t;
  typedef std::string                       string_t;
  typedef std::map <std::string, real_t>    map_t;

  // triangle with vertices in counter-clockwise order
  typedef struct {
      real_t v1[3];   // coords of v1
      real_t v2[3];   // coords of v2
      real_t v3[3];   // coords of v3
  }                                         triangle_t;

  #ifdef USE_GPU

    typedef std::vector<cucomplex_t>        cucomplex_vec_t;

  #elif defined INTEL_AVX

    // AVX vector types
    #ifdef DOUBLEP  // double precision

      typedef __m256d                       avx_m256_t;
      typedef struct {  // complex
        __m256d real;
        __m256d imag;
      }                                     avx_m256c_t;

    #else           // single precision

      typedef __m256                        avx_m256_t;
      typedef struct {  // complex
        __m256 real;
        __m256 imag;
      }                                     avx_m256c_t;

    #endif  // DOUBLEP

    typedef struct {  // triangle structure of arrays
      // in DP, each object defines 4 triangles = 4 * 3 * 3 = 36 reals = 288 bytes
      // in SP, each object defines 8 triangles = 8 * 3 * 3 = 72 reals = 288 bytes
      // each vertex vector stores the x, y, z coords of each of the triangles
      // there are total of 9 avx vectors
      avx_m256_t a[3];    // coords of vertex a of each triangle
      avx_m256_t b[3];    // coords of vertex b of each trinagle
      avx_m256_t c[3];    // coords of vertex c of each triangle
    }                                       avx_triangle_t;

  #endif    // INTEL_AVX

} // namespace


#endif /* __HIG_TYPEDEFS_HPP__ */
