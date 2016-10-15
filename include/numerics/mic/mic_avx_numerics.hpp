/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: mic_avx_numerics.hpp
 *  Created: Apr 23, 2013
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

#ifdef USE_MIC
#ifdef __MIC__

#ifndef __MIC_AVX_NUMERICS_HPP__
#define __MIC_AVX_NUMERICS_HPP__

#include <immintrin.h>
//#include <amdlibm.h>		// for exp, sin, cos
/*#ifdef __cplusplus
extern "C" {
#include <mkl.h>
}
#endif*/

#include <common/typedefs.hpp>

namespace hig {

	// ////////////////////////////////////////////////////
	// intrinsic wrapper naming (only for floating-point)
	// mic_xxx_abc  => a = r|c, b = p|s, c = s|d
	// mic_xxx_abcd => a = r|c, b = r|c, c = p|s, d = s|d
	// r = real,             c = complex,
	// p = packed (vector),  s = scalar,
	// s = single-precision, d = double-precision
	// ////////////////////////////////////////////////////

#include "mic_avx_mathfun.hpp"


#pragma offload_attribute(push, target(mic))

	/**
	 * load store
	 */

	static inline mic_m512_t mic_load_rps(real_t* p) {
		return _mm512_load_ps(p);
	} // mic_load_rps()


	/**
	 * set
	 */
	
	static inline mic_m512c_t mic_setzero_cps() {
		mic_m512c_t vec;
		vec.xvec = _mm512_setzero_ps();
		vec.yvec = _mm512_setzero_ps();
		return vec;
	} // mic_setzero_cps()

	static inline mic_m512_t mic_set1_rps(real_t a) {
		return _mm512_set1_ps(a);
	} // mic_load1_rps()

	static inline mic_m512c_t mic_set1_cps(scomplex_t a) {
		mic_m512c_t vec;
		vec.xvec = _mm512_set1_ps(a.x);
		vec.yvec = _mm512_set1_ps(a.y);
		return vec;
	} // mic_set1_cps()


	/**
	 * addition
	 */

	static inline mic_m512_t mic_add_rrps(mic_m512_t a, mic_m512_t b) {
		return _mm512_add_ps(a, b);
	} // mic_add_rrps()

	static inline mic_m512c_t mic_add_rcps(mic_m512_t a, mic_m512c_t b) {
		mic_m512c_t vec;
		vec.xvec = _mm512_add_ps(a, b.xvec);
		vec.yvec = b.yvec;
		return vec;
	} // mic_add_rcps()

	static inline mic_m512c_t mic_add_ccps(mic_m512c_t a, mic_m512c_t b) {
		mic_m512c_t vec;
		vec.xvec = _mm512_add_ps(a.xvec, b.xvec);
		vec.yvec = _mm512_add_ps(a.yvec, b.yvec);
		return vec;
	} // mic_add_ccps()


	/**
	 * multiplication
	 */

	static inline mic_m512_t mic_mul_rrps(mic_m512_t a, mic_m512_t b) {
		return _mm512_mul_ps(a, b);
	} // mic_mul_rrps()

	static inline mic_m512c_t mic_mul_crps(mic_m512c_t a, mic_m512_t b) {
		mic_m512c_t vec;
		vec.xvec = _mm512_mul_ps(a.xvec, b);
		vec.yvec = _mm512_mul_ps(a.yvec, b);
		return vec;
	} // mic_mul_cpps()

	static inline mic_m512c_t mic_mul_ccps(mic_m512c_t a, mic_m512c_t b) {
		mic_m512c_t vec;
		vec.xvec = _mm512_sub_ps(_mm512_mul_ps(a.xvec, b.xvec), _mm512_mul_ps(a.yvec, b.yvec));
		vec.yvec = _mm512_add_ps(_mm512_mul_ps(a.xvec, b.yvec), _mm512_mul_ps(a.yvec, b.xvec));
		return vec;
	} // mic_mul_ccps()


	/**
	 * reciprocal
	 */

	static inline mic_m512c_t mic_rcp_cps(mic_m512c_t a) {
		mic_m512_t temp1 = _mm512_add_ps(_mm512_mul_ps(a.xvec, a.xvec), _mm512_mul_ps(a.yvec, a.yvec));
		mic_m512_t temp2 = _mm512_rcp23_ps(temp1);
		mic_m512c_t vec;
		//__m512 neg_mask = _mm512_castsi512_ps(_mm512_set1_epi32(0x80000000));
		vec.xvec = _mm512_mul_ps(a.xvec, temp2);
		//vec.yvec = _mm512_xor_ps(_mm512_mul_ps(a.yvec, temp2), neg_mask);
		mic_m512_t zero = _mm512_setzero_ps();
		vec.yvec = _mm512_sub_ps(zero, _mm512_mul_ps(a.yvec, temp2));
		return vec;
	} // mic_rcp_cps()


	/**
	 * reduction
	 */

	static inline scomplex_t mic_reduce_add_cps(mic_m512c_t v) {
		scomplex_t temp;
		temp.x = _mm512_reduce_add_ps(v.xvec);
		temp.y = _mm512_reduce_add_ps(v.yvec);
		return temp;
	} // mic_reduce_add_cps()


	/**
	 * exponential
	 */

	static inline mic_m512_t mic_exp_rps(mic_m512_t a) {
		// svml:
		//return _mm512_exp_ps(a);
		// custom:
		return mic_exp_ps(a);
	} // mic_exp_rps()


	/**
	 * trigonometry
	 */

	static inline mic_m512_t mic_sin_rps(mic_m512_t a) {
		// svml:
		//return _mm512_sin_ps(a);
		// custom:
		return mic_sin_ps(a);
	} // mic_sin_rps()

	static inline mic_m512_t mic_cos_rps(mic_m512_t a) {
		// svml:
		//return _mm512_cos_ps(a);
		// custom:
		return mic_cos_ps(a);
	} // mic_cos_rps()

	static inline void mic_sincos_rps(mic_m512_t a, mic_m512_t* s, mic_m512_t* c) {
		// svml:
		//*s = _mm512_sincos_ps(c, a);
		//*s = _mm512_sin_ps(a);
		//*c = _mm512_cos_ps(a);
		// custom:
		return mic_sincos_ps(a, s, c);
	} // mic_sincos_rps()


#pragma offload_attribute(pop)

} // namespace hig

#endif // __MIC_AVX_NUMERICS_HPP__

#endif // __MIC__
#endif // USE_MIC
