/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sse3_numerics.hpp
 *  Created: Apr 19, 2013
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

#ifdef __SSE3__

#ifndef __SSE3_NUMERICS_HPP__
#define __SSE3_NUMERICS_HPP__

#include <pmmintrin.h>
//#include <amdlibm.h>		// for exp, sin, cos
/*#ifdef __cplusplus
extern "C" {
#include <mkl.h>
}
#endif*/

#include <common/typedefs.hpp>
#include <numerics/cpu/sse_mathfun.h>

namespace hig {

	// ////////////////////////////////////////////////////
	// intrinsic wrapper naming (only for floating-point)
	// _mm_xxx_abc  => a = r|c, b = p|s, c = s|d
	// _mm_xxx_abcd => a = r|c, b = r|c, c = p|s, d = s|d
	// r = real,             c = complex,
	// p = packed (vector),  s = scalar,
	// s = single-precision, d = double-precision
	// ////////////////////////////////////////////////////

	/**
	 * load store
	 */

	static inline sse_m128_t sse_load_rps(real_t* p) {
		return _mm_load_ps(p);
	} // sse_load_rps()

	static inline void sse_addstore_css(complex_t* p, sse_m128c_t v) {
		real_t real = _mm_cvtss_f32(v.xvec);
		real_t imag = _mm_cvtss_f32(v.yvec);
		//real_t real, imag;
		//_mm_store_ss(&real, v.xvec);
		//_mm_store_ss(&imag, v.yvec);
		(*p) += complex_t(real, imag);
	} // sse_store_css()


	/**
	 * set
	 */
	
	static inline sse_m128c_t sse_setzero_cps() {
		sse_m128c_t vec;
		vec.xvec = _mm_setzero_ps();
		vec.yvec = _mm_setzero_ps();
		return vec;
	} // sse_setzero_cps()

	static inline sse_m128_t sse_set1_rps(real_t a) {
		return _mm_set1_ps(a);
	} // sse_load1_rps()

	static inline sse_m128c_t sse_set1_cps(complex_t a) {
		sse_m128c_t vec;
		vec.xvec = _mm_set1_ps(a.real());
		vec.yvec = _mm_set1_ps(a.imag());
		return vec;
	} // sse_set1_cps()


	/**
	 * addition
	 */

	static inline sse_m128_t sse_add_rrps(sse_m128_t a, sse_m128_t b) {
		return _mm_add_ps(a, b);
	} // sse_add_rrps()

	static inline sse_m128c_t sse_add_rcps(sse_m128_t a, sse_m128c_t b) {
		sse_m128c_t vec;
		vec.xvec = _mm_add_ps(a, b.xvec);
		vec.yvec = b.yvec;
		return vec;
	} // sse_add_rcps()

	static inline sse_m128c_t sse_add_ccps(sse_m128c_t a, sse_m128c_t b) {
		sse_m128c_t vec;
		vec.xvec = _mm_add_ps(a.xvec, b.xvec);
		vec.yvec = _mm_add_ps(a.yvec, b.yvec);
		return vec;
	} // sse_add_ccps()

	static inline sse_m128c_t sse_hadd_ccps(sse_m128c_t a, sse_m128c_t b) {
		sse_m128c_t vec;
		vec.xvec = _mm_hadd_ps(a.xvec, b.xvec);
		vec.yvec = _mm_hadd_ps(a.yvec, b.yvec);
		return vec;
	} // sse_hadd_ccps()


	/**
	 * multiplication
	 */

	static inline sse_m128_t sse_mul_rrps(sse_m128_t a, sse_m128_t b) {
		return _mm_mul_ps(a, b);
	} // sse_mul_rrps()

	static inline sse_m128c_t sse_mul_crps(sse_m128c_t a, sse_m128_t b) {
		sse_m128c_t vec;
		vec.xvec = _mm_mul_ps(a.xvec, b);
		vec.yvec = _mm_mul_ps(a.yvec, b);
		return vec;
	} // sse_mul_cpps()

	static inline sse_m128c_t sse_mul_ccps(sse_m128c_t a, sse_m128c_t b) {
		sse_m128c_t vec;
		vec.xvec = _mm_sub_ps(_mm_mul_ps(a.xvec, b.xvec), _mm_mul_ps(a.yvec, b.yvec));
		vec.yvec = _mm_add_ps(_mm_mul_ps(a.xvec, b.yvec), _mm_mul_ps(a.yvec, b.xvec));
		return vec;
	} // sse_mul_ccps()


	/**
	 * reciprocal
	 */

	static inline sse_m128c_t sse_rcp_cps(sse_m128c_t a) {
		sse_m128_t temp1 = _mm_add_ps(_mm_mul_ps(a.xvec, a.xvec), _mm_mul_ps(a.yvec, a.yvec));
		sse_m128_t temp2 = _mm_rcp_ps(temp1);
		sse_m128c_t vec;
		__m128 neg_mask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
		vec.xvec = _mm_mul_ps(a.xvec, temp2);
		vec.yvec = _mm_xor_ps(_mm_mul_ps(a.yvec, temp2), neg_mask);
		return vec;
	} // sse_rcp_cps()


	/**
	 * exponential
	 */

	static inline sse_m128_t sse_exp_rps(sse_m128_t a) {
		//return amd_vrs4_expf(a);
		//return sse_mathfun::exp_ps(a);
		return sse_mathfun::newexp_ps(a);
		/*float* mivec = (float*) _mm_malloc(4 * sizeof(float), 16);
		float* movec = (float*) _mm_malloc(4 * sizeof(float), 16);
		_mm_store_ps(mivec, a);
		vsExp(4, mivec, movec);
		sse_m128_t ret_vec = _mm_load_ps(movec);
		_mm_free(movec);
		_mm_free(mivec);
		return ret_vec;*/
	} // sse_exp_rps()

	static inline void sse_exp_rps_dual(sse_m128_t a1, sse_m128_t a2, sse_m128_t* exp1, sse_m128_t* exp2) {
		sse_mathfun::newexp_ps_dual(a1, a2, exp1, exp2);
	} // sse_exp_rps()


	/**
	 * trigonometry
	 */

	static inline sse_m128_t sse_sin_rps(sse_m128_t a) {
		//return amd_vrs4_sinf(a);
		return sse_mathfun::newsin_ps(a);
	} // sse_exp_rps()

	static inline sse_m128_t sse_cos_rps(sse_m128_t a) {
		//return amd_vrs4_cosf(a);
		return sse_mathfun::newcos_ps(a);
	} // sse_exp_rps()

	static inline void sse_sincos_rps(sse_m128_t a, sse_m128_t* s, sse_m128_t* c) {
		sse_mathfun::newsincos_ps(a, s, c);
	} // sse_exp_rps()

	static inline void sse_sincos_rps_dual(sse_m128_t a1, sse_m128_t a2, sse_m128_t* s1, sse_m128_t* s2,
											sse_m128_t* c1, sse_m128_t* c2) {
		sse_mathfun::newsincos_ps_dual(a1, a2, s1, s2, c1, c2);
	} // sse_exp_rps()

} // namespace hig

#endif // __SSE3_NUMERICS_HPP__

#endif // __SSE3__
