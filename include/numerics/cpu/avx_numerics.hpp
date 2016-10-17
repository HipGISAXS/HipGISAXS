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

#ifdef INTEL_SB_AVX

#ifndef __AVX_NUMERICS_HPP__
#define __AVX_NUMERICS_HPP__

#include <immintrin.h>

#include <common/typedefs.hpp>
#include <numerics/cpu/avx_mathfun.hpp>

namespace hig {

	// ////////////////////////////////////////////////////
	// intrinsic wrapper naming (only for floating-point)
	// _mm256_xxx_abc  => a = r|c, b = p|s, c = s|d
	// _mm256_xxx_abcd => a = r|c, b = r|c, c = p|s, d = s|d
	// r = real,             c = complex,
	// p = packed (vector),  s = scalar,
	// s = single-precision, d = double-precision
	// ////////////////////////////////////////////////////

	/**
	 * load store
	 */

	static inline avx_m256_t avx_load_rps(real_t* p) {
		return _mm256_load_ps(p);
	} // avx_load_rps()

	static inline void avx_addstore_css(complex_t* p, avx_m256c_t v) {
		real_t real = _mm256_cvtss_f32(v.xvec);
		real_t imag = _mm256_cvtss_f32(v.yvec);
		(*p) += complex_t(real, imag);
	} // avx_store_css()


	/**
	 * set
	 */
	
	static inline avx_m256c_t avx_setzero_cps() {
		avx_m256c_t vec;
		vec.xvec = _mm256_setzero_ps();
		vec.yvec = _mm256_setzero_ps();
		return vec;
	} // avx_setzero_cps()

	static inline avx_m256_t avx_set1_rps(real_t a) {
		return _mm256_set1_ps(a);
	} // avx_load1_rps()

	static inline avx_m256c_t avx_set1_cps(complex_t a) {
		avx_m256c_t vec;
		vec.xvec = _mm256_set1_ps(a.real());
		vec.yvec = _mm256_set1_ps(a.imag());
		return vec;
	} // avx_set1_cps()


	/**
	 * addition
	 */

	static inline avx_m256_t avx_add_rrps(avx_m256_t a, avx_m256_t b) {
		return _mm256_add_ps(a, b);
	} // avx_add_rrps()

	static inline avx_m256c_t avx_add_rcps(avx_m256_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_add_ps(a, b.xvec);
		vec.yvec = b.yvec;
		return vec;
	} // avx_add_rcps()

	static inline avx_m256c_t avx_add_ccps(avx_m256c_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_add_ps(a.xvec, b.xvec);
		vec.yvec = _mm256_add_ps(a.yvec, b.yvec);
		return vec;
	} // avx_add_ccps()

	static inline avx_m256c_t avx_hadd_ccps(avx_m256c_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_hadd_ps(a.xvec, b.xvec);
		vec.yvec = _mm256_hadd_ps(a.yvec, b.yvec);
		return vec;
	} // avx_hadd_ccps()


	/**
	 * multiplication
	 */

	static inline avx_m256_t avx_mul_rrps(avx_m256_t a, avx_m256_t b) {
		return _mm256_mul_ps(a, b);
	} // avx_mul_rrps()

	static inline avx_m256c_t avx_mul_crps(avx_m256c_t a, avx_m256_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_mul_ps(a.xvec, b);
		vec.yvec = _mm256_mul_ps(a.yvec, b);
		return vec;
	} // avx_mul_cpps()

	static inline avx_m256c_t avx_mul_ccps(avx_m256c_t a, avx_m256c_t b) {
		avx_m256c_t vec;
		vec.xvec = _mm256_sub_ps(_mm256_mul_ps(a.xvec, b.xvec), _mm256_mul_ps(a.yvec, b.yvec));
		vec.yvec = _mm256_add_ps(_mm256_mul_ps(a.xvec, b.yvec), _mm256_mul_ps(a.yvec, b.xvec));
		return vec;
	} // avx_mul_ccps()


	/**
	 * reciprocal
	 */

	static inline avx_m256c_t avx_rcp_cps(avx_m256c_t a) {
		avx_m256_t temp1 = _mm256_add_ps(_mm256_mul_ps(a.xvec, a.xvec), _mm256_mul_ps(a.yvec, a.yvec));
		avx_m256_t temp2 = _mm256_rcp_ps(temp1);
		avx_m256c_t vec;
		avx_m256_t neg_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000));
		vec.xvec = _mm256_mul_ps(a.xvec, temp2);
		vec.yvec = _mm256_xor_ps(_mm256_mul_ps(a.yvec, temp2), neg_mask);
		return vec;
	} // avx_rcp_cps()


	/**
	 * exponential
	 */

	static inline avx_m256_t avx_exp_rps(avx_m256_t a) {
		return avx_mathfun::newexp_ps(a);
	} // avx_exp_rps()

	static inline void avx_exp_rps_dual(avx_m256_t a1, avx_m256_t a2, avx_m256_t* exp1, avx_m256_t* exp2) {
		avx_mathfun::newexp_ps_dual(a1, a2, exp1, exp2);
	} // avx_exp_rps()


	/**
	 * trigonometry
	 */

	static inline avx_m256_t avx_sin_rps(avx_m256_t a) {
		//return amd_vrs4_sinf(a);
		return avx_mathfun::newsin_ps(a);
	} // avx_exp_rps()

	static inline avx_m256_t avx_cos_rps(avx_m256_t a) {
		//return amd_vrs4_cosf(a);
		return avx_mathfun::newcos_ps(a);
	} // avx_exp_rps()

	static inline void avx_sincos_rps(avx_m256_t a, avx_m256_t* s, avx_m256_t* c) {
		avx_mathfun::newsincos_ps(a, s, c);
		//*s = avx_mathfun::newsin_ps(a);
		//*c = avx_mathfun::newcos_ps(a);
	} // avx_exp_rps()

	static inline void avx_sincos_rps_dual(avx_m256_t a1, avx_m256_t a2, avx_m256_t* s1, avx_m256_t* s2,
											avx_m256_t* c1, avx_m256_t* c2) {
		avx_mathfun::newsincos_ps_dual(a1, a2, s1, s2, c1, c2);
	} // avx_exp_rps()

} // namespace hig

#endif // __AVX_NUMERICS_HPP__

#endif // INTEL_SB_AVX
