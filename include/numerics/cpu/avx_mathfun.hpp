/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: avx_mathfun.hpp
 *  Created: May 03, 2013
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

#ifndef __AVX_MATHFUN_HPP__
#define __AVX_MATHFUN_HPP__

# include <immintrin.h>

namespace avx_mathfun {

	// TODO: clean this file ... very messy and yuckky ...

typedef __m256 avx_m256_t;
typedef __m256i avx_m256i_t;

static const avx_m256_t _ps_sign_mask = _mm256_set1_ps(0x80000000);
static const avx_m256_t _ps_inv_sign_mask = _mm256_set1_ps(~0x80000000);

static const avx_m256i_t _pi32_1 = _mm256_set1_epi32(1);
static const avx_m256i_t _pi32_inv1 = _mm256_set1_epi32(~1);
static const avx_m256i_t _pi32_2 = _mm256_set1_epi32(2);
static const avx_m256i_t _pi32_4 = _mm256_set1_epi32(4);
static const avx_m256i_t _pi32_0x7f = _mm256_set1_epi32(0x7f);

static const avx_m256_t _ps_0 = _mm256_set1_ps(0.0f);
static const avx_m256_t _ps_0p5 = _mm256_set1_ps(0.5f);
static const avx_m256_t _ps_1 = _mm256_set1_ps(1.0f);
static const avx_m256_t _ps_2 = _mm256_set1_ps(2.0f);

static const avx_m256_t _ps_exp_hi = _mm256_set1_ps(88.3762626647949f);
static const avx_m256_t _ps_exp_lo = _mm256_set1_ps(-88.3762626647949f);

static const avx_m256_t _ps_cephes_LOG2EF = _mm256_set1_ps(1.44269504088896341f);
static const avx_m256_t _ps_cephes_exp_C1 = _mm256_set1_ps(0.693359375f);
static const avx_m256_t _ps_cephes_exp_C2 = _mm256_set1_ps(-2.12194440e-4f);
static const avx_m256_t _ps_cephes_exp_C12 = _mm256_set1_ps(0.69314718056f);

static const avx_m256_t _ps_cephes_exp_p0 = _mm256_set1_ps(1.9875691500E-4f);
static const avx_m256_t _ps_cephes_exp_p1 = _mm256_set1_ps(1.3981999507E-3f);
static const avx_m256_t _ps_cephes_exp_p2 = _mm256_set1_ps(8.3334519073E-3f);
static const avx_m256_t _ps_cephes_exp_p3 = _mm256_set1_ps(4.1665795894E-2f);
static const avx_m256_t _ps_cephes_exp_p4 = _mm256_set1_ps(1.6666665459E-1f);
static const avx_m256_t _ps_cephes_exp_p5 = _mm256_set1_ps(5.0000001201E-1f);

static const avx_m256_t _ps_minus_cephes_DP1 = _mm256_set1_ps(-0.78515625f);
static const avx_m256_t _ps_minus_cephes_DP2 = _mm256_set1_ps(-2.4187564849853515625e-4f);
static const avx_m256_t _ps_minus_cephes_DP3 = _mm256_set1_ps(-3.77489497744594108e-8f);
static const avx_m256_t _ps_minus_cephes_DP123 = _mm256_set1_ps(-0.7853981633974483096156608f);
static const avx_m256_t _ps_sincof_p0 = _mm256_set1_ps(-1.9515295891E-4f);
static const avx_m256_t _ps_sincof_p1 = _mm256_set1_ps( 8.3321608736E-3f);
static const avx_m256_t _ps_sincof_p2 = _mm256_set1_ps(-1.6666654611E-1f);
static const avx_m256_t _ps_coscof_p0 = _mm256_set1_ps( 2.443315711809948E-005f);
static const avx_m256_t _ps_coscof_p1 = _mm256_set1_ps(-1.388731625493765E-003f);
static const avx_m256_t _ps_coscof_p2 = _mm256_set1_ps( 4.166664568298827E-002f);
static const avx_m256_t _ps_cephes_FOPI = _mm256_set1_ps(1.27323954473516f);


inline avx_m256_t newexp_ps(avx_m256_t x) {
	avx_m256_t one = _ps_1;
	avx_m256_t zero = _ps_0;

	x = _mm256_min_ps(x, _ps_exp_hi);
	x = _mm256_max_ps(x, _ps_exp_lo);

	avx_m256_t temp_2 = _mm256_mul_ps(x, _ps_cephes_LOG2EF);
	temp_2 = _mm256_add_ps(temp_2, _ps_0p5);

	avx_m256i_t emm0 = _mm256_cvttps_epi32(temp_2);
	avx_m256_t temp_1 = _mm256_cvtepi32_ps(emm0);
	avx_m256_t temp_3 = _mm256_sub_ps(temp_1, temp_2);
	avx_m256_t mask = _mm256_cmp_ps(temp_3, zero, _CMP_GT_OQ);

	mask = _mm256_and_ps(mask, one);
	temp_2 = _mm256_sub_ps(temp_1, mask);
	emm0 = _mm256_cvttps_epi32(temp_2);

	temp_1 = _mm256_mul_ps(temp_2, _ps_cephes_exp_C12);
	x = _mm256_sub_ps(x, temp_1);

	avx_m256_t x2 = _mm256_mul_ps(x, x);
	avx_m256_t x3 = _mm256_mul_ps(x2, x);
	avx_m256_t x4 = _mm256_mul_ps(x2, x2);
 
	temp_1 = _mm256_add_ps(x, one);
	temp_2 = _mm256_mul_ps(x2, _ps_cephes_exp_p5);
	temp_3 = _mm256_mul_ps(x3, _ps_cephes_exp_p4);
	temp_1 = _mm256_add_ps(temp_1, temp_2);

	temp_2 = _mm256_mul_ps(x3, _ps_cephes_exp_p0);

	temp_1 = _mm256_add_ps(temp_1, temp_3);

	avx_m256_t temp_4 = _mm256_mul_ps(x, _ps_cephes_exp_p2);
	temp_3 = _mm256_mul_ps(x2, _ps_cephes_exp_p1);

	emm0 = _mm256_castps_si256(_mm256_add_ps(_mm256_castsi256_ps(emm0), _mm256_castsi256_ps(_pi32_0x7f)));

	temp_2 = _mm256_add_ps(temp_2, temp_3);
	temp_3 = _mm256_add_ps(temp_3, temp_4);

	//emm0 = _mm256_slli_epi32(emm0, 23);
	// convert emm0 into two 128-bit integer vectors
	// perform shift on both vectors
	// combine both vectors into 256-bit emm0
	__m128i emm0hi = _mm256_extractf128_si256(emm0, 0);
	__m128i emm0lo = _mm256_extractf128_si256(emm0, 1);
	emm0hi = _mm_slli_epi32(emm0hi, 23);
	emm0lo = _mm_slli_epi32(emm0lo, 23);
	emm0 = _mm256_insertf128_si256(emm0, emm0hi, 0);
	emm0 = _mm256_insertf128_si256(emm0, emm0lo, 1);

	avx_m256_t pow2n = _mm256_castsi256_ps(emm0);

	temp_2 = _mm256_add_ps(temp_2, temp_3);
	temp_2 = _mm256_mul_ps(temp_2, x4);

	avx_m256_t y = _mm256_add_ps(temp_1, temp_2);

	y = _mm256_mul_ps(y, pow2n);
	return y;
} // newexp_ps()


inline void newexp_ps_dual(avx_m256_t x1, avx_m256_t x2, avx_m256_t* exp1, avx_m256_t* exp2) {
	avx_m256_t one = _ps_1;
	avx_m256_t zero = _ps_0;

	x1 = _mm256_min_ps(x1, _ps_exp_hi);
	x2 = _mm256_min_ps(x2, _ps_exp_hi);
	x1 = _mm256_max_ps(x1, _ps_exp_lo);
	x2 = _mm256_max_ps(x2, _ps_exp_lo);

	avx_m256_t temp_21 = _mm256_mul_ps(x1, _ps_cephes_LOG2EF);
	avx_m256_t temp_22 = _mm256_mul_ps(x2, _ps_cephes_LOG2EF);
	temp_21 = _mm256_add_ps(temp_21, _ps_0p5);
	temp_22 = _mm256_add_ps(temp_22, _ps_0p5);

	avx_m256i_t emm01 = _mm256_cvttps_epi32(temp_21);
	avx_m256i_t emm02 = _mm256_cvttps_epi32(temp_22);
	avx_m256_t temp_11 = _mm256_cvtepi32_ps(emm01);
	avx_m256_t temp_12 = _mm256_cvtepi32_ps(emm02);
	avx_m256_t temp_31 = _mm256_sub_ps(temp_11, temp_21);
	avx_m256_t temp_32 = _mm256_sub_ps(temp_12, temp_22);
	avx_m256_t mask1 = _mm256_cmp_ps(temp_31, zero, _CMP_GT_OQ);
	avx_m256_t mask2 = _mm256_cmp_ps(temp_32, zero, _CMP_GT_OQ);

	mask1 = _mm256_and_ps(mask1, one);
	mask2 = _mm256_and_ps(mask2, one);
	temp_21 = _mm256_sub_ps(temp_11, mask1);
	temp_22 = _mm256_sub_ps(temp_12, mask2);
	emm01 = _mm256_cvttps_epi32(temp_21);
	emm02 = _mm256_cvttps_epi32(temp_22);

	temp_11 = _mm256_mul_ps(temp_21, _ps_cephes_exp_C12);
	temp_12 = _mm256_mul_ps(temp_22, _ps_cephes_exp_C12);
	x1 = _mm256_sub_ps(x1, temp_11);
	x2 = _mm256_sub_ps(x2, temp_12);

	avx_m256_t x21 = _mm256_mul_ps(x1, x1);
	avx_m256_t x22 = _mm256_mul_ps(x2, x2);
	avx_m256_t x31 = _mm256_mul_ps(x21, x1);
	avx_m256_t x32 = _mm256_mul_ps(x22, x2);
	avx_m256_t x41 = _mm256_mul_ps(x21, x21);
	avx_m256_t x42 = _mm256_mul_ps(x22, x22);
 
	temp_11 = _mm256_add_ps(x1, one);
	temp_12 = _mm256_add_ps(x2, one);
	temp_21 = _mm256_mul_ps(x21, _ps_cephes_exp_p5);
	temp_22 = _mm256_mul_ps(x22, _ps_cephes_exp_p5);
	temp_31 = _mm256_mul_ps(x31, _ps_cephes_exp_p4);
	temp_32 = _mm256_mul_ps(x32, _ps_cephes_exp_p4);
	temp_11 = _mm256_add_ps(temp_11, temp_21);
	temp_12 = _mm256_add_ps(temp_12, temp_22);

	temp_21 = _mm256_mul_ps(x31, _ps_cephes_exp_p0);
	temp_22 = _mm256_mul_ps(x32, _ps_cephes_exp_p0);

	temp_11 = _mm256_add_ps(temp_11, temp_31);
	temp_12 = _mm256_add_ps(temp_12, temp_32);

	avx_m256_t temp_41 = _mm256_mul_ps(x1, _ps_cephes_exp_p2);
	avx_m256_t temp_42 = _mm256_mul_ps(x2, _ps_cephes_exp_p2);
	temp_31 = _mm256_mul_ps(x21, _ps_cephes_exp_p1);
	temp_32 = _mm256_mul_ps(x22, _ps_cephes_exp_p1);

	//emm01 = _mm256_add_epi32(emm01, _pi32_0x7f);
	//emm02 = _mm256_add_epi32(emm02, _pi32_0x7f);
	emm01 = _mm256_castps_si256(_mm256_add_ps(_mm256_castsi256_ps(emm01), _mm256_castsi256_ps(_pi32_0x7f)));
	emm02 = _mm256_castps_si256(_mm256_add_ps(_mm256_castsi256_ps(emm02), _mm256_castsi256_ps(_pi32_0x7f)));

	temp_21 = _mm256_add_ps(temp_21, temp_31);
	temp_22 = _mm256_add_ps(temp_22, temp_32);
	temp_31 = _mm256_add_ps(temp_31, temp_41);
	temp_32 = _mm256_add_ps(temp_32, temp_42);

	//emm01 = _mm256_slli_epi32(emm01, 23);
	__m128i emm0hi1 = _mm256_extractf128_si256(emm01, 0);
	__m128i emm0lo1 = _mm256_extractf128_si256(emm01, 1);
	emm0hi1 = _mm_slli_epi32(emm0hi1, 23);
	emm0lo1 = _mm_slli_epi32(emm0lo1, 23);
	emm01 = _mm256_insertf128_si256(emm01, emm0hi1, 0);
	emm01 = _mm256_insertf128_si256(emm01, emm0lo1, 1);

	//emm02 = _mm256_slli_epi32(emm02, 23);
	__m128i emm0hi2 = _mm256_extractf128_si256(emm02, 0);
	__m128i emm0lo2 = _mm256_extractf128_si256(emm02, 1);
	emm0hi2 = _mm_slli_epi32(emm0hi2, 23);
	emm0lo2 = _mm_slli_epi32(emm0lo2, 23);
	emm02 = _mm256_insertf128_si256(emm02, emm0hi2, 0);
	emm02 = _mm256_insertf128_si256(emm02, emm0lo2, 1);

	avx_m256_t pow2n1 = _mm256_castsi256_ps(emm01);
	avx_m256_t pow2n2 = _mm256_castsi256_ps(emm02);

	temp_21 = _mm256_add_ps(temp_21, temp_31);
	temp_22 = _mm256_add_ps(temp_22, temp_32);
	temp_21 = _mm256_mul_ps(temp_21, x41);
	temp_22 = _mm256_mul_ps(temp_22, x42);

	avx_m256_t y1 = _mm256_add_ps(temp_11, temp_21);
	avx_m256_t y2 = _mm256_add_ps(temp_12, temp_22);

	*exp1 = _mm256_mul_ps(y1, pow2n1);
	*exp2 = _mm256_mul_ps(y2, pow2n2);
} // newexp_ps_dual()



inline avx_m256_t newsin_ps(avx_m256_t x) {
	avx_m256_t sign_bit = _mm256_and_ps(x, _ps_sign_mask);
	x = _mm256_and_ps(x, _ps_inv_sign_mask);
	
	avx_m256_t y = _mm256_mul_ps(x, _ps_cephes_FOPI);

	avx_m256i_t emm2 = _mm256_cvttps_epi32(y);
	emm2 = _mm256_add_epi32(emm2, _pi32_1);
	emm2 = _mm256_and_si256(emm2, _pi32_inv1);
	y = _mm256_cvtepi32_ps(emm2);

	avx_m256i_t emm0 = _mm256_and_si256(emm2, _pi32_4);
	emm0 = _mm256_slli_epi32(emm0, 29);

	emm2 = _mm256_and_si256(emm2, _pi32_2);
	emm2 = _mm256_cmpeq_epi32(emm2, _mm256_setzero_si256());
	
	avx_m256_t swap_sign_bit = _mm256_castsi256_ps(emm0);
	avx_m256_t poly_mask = _mm256_castsi256_ps(emm2);
	sign_bit = _mm256_xor_ps(sign_bit, swap_sign_bit);
	
	avx_m256_t temp = _ps_minus_cephes_DP123;
	temp = _mm256_mul_ps(y, temp);
	x = _mm256_add_ps(x, temp);

	avx_m256_t x2 = _mm256_mul_ps(x, x);
	avx_m256_t x3 = _mm256_mul_ps(x2, x);
	avx_m256_t x4 = _mm256_mul_ps(x2, x2);

	y = _ps_coscof_p0;
	avx_m256_t y2 = _ps_sincof_p0;
	y = _mm256_mul_ps(y, x2);
	y2 = _mm256_mul_ps(y2, x2);
	y = _mm256_add_ps(y, _ps_coscof_p1);
	y2 = _mm256_add_ps(y2, _ps_sincof_p1);
	y = _mm256_mul_ps(y, x2);
	y2 = _mm256_mul_ps(y2, x2);
	y = _mm256_add_ps(y, _ps_coscof_p2);
	y2 = _mm256_add_ps(y2, _ps_sincof_p2);
	y = _mm256_mul_ps(y, x4);
	y2 = _mm256_mul_ps(y2, x3);
	temp = _mm256_mul_ps(x2, _ps_0p5);
	temp = _mm256_sub_ps(temp, _ps_1);
	y = _mm256_sub_ps(y, temp);
	y2 = _mm256_add_ps(y2, x);

	y = _mm256_andnot_ps(poly_mask, y);
	y2 = _mm256_and_ps(poly_mask, y2);
	y = _mm256_add_ps(y, y2);

	y = _mm256_xor_ps(y, sign_bit);

	return y;
} // newsin_ps()


inline avx_m256_t newcos_ps(avx_m256_t x) {
	x = _mm256_and_ps(x, _ps_inv_sign_mask);
	avx_m256_t y = _mm256_mul_ps(x, _ps_cephes_FOPI);

	avx_m256i_t emm2 = _mm256_cvttps_epi32(y);
	emm2 = _mm256_add_epi32(emm2, _pi32_1);
	emm2 = _mm256_and_si256(emm2, _pi32_inv1);
	y = _mm256_cvtepi32_ps(emm2);

	emm2 = _mm256_sub_epi32(emm2, _pi32_2);

	avx_m256i_t emm0 = _mm256_andnot_si256(emm2, _pi32_4);
	emm0 = _mm256_slli_epi32(emm0, 29);

	emm2 = _mm256_and_si256(emm2, _pi32_2);
	emm2 = _mm256_cmpeq_epi32(emm2, _mm256_setzero_si256());
	
	avx_m256_t sign_bit = _mm256_castsi256_ps(emm0);
	avx_m256_t poly_mask = _mm256_castsi256_ps(emm2);
	
	avx_m256_t temp = _ps_minus_cephes_DP123;
	temp = _mm256_mul_ps(y, temp);
	x = _mm256_add_ps(x, temp);

	avx_m256_t x2 = _mm256_mul_ps(x, x);
	avx_m256_t x3 = _mm256_mul_ps(x2, x);
	avx_m256_t x4 = _mm256_mul_ps(x2, x2);

	y = _ps_coscof_p0;
	avx_m256_t y2 = _ps_sincof_p0;
	y = _mm256_mul_ps(y, x2);
	y2 = _mm256_mul_ps(y2, x2);
	y = _mm256_add_ps(y, _ps_coscof_p1);
	y2 = _mm256_add_ps(y2, _ps_sincof_p1);
	y = _mm256_mul_ps(y, x2);
	y2 = _mm256_mul_ps(y2, x2);
	y = _mm256_add_ps(y, _ps_coscof_p2);
	y2 = _mm256_add_ps(y2, _ps_sincof_p2);
	y = _mm256_mul_ps(y, x4);
	y2 = _mm256_mul_ps(y2, x3);
	temp = _mm256_mul_ps(x2, _ps_0p5);
	temp = _mm256_sub_ps(temp, _ps_1);
	y = _mm256_sub_ps(y, temp);
	y2 = _mm256_add_ps(y2, x);

	y = _mm256_andnot_ps(poly_mask, y);
	y2 = _mm256_and_ps(poly_mask, y2);
	y = _mm256_add_ps(y, y2);

	y = _mm256_xor_ps(y, sign_bit);

	return y;
} // newcos_ps()


inline void newsincos_ps(avx_m256_t x, avx_m256_t *s, avx_m256_t *c) {
	avx_m256_t sign_bit = _mm256_and_ps(x, _ps_sign_mask);
	x = _mm256_and_ps(x, _ps_inv_sign_mask);

	avx_m256_t y = _mm256_mul_ps(x, _ps_cephes_FOPI);

	//avx_m256i_t emm2 = _mm256_cvttps_epi32(y);
	//emm2 = _mm256_add_epi32(emm2, _pi32_1);
	avx_m256i_t emm2 = _mm256_cvttps_epi32(_mm256_add_ps(y, _ps_1));

	//emm2 = _mm256_and_si256(emm2, _pi32_inv1);
	emm2 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm2), _mm256_castsi256_ps(_pi32_inv1)));

	y = _mm256_cvtepi32_ps(emm2);

	//avx_m256i_t cos_emm2 = _mm256_sub_epi32(emm2, _pi32_2);
	avx_m256i_t cos_emm2 = _mm256_cvtps_epi32(_mm256_sub_ps(_mm256_cvtepi32_ps(emm2), _ps_2));

	//avx_m256i_t emm0 = _mm256_and_si256(emm2, _pi32_4);
	avx_m256i_t emm0 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm2),
											_mm256_castsi256_ps(_pi32_4)));

	//avx_m256i_t cos_emm0 = _mm256_andnot_si256(cos_emm2, _pi32_4);
	avx_m256i_t cos_emm0 = _mm256_castps_si256(_mm256_andnot_ps(_mm256_castsi256_ps(cos_emm2),
											_mm256_castsi256_ps(_pi32_4)));

	//emm0 = _mm256_slli_epi32(emm0, 29);
	__m128i emm0hi = _mm256_extractf128_si256(emm0, 0);
	__m128i emm0lo = _mm256_extractf128_si256(emm0, 1);
	emm0hi = _mm_slli_epi32(emm0hi, 29);
	emm0lo = _mm_slli_epi32(emm0lo, 29);
	emm0 = _mm256_insertf128_si256(emm0, emm0hi, 0);
	emm0 = _mm256_insertf128_si256(emm0, emm0lo, 1);

	//cos_emm0 = _mm256_slli_epi32(cos_emm0, 29);
	__m128i cos_emm0hi = _mm256_extractf128_si256(cos_emm0, 0);
	__m128i cos_emm0lo = _mm256_extractf128_si256(cos_emm0, 1);
	cos_emm0hi = _mm_slli_epi32(cos_emm0hi, 29);
	cos_emm0lo = _mm_slli_epi32(cos_emm0lo, 29);
	cos_emm0 = _mm256_insertf128_si256(cos_emm0, cos_emm0hi, 0);
	cos_emm0 = _mm256_insertf128_si256(cos_emm0, cos_emm0lo, 1);

	//emm2 = _mm256_and_si256(emm2, _pi32_2);
	emm2 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm2),
											_mm256_castsi256_ps(_pi32_2)));

	//cos_emm2 = _mm256_and_si256(cos_emm2, _pi32_2);
	cos_emm2 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(cos_emm2),
											_mm256_castsi256_ps(_pi32_2)));

	//emm2 = _mm256_cmpeq_epi32(emm2, _mm256_setzero_si256());
	emm2 = _mm256_castps_si256(_mm256_cmp_ps(_mm256_castsi256_ps(emm2), _mm256_setzero_ps(), _CMP_EQ_UQ));

	//cos_emm2 = _mm256_cmpeq_epi32(cos_emm2, _mm256_setzero_si256());
	cos_emm2 = _mm256_castps_si256(_mm256_cmp_ps(_mm256_castsi256_ps(cos_emm2), _mm256_setzero_ps(), _CMP_EQ_UQ));

	avx_m256_t emm0f = _mm256_castsi256_ps(emm0);
	avx_m256_t emm2f = _mm256_castsi256_ps(emm2);
	avx_m256_t cos_emm0f = _mm256_castsi256_ps(cos_emm0);
	avx_m256_t cos_emm2f = _mm256_castsi256_ps(cos_emm2);

	sign_bit = _mm256_xor_ps(sign_bit, emm0f);

	avx_m256_t temp_2 = _ps_minus_cephes_DP123;
	temp_2 = _mm256_mul_ps(y, temp_2);
	x = _mm256_add_ps(x, temp_2);

	avx_m256_t x2 = _mm256_mul_ps(x, x);
	avx_m256_t x3 = _mm256_mul_ps(x2, x);
	avx_m256_t x4 = _mm256_mul_ps(x2, x2);

	y = _ps_coscof_p0;
	avx_m256_t y2 = _ps_sincof_p0;
	y = _mm256_mul_ps(y, x2);
	y2 = _mm256_mul_ps(y2, x2);
	y = _mm256_add_ps(y, _ps_coscof_p1);
	y2 = _mm256_add_ps(y2, _ps_sincof_p1);
	y = _mm256_mul_ps(y, x2);
	y2 = _mm256_mul_ps(y2, x2);
	y = _mm256_add_ps(y, _ps_coscof_p2);
	y2 = _mm256_add_ps(y2, _ps_sincof_p2);
	y = _mm256_mul_ps(y, x4);
	y2 = _mm256_mul_ps(y2, x3);
	temp_2 = _mm256_mul_ps(x2, _ps_0p5);
	y2 = _mm256_add_ps(y2, x);
	temp_2 = _mm256_sub_ps(temp_2, _ps_1);
	y = _mm256_sub_ps(y, temp_2);

	avx_m256_t cos_y = y;
	avx_m256_t cos_y2 = y2;
	y = _mm256_andnot_ps(emm2f, y);
	cos_y = _mm256_andnot_ps(cos_emm2f, cos_y);
	y2 = _mm256_and_ps(emm2f, y2);
	cos_y2 = _mm256_and_ps(cos_emm2f, cos_y2);
	y = _mm256_add_ps(y, y2);
	cos_y = _mm256_add_ps(cos_y, cos_y2);

	*s = _mm256_xor_ps(y, sign_bit);
	*c = _mm256_xor_ps(cos_y, cos_emm0f);
} // newsincos_ps()


inline void newsincos_ps_dual(avx_m256_t x1, avx_m256_t x2, avx_m256_t *s1, avx_m256_t *s2,
						avx_m256_t *c1, avx_m256_t *c2) {
	avx_m256_t tempa = _ps_sign_mask;
	avx_m256_t tempb = _ps_inv_sign_mask;
	avx_m256_t sign_bit1 = _mm256_and_ps(x1, tempa);
	avx_m256_t sign_bit2 = _mm256_and_ps(x2, tempa);
	x1 = _mm256_and_ps(x1, tempb);
	x2 = _mm256_and_ps(x2, tempb);

	tempa = _ps_cephes_FOPI;
	avx_m256_t y1 = _mm256_mul_ps(x1, tempa);
	avx_m256_t y2 = _mm256_mul_ps(x2, tempa);

	//avx_m256i_t emm21 = _mm256_cvttps_epi32(y1);
	//avx_m256i_t emm22 = _mm256_cvttps_epi32(y2);
	//emm21 = _mm256_add_epi32(emm21, _pi32_1);
	//emm22 = _mm256_add_epi32(emm22, _pi32_1);
	avx_m256i_t emm21 = _mm256_cvttps_epi32(_mm256_add_ps(y1, _ps_1));
	avx_m256i_t emm22 = _mm256_cvttps_epi32(_mm256_add_ps(y2, _ps_1));

	//emm21 = _mm256_and_si256(emm21, _pi32_inv1);
	//emm22 = _mm256_and_si256(emm22, _pi32_inv1);
	emm21 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm21), _mm256_castsi256_ps(_pi32_inv1)));
	emm22 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm22), _mm256_castsi256_ps(_pi32_inv1)));

	y1 = _mm256_cvtepi32_ps(emm21);
	y2 = _mm256_cvtepi32_ps(emm22);

	//avx_m256i_t tempia = _pi32_2;
	//avx_m256i_t cos_emm21 = _mm256_sub_epi32(emm21, tempia);
	//avx_m256i_t cos_emm22 = _mm256_sub_epi32(emm22, tempia);
	avx_m256i_t cos_emm21 = _mm256_cvtps_epi32(_mm256_sub_ps(_mm256_cvtepi32_ps(emm21), _ps_2));
	avx_m256i_t cos_emm22 = _mm256_cvtps_epi32(_mm256_sub_ps(_mm256_cvtepi32_ps(emm22), _ps_2));

	//avx_m256i_t tempib = _pi32_4;
	//avx_m256i_t emm01 = _mm256_and_si256(emm21, tempib);
	//avx_m256i_t emm02 = _mm256_and_si256(emm22, tempib);
	avx_m256i_t emm01 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm21),
											_mm256_castsi256_ps(_pi32_4)));
	avx_m256i_t emm02 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm22),
											_mm256_castsi256_ps(_pi32_4)));

	//avx_m256i_t cos_emm01 = _mm256_andnot_si256(cos_emm21, tempib);
	//avx_m256i_t cos_emm02 = _mm256_andnot_si256(cos_emm22, tempib);
	avx_m256i_t cos_emm01 = _mm256_castps_si256(_mm256_andnot_ps(_mm256_castsi256_ps(cos_emm21),
											_mm256_castsi256_ps(_pi32_4)));
	avx_m256i_t cos_emm02 = _mm256_castps_si256(_mm256_andnot_ps(_mm256_castsi256_ps(cos_emm22),
											_mm256_castsi256_ps(_pi32_4)));

	//emm01 = _mm256_slli_epi32(emm01, 29);
	__m128i emm0hi1 = _mm256_extractf128_si256(emm01, 0);
	__m128i emm0lo1 = _mm256_extractf128_si256(emm01, 1);
	emm0hi1 = _mm_slli_epi32(emm0hi1, 29);
	emm0lo1 = _mm_slli_epi32(emm0lo1, 29);
	emm01 = _mm256_insertf128_si256(emm01, emm0hi1, 0);
	emm01 = _mm256_insertf128_si256(emm01, emm0lo1, 1);

	//emm02 = _mm256_slli_epi32(emm02, 29);
	__m128i emm0hi2 = _mm256_extractf128_si256(emm02, 0);
	__m128i emm0lo2 = _mm256_extractf128_si256(emm02, 1);
	emm0hi2 = _mm_slli_epi32(emm0hi2, 29);
	emm0lo2 = _mm_slli_epi32(emm0lo2, 29);
	emm02 = _mm256_insertf128_si256(emm02, emm0hi1, 0);
	emm02 = _mm256_insertf128_si256(emm02, emm0lo1, 1);

	//cos_emm01 = _mm256_slli_epi32(cos_emm01, 29);
	__m128i cos_emm0hi1 = _mm256_extractf128_si256(cos_emm01, 0);
	__m128i cos_emm0lo1 = _mm256_extractf128_si256(cos_emm01, 1);
	cos_emm0hi1 = _mm_slli_epi32(cos_emm0hi1, 29);
	cos_emm0lo1 = _mm_slli_epi32(cos_emm0lo1, 29);
	cos_emm01 = _mm256_insertf128_si256(cos_emm01, cos_emm0hi1, 0);
	cos_emm01 = _mm256_insertf128_si256(cos_emm01, cos_emm0lo1, 1);

	//cos_emm02 = _mm256_slli_epi32(cos_emm02, 29);
	__m128i cos_emm0hi2 = _mm256_extractf128_si256(cos_emm02, 0);
	__m128i cos_emm0lo2 = _mm256_extractf128_si256(cos_emm02, 1);
	cos_emm0hi2 = _mm_slli_epi32(cos_emm0hi2, 29);
	cos_emm0lo2 = _mm_slli_epi32(cos_emm0lo2, 29);
	cos_emm02 = _mm256_insertf128_si256(cos_emm02, cos_emm0hi2, 0);
	cos_emm02 = _mm256_insertf128_si256(cos_emm02, cos_emm0lo2, 1);

	//tempia = _pi32_2;
	//tempib = _mm256_setzero_si256();
	//emm21 = _mm256_and_si256(emm21, tempia);
	//emm22 = _mm256_and_si256(emm22, tempia);
	emm21 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm21),
											_mm256_castsi256_ps(_pi32_2)));
	emm22 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(emm22),
											_mm256_castsi256_ps(_pi32_2)));

	//cos_emm21 = _mm256_and_si256(cos_emm21, tempia);
	//cos_emm22 = _mm256_and_si256(cos_emm22, tempia);
	cos_emm21 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(cos_emm21),
											_mm256_castsi256_ps(_pi32_2)));
	cos_emm22 = _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(cos_emm22),
											_mm256_castsi256_ps(_pi32_2)));

	//emm21 = _mm256_cmpeq_epi32(emm21, tempib);
	//emm22 = _mm256_cmpeq_epi32(emm22, tempib);
	emm21 = _mm256_castps_si256(_mm256_cmp_ps(_mm256_castsi256_ps(emm21), _mm256_setzero_ps(), _CMP_EQ_UQ));
	emm22 = _mm256_castps_si256(_mm256_cmp_ps(_mm256_castsi256_ps(emm22), _mm256_setzero_ps(), _CMP_EQ_UQ));

	//cos_emm21 = _mm256_cmpeq_epi32(cos_emm21, tempib);
	//cos_emm22 = _mm256_cmpeq_epi32(cos_emm22, tempib);
	cos_emm21 = _mm256_castps_si256(_mm256_cmp_ps(_mm256_castsi256_ps(cos_emm21), _mm256_setzero_ps(), _CMP_EQ_UQ));
	cos_emm22 = _mm256_castps_si256(_mm256_cmp_ps(_mm256_castsi256_ps(cos_emm22), _mm256_setzero_ps(), _CMP_EQ_UQ));
	
	avx_m256_t emm0f1 = _mm256_castsi256_ps(emm01);
	avx_m256_t emm0f2 = _mm256_castsi256_ps(emm02);
	avx_m256_t emm2f1 = _mm256_castsi256_ps(emm21);
	avx_m256_t emm2f2 = _mm256_castsi256_ps(emm22);
	avx_m256_t cos_emm0f1 = _mm256_castsi256_ps(cos_emm01);
	avx_m256_t cos_emm0f2 = _mm256_castsi256_ps(cos_emm02);
	avx_m256_t cos_emm2f1 = _mm256_castsi256_ps(cos_emm21);
	avx_m256_t cos_emm2f2 = _mm256_castsi256_ps(cos_emm22);

	sign_bit1 = _mm256_xor_ps(sign_bit1, emm0f1);
	sign_bit2 = _mm256_xor_ps(sign_bit2, emm0f2);

	tempa = _ps_minus_cephes_DP123;
	tempb = _mm256_mul_ps(y2, tempa);
	tempa = _mm256_mul_ps(y1, tempa);
	x2 = _mm256_add_ps(x2, tempb);
	x1 = _mm256_add_ps(x1, tempa);

	avx_m256_t x21 = _mm256_mul_ps(x1, x1);
	avx_m256_t x22 = _mm256_mul_ps(x2, x2);
	avx_m256_t x31 = _mm256_mul_ps(x21, x1);
	avx_m256_t x32 = _mm256_mul_ps(x22, x2);
	avx_m256_t x41 = _mm256_mul_ps(x21, x21);
	avx_m256_t x42 = _mm256_mul_ps(x22, x22);

	tempa = _ps_coscof_p0;
	tempb = _ps_sincof_p0;

	y1 = _mm256_mul_ps(x21, tempa);
	y2 = _mm256_mul_ps(x22, tempa);
	avx_m256_t y21 = _mm256_mul_ps(x21, tempb);
	avx_m256_t y22 = _mm256_mul_ps(x22, tempb);
	tempa = _ps_coscof_p1;
	tempb = _ps_sincof_p1;
	y1 = _mm256_add_ps(y1, tempa);
	y2 = _mm256_add_ps(y2, tempa);
	y21 = _mm256_add_ps(y21, tempb);
	y22 = _mm256_add_ps(y22, tempb);
	y1 = _mm256_mul_ps(y1, x21);
	y2 = _mm256_mul_ps(y2, x22);
	y21 = _mm256_mul_ps(y21, x21);
	y22 = _mm256_mul_ps(y22, x22);
	tempa = _ps_coscof_p2;
	tempb = _ps_sincof_p2;
	y1 = _mm256_add_ps(y1, tempa);
	y2 = _mm256_add_ps(y2, tempa);
	y21 = _mm256_add_ps(y21, tempb);
	y22 = _mm256_add_ps(y22, tempb);
	y1 = _mm256_mul_ps(y1, x41);
	y2 = _mm256_mul_ps(y2, x42);
	y21 = _mm256_mul_ps(y21, x31);
	y22 = _mm256_mul_ps(y22, x32);
	tempa = _ps_0p5;
	tempb = _ps_1;
	avx_m256_t temp_21 = _mm256_mul_ps(x21, tempa);
	avx_m256_t temp_22 = _mm256_mul_ps(x22, tempa);
	y21 = _mm256_add_ps(y21, x1);
	y22 = _mm256_add_ps(y22, x2);
	temp_21 = _mm256_sub_ps(temp_21, tempb);
	temp_22 = _mm256_sub_ps(temp_22, tempb);
	y1 = _mm256_sub_ps(y1, temp_21);
	y2 = _mm256_sub_ps(y2, temp_22);

	avx_m256_t cos_y1 = y1;
	avx_m256_t cos_y2 = y2;
	avx_m256_t cos_y21 = y21;
	avx_m256_t cos_y22 = y22;
	y1 = _mm256_andnot_ps(emm2f1, y1);
	y2 = _mm256_andnot_ps(emm2f2, y2);
	cos_y1 = _mm256_andnot_ps(cos_emm2f1, cos_y1);
	cos_y2 = _mm256_andnot_ps(cos_emm2f2, cos_y2);
	y21 = _mm256_and_ps(emm2f1, y21);
	y22 = _mm256_and_ps(emm2f2, y22);
	cos_y21 = _mm256_and_ps(cos_emm2f1, cos_y21);
	cos_y22 = _mm256_and_ps(cos_emm2f2, cos_y22);
	y1 = _mm256_add_ps(y1, y21);
	y2 = _mm256_add_ps(y2, y22);
	cos_y1 = _mm256_add_ps(cos_y1, cos_y21);
	cos_y2 = _mm256_add_ps(cos_y2, cos_y22);

	*s1 = _mm256_xor_ps(y1, sign_bit1);
	*s2 = _mm256_xor_ps(y2, sign_bit2);
	*c1 = _mm256_xor_ps(cos_y1, cos_emm0f1);
	*c2 = _mm256_xor_ps(cos_y2, cos_emm0f2);
} // newsincos_ps_dual()

} // namespace avx_mathfun

#endif // __AVX_MATHFUN_HPP__
