/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sse_mathfun.h
 *  Created: Apr 20, 2013
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

# include <emmintrin.h>

namespace sse_mathfun {

	// TODO: clean this file ... very messy and yuckky ...

# define ALIGN16_BEG
# define ALIGN16_END __attribute__((aligned(16)))

typedef __m128 sse_m128_t;

#define _PS_CONST(Name, Val)                                            \
	static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)                                            \
	static const ALIGN16_BEG int _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PS_CONST_TYPE(Name, Type, Val)                                 \
	static const ALIGN16_BEG Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

_PS_CONST(1  , 1.0f);
_PS_CONST(0  , 0.0f);
_PS_CONST(0p5, 0.5f);

_PS_CONST_TYPE(sign_mask, int, (int)0x80000000);
_PS_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);

_PS_CONST(exp_hi,	88.3762626647949f);
_PS_CONST(exp_lo,	-88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341);
_PS_CONST(cephes_exp_C1, 0.693359375);
_PS_CONST(cephes_exp_C2, -2.12194440e-4);
_PS_CONST(cephes_exp_C12, 0.69314718056);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1);


inline sse_m128_t newexp_ps(sse_m128_t x) {
	sse_m128_t one = *(sse_m128_t*)_ps_1;
	sse_m128_t zero = *(sse_m128_t*)_ps_0;

	x = _mm_min_ps(x, *(sse_m128_t*)_ps_exp_hi);
	x = _mm_max_ps(x, *(sse_m128_t*)_ps_exp_lo);

	sse_m128_t temp_2 = _mm_mul_ps(x, *(sse_m128_t*)_ps_cephes_LOG2EF);
	temp_2 = _mm_add_ps(temp_2, *(sse_m128_t*)_ps_0p5);

	__m128i emm0 = _mm_cvttps_epi32(temp_2);
	sse_m128_t temp_1 = _mm_cvtepi32_ps(emm0);
	sse_m128_t temp_3 = _mm_sub_ps(temp_1, temp_2);
	sse_m128_t mask = _mm_cmpgt_ps(temp_3, zero);

	mask = _mm_and_ps(mask, one);
	temp_2 = _mm_sub_ps(temp_1, mask);
	emm0 = _mm_cvttps_epi32(temp_2);

	temp_1 = _mm_mul_ps(temp_2, *(sse_m128_t*)_ps_cephes_exp_C12);
	x = _mm_sub_ps(x, temp_1);

	sse_m128_t x2 = _mm_mul_ps(x, x);
	sse_m128_t x3 = _mm_mul_ps(x2, x);
	sse_m128_t x4 = _mm_mul_ps(x2, x2);
 
	temp_1 = _mm_add_ps(x, one);
	temp_2 = _mm_mul_ps(x2, *(sse_m128_t*)_ps_cephes_exp_p5);
	temp_3 = _mm_mul_ps(x3, *(sse_m128_t*)_ps_cephes_exp_p4);
	temp_1 = _mm_add_ps(temp_1, temp_2);

	temp_2 = _mm_mul_ps(x3, *(sse_m128_t*)_ps_cephes_exp_p0);

	temp_1 = _mm_add_ps(temp_1, temp_3);

	sse_m128_t temp_4 = _mm_mul_ps(x, *(sse_m128_t*)_ps_cephes_exp_p2);
	temp_3 = _mm_mul_ps(x2, *(sse_m128_t*)_ps_cephes_exp_p1);

	emm0 = _mm_add_epi32(emm0, *(__m128i*)_pi32_0x7f);

	temp_2 = _mm_add_ps(temp_2, temp_3);
	temp_3 = _mm_add_ps(temp_3, temp_4);

	emm0 = _mm_slli_epi32(emm0, 23);
	sse_m128_t pow2n = _mm_castsi128_ps(emm0);

	temp_2 = _mm_add_ps(temp_2, temp_3);
	temp_2 = _mm_mul_ps(temp_2, x4);

	sse_m128_t y = _mm_add_ps(temp_1, temp_2);

	y = _mm_mul_ps(y, pow2n);
	return y;
} // newexp_ps()


inline void newexp_ps_dual(sse_m128_t x1, sse_m128_t x2, sse_m128_t* exp1, sse_m128_t* exp2) {
	sse_m128_t one = *(sse_m128_t*)_ps_1;
	sse_m128_t zero = *(sse_m128_t*)_ps_0;

	x1 = _mm_min_ps(x1, *(sse_m128_t*)_ps_exp_hi);
	x2 = _mm_min_ps(x2, *(sse_m128_t*)_ps_exp_hi);
	x1 = _mm_max_ps(x1, *(sse_m128_t*)_ps_exp_lo);
	x2 = _mm_max_ps(x2, *(sse_m128_t*)_ps_exp_lo);

	sse_m128_t temp_21 = _mm_mul_ps(x1, *(sse_m128_t*)_ps_cephes_LOG2EF);
	sse_m128_t temp_22 = _mm_mul_ps(x2, *(sse_m128_t*)_ps_cephes_LOG2EF);
	temp_21 = _mm_add_ps(temp_21, *(sse_m128_t*)_ps_0p5);
	temp_22 = _mm_add_ps(temp_22, *(sse_m128_t*)_ps_0p5);

	__m128i emm01 = _mm_cvttps_epi32(temp_21);
	__m128i emm02 = _mm_cvttps_epi32(temp_22);
	sse_m128_t temp_11 = _mm_cvtepi32_ps(emm01);
	sse_m128_t temp_12 = _mm_cvtepi32_ps(emm02);
	sse_m128_t temp_31 = _mm_sub_ps(temp_11, temp_21);
	sse_m128_t temp_32 = _mm_sub_ps(temp_12, temp_22);
	sse_m128_t mask1 = _mm_cmpgt_ps(temp_31, zero);
	sse_m128_t mask2 = _mm_cmpgt_ps(temp_32, zero);

	mask1 = _mm_and_ps(mask1, one);
	mask2 = _mm_and_ps(mask2, one);
	temp_21 = _mm_sub_ps(temp_11, mask1);
	temp_22 = _mm_sub_ps(temp_12, mask2);
	emm01 = _mm_cvttps_epi32(temp_21);
	emm02 = _mm_cvttps_epi32(temp_22);

	temp_11 = _mm_mul_ps(temp_21, *(sse_m128_t*)_ps_cephes_exp_C12);
	temp_12 = _mm_mul_ps(temp_22, *(sse_m128_t*)_ps_cephes_exp_C12);
	x1 = _mm_sub_ps(x1, temp_11);
	x2 = _mm_sub_ps(x2, temp_12);

	sse_m128_t x21 = _mm_mul_ps(x1, x1);
	sse_m128_t x22 = _mm_mul_ps(x2, x2);
	sse_m128_t x31 = _mm_mul_ps(x21, x1);
	sse_m128_t x32 = _mm_mul_ps(x22, x2);
	sse_m128_t x41 = _mm_mul_ps(x21, x21);
	sse_m128_t x42 = _mm_mul_ps(x22, x22);
 
	temp_11 = _mm_add_ps(x1, one);
	temp_12 = _mm_add_ps(x2, one);
	temp_21 = _mm_mul_ps(x21, *(sse_m128_t*)_ps_cephes_exp_p5);
	temp_22 = _mm_mul_ps(x22, *(sse_m128_t*)_ps_cephes_exp_p5);
	temp_31 = _mm_mul_ps(x31, *(sse_m128_t*)_ps_cephes_exp_p4);
	temp_32 = _mm_mul_ps(x32, *(sse_m128_t*)_ps_cephes_exp_p4);
	temp_11 = _mm_add_ps(temp_11, temp_21);
	temp_12 = _mm_add_ps(temp_12, temp_22);

	temp_21 = _mm_mul_ps(x31, *(sse_m128_t*)_ps_cephes_exp_p0);
	temp_22 = _mm_mul_ps(x32, *(sse_m128_t*)_ps_cephes_exp_p0);

	temp_11 = _mm_add_ps(temp_11, temp_31);
	temp_12 = _mm_add_ps(temp_12, temp_32);

	sse_m128_t temp_41 = _mm_mul_ps(x1, *(sse_m128_t*)_ps_cephes_exp_p2);
	sse_m128_t temp_42 = _mm_mul_ps(x2, *(sse_m128_t*)_ps_cephes_exp_p2);
	temp_31 = _mm_mul_ps(x21, *(sse_m128_t*)_ps_cephes_exp_p1);
	temp_32 = _mm_mul_ps(x22, *(sse_m128_t*)_ps_cephes_exp_p1);

	emm01 = _mm_add_epi32(emm01, *(__m128i*)_pi32_0x7f);
	emm02 = _mm_add_epi32(emm02, *(__m128i*)_pi32_0x7f);

	temp_21 = _mm_add_ps(temp_21, temp_31);
	temp_22 = _mm_add_ps(temp_22, temp_32);
	temp_31 = _mm_add_ps(temp_31, temp_41);
	temp_32 = _mm_add_ps(temp_32, temp_42);

	emm01 = _mm_slli_epi32(emm01, 23);
	emm02 = _mm_slli_epi32(emm02, 23);
	sse_m128_t pow2n1 = _mm_castsi128_ps(emm01);
	sse_m128_t pow2n2 = _mm_castsi128_ps(emm02);

	temp_21 = _mm_add_ps(temp_21, temp_31);
	temp_22 = _mm_add_ps(temp_22, temp_32);
	temp_21 = _mm_mul_ps(temp_21, x41);
	temp_22 = _mm_mul_ps(temp_22, x42);

	sse_m128_t y1 = _mm_add_ps(temp_11, temp_21);
	sse_m128_t y2 = _mm_add_ps(temp_12, temp_22);

	*exp1 = _mm_mul_ps(y1, pow2n1);
	*exp2 = _mm_mul_ps(y2, pow2n2);
} // newexp_ps_dual()


_PS_CONST(minus_cephes_DP1, -0.78515625);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8);
_PS_CONST(minus_cephes_DP123, -0.7853981633974483096156608);
_PS_CONST(sincof_p0, -1.9515295891E-4);
_PS_CONST(sincof_p1,  8.3321608736E-3);
_PS_CONST(sincof_p2, -1.6666654611E-1);
_PS_CONST(coscof_p0,  2.443315711809948E-005);
_PS_CONST(coscof_p1, -1.388731625493765E-003);
_PS_CONST(coscof_p2,  4.166664568298827E-002);
_PS_CONST(cephes_FOPI, 1.27323954473516); // 4 / M_PI


inline sse_m128_t newsin_ps(sse_m128_t x) {
	sse_m128_t sign_bit = _mm_and_ps(x, *(sse_m128_t*)_ps_sign_mask);
	x = _mm_and_ps(x, *(sse_m128_t*)_ps_inv_sign_mask);
	
	sse_m128_t y = _mm_mul_ps(x, *(sse_m128_t*)_ps_cephes_FOPI);

	__m128i emm2 = _mm_cvttps_epi32(y);
	emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1);
	emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);
	y = _mm_cvtepi32_ps(emm2);

	__m128i emm0 = _mm_and_si128(emm2, *(__m128i*)_pi32_4);
	emm0 = _mm_slli_epi32(emm0, 29);

	emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_2);
	emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
	
	sse_m128_t swap_sign_bit = _mm_castsi128_ps(emm0);
	sse_m128_t poly_mask = _mm_castsi128_ps(emm2);
	sign_bit = _mm_xor_ps(sign_bit, swap_sign_bit);
	
	sse_m128_t temp = *(sse_m128_t*)_ps_minus_cephes_DP123;
	temp = _mm_mul_ps(y, temp);
	x = _mm_add_ps(x, temp);

	sse_m128_t x2 = _mm_mul_ps(x, x);
	sse_m128_t x3 = _mm_mul_ps(x2, x);
	sse_m128_t x4 = _mm_mul_ps(x2, x2);

	y = *(sse_m128_t*)_ps_coscof_p0;
	sse_m128_t y2 = *(sse_m128_t*)_ps_sincof_p0;
	y = _mm_mul_ps(y, x2);
	y2 = _mm_mul_ps(y2, x2);
	y = _mm_add_ps(y, *(sse_m128_t*)_ps_coscof_p1);
	y2 = _mm_add_ps(y2, *(sse_m128_t*)_ps_sincof_p1);
	y = _mm_mul_ps(y, x2);
	y2 = _mm_mul_ps(y2, x2);
	y = _mm_add_ps(y, *(sse_m128_t*)_ps_coscof_p2);
	y2 = _mm_add_ps(y2, *(sse_m128_t*)_ps_sincof_p2);
	y = _mm_mul_ps(y, x4);
	y2 = _mm_mul_ps(y2, x3);
	temp = _mm_mul_ps(x2, *(sse_m128_t*)_ps_0p5);
	temp = _mm_sub_ps(temp, *(sse_m128_t*)_ps_1);
	y = _mm_sub_ps(y, temp);
	y2 = _mm_add_ps(y2, x);

	y = _mm_andnot_ps(poly_mask, y);
	y2 = _mm_and_ps(poly_mask, y2);
	y = _mm_add_ps(y, y2);

	y = _mm_xor_ps(y, sign_bit);

	return y;
} // newsin_ps()


inline sse_m128_t newcos_ps(sse_m128_t x) {
	x = _mm_and_ps(x, *(sse_m128_t*)_ps_inv_sign_mask);
	sse_m128_t y = _mm_mul_ps(x, *(sse_m128_t*)_ps_cephes_FOPI);

	__m128i emm2 = _mm_cvttps_epi32(y);
	emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1);
	emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);
	y = _mm_cvtepi32_ps(emm2);

	emm2 = _mm_sub_epi32(emm2, *(__m128i*)_pi32_2);

	__m128i emm0 = _mm_andnot_si128(emm2, *(__m128i*)_pi32_4);
	emm0 = _mm_slli_epi32(emm0, 29);

	emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_2);
	emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
	
	sse_m128_t sign_bit = _mm_castsi128_ps(emm0);
	sse_m128_t poly_mask = _mm_castsi128_ps(emm2);
	
	sse_m128_t temp = *(sse_m128_t*)_ps_minus_cephes_DP123;
	temp = _mm_mul_ps(y, temp);
	x = _mm_add_ps(x, temp);

	sse_m128_t x2 = _mm_mul_ps(x, x);
	sse_m128_t x3 = _mm_mul_ps(x2, x);
	sse_m128_t x4 = _mm_mul_ps(x2, x2);

	y = *(sse_m128_t*)_ps_coscof_p0;
	sse_m128_t y2 = *(sse_m128_t*)_ps_sincof_p0;
	y = _mm_mul_ps(y, x2);
	y2 = _mm_mul_ps(y2, x2);
	y = _mm_add_ps(y, *(sse_m128_t*)_ps_coscof_p1);
	y2 = _mm_add_ps(y2, *(sse_m128_t*)_ps_sincof_p1);
	y = _mm_mul_ps(y, x2);
	y2 = _mm_mul_ps(y2, x2);
	y = _mm_add_ps(y, *(sse_m128_t*)_ps_coscof_p2);
	y2 = _mm_add_ps(y2, *(sse_m128_t*)_ps_sincof_p2);
	y = _mm_mul_ps(y, x4);
	y2 = _mm_mul_ps(y2, x3);
	temp = _mm_mul_ps(x2, *(sse_m128_t*)_ps_0p5);
	temp = _mm_sub_ps(temp, *(sse_m128_t*)_ps_1);
	y = _mm_sub_ps(y, temp);
	y2 = _mm_add_ps(y2, x);

	y = _mm_andnot_ps(poly_mask, y);
	y2 = _mm_and_ps(poly_mask, y2);
	y = _mm_add_ps(y, y2);

	y = _mm_xor_ps(y, sign_bit);

	return y;
} // newcos_ps()


inline void newsincos_ps(sse_m128_t x, sse_m128_t *s, sse_m128_t *c) {
	sse_m128_t sign_bit = _mm_and_ps(x, *(sse_m128_t*)_ps_sign_mask);
	x = _mm_and_ps(x, *(sse_m128_t*)_ps_inv_sign_mask);

	sse_m128_t y = _mm_mul_ps(x, *(sse_m128_t*)_ps_cephes_FOPI);

	__m128i emm2 = _mm_cvttps_epi32(y);
	emm2 = _mm_add_epi32(emm2, *(__m128i*)_pi32_1);
	emm2 = _mm_and_si128(emm2, *(__m128i*)_pi32_inv1);
	y = _mm_cvtepi32_ps(emm2);

	__m128i cos_emm2 = _mm_sub_epi32(emm2, *(__m128i*)_pi32_2);

	__m128i temp_1 = *(__m128i*)_pi32_4;
	__m128i emm0 = _mm_and_si128(emm2, temp_1);
	__m128i cos_emm0 = _mm_andnot_si128(cos_emm2, temp_1);
	emm0 = _mm_slli_epi32(emm0, 29);
	cos_emm0 = _mm_slli_epi32(cos_emm0, 29);

	temp_1 = *(__m128i*)_pi32_2;
	emm2 = _mm_and_si128(emm2, temp_1);
	cos_emm2 = _mm_and_si128(cos_emm2, temp_1);
	emm2 = _mm_cmpeq_epi32(emm2, _mm_setzero_si128());
	cos_emm2 = _mm_cmpeq_epi32(cos_emm2, _mm_setzero_si128());
	
	sse_m128_t emm0f = _mm_castsi128_ps(emm0);
	sse_m128_t emm2f = _mm_castsi128_ps(emm2);
	sse_m128_t cos_emm0f = _mm_castsi128_ps(cos_emm0);
	sse_m128_t cos_emm2f = _mm_castsi128_ps(cos_emm2);

	sign_bit = _mm_xor_ps(sign_bit, emm0f);

	sse_m128_t temp_2 = *(sse_m128_t*)_ps_minus_cephes_DP123;
	temp_2 = _mm_mul_ps(y, temp_2);
	x = _mm_add_ps(x, temp_2);

	sse_m128_t x2 = _mm_mul_ps(x, x);
	sse_m128_t x3 = _mm_mul_ps(x2, x);
	sse_m128_t x4 = _mm_mul_ps(x2, x2);

	y = *(sse_m128_t*)_ps_coscof_p0;
	sse_m128_t y2 = *(sse_m128_t*)_ps_sincof_p0;
	y = _mm_mul_ps(y, x2);
	y2 = _mm_mul_ps(y2, x2);
	y = _mm_add_ps(y, *(sse_m128_t*)_ps_coscof_p1);
	y2 = _mm_add_ps(y2, *(sse_m128_t*)_ps_sincof_p1);
	y = _mm_mul_ps(y, x2);
	y2 = _mm_mul_ps(y2, x2);
	y = _mm_add_ps(y, *(sse_m128_t*)_ps_coscof_p2);
	y2 = _mm_add_ps(y2, *(sse_m128_t*)_ps_sincof_p2);
	y = _mm_mul_ps(y, x4);
	y2 = _mm_mul_ps(y2, x3);
	temp_2 = _mm_mul_ps(x2, *(sse_m128_t*)_ps_0p5);
	y2 = _mm_add_ps(y2, x);
	temp_2 = _mm_sub_ps(temp_2, *(sse_m128_t*)_ps_1);
	y = _mm_sub_ps(y, temp_2);

	sse_m128_t cos_y = y;
	sse_m128_t cos_y2 = y2;
	y = _mm_andnot_ps(emm2f, y);
	cos_y = _mm_andnot_ps(cos_emm2f, cos_y);
	y2 = _mm_and_ps(emm2f, y2);
	cos_y2 = _mm_and_ps(cos_emm2f, cos_y2);
	y = _mm_add_ps(y, y2);
	cos_y = _mm_add_ps(cos_y, cos_y2);

	*s = _mm_xor_ps(y, sign_bit);
	*c = _mm_xor_ps(cos_y, cos_emm0f);
} // newsincos_ps()


inline void newsincos_ps_dual(sse_m128_t x1, sse_m128_t x2, sse_m128_t *s1, sse_m128_t *s2,
						sse_m128_t *c1, sse_m128_t *c2) {
	sse_m128_t tempa = *(sse_m128_t*)_ps_sign_mask;
	sse_m128_t tempb = *(sse_m128_t*)_ps_inv_sign_mask;
	sse_m128_t sign_bit1 = _mm_and_ps(x1, tempa);
	sse_m128_t sign_bit2 = _mm_and_ps(x2, tempa);
	x1 = _mm_and_ps(x1, tempb);
	x2 = _mm_and_ps(x2, tempb);

	tempa = *(sse_m128_t*)_ps_cephes_FOPI;
	sse_m128_t y1 = _mm_mul_ps(x1, tempa);
	sse_m128_t y2 = _mm_mul_ps(x2, tempa);

	__m128i emm21 = _mm_cvttps_epi32(y1);
	__m128i emm22 = _mm_cvttps_epi32(y2);
	__m128i tempia = *(__m128i*)_pi32_1;
	__m128i tempib = *(__m128i*)_pi32_inv1;
	emm21 = _mm_add_epi32(emm21, tempia);
	emm22 = _mm_add_epi32(emm22, tempia);
	emm21 = _mm_and_si128(emm21, tempib);
	emm22 = _mm_and_si128(emm22, tempib);
	y1 = _mm_cvtepi32_ps(emm21);
	y2 = _mm_cvtepi32_ps(emm22);

	tempia = *(__m128i*)_pi32_2;
	__m128i cos_emm21 = _mm_sub_epi32(emm21, tempia);
	__m128i cos_emm22 = _mm_sub_epi32(emm22, tempia);

	tempib = *(__m128i*)_pi32_4;
	__m128i emm01 = _mm_and_si128(emm21, tempib);
	__m128i emm02 = _mm_and_si128(emm22, tempib);
	__m128i cos_emm01 = _mm_andnot_si128(cos_emm21, tempib);
	__m128i cos_emm02 = _mm_andnot_si128(cos_emm22, tempib);
	emm01 = _mm_slli_epi32(emm01, 29);
	emm02 = _mm_slli_epi32(emm02, 29);
	cos_emm01 = _mm_slli_epi32(cos_emm01, 29);
	cos_emm02 = _mm_slli_epi32(cos_emm02, 29);

	tempia = *(__m128i*)_pi32_2;
	tempib = _mm_setzero_si128();
	emm21 = _mm_and_si128(emm21, tempia);
	emm22 = _mm_and_si128(emm22, tempia);
	cos_emm21 = _mm_and_si128(cos_emm21, tempia);
	cos_emm22 = _mm_and_si128(cos_emm22, tempia);
	emm21 = _mm_cmpeq_epi32(emm21, tempib);
	emm22 = _mm_cmpeq_epi32(emm22, tempib);
	cos_emm21 = _mm_cmpeq_epi32(cos_emm21, tempib);
	cos_emm22 = _mm_cmpeq_epi32(cos_emm22, tempib);
	
	sse_m128_t emm0f1 = _mm_castsi128_ps(emm01);
	sse_m128_t emm0f2 = _mm_castsi128_ps(emm02);
	sse_m128_t emm2f1 = _mm_castsi128_ps(emm21);
	sse_m128_t emm2f2 = _mm_castsi128_ps(emm22);
	sse_m128_t cos_emm0f1 = _mm_castsi128_ps(cos_emm01);
	sse_m128_t cos_emm0f2 = _mm_castsi128_ps(cos_emm02);
	sse_m128_t cos_emm2f1 = _mm_castsi128_ps(cos_emm21);
	sse_m128_t cos_emm2f2 = _mm_castsi128_ps(cos_emm22);

	sign_bit1 = _mm_xor_ps(sign_bit1, emm0f1);
	sign_bit2 = _mm_xor_ps(sign_bit2, emm0f2);

	tempa = *(sse_m128_t*)_ps_minus_cephes_DP123;
	tempb = _mm_mul_ps(y2, tempa);
	tempa = _mm_mul_ps(y1, tempa);
	x2 = _mm_add_ps(x2, tempb);
	x1 = _mm_add_ps(x1, tempa);

	sse_m128_t x21 = _mm_mul_ps(x1, x1);
	sse_m128_t x22 = _mm_mul_ps(x2, x2);
	sse_m128_t x31 = _mm_mul_ps(x21, x1);
	sse_m128_t x32 = _mm_mul_ps(x22, x2);
	sse_m128_t x41 = _mm_mul_ps(x21, x21);
	sse_m128_t x42 = _mm_mul_ps(x22, x22);

	tempa = *(sse_m128_t*)_ps_coscof_p0;
	tempb = *(sse_m128_t*)_ps_sincof_p0;

	y1 = _mm_mul_ps(x21, tempa);
	y2 = _mm_mul_ps(x22, tempa);
	sse_m128_t y21 = _mm_mul_ps(x21, tempb);
	sse_m128_t y22 = _mm_mul_ps(x22, tempb);
	tempa = *(sse_m128_t*)_ps_coscof_p1;
	tempb = *(sse_m128_t*)_ps_sincof_p1;
	y1 = _mm_add_ps(y1, tempa);
	y2 = _mm_add_ps(y2, tempa);
	y21 = _mm_add_ps(y21, tempb);
	y22 = _mm_add_ps(y22, tempb);
	y1 = _mm_mul_ps(y1, x21);
	y2 = _mm_mul_ps(y2, x22);
	y21 = _mm_mul_ps(y21, x21);
	y22 = _mm_mul_ps(y22, x22);
	tempa = *(sse_m128_t*)_ps_coscof_p2;
	tempb = *(sse_m128_t*)_ps_sincof_p2;
	y1 = _mm_add_ps(y1, tempa);
	y2 = _mm_add_ps(y2, tempa);
	y21 = _mm_add_ps(y21, tempb);
	y22 = _mm_add_ps(y22, tempb);
	y1 = _mm_mul_ps(y1, x41);
	y2 = _mm_mul_ps(y2, x42);
	y21 = _mm_mul_ps(y21, x31);
	y22 = _mm_mul_ps(y22, x32);
	tempa = *(sse_m128_t*)_ps_0p5;
	tempb = *(sse_m128_t*)_ps_1;
	sse_m128_t temp_21 = _mm_mul_ps(x21, tempa);
	sse_m128_t temp_22 = _mm_mul_ps(x22, tempa);
	y21 = _mm_add_ps(y21, x1);
	y22 = _mm_add_ps(y22, x2);
	temp_21 = _mm_sub_ps(temp_21, tempb);
	temp_22 = _mm_sub_ps(temp_22, tempb);
	y1 = _mm_sub_ps(y1, temp_21);
	y2 = _mm_sub_ps(y2, temp_22);

	sse_m128_t cos_y1 = y1;
	sse_m128_t cos_y2 = y2;
	sse_m128_t cos_y21 = y21;
	sse_m128_t cos_y22 = y22;
	y1 = _mm_andnot_ps(emm2f1, y1);
	y2 = _mm_andnot_ps(emm2f2, y2);
	cos_y1 = _mm_andnot_ps(cos_emm2f1, cos_y1);
	cos_y2 = _mm_andnot_ps(cos_emm2f2, cos_y2);
	y21 = _mm_and_ps(emm2f1, y21);
	y22 = _mm_and_ps(emm2f2, y22);
	cos_y21 = _mm_and_ps(cos_emm2f1, cos_y21);
	cos_y22 = _mm_and_ps(cos_emm2f2, cos_y22);
	y1 = _mm_add_ps(y1, y21);
	y2 = _mm_add_ps(y2, y22);
	cos_y1 = _mm_add_ps(cos_y1, cos_y21);
	cos_y2 = _mm_add_ps(cos_y2, cos_y22);

	*s1 = _mm_xor_ps(y1, sign_bit1);
	*s2 = _mm_xor_ps(y2, sign_bit2);
	*c1 = _mm_xor_ps(cos_y1, cos_emm0f1);
	*c2 = _mm_xor_ps(cos_y2, cos_emm0f2);
} // newsincos_ps_dual()

}
