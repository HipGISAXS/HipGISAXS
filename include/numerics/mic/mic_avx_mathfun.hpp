/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: mic_avx_mathfun.hpp
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

#pragma offload_attribute(push, target(mic))

static const __m512i _pi32_sign_mask = _mm512_set1_epi32(0x80000000);
static const __m512i _pi32_inv_sign_mask = _mm512_set1_epi32(~0x80000000);
static const __m512i _pi32_0 = _mm512_set1_epi32(0);
static const __m512i _pi32_1 = _mm512_set1_epi32(1);
static const __m512i _pi32_2 = _mm512_set1_epi32(2);
static const __m512i _pi32_4 = _mm512_set1_epi32(4);
static const __m512i _pi32_inv1 = _mm512_set1_epi32(~1);
static const __m512i _pi32_0x7f = _mm512_set1_epi32(0x7f);
static const __m512i _pi32_ffff = _mm512_set1_epi32(0xffffffff);

static const mic_m512_t _ps_1 = _mm512_set1_ps(1.0f);
static const mic_m512_t _ps_0point5 = _mm512_set1_ps(0.5f);
static const mic_m512_t _ps_0 = _mm512_set1_ps(0.0f);
static const mic_m512_t _ps_exp_hi = _mm512_set1_ps(88.3762626647949f);
static const mic_m512_t _ps_exp_lo = _mm512_set1_ps(-88.3762626647949f);
static const mic_m512_t _ps_cephes_LOG2EF = _mm512_set1_ps(1.44269504088896341f);
static const mic_m512_t _ps_cephes_exp_C12 = _mm512_set1_ps(0.69314718056f);
static const mic_m512_t _ps_cephes_exp_p0 = _mm512_set1_ps(1.9875691500E-4f);
static const mic_m512_t _ps_cephes_exp_p1 = _mm512_set1_ps(1.3981999507E-3f);
static const mic_m512_t _ps_cephes_exp_p2 = _mm512_set1_ps(8.3334519073E-3f);
static const mic_m512_t _ps_cephes_exp_p4 = _mm512_set1_ps(1.6666665459E-1f);
static const mic_m512_t _ps_cephes_exp_p5 = _mm512_set1_ps(5.0000001201E-1f);
static const mic_m512_t _ps_minus_cephes_DP1 = _mm512_set1_ps(-0.78515625f);
static const mic_m512_t _ps_minus_cephes_DP2 = _mm512_set1_ps(-2.4187564849853515625e-4f);
static const mic_m512_t _ps_minus_cephes_DP3 = _mm512_set1_ps(-3.77489497744594108e-8f);
static const mic_m512_t _ps_minus_cephes_DP123 = _mm512_set1_ps(-0.7853981633974483096156608f);
static const mic_m512_t _ps_sincof_p0 = _mm512_set1_ps(-1.9515295891E-4f);
static const mic_m512_t _ps_sincof_p1 = _mm512_set1_ps(8.3321608736E-3f);
static const mic_m512_t _ps_sincof_p2 = _mm512_set1_ps(-1.6666654611E-1f);
static const mic_m512_t _ps_coscof_p0 = _mm512_set1_ps(2.443315711809948E-005f);
static const mic_m512_t _ps_coscof_p1 = _mm512_set1_ps(-1.388731625493765E-003f);
static const mic_m512_t _ps_coscof_p2 = _mm512_set1_ps(4.166664568298827E-002f);
static const mic_m512_t _ps_cephes_FOPI = _mm512_set1_ps(1.27323954473516f);


// exp()
inline mic_m512_t mic_exp_ps(mic_m512_t x) {
	x = _mm512_min_ps(x, _ps_exp_hi);
	x = _mm512_max_ps(x, _ps_exp_lo);

	mic_m512_t temp_2 = _mm512_fmadd_ps(x, _ps_cephes_LOG2EF, _ps_0point5);

	mic_m512_t temp_1 = _mm512_round_ps(temp_2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	mic_m512_t temp_3 = _mm512_sub_ps(temp_1, temp_2);
	__mmask16 mask = _mm512_cmp_ps_mask(temp_3, _ps_0, _MM_CMPINT_GT);

	temp_2 = _mm512_mask_sub_ps(temp_1, mask, temp_1, _ps_1);
	__m512i emm0 = _mm512_cvtfxpnt_round_adjustps_epi32(temp_2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	x = _mm512_fnmadd_ps(temp_2, _ps_cephes_exp_C12, x);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);
 
	temp_1 = _mm512_add_ps(x, _ps_1);
	temp_1 = _mm512_fmadd_ps(x2, _ps_cephes_exp_p5, temp_1);
	temp_1 = _mm512_fmadd_ps(x3, _ps_cephes_exp_p4, temp_1);

	temp_2 = _mm512_mul_ps(x3, _ps_cephes_exp_p0);
	temp_3 = _mm512_mul_ps(x2, _ps_cephes_exp_p1);

	mic_m512_t temp_4 = _mm512_mul_ps(x, _ps_cephes_exp_p2);

	emm0 = _mm512_add_epi32(emm0, _pi32_0x7f);

	temp_2 = _mm512_add_ps(temp_2, temp_3);
	temp_3 = _mm512_add_ps(temp_3, temp_4);
	temp_2 = _mm512_add_ps(temp_2, temp_3);

	emm0 = _mm512_slli_epi32(emm0, 23);
	mic_m512_t pow2n = _mm512_castsi512_ps(emm0);

	temp_2 = _mm512_mul_ps(temp_2, x4);

	mic_m512_t y = _mm512_add_ps(temp_1, temp_2);

	y = _mm512_mul_ps(y, pow2n);
	return y;
} // newexp_ps()



// sin()
static inline mic_m512_t mic_sin_ps(mic_m512_t x) {
	__m512i sign_bit;
	sign_bit = _mm512_and_epi32(_mm512_castps_si512(x), _pi32_sign_mask);
	x = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x), _pi32_inv_sign_mask));
	
	mic_m512_t y = _mm512_mul_ps(x, _ps_cephes_FOPI);

	__m512i emm2 = _mm512_cvtfxpnt_round_adjustps_epi32(y, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	emm2 = _mm512_add_epi32(emm2, _pi32_1);
	emm2 = _mm512_and_epi32(emm2, _pi32_inv1);
	y = _mm512_cvtfxpnt_round_adjustepu32_ps(emm2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	__m512i emm0 = _mm512_and_epi32(emm2, _pi32_4);
	emm0 = _mm512_slli_epi32(emm0, 29);

	emm2 = _mm512_and_epi32(emm2, _pi32_2);
	__mmask16 mask = _mm512_cmp_epi32_mask(emm2, _pi32_0, _MM_CMPINT_EQ);
	emm2 = _mm512_mask_add_epi32(_pi32_0, mask, _pi32_ffff, _pi32_0);
	
	sign_bit = _mm512_xor_epi32(sign_bit, emm0);
	
	mic_m512_t temp = _ps_minus_cephes_DP123;
	temp = _mm512_mul_ps(y, temp);
	x = _mm512_add_ps(x, temp);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);

	y = _mm512_mul_ps(_ps_coscof_p0, x2);
	mic_m512_t y2 = _mm512_mul_ps(_ps_sincof_p0, x2);
	y = _mm512_add_ps(y, _ps_coscof_p1);
	y2 = _mm512_add_ps(y2, _ps_sincof_p1);
	y = _mm512_mul_ps(y, x2);
	y2 = _mm512_mul_ps(y2, x2);
	y = _mm512_add_ps(y, _ps_coscof_p2);
	y2 = _mm512_add_ps(y2, _ps_sincof_p2);
	y = _mm512_mul_ps(y, x4);
	y2 = _mm512_mul_ps(y2, x3);
	temp = _mm512_mul_ps(x2, _ps_0point5);
	temp = _mm512_sub_ps(temp, _ps_1);
	y = _mm512_sub_ps(y, temp);
	y2 = _mm512_add_ps(y2, x);

	y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(y)));
	y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(y2)));

	y = _mm512_add_ps(y, y2);
	y = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(y), sign_bit));

	return y;
} // sin_ps()


static inline mic_m512_t mic_cos_ps(mic_m512_t x) {
	x = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x), _pi32_inv_sign_mask));

	mic_m512_t y = _mm512_mul_ps(x, _ps_cephes_FOPI);

	__m512i emm2 = _mm512_cvtfxpnt_round_adjustps_epi32(y, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	emm2 = _mm512_add_epi32(emm2, _pi32_1);
	emm2 = _mm512_and_epi32(emm2, _pi32_inv1);
	y = _mm512_cvtfxpnt_round_adjustepu32_ps(emm2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	emm2 = _mm512_sub_epi32(emm2, _pi32_2);

	__m512i emm0 = _mm512_andnot_epi32(emm2, _pi32_4);
	emm0 = _mm512_slli_epi32(emm0, 29);

	emm2 = _mm512_and_epi32(emm2, _pi32_2);
	__mmask16 mask = _mm512_cmp_epi32_mask(emm2, _pi32_0, _MM_CMPINT_EQ);
	emm2 = _mm512_mask_add_epi32(_pi32_0, mask, _pi32_ffff, _pi32_0);
	
	mic_m512_t temp = _ps_minus_cephes_DP123;
	temp = _mm512_mul_ps(y, temp);
	x = _mm512_add_ps(x, temp);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);

	y = _mm512_mul_ps(_ps_coscof_p0, x2);
	mic_m512_t y2 = _mm512_mul_ps(_ps_sincof_p0, x2);
	y = _mm512_add_ps(y, _ps_coscof_p1);
	y2 = _mm512_add_ps(y2, _ps_sincof_p1);
	y = _mm512_mul_ps(y, x2);
	y2 = _mm512_mul_ps(y2, x2);
	y = _mm512_add_ps(y, _ps_coscof_p2);
	y2 = _mm512_add_ps(y2, _ps_sincof_p2);
	y = _mm512_mul_ps(y, x4);
	y2 = _mm512_mul_ps(y2, x3);
	temp = _mm512_mul_ps(x2, _ps_0point5);
	temp = _mm512_sub_ps(temp, _ps_1);
	y = _mm512_sub_ps(y, temp);
	y2 = _mm512_add_ps(y2, x);

	y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(y)));
	y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(y2)));

	y = _mm512_add_ps(y, y2);
	y = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(y), emm0));

	return y;
} // cos_ps()


inline void mic_sincos_ps(mic_m512_t x, mic_m512_t *s, mic_m512_t *c) {
	__m512i sign_bit = _mm512_and_epi32(_mm512_castps_si512(x), _pi32_sign_mask);
	x = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x), _pi32_inv_sign_mask));

	mic_m512_t y = _mm512_mul_ps(x, _ps_cephes_FOPI);

	__m512i emm2 = _mm512_cvtfxpnt_round_adjustps_epi32(y, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	emm2 = _mm512_add_epi32(emm2, _pi32_1);
	emm2 = _mm512_and_epi32(emm2, _pi32_inv1);
	y = _mm512_cvtfxpnt_round_adjustepu32_ps(emm2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	__m512i cos_emm2 = _mm512_sub_epi32(emm2, _pi32_2);

	__m512i emm0 = _mm512_and_epi32(emm2, _pi32_4);
	__m512i cos_emm0 = _mm512_andnot_epi32(cos_emm2, _pi32_4);
	emm0 = _mm512_slli_epi32(emm0, 29);
	cos_emm0 = _mm512_slli_epi32(cos_emm0, 29);
	sign_bit = _mm512_xor_epi32(sign_bit, emm0);

	emm2 = _mm512_and_epi32(emm2, _pi32_2);
	cos_emm2 = _mm512_and_epi32(cos_emm2, _pi32_2);
	__mmask16 mask = _mm512_cmp_epi32_mask(emm2, _pi32_0, _MM_CMPINT_EQ);
	emm2 = _mm512_mask_add_epi32(_pi32_0, mask, _pi32_ffff, _pi32_0);
	__mmask16 cos_mask = _mm512_cmp_epi32_mask(cos_emm2, _pi32_0, _MM_CMPINT_EQ);
	cos_emm2 = _mm512_mask_add_epi32(_pi32_0, cos_mask, _pi32_ffff, _pi32_0);
	
	x = _mm512_fmadd_ps(y, _ps_minus_cephes_DP123, x);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);

	y = _mm512_fmadd_ps(_ps_coscof_p0, x2, _ps_coscof_p1);
	y = _mm512_fmadd_ps(y, x2, _ps_coscof_p2);
	mic_m512_t temp_2 = _mm512_fmsub_ps(x2, _ps_0point5, _ps_1);
	y = _mm512_fmsub_ps(y, x4, temp_2);
	mic_m512_t y2 = _mm512_fmadd_ps(_ps_sincof_p0, x2, _ps_sincof_p1);
	y2 = _mm512_fmadd_ps(y2, x2, _ps_sincof_p2);
	y2 = _mm512_fmadd_ps(y2, x3, x);

	mic_m512_t cos_y = y;
	mic_m512_t cos_y2 = y2;

	y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(y)));
	cos_y = _mm512_castsi512_ps(_mm512_andnot_epi32(cos_emm2, _mm512_castps_si512(cos_y)));
	y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(y2)));
	cos_y2 = _mm512_castsi512_ps(_mm512_and_epi32(cos_emm2, _mm512_castps_si512(cos_y2)));

	y = _mm512_add_ps(y, y2);
	cos_y = _mm512_add_ps(cos_y, cos_y2);

	*s = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(y), sign_bit));
	*c = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(cos_y), cos_emm0));
} // sincos_ps()

#pragma offload_attribute(pop)
