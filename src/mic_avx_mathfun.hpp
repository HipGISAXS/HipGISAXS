# define ALIGN64_BEG
# define ALIGN64_END __attribute__((aligned(64)))

#define _PS_CONST(Name, Val)	\
	static const ALIGN64_BEG float _ps_##Name[16] ALIGN64_END = { Val, Val, Val, Val, \
 																Val, Val, Val, Val, \
 																Val, Val, Val, Val, \
 																Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)	\
	static const ALIGN64_BEG int _pi32_##Name[16] ALIGN64_END = { Val, Val, Val, Val, \
 																Val, Val, Val, Val, \
 																Val, Val, Val, Val, \
 																Val, Val, Val, Val }
#define _PS_CONST_TYPE(Name, Type, Val)	\
	static const ALIGN64_BEG Type _ps_##Name[16] ALIGN64_END  = { Val, Val, Val, Val, \
 																Val, Val, Val, Val, \
 																Val, Val, Val, Val, \
 																Val, Val, Val, Val }

#pragma offload_attribute(push, target(mic))

_PS_CONST(0	, 0.0f);
_PS_CONST(1	, 1.0f);
_PS_CONST(0p5, 0.5f);

_PS_CONST_TYPE(min_norm_pos, int, 0x00800000);
_PS_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS_CONST_TYPE(sign_mask_ps, int, (int)0x80000000);
_PS_CONST_TYPE(inv_sign_mask_ps, int, ~0x80000000);
_PI32_CONST(sign_mask, 0x80000000);
_PI32_CONST(inv_sign_mask, ~0x80000000);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);
_PI32_CONST(0, 0);
_PI32_CONST(ffff, 0xffffffff);

_PS_CONST(cephes_SQRTHF, 0.707106781186547524);
_PS_CONST(cephes_log_p0, 7.0376836292E-2);
_PS_CONST(cephes_log_p1, - 1.1514610310E-1);
_PS_CONST(cephes_log_p2, 1.1676998740E-1);
_PS_CONST(cephes_log_p3, - 1.2420140846E-1);
_PS_CONST(cephes_log_p4, + 1.4249322787E-1);
_PS_CONST(cephes_log_p5, - 1.6668057665E-1);
_PS_CONST(cephes_log_p6, + 2.0000714765E-1);
_PS_CONST(cephes_log_p7, - 2.4999993993E-1);
_PS_CONST(cephes_log_p8, + 3.3333331174E-1);
_PS_CONST(cephes_log_q1, -2.12194440e-4);
_PS_CONST(cephes_log_q2, 0.693359375);

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


// exp()
static inline mic_m512_t mic_exp_ps(mic_m512_t x) {
	mic_m512_t zero = *(mic_m512_t*)_ps_0;
	mic_m512_t one = *(mic_m512_t*)_ps_1;

	x = _mm512_min_ps(x, *(mic_m512_t*)_ps_exp_hi);
	x = _mm512_max_ps(x, *(mic_m512_t*)_ps_exp_lo);

	//mic_m512_t temp_2 = _mm512_mul_ps(x, *(mic_m512_t*)_ps_cephes_LOG2EF);
	//temp_2 = _mm512_add_ps(temp_2, *(mic_m512_t*)_ps_0p5);
	mic_m512_t temp_2 = _mm512_fmadd_ps(x, *(mic_m512_t*)_ps_cephes_LOG2EF, *(mic_m512_t*)_ps_0p5);

	mic_m512_t temp_1 = _mm512_round_ps(temp_2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	mic_m512_t temp_3 = _mm512_sub_ps(temp_1, temp_2);
	__mmask16 mask = _mm512_cmp_ps_mask(temp_3, zero, _MM_CMPINT_GT);

	temp_2 = _mm512_mask_sub_ps(temp_1, mask, temp_1, one);
	__m512i emm0 = _mm512_cvtfxpnt_round_adjustps_epi32(temp_2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	//temp_1 = _mm512_mul_ps(temp_2, *(mic_m512_t*)_ps_cephes_exp_C12);
	//x = _mm512_sub_ps(x, temp_1);
	x = _mm512_fnmadd_ps(temp_2, *(mic_m512_t*)_ps_cephes_exp_C12, x);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);
 
	//temp_2 = _mm512_mul_ps(x2, *(mic_m512_t*)_ps_cephes_exp_p5);
	//temp_1 = _mm512_add_ps(temp_1, temp_2);
	//temp_3 = _mm512_mul_ps(x3, *(mic_m512_t*)_ps_cephes_exp_p4);
	//temp_1 = _mm512_add_ps(temp_1, temp_3);

	temp_1 = _mm512_add_ps(x, one);
	temp_1 = _mm512_fmadd_ps(x2, *(mic_m512_t*)_ps_cephes_exp_p5, temp_1);
	temp_1 = _mm512_fmadd_ps(x3, *(mic_m512_t*)_ps_cephes_exp_p4, temp_1);

	temp_2 = _mm512_mul_ps(x3, *(mic_m512_t*)_ps_cephes_exp_p0);
	temp_3 = _mm512_mul_ps(x2, *(mic_m512_t*)_ps_cephes_exp_p1);

	mic_m512_t temp_4 = _mm512_mul_ps(x, *(mic_m512_t*)_ps_cephes_exp_p2);

	emm0 = _mm512_add_epi32(emm0, *(__m512i*)_pi32_0x7f);

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


_PS_CONST(minus_cephes_DP1, -0.78515625);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8);
_PS_CONST(minus_cephes_DP123, -0.7853981633974483096156608);
_PS_CONST(sincof_p0, -1.9515295891E-4);
_PS_CONST(sincof_p1,	8.3321608736E-3);
_PS_CONST(sincof_p2, -1.6666654611E-1);
_PS_CONST(coscof_p0,	2.443315711809948E-005);
_PS_CONST(coscof_p1, -1.388731625493765E-003);
_PS_CONST(coscof_p2,	4.166664568298827E-002);
_PS_CONST(cephes_FOPI, 1.27323954473516); // 4 / M_PI


// sin()
static inline mic_m512_t mic_sin_ps(mic_m512_t x) {
	__m512i sign_bit;
	sign_bit = _mm512_and_epi32(_mm512_castps_si512(x), *(__m512i*)_pi32_sign_mask);
	x = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x), *(__m512i*)_pi32_inv_sign_mask));
	
	mic_m512_t y = _mm512_mul_ps(x, *(mic_m512_t*)_ps_cephes_FOPI);

	__m512i emm2 = _mm512_cvtfxpnt_round_adjustps_epi32(y, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	emm2 = _mm512_add_epi32(emm2, *(__m512i*)_pi32_1);
	emm2 = _mm512_and_epi32(emm2, *(__m512i*)_pi32_inv1);
	y = _mm512_cvtfxpnt_round_adjustepu32_ps(emm2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	__m512i emm0 = _mm512_and_epi32(emm2, *(__m512i*)_pi32_4);
	emm0 = _mm512_slli_epi32(emm0, 29);

	emm2 = _mm512_and_epi32(emm2, *(__m512i*)_pi32_2);
	__mmask16 mask = _mm512_cmp_epi32_mask(emm2, *(__m512i*)_pi32_0, _MM_CMPINT_EQ);
	emm2 = _mm512_mask_add_epi32(*(__m512i*)_pi32_0, mask, *(__m512i*)_pi32_ffff, *(__m512i*)_pi32_0);
	
	sign_bit = _mm512_xor_epi32(sign_bit, emm0);
	
	mic_m512_t temp = *(mic_m512_t*)_ps_minus_cephes_DP123;
	temp = _mm512_mul_ps(y, temp);
	x = _mm512_add_ps(x, temp);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);

	y = *(mic_m512_t*)_ps_coscof_p0;
	mic_m512_t y2 = *(mic_m512_t*)_ps_sincof_p0;
	y = _mm512_mul_ps(y, x2);
	y2 = _mm512_mul_ps(y2, x2);
	y = _mm512_add_ps(y, *(mic_m512_t*)_ps_coscof_p1);
	y2 = _mm512_add_ps(y2, *(mic_m512_t*)_ps_sincof_p1);
	y = _mm512_mul_ps(y, x2);
	y2 = _mm512_mul_ps(y2, x2);
	y = _mm512_add_ps(y, *(mic_m512_t*)_ps_coscof_p2);
	y2 = _mm512_add_ps(y2, *(mic_m512_t*)_ps_sincof_p2);
	y = _mm512_mul_ps(y, x4);
	y2 = _mm512_mul_ps(y2, x3);
	temp = _mm512_mul_ps(x2, *(mic_m512_t*)_ps_0p5);
	temp = _mm512_sub_ps(temp, *(mic_m512_t*)_ps_1);
	y = _mm512_sub_ps(y, temp);
	y2 = _mm512_add_ps(y2, x);

	y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(y)));
	y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(y2)));

	y = _mm512_add_ps(y, y2);
	y = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(y), sign_bit));

	return y;
} // sin_ps()


static inline mic_m512_t mic_cos_ps(mic_m512_t x) {
	x = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x), *(__m512i*)_pi32_inv_sign_mask));

	mic_m512_t y = _mm512_mul_ps(x, *(mic_m512_t*)_ps_cephes_FOPI);

	__m512i emm2 = _mm512_cvtfxpnt_round_adjustps_epi32(y, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	emm2 = _mm512_add_epi32(emm2, *(__m512i*)_pi32_1);
	emm2 = _mm512_and_epi32(emm2, *(__m512i*)_pi32_inv1);
	y = _mm512_cvtfxpnt_round_adjustepu32_ps(emm2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	emm2 = _mm512_sub_epi32(emm2, *(__m512i*)_pi32_2);

	__m512i emm0 = _mm512_andnot_epi32(emm2, *(__m512i*)_pi32_4);
	emm0 = _mm512_slli_epi32(emm0, 29);

	emm2 = _mm512_and_epi32(emm2, *(__m512i*)_pi32_2);
	__mmask16 mask = _mm512_cmp_epi32_mask(emm2, *(__m512i*)_pi32_0, _MM_CMPINT_EQ);
	emm2 = _mm512_mask_add_epi32(*(__m512i*)_pi32_0, mask, *(__m512i*)_pi32_ffff, *(__m512i*)_pi32_0);
	
	mic_m512_t temp = *(mic_m512_t*)_ps_minus_cephes_DP123;
	temp = _mm512_mul_ps(y, temp);
	x = _mm512_add_ps(x, temp);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);

	y = *(mic_m512_t*)_ps_coscof_p0;
	mic_m512_t y2 = *(mic_m512_t*)_ps_sincof_p0;
	y = _mm512_mul_ps(y, x2);
	y2 = _mm512_mul_ps(y2, x2);
	y = _mm512_add_ps(y, *(mic_m512_t*)_ps_coscof_p1);
	y2 = _mm512_add_ps(y2, *(mic_m512_t*)_ps_sincof_p1);
	y = _mm512_mul_ps(y, x2);
	y2 = _mm512_mul_ps(y2, x2);
	y = _mm512_add_ps(y, *(mic_m512_t*)_ps_coscof_p2);
	y2 = _mm512_add_ps(y2, *(mic_m512_t*)_ps_sincof_p2);
	y = _mm512_mul_ps(y, x4);
	y2 = _mm512_mul_ps(y2, x3);
	temp = _mm512_mul_ps(x2, *(mic_m512_t*)_ps_0p5);
	temp = _mm512_sub_ps(temp, *(mic_m512_t*)_ps_1);
	y = _mm512_sub_ps(y, temp);
	y2 = _mm512_add_ps(y2, x);

	y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(y)));
	y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(y2)));

	y = _mm512_add_ps(y, y2);
	y = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(y), emm0));

	return y;
} // cos_ps()


// FIXME: this is not working correctly ... individual sin and cos are working good ...
static inline void mic_sincos_ps(mic_m512_t x, mic_m512_t *s, mic_m512_t *c) {
	__m512i sign_bit = _mm512_and_epi32(_mm512_castps_si512(x), *(__m512i*)_pi32_sign_mask);
	x = _mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(x), *(__m512i*)_pi32_inv_sign_mask));

	mic_m512_t y = _mm512_mul_ps(x, *(mic_m512_t*)_ps_cephes_FOPI);

	__m512i emm2 = _mm512_cvtfxpnt_round_adjustps_epi32(y, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);
	emm2 = _mm512_add_epi32(emm2, *(__m512i*)_pi32_1);
	emm2 = _mm512_and_epi32(emm2, *(__m512i*)_pi32_inv1);
	y = _mm512_cvtfxpnt_round_adjustepu32_ps(emm2, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_NONE);

	__m512i cos_emm2 = _mm512_sub_epi32(emm2, *(__m512i*)_pi32_2);

	__m512i temp_1 = *(__m512i*)_pi32_4;
	__m512i emm0 = _mm512_and_epi32(emm2, temp_1);
	__m512i cos_emm0 = _mm512_andnot_epi32(cos_emm2, temp_1);
	emm0 = _mm512_slli_epi32(emm0, 29);
	cos_emm0 = _mm512_slli_epi32(cos_emm0, 29);

	temp_1 = *(__m512i*)_pi32_2;
	emm2 = _mm512_and_epi32(emm2, temp_1);
	cos_emm2 = _mm512_and_epi32(cos_emm2, temp_1);
	__mmask16 mask = _mm512_cmp_epi32_mask(emm2, *(__m512i*)_pi32_0, _MM_CMPINT_EQ);
	emm2 = _mm512_mask_add_epi32(*(__m512i*)_pi32_0, mask, *(__m512i*)_pi32_ffff, *(__m512i*)_pi32_0);
	__mmask16 cos_mask = _mm512_cmp_epi32_mask(cos_emm2, *(__m512i*)_pi32_0, _MM_CMPINT_EQ);
	cos_emm2 = _mm512_mask_add_epi32(*(__m512i*)_pi32_0, cos_mask, *(__m512i*)_pi32_ffff, *(__m512i*)_pi32_0);
	
	sign_bit = _mm512_xor_epi32(sign_bit, emm0);

	x = _mm512_fmadd_ps(y, *(mic_m512_t*)_ps_minus_cephes_DP123, x);

	mic_m512_t x2 = _mm512_mul_ps(x, x);
	mic_m512_t x3 = _mm512_mul_ps(x2, x);
	mic_m512_t x4 = _mm512_mul_ps(x2, x2);

	//y = *(mic_m512_t*)_ps_coscof_p0;
	//y = _mm512_mul_ps(y, x2);
	//y = _mm512_add_ps(y, *(mic_m512_t*)_ps_coscof_p1);
	//y = _mm512_mul_ps(y, x2);
	//y = _mm512_add_ps(y, *(mic_m512_t*)_ps_coscof_p2);
	//temp_2 = _mm512_mul_ps(x2, *(mic_m512_t*)_ps_0p5);
	//temp_2 = _mm512_sub_ps(temp_2, *(mic_m512_t*)_ps_1);
	//y = _mm512_mul_ps(y, x4);
	//y = _mm512_sub_ps(y, temp_2);

	//mic_m512_t y2 = *(mic_m512_t*)_ps_sincof_p0;
	//y2 = _mm512_mul_ps(y2, x2);
	//y2 = _mm512_add_ps(y2, *(mic_m512_t*)_ps_sincof_p1);
	//y2 = _mm512_mul_ps(y2, x2);
	//y2 = _mm512_add_ps(y2, *(mic_m512_t*)_ps_sincof_p2);
	//y2 = _mm512_mul_ps(y2, x3);
	//y2 = _mm512_add_ps(y2, x);

	y = _mm512_fmadd_ps(*(mic_m512_t*)_ps_coscof_p0, x2, *(mic_m512_t*)_ps_coscof_p1);
	mic_m512_t y2 = _mm512_fmadd_ps(*(mic_m512_t*)_ps_sincof_p0, x2, *(mic_m512_t*)_ps_sincof_p1);
	y = _mm512_fmadd_ps(y, x2, *(mic_m512_t*)_ps_coscof_p2);
	y2 = _mm512_fmadd_ps(y2, x2, *(mic_m512_t*)_ps_sincof_p2);
	mic_m512_t temp_2 = _mm512_fmsub_ps(x2, *(mic_m512_t*)_ps_0p5, *(mic_m512_t*)_ps_1);
	y2 = _mm512_fmadd_ps(y2, x3, x);
	y = _mm512_fmsub_ps(y, x4, temp_2);

	mic_m512_t cos_y = y;
	mic_m512_t cos_y2 = y2;

	y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(y)));
	cos_y = _mm512_castsi512_ps(_mm512_andnot_epi32(emm2, _mm512_castps_si512(cos_y)));
	y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(y2)));
	cos_y2 = _mm512_castsi512_ps(_mm512_and_epi32(emm2, _mm512_castps_si512(cos_y2)));

	y = _mm512_add_ps(y, y2);
	cos_y = _mm512_add_ps(cos_y, cos_y2);

	*s = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(y), sign_bit));
	*c = _mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(cos_y), cos_emm0));
} // sincos_ps()

#pragma offload_attribute(pop)
