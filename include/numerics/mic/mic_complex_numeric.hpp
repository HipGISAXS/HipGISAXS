/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: mic_complex_numeric.hpp
 *  Created: Apr 02, 2013
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

#ifndef __MIC_COMPLEX_NUMERIC_HPP__
#define __MIC_COMPLEX_NUMERIC_HPP__

#include <cmath>		// for real valued functions


namespace hig {

	// construction

	__attribute__((target(mic:0)))
	static inline float2_t make_sC(float r, float i) {
		float2_t temp;
		temp.x = r; temp.y = i;
		return temp;
	} // make_sC()

	__attribute__((target(mic:0)))
	static inline double2_t make_sC(double r, double i) {
		double2_t temp;
		temp.x = r; temp.y = i;
		return temp;
	} // make_sC()


	// addition

	__attribute__((target(mic:0)))
	static inline float2_t operator+(float2_t a, float2_t b) {
		return make_sC(a.x + b.x, a.y + b.y);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline float2_t operator+(float2_t a, float b) {
		return make_sC(a.x + b, a.y);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline float2_t operator+(float a, float2_t b) {
		return make_sC(a + b.x, b.y);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator+(double2_t a, double2_t b) {
		return make_sC(a.x + b.x, a.y + b.y);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator+(double2_t a, double b) {
		return make_sC(a.x + b, a.y);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator+(double a, double2_t b) {
		return make_sC(a + b.x, b.y);
	} // operator+()


	// subtraction

	__attribute__((target(mic:0)))
	static inline float2_t operator-(float2_t a, float2_t b) {
		return make_sC(a.x - b.x, a.y - b.y);
	} // operator-()

	__attribute__((target(mic:0)))
	static inline float2_t operator-(float2_t a, float b) {
		return make_sC(a.x - b, a.y);
	} // operator-()

	__attribute__((target(mic:0)))
	static inline float2_t operator-(float a, float2_t b) {
		return make_sC(a - b.x, - b.y);
	} // operator-()

	__attribute__((target(mic:0)))
	static inline double2_t operator-(double2_t a, double2_t b) {
		return make_sC(a.x - b.x, a.y - b.y);
	} // operator-()

	__attribute__((target(mic:0)))
	static inline double2_t operator-(double2_t a, double b) {
		return make_sC(a.x - b, a.y);
	} // operator-()

	__attribute__((target(mic:0)))
	static inline double2_t operator-(double a, double2_t b) {
		return make_sC(a - b.x, - b.y);
	} // operator-()


	// multiplication

	__attribute__((target(mic:0)))
	static inline float2_t operator*(float2_t a, float2_t b) {
		return make_sC(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline float2_t operator*(float2_t a, float b) {
		return make_sC(a.x * b, a.y * b);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline float2_t operator*(float a, float2_t b) {
		return make_sC(a * b.x, a * b.y);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator*(double2_t a, double2_t b) {
		return make_sC(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator*(double2_t a, double b) {
		return make_sC(a.x * b, a.y * b);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator*(double a, double2_t b) {
		return make_sC(a * b.x, a * b.y);
	} // operator+()


	// division

	__attribute__((target(mic:0)))
	static inline float2_t operator/(float2_t a, float2_t b) {
		float s = fabsf(b.x) + fabsf(b.x);
		float oos = 1.0f / s;
		float ars = a.x * oos;
		float ais = a.y * oos;
		float brs = b.x * oos;
		float bis = b.y * oos;
		s = (brs * brs) + (bis * bis);
		oos = 1.0f / s;
		float2_t quot = make_sC(((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos);
		return quot;
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator/(double2_t a, double2_t b) {
		double s = fabs(b.x) + fabs(b.x);
		double oos = 1.0 / s;
		double ars = a.x * oos;
		double ais = a.y * oos;
		double brs = b.x * oos;
		double bis = b.y * oos;
		s = (brs * brs) + (bis * bis);
		oos = 1.0 / s;
		double2_t quot = make_sC(((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos);
		return quot;
	} // operator+()

	__attribute__((target(mic:0)))
	static inline float2_t operator/(float2_t a, float b) {
		return make_sC(a.x / b, a.y / b);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator/(double2_t a, double b) {
		return make_sC(a.x / b, a.y / b);
	} // operator+()

	__attribute__((target(mic:0)))
	static inline float2_t operator/(float a, float2_t b) {
		return make_sC(a, 0.0f) / b;
	} // operator+()

	__attribute__((target(mic:0)))
	static inline double2_t operator/(double a, double2_t b) {
		return make_sC(a, 0.0) / b;
	} // operator+()


	// square root

	__attribute__((target(mic:0)))
	static inline float2_t sCsqrt(float2_t z) {
		float x = z.x;
		float y = z.y;
		if(x == 0) {
			float t = sqrtf(fabsf(y) / 2.0f);
			return make_sC(t, y < 0 ? -t : t);
		} else if(y == 0) {
			float t = sqrtf(fabsf(x));
			return x < 0 ? make_sC(0.0f, t) : make_sC(t, 0.0f);
		} else {
			float t = sqrtf(2.0f * (sqrtf(x * x + y * y) + fabsf(x)));
			float u = t / 2.0f;
			return x > 0 ? make_sC(u, y / t) : make_sC(fabsf(y) / t, y < 0 ? -u : u);
		} // if-else
	} // sCsqrt()

	__attribute__((target(mic:0)))
	static inline double2_t sCsqrt(double2_t z) {
		double x = z.x;
		double y = z.y;
		if(x == 0) {
			double t = sqrt(fabs(y) / 2.0);
			return make_sC(t, y < 0 ? -t : t);
		} else if (y == 0) {
			return x < 0 ? make_sC(0.0, sqrt(fabs(x))) : make_sC(sqrt(x), 0.0);
		} else {
			double t = sqrt(2.0 * (sqrt(x * x + y * y) + fabs(x)));
			double u = t / 2.0;
			return x > 0 ? make_sC(u, y / t) : make_sC(fabs(y) / t, y < 0 ? -u : u);
		} // if-else
	} // sCsqrt()


	// eucledian norm

	/*static inline float2_t sCnorm3(
				float2_t a, float2_t b, float2_t c) {
		float2_t a2 = sCmulf(a, a);
		float2_t b2 = sCmulf(b, b);
		float2_t c2 = sCmulf(c, c);
		//float2_t n = sCsqrt(sCaddf(sCaddf(a2, b2), c2));
		float2_t n = sCsqrt(make_sC(a2.x + b2.x + c2.x, a2.y + b2.y + c2.y));
		return n;
	} // sCnorm3()

	static inline double2_t sCnorm3(
				double2_t a, double2_t b, double2_t c) {
		double2_t a2 = sCmul(a, a);
		double2_t b2 = sCmul(b, b);
		double2_t c2 = sCmul(c, c);
		//double2_t n = sCsqrt(sCadd(sCadd(a2, b2), c2));
		double2_t n = sCsqrt(make_sC(a2.x + b2.x + c2.x, a2.y + b2.y + c2.y));
		return n;
	} // sCnorm3()

	static inline float2_t sCnorm3(float a, float b, float2_t c) {
		return sCnorm3(make_float2_t(a, 0), make_float2_t(b, 0), c);
	} // sCnorm3()

	static inline double2_t sCnorm3(double a, double b, double2_t c) {
		return sCnorm3(make_double2_t(a, 0), make_double2_t(b, 0), c);
	} // sCnorm3()*/


	// e^z = e^z.x (cos(z.y) + isin(z.y))
	
/*	static inline float2_t sCexp(float2_t z) {
		float temp1 = cosf(z.y);
		float temp2 = sinf(z.y);
		float temp3 = expf(z.x);
		return make_float2_t(temp1 * temp3, temp2 * temp3);
	} // sCexp()

	static inline double2_t sCexp(double2_t z) {
		double temp1 = cos(z.y);
		double temp2 = sin(z.y);
		double temp3 = exp(z.x);
		return make_double2_t(temp1 * temp3, temp2 * temp3);
	} // sCexp()


	// e^if = cos(f) + isin(f)
	
	static inline float2_t sCexpi(float f) {
		return make_float2_t(cosf(f), sinf(f));
	} // sCexpi()

	static inline double2_t sCexpi(double f) {
		return make_double2_t(cos(f), sin(f));
	} // sCexpi()


	// e^iz
	
	static inline float2_t sCexpi(float2_t z) {
		return sCexp(make_float2_t(-z.y, z.x));
	} // sCexpi()

	static inline double2_t sCexpi(double2_t z) {
		return sCexp(make_double2_t(-z.y, z.x));
	} // sCexpi()
*/

	// sine cosine

	__attribute__((target(mic:0)))
	static inline float2_t sCsin(float2_t z) {
		float x = z.x;
		float y = z.y;
		return make_sC(sinf(x) * coshf(y), cosf(x) * sinhf(y));
	} // sCsin()

	__attribute__((target(mic:0)))
	static inline double2_t sCsin(double2_t z) {
		double x = z.x;
		double y = z.y;
		return make_sC(sin(x) * cosh(y), cos(x) * sinh(y));
	} // sCsin()

	__attribute__((target(mic:0)))
	static inline float2_t sCcos(float2_t z) {
		float x = z.x;
		float y = z.y;
		return make_sC(cosf(x) * coshf(y), -sinf(x) * sinhf(y));
	} // sCsin()

	__attribute__((target(mic:0)))
	static inline double2_t sCcos(double2_t z) {
		double x = z.x;
		double y = z.y;
		return make_sC(cos(x) * cosh(y), -sin(x) * sinh(y));
	} // sCsin()


	// other trignometrics

/*	static inline float2_t sCsinc(float2_t value) {
		float2_t temp;
		if(fabsf(value.x) < 1e-14 && fabsf(value.y) < 1e-14) temp = make_float2_t(1.0, 0.0);
		else temp = sCsin(value) / value;
		return temp;
	} // sCsinc()

	static inline double2_t sCsinc(double2_t value) {
		double2_t temp;
		if(fabs(value.x) < 1e-14 && fabs(value.y) < 1e-14) temp = make_double2_t(1.0, 0.0);
		else temp = sCsin(value) / value;
		return temp;
	} // sCsinc()


	// check if zero

	static inline bool sCiszero(float2_t value) {
		return (fabsf(value.x) < 1e-14 && fabsf(value.y) < 1e-14);
	} // sCiszero()

	static inline bool sCiszero(double2_t value) {
		return (fabs(value.x) < 1e-14 && fabs(value.y) < 1e-14);
	} // sCiszero()
*/

	/**
	 * Atomics
	 */

	// atomic add for double precision

/*	static inline double atomicAdd(double* address, double val) {
		unsigned long long int* ull_address = (unsigned long long int*) address;
		unsigned long long int old = *ull_address;
		unsigned long long int assumed;
		do {
			assumed = old;
			old = atomicCAS(ull_address, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
		} while (assumed != old);
		return __longlong_as_double(old);
	} // atomicAdd()
*/

	/**
	 * Others
	 */

/*	static inline double cugamma(double x) {
		double coef[6];
		coef[0] = 76.18009173;
		coef[1] = -86.50532033;
		coef[2] = 24.01409822;
		coef[3] = -1.231739516;
		coef[4] = 0.120858003e-2;
		coef[5] = -0.536382e-5;
		double stp = 2.50662827465;
		double temp = x + 5.5;
		temp = (x + 0.5) * log(temp) - temp;
		double ser = 1.0;
		for(int j = 0; j < 6; ++ j) {
			x += 1.0;
			ser += coef[j] / x;
		} // for
		return exp(temp + log(stp * ser));
	} // cugamma()


	const int MAXK = 20;

	// bessel

	static inline double2_t sCcbessj(double2_t zz, int order) {
		double2_t z = zz;
		double2_t temp_z = z / 2.0;
		double2_t temp1 = sCpow(temp_z, order);
		double2_t z2 = -1.0 * temp_z * temp_z;
		double2_t sum = make_double2_t(0.0, 0.0);
		double factorial_k = 1.0;
		double2_t pow_z2_k = make_double2_t(1.0, 0.0);
		for(unsigned int k = 0; k <= MAXK; ++ k) {
			if(k == 0) factorial_k = 1.0;
			else factorial_k *= k;
			double2_t temp2 = pow_z2_k / (factorial_k * cugamma((double) order + k));
			sum = sum + temp2;
			pow_z2_k = pow_z2_k * z2;
		} // for
		temp1 = temp1 * sum;
		return temp1;
	} // sCcbessj()

	static inline float2_t sCcbessj(float2_t zz, int order) {
		double2_t temp = sCcbessj(make_double2_t((double) zz.x, (double) zz.y), order);
		return make_float2_t((float) temp.x, (float) temp.y);
	} // sCcbessj()*/
	
} // namespace

#endif // _CU_COMPLEX_NUMERIC_CUH_
