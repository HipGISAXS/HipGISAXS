/***
  *  $Id$
  *
  *  Project:
  *
  *  File: cu_complex_numeric.cuh
  *  Created: Oct 17, 2012
  *  Modified: Sat 23 Feb 2013 01:11:34 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _CU_COMPLEX_NUMERIC_CUH_
#define _CU_COMPLEX_NUMERIC_CUH_

#include <cuComplex.h>
#include <cmath>		// for real valued functions

// already defined:
// cuCadd, cuCsub, cuCmul, cuCdiv, 

namespace hig {

	__device__ static __inline__ cuFloatComplex make_cuC(float r, float i) {
		return make_cuFloatComplex(r, i);
	} // make_cuC()

	__device__ static __inline__ cuDoubleComplex make_cuC(double r, double i) {
		return make_cuDoubleComplex(r, i);
	} // make_cuC()


	// addition

	__device__ static __inline__ cuFloatComplex operator+(cuFloatComplex a, cuFloatComplex b) {
		return make_cuFloatComplex(a.x + b.x, a.y + b.y);
	} // operator+()

	__device__ static __inline__ cuFloatComplex operator+(cuFloatComplex a, float b) {
		return make_cuFloatComplex(a.x + b, a.y);
	} // operator+()

	__device__ static __inline__ cuFloatComplex operator+(float a, cuFloatComplex b) {
		return make_cuFloatComplex(a + b.x, b.y);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator+(cuDoubleComplex a, cuDoubleComplex b) {
		return make_cuDoubleComplex(a.x + b.x, a.y + b.y);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator+(cuDoubleComplex a, double b) {
		return make_cuDoubleComplex(a.x + b, a.y);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator+(double a, cuDoubleComplex b) {
		return make_cuDoubleComplex(a + b.x, b.y);
	} // operator+()


	// subtraction

	__device__ static __inline__ cuFloatComplex operator-(cuFloatComplex a, cuFloatComplex b) {
		return make_cuFloatComplex(a.x - b.x, a.y - b.y);
	} // operator-()

	__device__ static __inline__ cuFloatComplex operator-(cuFloatComplex a, float b) {
		return make_cuFloatComplex(a.x - b, a.y);
	} // operator-()

	__device__ static __inline__ cuFloatComplex operator-(float a, cuFloatComplex b) {
		return make_cuFloatComplex(a - b.x, - b.y);
	} // operator-()

	__device__ static __inline__ cuDoubleComplex operator-(cuDoubleComplex a, cuDoubleComplex b) {
		return make_cuDoubleComplex(a.x - b.x, a.y - b.y);
	} // operator-()

	__device__ static __inline__ cuDoubleComplex operator-(cuDoubleComplex a, double b) {
		return make_cuDoubleComplex(a.x - b, a.y);
	} // operator-()

	__device__ static __inline__ cuDoubleComplex operator-(double a, cuDoubleComplex b) {
		return make_cuDoubleComplex(a - b.x, - b.y);
	} // operator-()


	// multiplication

	__device__ static __inline__ cuFloatComplex operator*(cuFloatComplex a, cuFloatComplex b) {
		return make_cuFloatComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
	} // operator+()

	__device__ static __inline__ cuFloatComplex operator*(cuFloatComplex a, float b) {
		return make_cuFloatComplex(a.x * b, a.y);
	} // operator+()

	__device__ static __inline__ cuFloatComplex operator*(float a, cuFloatComplex b) {
		return make_cuFloatComplex(a * b.x, b.y);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b) {
		return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator*(cuDoubleComplex a, double b) {
		return make_cuDoubleComplex(a.x * b, a.y);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator*(double a, cuDoubleComplex b) {
		return make_cuDoubleComplex(a * b.x, b.y);
	} // operator+()


	// division

	__device__ static __inline__ cuFloatComplex operator/(cuFloatComplex a, cuFloatComplex b) {
		return cuCdivf(a, b);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator/(cuDoubleComplex a, cuDoubleComplex b) {
		return cuCdiv(a, b);
	} // operator+()

	__device__ static __inline__ cuFloatComplex operator/(cuFloatComplex a, float b) {
		return make_cuFloatComplex(a.x / b, a.y / b);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator/(cuDoubleComplex a, double b) {
		return make_cuDoubleComplex(a.x / b, a.y / b);
	} // operator+()

	__device__ static __inline__ cuFloatComplex operator/(float a, cuFloatComplex b) {
		return cuCdivf(make_cuFloatComplex(a, 0.0f), b);
	} // operator+()

	__device__ static __inline__ cuDoubleComplex operator/(double a, cuDoubleComplex b) {
		return cuCdiv(make_cuDoubleComplex(a, 0.0), b);
	} // operator+()


/*	__device__ static __inline__ float cuCabs(cuFloatComplex z) {
		float x = z.x;
		float y = z.y;
		return sqrtf(x * x + y * y);
	} // cuCabs()
	__device__ static __inline__ double cuCabs(cuDoubleComplex z) {
		double x = z.x;
		double y = z.y;
		return sqrt(x * x + y * y);
	} // cuCabs()
*/
	__device__ static __inline__ cuFloatComplex cuCsqrt(cuFloatComplex z) {
		float x = z.x;
		float y = z.y;
		if (x == 0) {
			float t = sqrtf(fabsf(y) / 2);
			return make_cuC(t, y < 0 ? -t : t);
		} else if (y == 0) {
			float t = sqrtf(fabsf(x));
			return x < 0 ? make_cuC(0, t) : make_cuC(t, 0);
		} else {
			//float t = sqrtf(2 * (cuCabsf(z) + fabsf(x)));
			float t = sqrtf(2 * (sqrtf(x * x + y * y) + fabsf(x)));
			float u = t / 2;
			return x > 0 ? make_cuC(u, y / t) : make_cuC(fabsf(y) / t, y < 0 ? -u : u);
		} // if-else
	} // cuCsqrt()

	__device__ static __inline__ cuDoubleComplex cuCsqrt(cuDoubleComplex z) {
		double x = z.x;
		double y = z.y;
		if (x == 0) {
			double t = sqrt(fabs(y) / 2);
			return make_cuC(t, y < 0 ? -t : t);
		} else if (y == 0) {
			return x < 0 ? make_cuC(0, sqrt(fabs(x))) : make_cuC(sqrt(x), 0);
		} else {
			//double t = sqrt(2 * (cuCabs(z) + fabs(x)));
			double t = sqrt(2 * (sqrt(x * x + y * y) + fabs(x)));
			double u = t / 2;
			return x > 0 ? make_cuC(u, y / t) : make_cuC(fabs(y) / t, y < 0 ? -u : u);
		} // if-else
	} // cuCsqrt()

	__device__ static __inline__ cuFloatComplex cuCpow(cuFloatComplex a, int p) {
		// ...
		return a;
	} // cuCpow()

	__device__ static __inline__ cuDoubleComplex cuCpow(cuDoubleComplex a, int p) {
		// ...
		return a;
	} // cuCpow()

	// eucledian norm
	__device__ static __inline__ cuFloatComplex cuCnorm3(
				cuFloatComplex a, cuFloatComplex b, cuFloatComplex c) {
		cuFloatComplex a2 = cuCmulf(a, a);
		cuFloatComplex b2 = cuCmulf(b, b);
		cuFloatComplex c2 = cuCmulf(c, c);
		//cuFloatComplex n = cuCsqrt(cuCaddf(cuCaddf(a2, b2), c2));
		cuFloatComplex n = cuCsqrt(make_cuC(a2.x + b2.x + c2.x, a2.y + b2.y + c2.y));
		return n;
	} // cuCnorm3()

	__device__ static __inline__ cuDoubleComplex cuCnorm3(
				cuDoubleComplex a, cuDoubleComplex b, cuDoubleComplex c) {
		cuDoubleComplex a2 = cuCmul(a, a);
		cuDoubleComplex b2 = cuCmul(b, b);
		cuDoubleComplex c2 = cuCmul(c, c);
		//cuDoubleComplex n = cuCsqrt(cuCadd(cuCadd(a2, b2), c2));
		cuDoubleComplex n = cuCsqrt(make_cuC(a2.x + b2.x + c2.x, a2.y + b2.y + c2.y));
		return n;
	} // cuCnorm3()

	__device__ static __inline__ cuFloatComplex cuCnorm3(float a, float b, cuFloatComplex c) {
		return cuCnorm3(make_cuFloatComplex(a, 0), make_cuFloatComplex(b, 0), c);
	} // cuCnorm3()

	__device__ static __inline__ cuDoubleComplex cuCnorm3(double a, double b, cuDoubleComplex c) {
		return cuCnorm3(make_cuDoubleComplex(a, 0), make_cuDoubleComplex(b, 0), c);
	} // cuCnorm3()

	// e^z = e^z.x (cos(z.y) + isin(z.y))
	__device__ static __inline__ cuFloatComplex cuCexp(cuFloatComplex z) {
		float temp1 = cosf(z.y);
		float temp2 = sinf(z.y);
		float temp3 = expf(z.x);
		return make_cuFloatComplex(temp1 * temp3, temp2 * temp3);
	} // cuCexp()

	__device__ static __inline__ cuDoubleComplex cuCexp(cuDoubleComplex z) {
		double temp1 = cos(z.y);
		double temp2 = sin(z.y);
		double temp3 = exp(z.x);
		return make_cuDoubleComplex(temp1 * temp3, temp2 * temp3);
	} // cuCexp()

	// e^if = cos(f) + isin(f)
	__device__ static __inline__ cuFloatComplex cuCexpi(float f) {
		return make_cuFloatComplex(cosf(f), sinf(f));
	} // cuCexpi()

	__device__ static __inline__ cuDoubleComplex cuCexpi(double f) {
		return make_cuDoubleComplex(cos(f), sin(f));
	} // cuCexpi()

	// e^iz
	__device__ static __inline__ cuFloatComplex cuCexpi(cuFloatComplex z) {
		return cuCexp(make_cuFloatComplex(-z.y, z.x));
	} // cuCexpi()

	__device__ static __inline__ cuDoubleComplex cuCexpi(cuDoubleComplex z) {
		return cuCexp(make_cuDoubleComplex(-z.y, z.x));
	} // cuCexpi()

	__device__ static __inline__ cuFloatComplex cuCsin(cuFloatComplex z) {
		float x = z.x;
		float y = z.y;
		return make_cuC(sinf(x) * coshf(y), cosf(x) * sinhf(y));
	} // cuCsin()

	__device__ static __inline__ cuDoubleComplex cuCsin(cuDoubleComplex z) {
		double x = z.x;
		double y = z.y;
		return make_cuC(sin(x) * cosh(y), cos(x) * sinh(y));
	} // cuCsin()

	__device__ static __inline__ cuFloatComplex cuCcos(cuFloatComplex z) {
		float x = z.x;
		float y = z.y;
		return make_cuC(cosf(x) * coshf(y), -sinf(x) * sinhf(y));
	} // cuCsin()

	__device__ static __inline__ cuDoubleComplex cuCcos(cuDoubleComplex z) {
		double x = z.x;
		double y = z.y;
		return make_cuC(cos(x) * cosh(y), -sin(x) * sinh(y));
	} // cuCsin()


	__device__ static __inline__ cuFloatComplex cuCsinc(cuFloatComplex value) {
		cuFloatComplex temp;
		if(fabsf(value.x) < 1e-14 && fabsf(value.y) < 1e-14) temp = make_cuFloatComplex(1.0, 0.0);
		else temp = cuCsin(value) / value;
		return temp;
	} // cuCsinc()

	__device__ static __inline__ cuDoubleComplex cuCsinc(cuDoubleComplex value) {
		cuDoubleComplex temp;
		if(fabs(value.x) < 1e-14 && fabs(value.y) < 1e-14) temp = make_cuDoubleComplex(1.0, 0.0);
		else temp = cuCsin(value) / value;
		return temp;
	} // cuCsinc()

	__device__ static __inline__ bool cuCiszero(cuFloatComplex value) {
		return (fabsf(value.x) < 1e-14 && fabsf(value.y) < 1e-14);
	} // cuCiszero()

	__device__ static __inline__ bool cuCiszero(cuDoubleComplex value) {
		return (fabs(value.x) < 1e-14 && fabs(value.y) < 1e-14);
	} // cuCiszero()


	/**
	 * Atomics
	 */

	// atomic add for double precision
	__device__ static __inline__ double atomicAdd(double* address, double val) {
		unsigned long long int* ull_address = (unsigned long long int*) address;
		unsigned long long int old = *ull_address;
		unsigned long long int assumed;
		do {
			assumed = old;
			old = atomicCAS(ull_address, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
		} while (assumed != old);
		return __longlong_as_double(old);
	} // atomicAdd()


	/**
	 * Others
	 */

	__device__ static __inline__ double cugamma(double x) {
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
	__device__ static __inline__ cuDoubleComplex cuCcbessj(cuDoubleComplex zz, int order) {
		cuDoubleComplex z = zz;
		cuDoubleComplex temp_z = z / 2.0;
		cuDoubleComplex temp1 = cuCpow(temp_z, order);
		cuDoubleComplex z2 = -1.0 * temp_z * temp_z;
		cuDoubleComplex sum = make_cuDoubleComplex(0.0, 0.0);
		double factorial_k = 1.0;
		cuDoubleComplex pow_z2_k = make_cuDoubleComplex(1.0, 0.0);
		for(unsigned int k = 0; k <= MAXK; ++ k) {
			if(k == 0) factorial_k = 1.0;
			else factorial_k *= k;
			cuDoubleComplex temp2 = pow_z2_k / (factorial_k * cugamma((double) order + k));
			sum = sum + temp2;
			pow_z2_k = pow_z2_k * z2;
		} // for
		temp1 = temp1 * sum;
		return temp1;
	} // cuCcbessj()

	__device__ static __inline__ cuFloatComplex cuCcbessj(cuFloatComplex zz, int order) {
		cuDoubleComplex temp = cuCcbessj(make_cuDoubleComplex((double) zz.x, (double) zz.y), order);
		return make_cuFloatComplex((float) temp.x, (float) temp.y);
	} // cuCcbessj()
	
} // namespace

#endif // _CU_COMPLEX_NUMERIC_CUH_
