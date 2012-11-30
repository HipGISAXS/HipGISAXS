/***
  *  $Id$
  *
  *  Project:
  *
  *  File: cu_complex_numeric.cuh
  *  Created: Oct 17, 2012
  *  Modified: Thu 29 Nov 2012 03:21:39 PM PST
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

	/**
	 * Atomics
	 */

	// atomic add for double precision
	__device__ double atomicAdd(double* address, double val) {
		unsigned long long int* ull_address = (unsigned long long int*) address;
		unsigned long long int old = *ull_address;
		unsigned long long int assumed;
		do {
			assumed = old;
			old = atomicCAS(ull_address, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
		} while (assumed != old);
		return __longlong_as_double(old);
	} // atomicAdd()


} // namespace

#endif // _CU_COMPLEX_NUMERIC_CUH_
