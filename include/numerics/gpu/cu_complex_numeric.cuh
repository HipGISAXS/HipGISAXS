/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: cu_complex_numeric.cuh
 *  Created: Oct 17, 2012
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

#ifndef _CU_COMPLEX_NUMERIC_CUH_
#define _CU_COMPLEX_NUMERIC_CUH_

#include <cuComplex.h>
#include <cmath>    // for real valued functions

// already defined:
// cuCadd, cuCsub, cuCmul, cuCdiv, 

namespace hig {

  __device__ static __inline__ float cuRtan(float a){
    return tanf(a);
  }

  __device__ static __inline__ double cuRtan(double a){
    return tan(a);
  }

  __device__ static __inline__ cuFloatComplex make_cuC(float r, float i) {
    return make_cuFloatComplex(r, i);
  } // make_cuC()

  __device__ static __inline__ cuDoubleComplex make_cuC(double r, double i) {
    return make_cuDoubleComplex(r, i);
  } // make_cuC()

  // get the real part from float type complex
  __device__ static __inline__ float cu_real (cuFloatComplex x) {
    return cuCrealf (x);
  }

  // get the real part from double type complex
  __device__ static __inline__ double cu_real (cuDoubleComplex x) {
    return cuCreal (x);
  }

  // get the imaginary part from float type complex
  __device__ static __inline__ float cu_imag (cuFloatComplex x) {
    return cuCimagf (x);
  }

  // get the imaginary part from double type complex
  __device__ static __inline__ double cu_imag (cuDoubleComplex x) {
    return cuCimag (x);
  }


  // single-precision conjugate
  __device__ static __inline__ cuFloatComplex cuCconj(cuFloatComplex x){
    return make_cuFloatComplex(x.x, -x.y);
  }

  // double-precision conjugate
  __device__ static __inline__ cuDoubleComplex cuCconj(cuDoubleComplex x){
    return make_cuDoubleComplex(x.x, -x.y);
  }

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
  } // operator*()

  __device__ static __inline__ cuFloatComplex operator*(cuFloatComplex a, float b) {
    return make_cuFloatComplex(a.x * b, a.y * b);
  } // operator*()

  __device__ static __inline__ cuFloatComplex operator*(float a, cuFloatComplex b) {
    return make_cuFloatComplex(a * b.x, a * b.y);
  } // operator*()

  __device__ static __inline__ cuDoubleComplex operator*(cuDoubleComplex a, cuDoubleComplex b) {
    return make_cuDoubleComplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
  } // operator*()

  __device__ static __inline__ cuDoubleComplex operator*(cuDoubleComplex a, double b) {
    return make_cuDoubleComplex(a.x * b, a.y * b);
  } // operator*()

  __device__ static __inline__ cuDoubleComplex operator*(double a, cuDoubleComplex b) {
    return make_cuDoubleComplex(a * b.x, a * b.y);
  } // operator*()


  // division

  __device__ static __inline__ cuFloatComplex operator/(cuFloatComplex a, cuFloatComplex b) {
    return cuCdivf(a, b);
  } // operator/()

  __device__ static __inline__ cuDoubleComplex operator/(cuDoubleComplex a, cuDoubleComplex b) {
    return cuCdiv(a, b);
  } // operator/()

  __device__ static __inline__ cuFloatComplex operator/(cuFloatComplex a, float b) {
    return make_cuFloatComplex(a.x / b, a.y / b);
  } // operator/()

  __device__ static __inline__ cuDoubleComplex operator/(cuDoubleComplex a, double b) {
    return make_cuDoubleComplex(a.x / b, a.y / b);
  } // operator/()

  __device__ static __inline__ cuFloatComplex operator/(float a, cuFloatComplex b) {
    //return cuCdivf(make_cuFloatComplex(a, 0.0f), b);
    float den = b.x * b.x + b.y * b.y;
    float den_inv = 1.0f / den;
    return make_cuFloatComplex((a * b.x) * den_inv, (a * b.y) * den_inv);
  } // operator/()

  __device__ static __inline__ cuFloatComplex cuCrcpf(cuFloatComplex b) {
    float den = __fadd_ru(__fmul_ru(b.x, b.x), __fmul_ru(b.y, b.y));
    float den_inv = __frcp_rn(den);
    return make_cuFloatComplex(__fmul_ru(b.x, den_inv), __fmul_ru(b.y, den_inv));
  } // operator/()

  __device__ static __inline__ cuDoubleComplex operator/(double a, cuDoubleComplex b) {
    return cuCdiv(make_cuDoubleComplex(a, 0.0), b);
  } // operator/()


  __device__ static __inline__ float cuC_abs(cuFloatComplex a) {
    return sqrtf(a.x * a.x + a.y * a.y);
  } // cuC_abs()

  __device__ static __inline__ double cuC_abs(cuDoubleComplex a) {
    return sqrt(a.x * a.x + a.y * a.y);
  } // cuC_abs()

  __device__ static __inline__ float cuC_norm(cuFloatComplex a){
    return (a.x * a.x + a.y * a.y);
  } // cuC_norm()

  __device__ static __inline__ double cuC_norm(cuDoubleComplex a){
    return (a.x * a.x + a.y * a.y);
  } // cuC_norm()

  // rotate the Q-vector float
  __device__ static __inline__ void rotate_q(float * rot, float qx, float qy, cuFloatComplex qz,
          cuFloatComplex & mqx, cuFloatComplex & mqy, cuFloatComplex & mqz) {
      mqx.x = rot[0] * qx + rot[1] * qy + rot[2] * qz.x; mqx.y = rot[2] * qz.y;
      mqy.x = rot[3] * qx + rot[4] * qy + rot[5] * qz.x; mqy.y = rot[5] * qz.y;
      mqz.x = rot[6] * qx + rot[7] * qy + rot[8] * qz.x; mqz.y = rot[8] * qz.y;
  }

  // rotate the Q-vector double
  __device__ static __inline__ void rotate_q(double * rot, double qx, double qy, 
          cuDoubleComplex qz, cuDoubleComplex & mqx, cuDoubleComplex & mqy, cuDoubleComplex & mqz){

      mqx.x = rot[0] * qx + rot[1] * qy + rot[2] * qz.x; mqx.y = rot[2] * qz.y;
      mqy.x = rot[3] * qx + rot[4] * qy + rot[5] * qz.x; mqy.y = rot[5] * qz.y;
      mqz.x = rot[6] * qx + rot[7] * qy + rot[8] * qz.x; mqz.y = rot[8] * qz.y;
  }

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

  __device__ static __inline__ cuFloatComplex cuCsqr(cuFloatComplex a) {
    return make_cuFloatComplex(a.x * a.x - a.y * a.y, 2.0f * a.x * a.y);
  } // cuCsqr()

  __device__ static __inline__ cuDoubleComplex cuCsqr(cuDoubleComplex a) {
    return make_cuDoubleComplex(a.x * a.x - a.y * a.y, 2.0 * a.x * a.y);
  } // cuCsqr()

  // eucledian norm
  __device__ static __inline__ float cuCnorm3(cuFloatComplex * a){
    float res = 0.f;
    for (int i = 0; i < 3; i++) res += cuC_norm(a[i]);
    return res;
  } // cuCnorm3()

  __device__ static __inline__ double cuCnorm3(cuDoubleComplex * a){
    float res = 0.;
    for (int i = 0; i < 3; i++) res += cuC_norm(a[i]);
    return res;
  }
  __device__ static __inline__ float cuCnorm3(cuFloatComplex a, cuFloatComplex b, cuFloatComplex c) {
    float a2 = cuC_norm(a);
    float b2 = cuC_norm(b);
    float c2 = cuC_norm(c);
    return (a2 + b2 + c2);
  } // cuCnorm3()

  __device__ static __inline__ double cuCnorm3(cuDoubleComplex a, cuDoubleComplex b, cuDoubleComplex c) {
    double a2 = cuC_norm(a);
    double b2 = cuC_norm(b);
    double c2 = cuC_norm(c);
    return (a2 + b2 + c2);
  } // cuCnorm3()

  __device__ static __inline__ float cuCnorm3(float a, float b, cuFloatComplex c) {
    float c2 = cuC_norm(c);
    return (a * a + b * b + c2);
  } // cuCnorm3()

  __device__ static __inline__ double cuCnorm3(double a, double b, cuDoubleComplex c) {
    double c2 = cuC_norm(c);
    return (a * a + b * b + c2);
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

  __device__ static __inline__ float cuCabs(cuFloatComplex z) {
    return hypotf(z.x, z.y);
  } // cuCabs()

//  __device__ static __inline__ double cuCabs(cuDoubleComplex z) {
//    return hypot(z.x, z.y);
//  } // cuCabs()

  __device__ static __inline__ float cuCarg(cuFloatComplex z) {
    return atan2f(z.y, z.x);
  } // cuCarg()

  __device__ static __inline__ double cuCarg(cuDoubleComplex z) {
    return atan2(z.y, z.x);
  } // cuCarg()

  __device__ static __inline__ cuFloatComplex cuClog(cuFloatComplex z) {
    return make_cuFloatComplex(log(cuCabs(z)), cuCarg(z));
  } // cuClog()

  __device__ static __inline__ cuDoubleComplex cuClog(cuDoubleComplex z) {
    return make_cuDoubleComplex(log(cuCabs(z)), cuCarg(z));
  } // cuClog()

  __device__ static __inline__ cuFloatComplex cuCpow(cuFloatComplex a, float p) {
    return cuCexp(cuClog(a) * p);;
  } // cuCpow()

  __device__ static __inline__ cuDoubleComplex cuCpow(cuDoubleComplex a, double p) {
    return cuCexp(cuClog(a) * p);;
  } // cuCpow()

  __device__ static __inline__ cuFloatComplex cuCpow(cuFloatComplex a, cuFloatComplex p) {
    return cuCexp(cuClog(a) * p);;
  } // cuCpow()

  __device__ static __inline__ cuDoubleComplex cuCpow(cuDoubleComplex a, cuDoubleComplex p) {
    return cuCexp(cuClog(a) * p);
  } // cuCpow()

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

  __device__ static __inline__ float cuCsinc(float v) {
      if (fabs(v) > 1.0E-05)
        return (sin(v)/v);
      else
        return (1. - v*v/6.f + v*v*v*v/120.f);
  }

  __device__ static __inline__ double cuCsinc(double v) {
      if (abs(v) > 1.0E-05)
        return (sin(v)/v);
      else
        return (1. - v*v/6. + v*v*v*v/120.);
  }

  __device__ static __inline__ cuFloatComplex cuCsinc(cuFloatComplex value) {
    cuFloatComplex temp;
    //if(fabsf(value.x) < 1e-14 && fabsf(value.y) < 1e-14) temp = make_cuFloatComplex(1.0, 0.0);
    if(fabsf(value.x) < 1e-5 && fabsf(value.y) < 1e-5) temp = make_cuFloatComplex(1.0, 0.0);
    else temp = cuCsin(value) / value;
    return temp;
  } // cuCsinc()

  __device__ static __inline__ cuDoubleComplex cuCsinc(cuDoubleComplex value) {
    cuDoubleComplex temp;
    //if(fabs(value.x) < 1e-14 && fabs(value.y) < 1e-14) temp = make_cuDoubleComplex(1.0, 0.0);
    if(fabs(value.x) < 1e-9 && fabs(value.y) < 1e-9) temp = make_cuDoubleComplex(1.0, 0.0);
    else temp = cuCsin(value) / value;
    return temp;
  } // cuCsinc()

  __device__ static __inline__ bool cuCiszero(cuFloatComplex value) {
    //return (fabsf(value.x) < 1e-14 && fabsf(value.y) < 1e-14);
    return (fabsf(value.x) < 1e-5 && fabsf(value.y) < 1e-5);
  } // cuCiszero()

  __device__ static __inline__ bool cuCiszero(cuDoubleComplex value) {
    //return (fabs(value.x) < 1e-14 && fabs(value.y) < 1e-14);
    return (fabs(value.x) < 1e-9 && fabs(value.y) < 1e-9);
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
/*  __device__ static __inline__ cuDoubleComplex cuCcbessj(cuDoubleComplex zz, int order) {
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
  } // cuCcbessj() */
  
  // temporary fix. assuming imaginary component is 0
  __device__ static __inline__ cuDoubleComplex cuCcbessj(cuDoubleComplex zz, int order) {
    return make_cuDoubleComplex(j1(cu_real(zz)), 0.0);
  } // cuCcbessj()
  
  __device__ static __inline__ cuFloatComplex cuCcbessj(cuFloatComplex zz, int order) {
    return make_cuFloatComplex(j1f(cu_real(zz)), (float) 0.0);
  } // cuCcbessj()
  
} // namespace

#endif // _CU_COMPLEX_NUMERIC_CUH_
