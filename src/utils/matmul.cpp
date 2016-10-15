/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: matmul.cpp
 *  Created: Jun 25, 2012
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

#include <utils/matmul.hpp>

namespace hig {


  /**
   * matrix multiplication for two 3x3 matrices
   * operation is:
   * x1 x2 x3   a1 a2 a3   d1 d2 d3
   * y1 y2 y3 = b1 b2 b3 x e1 e2 e3
   * z1 z2 z3   c1 c2 c3   f1 f2 f3
   * use boost libs ... and make it general ...
   */
  bool mat_mul_3x3(vector3_t a, vector3_t b, vector3_t c,
          vector3_t d, vector3_t e, vector3_t f,
          vector3_t& x, vector3_t& y, vector3_t& z) {
    real_t *A = new (std::nothrow) real_t[9];
    real_t *B = new (std::nothrow) real_t[9];
    real_t *C = new (std::nothrow) real_t[9];

    A[0] = a[0]; A[1] = a[1]; A[2] = a[2];
    A[3] = b[0]; A[4] = b[1]; A[5] = b[2];
    A[6] = c[0]; A[7] = c[1]; A[8] = c[2];
    B[0] = d[0]; B[1] = d[1]; B[2] = d[2];
    B[3] = e[0]; B[4] = e[1]; B[5] = e[2];
    B[6] = f[0]; B[7] = f[1]; B[8] = f[2];

    for(int i = 0; i < 3; i ++) {
      for(int j = 0; j < 3; j ++) {
        C[3 * i + j] = 0.0;
        for(int k = 0; k < 3; k ++) {
          C[3 * i + j] += A[3 * i + k] * B[3 * k + j];
        } // for k
      } // for j
    } // for i

    x[0] = C[0]; x[1] = C[1]; x[2] = C[2];
    y[0] = C[3]; y[1] = C[4]; y[2] = C[5];
    z[0] = C[6]; z[1] = C[7]; z[2] = C[8];

    /*x[0] = a[0] * d[0] + a[1] * e[0] + a[2] * f[0];
    x[1] = a[0] * d[1] + a[1] * e[1] + a[2] * f[1];
    x[2] = a[0] * d[2] + a[1] * e[2] + a[2] * f[2];
    y[0] = b[0] * d[0] + b[1] * e[0] + b[2] * f[0];
    y[1] = b[0] * d[1] + b[1] * e[1] + b[2] * f[1];
    y[2] = b[0] * d[2] + b[1] * e[2] + b[2] * f[2];
    z[0] = c[0] * d[0] + c[1] * e[0] + c[2] * f[0];
    z[1] = c[0] * d[1] + c[1] * e[1] + c[2] * f[1];
    z[2] = c[0] * d[2] + c[1] * e[2] + c[2] * f[2];
*/
    delete[] C;
    delete[] B;
    delete[] A;
    return true;
  } // mat_mul_3x3()

  bool mat_mul_3x3(const real_vec_t a, const real_vec_t d, real_vec_t& x) {
    const real_t *A = &a[0], *B = &d[0];
    real_t *C = &x[0];
    return mat_mul_3x3(A, B, C);
  } // mat_mul_3x3()

  bool mat_mul_3x3(const real_vec_t a, const real_t* d, real_t*& x) {
    const real_t *A = &(a[0]), *B = d;
    real_t *C = x;
    return mat_mul_3x3(A, B, C);
  } // mat_mul_3x3()

  bool mat_mul_3x3(const real_t* A, const real_t* B, real_t*& C) {
    for(int i = 0; i < 3; i ++) {
      for(int j = 0; j < 3; j ++) {
        C[3 * i + j] = 0.0;
        for(int k = 0; k < 3; k ++) {
          C[3 * i + j] += A[3 * i + k] * B[3 * k + j];
        } // for k
      } // for j
    } // for i
  } // mat_mul_3x3()


  /**
   * matrix vector product for matrix of size 3x3 and vector of size 1x3
   * operation is:
   * x1   a1 a2 a3   d1
   * x2 = b1 b2 b3 x d2
   * x3   c1 c2 c3   d3
   * note: transpose of d is used
   * use boost libs ...
   */
  bool mat_mul_3x1(vector3_t a, vector3_t b, vector3_t c, vector3_t d, vector3_t& x) {
    x[0] = a[0] * d[0] + a[1] * d[1] + a[2] * d[2];
    x[1] = b[0] * d[0] + b[1] * d[1] + b[2] * d[2];
    x[2] = c[0] * d[0] + c[1] * d[1] + c[2] * d[2];

    return true;
  } // mat_mul_3x1()

  /**
   * scalar-matrix multiplication    WRONG ...
   */
  /*std::vector<complex_t>& mat_mul(real_t scalar,
                  //unsigned int x_size, unsigned int y_size, unsigned int z_size,
                  std::vector<complex_t>& matrix) {
    std::vector<complex_t> result;
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back((*i) * scalar);
    } // for
    return result;
  } // mat_mul()

  std::vector<complex_t>& mat_mul(complex_t scalar,
                  //unsigned int x_size, unsigned int y_size, unsigned int z_size,
                  std::vector<complex_t>& matrix) {
    std::vector<complex_t> result;
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back((*i) * scalar);
    } // for
    return result;
  } // mat_mul()*/

  /*std::vector<complex_t>& mat_mul(std::vector<complex_t>& matrix, real_t scalar) {
    return mat_mul(scalar, matrix);
  } // mat_mul()

  std::vector<complex_t>& mat_mul(std::vector<complex_t>& matrix, complex_t scalar) {
    return mat_mul(scalar, matrix);
  } // mat_mul()*/

  /**
   * scalar-matrix multiplication into result
   */
  bool mat_mul(real_t scalar, const complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(std::vector<complex_t>::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back((*i) * scalar);
    } // for
    return true;
  } // mat_mul()

  bool mat_mul(complex_t scalar, const complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back((*i) * scalar);
    } // for
    return true;
  } // mat_mul()

  bool mat_mul(const complex_vec_t& matrix, real_t scalar, complex_vec_t& result) {
    return mat_mul(scalar, matrix, result);
  } // mat_mul()

  bool mat_mul(const complex_vec_t& matrix, complex_t scalar, complex_vec_t& result) {
    return mat_mul(scalar, matrix, result);
  } // mat_mul()

  /**
   * in-place scalar-matrix multiplication
   */
  bool mat_mul_in(real_t scalar,  complex_vec_t& matrix) {
    for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      *i = (*i) * scalar;
    } // for
    return true;
  } // mat_mul()

  bool mat_mul_in(complex_t scalar, complex_vec_t& matrix) {
    for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      *i = (*i) * scalar;
    } // for
    return true;
  } // mat_mul()

  bool mat_mul_in(complex_vec_t& matrix, real_t scalar) {
    return mat_mul_in(scalar, matrix);
  } // mat_mul()

  bool mat_mul_in(complex_vec_t& matrix, complex_t scalar) {
    return mat_mul_in(scalar, matrix);
  } // mat_mul()

} // namespace hig
