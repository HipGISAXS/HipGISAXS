/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: utilities.cpp
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

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>
//#include <boost/math/special_functions/bessel.hpp>
//#include <pari/pari.h>  // for bessel functions

#include <common/constants.hpp>
#include <utils/utilities.hpp>
#include <numerics/numeric_utils.hpp>

namespace hig {


  /**
   * specialization for vector3_t
   */
  template <>
  vector3_t min <vector3_t> (vector3_t a, vector3_t b) {
    return vector3_t((a[0] < b[0] ? a[0] : b[0]), (a[1] < b[1] ? a[1] : b[1]),
            (a[2] < b[2] ? a[2] : b[2]));
  } // min <vector3_t>

  /**
   * specialization for vector3_t
   */
  template <>
  vector3_t max <vector3_t> (vector3_t a, vector3_t b) {
    return vector3_t((a[0] > b[0] ? a[0] : b[0]), (a[1] > b[1] ? a[1] : b[1]),
            (a[2] > b[2] ? a[2] : b[2]));
  } // max <vector3_t>


  /**
   * apply log10 to all elements of the 2D matrix
   */
  bool mat_log10_2d(unsigned int x_size, unsigned int y_size, real_t* &data) {
    if(data == NULL) {
      std::cerr << "error: data is null while calculating log10" << std::endl;
      return false;
    } // if
    for(unsigned int i = 0; i < x_size * y_size; ++ i) {
      if(data[i] <= 0) {
        if(data[i] == 0) {
          //std::cerr << "warning: matrix has a zero value. "
          //    << "cannot calculate logarithm. keeping zero."
          //    << std::endl;
          data[i] = 0.0;
          continue;
        } else {
          std::cerr << "error: matrix has a negative value. cannot calculate logarithm"
              << std::endl;
          return false;
        } // if-else
      } else
        data[i] = log10(data[i]);
    } // for
    return true;
  } // mat_log10()


  /**
   * compute the transpose of a matrix
   * use boost libs ...
   */
  bool transpose(unsigned int x_size, unsigned int y_size, const real_t *matrix, real_t* &transp) {
    if(matrix == NULL) {
      std::cerr << "error: matrix is NULL while tranposing" << std::endl;
      return false;
    } // if
    transp = new (std::nothrow) real_t[x_size * y_size];
    for(unsigned int y = 0; y < y_size; ++ y) {
      for(unsigned int x = 0; x < x_size; ++ x) {
        transp[y_size * x + y] = matrix[x_size * y + x];
      } // for x
    } // for y

    return true;
  } // transpose()


  /**
   * matrix multiplication for two 3x3 matrices
   * operation is:
   * x1 x2 x3   a1 a2 a3   d1 d2 d3
   * y1 y2 y3 = b1 b2 b3 x e1 e2 e3
   * z1 z2 z3   c1 c2 c3   f1 f2 f3
   * use boost libs ... and make it general ...
   */
/*  bool mat_mul_3x3(vector3_t a, vector3_t b, vector3_t c,
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
*/
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
/*    delete[] C;
    delete[] B;
    delete[] A;
    return true;
  } // mat_mul_3x3()

  bool mat_mul_3x3(real_vec_t a, real_vec_t d, real_vec_t& x) {
    real_t *A = &(a[0]), *B = &(d[0]), *C = &(x[0]);
    return mat_mul_3x3(A, B, C);
  } // mat_mul_3x3()

  bool mat_mul_3x3(real_vec_t a, real_t* d, real_t*& x) {
    real_t *A = &(a[0]), *B = d, *C = x;
    return mat_mul_3x3(A, B, C);
  } // mat_mul_3x3()

  bool mat_mul_3x3(real_t* A, real_t* B, real_t*& C) {
    for(int i = 0; i < 3; i ++) {
      for(int j = 0; j < 3; j ++) {
        C[3 * i + j] = 0.0;
        for(int k = 0; k < 3; k ++) {
          C[3 * i + j] += A[3 * i + k] * B[3 * k + j];
        } // for k
      } // for j
    } // for i
  } // mat_mul_3x3()
*/

  /**
   * matrix vector product for matrix of size 3x3 and vector of size 1x3
   * operation is:
   * x1   a1 a2 a3   d1
   * x2 = b1 b2 b3 x d2
   * x3   c1 c2 c3   d3
   * note: transpose of d is used
   * use boost libs ...
   */
/*  bool mat_mul_3x1(vector3_t a, vector3_t b, vector3_t c, vector3_t d, vector3_t& x) {
    x[0] = a[0] * d[0] + a[1] * d[1] + a[2] * d[2];
    x[1] = b[0] * d[0] + b[1] * d[1] + b[2] * d[2];
    x[2] = c[0] * d[0] + c[1] * d[1] + c[2] * d[2];

    return true;
  } // mat_mul_3x1()
*/

  /**
   * specialized floor function
   */
  vector3_t floor(vector3_t a) {
    return vector3_t(std::floor(a[0]), std::floor(a[1]), std::floor(a[2]));
  } // floor()


  /**
   * comparison of complex numbers
   */
  bool operator<(complex_t a, complex_t b) {
    if(std::real(a) < std::real(b)) return true;
    if(std::real(a) == std::real(b) && std::imag(a) < std::imag(b)) return true;
    //if(a.x < b.x) return true;
    //if(a.x == b.x && a.y < b.y) return true;
    return false;
  } // operator<()


  /**
   * arithmetic operators for complex types
   */

#if (!(defined(__ICC) || defined(__INTEL_COMPILER))) || \
    (defined(__ICC) || defined(__INTEL_COMPILER)) && \
    (__ICC < 1600 || __INTEL_COMPILER < 1600)
  complex_t operator*(complex_t c, real_t s) {
    return complex_t(c.real() * s, c.imag() * s);
  } // operator*()

  complex_t operator*(real_t s, complex_t c) {
    return complex_t(c.real() * s, c.imag() * s);
  } // operator*()
#endif

  std::complex<long double> operator*(std::complex<long double> c, long double s) {
    return std::complex<long double>(c.real() * s, c.imag() * s);
  } // operator*()

  #ifdef USE_GPU

    complex_t operator*(float2 c, float2 s) {
      return complex_t(c.x * s.x - c.y * s.y, c.x * s.y + c.y * s.x);
    } // operator*()

    complex_t operator*(cucomplex_t c, real_t s) {
      return complex_t(s * c.x, s * c.y);
    } // operator*()

    complex_t operator*(real_t s, cucomplex_t c) {
      return complex_t(s * c.x, s * c.y);
    } // operator*()

    complex_t operator*(complex_t s, cucomplex_t c) {
      return complex_t(s.real() * c.x - s.imag() * c.y, s.real() * c.y + s.imag() * c.x);
    } // operator*()

    complex_t operator*(cucomplex_t c, complex_t s) {
      return complex_t(s.real() * c.x - s.imag() * c.y, s.real() * c.y + s.imag() * c.x);
    } // operator*()

    complex_t operator+(complex_t s, cucomplex_t c) {
      return complex_t(s.real() + c.x, s.imag() + c.y);
    } // operator+()

    complex_t operator+(cucomplex_t c, complex_t s) {
      return complex_t(s.real() + c.x, s.imag() + c.y);
    } // operator+()

    complex_t operator+(cucomplex_t c, real_t s) {
      return complex_t(s + c.x, c.y);
    } // operator+()

    complex_t operator+(real_t s, cucomplex_t c) {
      return complex_t(s + c.x, c.y);
    } // operator+()

  #endif

  std::complex<long double> operator/(std::complex<long double> c, long double s) {
    return std::complex<long double>(c.real() / s, c.imag() / s);
  } // operator/()


  /**
   * returns element-by-element sum of two matrices  WRONG ...
   */
  /*std::vector<complex_t>& mat_add(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
                  std::vector<complex_t>& matrix1,
                  unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
                  std::vector<complex_t>& matrix2) {
    if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
        || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for addition operation" << std::endl;
      return matrix1;
    } // if
    std::vector<complex_t> result;
    std::vector<complex_t>::iterator i1 = matrix1.begin();
    std::vector<complex_t>::iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      result.push_back((*i1) + (*i2));
    } // for

    return result;
  } // mat_add()*/

  /**
   * constructs element-by-element sum of two matrices into result
   */
  bool mat_add(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
          const std::vector<complex_t>& matrix1,
          unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
          const std::vector<complex_t>& matrix2,
          std::vector<complex_t>& result) {
    if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
        || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for addition operation" << std::endl;
      return false;
    } // if
    result.clear();
    std::vector<complex_t>::const_iterator i1 = matrix1.begin();
    std::vector<complex_t>::const_iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      result.push_back((*i1) + (*i2));
    } // for
    return true;
  } // mat_add()

  /**
   * performs in-place element-by-element sum of two matrices into first matrix
   */
  bool mat_add_in(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
          std::vector<complex_t>& matrix1,
          unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
          std::vector<complex_t>& matrix2) {
    if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
        || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for addition operation" << std::endl;
      return false;
    } // if
    std::vector<complex_t>::iterator i1 = matrix1.begin();
    std::vector<complex_t>::iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      *i1 = (*i1) + (*i2);
    } // for
    return true;
  } // mat_add()


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
/*  bool mat_mul(real_t scalar, const complex_vec_t& matrix, complex_vec_t& result) {
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
*/
  /**
   * in-place scalar-matrix multiplication
   */
/*  bool mat_mul_in(real_t scalar,  complex_vec_t& matrix) {
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
*/

  /**
   * returns element-by-element product of two matrices    WRONG
   */
  /*std::vector<complex_t>& mat_dot_prod(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
                    std::vector<complex_t>& matrix1,
                    unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
                    std::vector<complex_t>& matrix2) {
    if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
        || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for dot product operation" << std::endl;
      return matrix1;
    } // if

    std::vector<complex_t> result;
    std::vector<complex_t>::iterator i1 = matrix1.begin();
    std::vector<complex_t>::iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      result.push_back((*i1) * (*i2));
    } // for

    return result;
  } // mat_dot_prod()*/

  /**
   * computes element-by-element product of two matrices into result
   */
  bool mat_dot_prod(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
            const complex_vec_t& matrix1,
            unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
            const complex_vec_t& matrix2,
            complex_vec_t& result) {
    if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
        || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for dot product operation" << std::endl;
      return false;
    } // if
    result.clear();
    complex_vec_t::const_iterator i1 = matrix1.begin();
    complex_vec_t::const_iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      result.push_back((*i1) * (*i2));
    } // for
    return true;
  } // mat_dot_prod()

  /**
   * performs in-place element-by-element product of two matrices
   */
  bool mat_dot_prod_in(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
              std::vector<complex_t>& matrix1,
              unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
              std::vector<complex_t>& matrix2) {
    if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
        || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for dot product operation" << std::endl;
      return false;
    } // if
    std::vector<complex_t>::iterator i1 = matrix1.begin();
    std::vector<complex_t>::iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      *i1 = (*i1) * (*i2);
    } // for
    return true;
  } // mat_dot_prod()


  /**
   * returns element-by-element division of two matrices (matrix1 / matrix2)  WRONG ...
   */
  /*std::vector<complex_t>& mat_dot_div(unsigned int nx1, unsigned int ny1, unsigned int nz1,
                    std::vector<complex_t>& matrix1,
                    unsigned int nx2, unsigned int ny2, unsigned int nz2,
                    std::vector<complex_t>& matrix2) {
    if(nx1 != nx2 || ny1 != ny2 || nz1 != nz2 || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for dot division operation"
            << std::endl;
      return matrix1;
    } // if
    std::vector<complex_t> result;
    std::vector<complex_t>::iterator i1 = matrix1.begin();
    std::vector<complex_t>::iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      result.push_back((*i1) / (*i2));
    } // for
    return result;
  } // mat_dot_div()*/

  /**
   * computes element-by-element division of two matrices (matrix1 / matrix2)  into result
   */
  bool mat_dot_div(unsigned int nx1, unsigned int ny1, unsigned int nz1,
            const complex_vec_t& matrix1,
            unsigned int nx2, unsigned int ny2, unsigned int nz2,
            const complex_vec_t& matrix2,
            complex_vec_t& result) {
    if(nx1 != nx2 || ny1 != ny2 || nz1 != nz2 || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for dot division operation"
            << std::endl;
      return false;
    } // if
    result.clear();
    complex_vec_t::const_iterator i1 = matrix1.begin();
    complex_vec_t::const_iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      result.push_back((*i1) / (*i2));
    } // for
    return true;
  } // mat_dot_div()

  /**
   * performs in-place element-by-element division of two matrices (matrix1/matrix2) into matrix1
   */
  bool mat_dot_div_in(unsigned int nx1, unsigned int ny1, unsigned int nz1,
            std::vector<complex_t>& matrix1,
            unsigned int nx2, unsigned int ny2, unsigned int nz2,
            std::vector<complex_t>& matrix2) {
    if(nx1 != nx2 || ny1 != ny2 || nz1 != nz2 || matrix1.size() != matrix2.size()) {
      std::cerr << "error: matrix sizes are not the same for dot division operation"
            << std::endl;
      return false;
    } // if
    std::vector<complex_t>::iterator i1 = matrix1.begin();
    std::vector<complex_t>::iterator i2 = matrix2.begin();
    for(; i1 != matrix1.end(); ++ i1, ++ i2) {
      *i1 = (*i1) / (*i2);
    } // for
    return true;
  } // mat_dot_div()


  /*std::vector<complex_t>& mat_sqr(//unsigned int nx, unsigned int ny, unsigned int nz,
                  std::vector<complex_t>& matrix) {
    std::vector<complex_t> result;
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back((*i) * (*i));
    } // for
    return result;
  } // mat_sqr()*/

  bool mat_sqr(const complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back((*i) * (*i));
    } // for
    return true;
  } // mat_sqr()

  bool mat_sqr_in(//unsigned int nx, unsigned int ny, unsigned int nz,
            std::vector<complex_t>& matrix) {
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      *i = (*i) * (*i);
    } // for
    return true;
  } // mat_sqr()


  /*std::vector<complex_t>& mat_sqrt(//unsigned int nx, unsigned int ny, unsigned int nz,
                  std::vector<complex_t>& matrix) {
    std::vector<complex_t> result;
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back(sqrt(*i));
    } // for
    return result;
  } // mat_sqrt()*/

  bool mat_sqrt(const complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back(sqrt(*i));
    } // for
    return true;
  } // mat_sqrt()

  bool mat_sqrt_in(//unsigned int nx, unsigned int ny, unsigned int nz,
            std::vector<complex_t>& matrix) {
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      *i = sqrt(*i);
    } // for
    return true;
  } // mat_sqrt()


  bool mat_exp(complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i)
      result.push_back(exp(*i));
    return true;
  } // mat_exp()


  bool mat_exp_in(complex_vec_t& matrix) {
    for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i) *i = exp(*i);
    return true;
  } // mat_exp_in()


  /*std::vector<complex_t>& mat_besselj(int j, unsigned int nx, unsigned int ny, unsigned int nz,
                    std::vector<complex_t>& matrix) {
    std::vector<complex_t> result;
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back(cbessj(*i, j));
    } // for
    return result;
  } // mat_besselj()*/

  bool mat_besselj(int j, unsigned int nx, unsigned int ny, unsigned int nz,
            const complex_vec_t& matrix, complex_vec_t& result) {
    result.clear();
    for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      result.push_back(cbessj(*i, j));
    } // for
    return true;
  } // mat_besselj()

  bool mat_besselj_in(int j, unsigned int nx, unsigned int ny, unsigned int nz,
                    std::vector<complex_t>& matrix) {
    for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
      //*i = besselj(j, *i);
      *i = cbessj(*i, j);
    } // for
    return true;
  } // mat_besselj()


  // pari
  /*complex_t besselj(int j, complex_t m) {
    pari_init(5000000, 100000);
    GEN inp_m;
    inp_m = cgetc(64);  // using double precision
    gel(inp_m, 1) = dbltor(m.real());
    gel(inp_m, 2) = dbltor(m.imag());
    GEN unit = cgetc(64);
    gel(unit, 1) = stoi(j);
    gel(unit, 2) = stoi(0);

    GEN out_res = cgetc(64);
    out_res = jbessel(unit, inp_m, 64);

    real_t out_real = rtodbl(gel(out_res, 1));
    real_t out_imag = rtodbl(gel(out_res, 2));

    pari_close();

    return complex_t(out_real, out_imag);
  } // besselj()*/


/*  extern "C" {
    void zbesj_(double* real, double* imag, double* fnu, int* kode, int* n,
          double** creal, double** cimag, int* nz, int* err);
  } // extern "C"

  complex_t besselj(int j, complex_t m) {
    // amos
    double real = m.real(), imag = m.imag(), fnu = 1;
    int n = 1, kode = 1;
    double *creal = new double[n];
    double *cimag = new double[n];
    int nz = 0, err = 0;

    zbesj_(&real, &imag, &fnu, &kode, &n, &creal, &cimag, &nz, &err);
    if(err != 0) {
      std::cerr << "error: something went wrong in computing bessel function of "
            << m << ". error code = " << err << std::endl;
      delete[] creal;
      delete[] cimag;
      return complex_t(0, 0);
    } // if
    complex_t result((real_t)creal[0], (real_t)cimag[0]);

    delete[] creal;
    delete[] cimag;
    
    return result; */
    /* // boost
    return boost::math::cyl_bessel_j(j, m);
  }*/

  // adding two data sets into one
  bool add_data_elements(real_t* &dst, const real_t* src1, const real_t* src2, int size) {
    if(dst == NULL || src1 == NULL || src2 == NULL) {
      std::cerr << "error: null pointer in adding data elements." << std::endl;
      return false;
    } // if
    #pragma omp parallel for
    for(int i = 0; i < size; ++ i) dst[i] = src1[i] + src2[i];
  } // add_data_elements()


  int count_naninfs(int nx, int ny, int nz, const complex_t* arr) {
    int count = 0;
    for(unsigned int i = 0; i < nx * ny * nz; ++ i) {
      if(!(boost::math::isfinite(arr[i].real()) && boost::math::isfinite(arr[i].imag()))) ++ count;
    } // for
    return count;
  } // count_naninfs()

  int count_naninfs(int nx, int ny, int nz, const std::vector<complex_t>& arr) {
    int count = 0;
    for(std::vector<complex_t>::const_iterator i = arr.begin(); i != arr.end(); ++ i) {
      if(!(boost::math::isfinite((*i).real()) && boost::math::isfinite((*i).imag()))) ++ count;
    } // for
    return count;
  } // count_naninfs()

  #ifdef USE_GPU
    int count_naninfs(int nx, int ny, int nz, const cucomplex_t* arr) {
      int count = 0;
      for(unsigned int i = 0; i < nx * ny * nz; ++ i) {
        if(!(boost::math::isfinite(arr[i].x) && boost::math::isfinite(arr[i].y))) ++ count;
      } // for
      return count;
    } // count_naninfs()
  #endif


  // compute integral of e^(ikx) between x1 and x2
  complex_t integral_e(real_t x1, real_t x2, complex_t k) {
    if(boost::math::fpclassify(k.real()) == FP_ZERO &&
        boost::math::fpclassify(k.imag()) == FP_ZERO) {
      return complex_t(x2 - x1, 0);
    } else {
      complex_t ik = complex_t(0.0, 1.0) * k;
      return (((real_t) 1.0 / ik) * (exp(ik * x2) - exp(ik * x1)));
    } // if-else
  } // integral_e

  // compute integral of (ax + b) e^(ikx) between x1 and x2
  complex_t integral_xe(real_t x1, real_t x2, real_t a, real_t b, complex_t k) {
    if(boost::math::fpclassify(k.real()) == FP_ZERO &&
        boost::math::fpclassify(k.imag()) == FP_ZERO) {
      return complex_t(a * (x2 * x2 - x1 * x1) / 2 + b * (x2 - x1), 0.0);
    } else {
      complex_t ik = complex_t(0.0, 1.0) * k;
      return (((real_t) 1.0 / ik) * ((a * x2 + b - a / ik) * exp(ik * x2) -
                      (a * x1 + b - a / ik) * exp(ik * x1)));
    } // if-else
  } // integral_xe()


  /**
   * some more operations on complex numbers
   */

  //inline real_t magnitude(complex_t z) {
  //  return sqrt(z.real() * z.real() + z.imag() * z.imag());
  //} // magnitude()

  inline bool normalize(complex_t* x, real_t* xn, int size) {
    if(x == NULL) return false;
    if(xn == NULL) xn = new (std::nothrow) real_t[size];

    real_t max = magnitude(x[0]);
    for(int i = 0; i < size; ++ i) {
      xn[i] = magnitude(x[i]);
      if(xn[i] > max) max = xn[i];
    } // for

    for(int i = 0; i < size; ++ i) xn[i] = xn[i] / max;

    return true;
  } // normalize()

  inline bool conjugate(complex_t* x, int size) {
    if(x == NULL) return false;
    for(int i = 0; i < size; ++ i) x[i] = complex_t(x[i].real(), -x[i].imag());
    return true;
  } // conjugate()


  inline real_t gaussian(real_t x, real_t y, real_t mux, real_t muy, real_t sig, bool isnorm) {
    real_t norm = -1.0;
    if(isnorm) norm = sig * sqrt(2 * PI_);
    return exp(-((x - mux) * (x - mux) + (y - muy) * (y - muy)) / (2 * sig * sig)) / norm;
  } // gaussian()

} // namespace hig
