/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: matmul.hpp
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

#ifndef __MATMUL_HPP__
#define __MATMUL_HPP__


#include <common/globals.hpp>
#include <common/typedefs.hpp>

namespace hig {

  extern complex_vec_t& mat_mul(real_t scalar, std::vector<complex_t>& matrix);
  extern complex_vec_t& mat_mul(complex_t scalar, std::vector<complex_t>& matrix);
  extern complex_vec_t& mat_mul(std::vector<complex_t>& matrix, real_t scalar);
  extern complex_vec_t& mat_mul(std::vector<complex_t>& matrix, complex_t scalar);
  extern bool mat_mul(real_t, const std::vector<complex_t>&, complex_vec_t&);
  extern bool mat_mul(complex_t, const std::vector<complex_t>&, complex_vec_t&);
  extern bool mat_mul(const std::vector<complex_t>&, real_t, complex_vec_t&);
  extern bool mat_mul(const std::vector<complex_t>&, complex_t, complex_vec_t&);
  extern bool mat_mul_in(real_t scalar, std::vector<complex_t>& matrix);
  extern bool mat_mul_in(complex_t scalar, std::vector<complex_t>& matrix);
  extern bool mat_mul_in(std::vector<complex_t>& matrix, real_t scalar);
  extern bool mat_mul_in(std::vector<complex_t>& matrix, complex_t scalar);

  /** matrix multiplication for two 3x3 matrices
   * operation is:
   * x1 x2 x3   a1 a2 a3   d1 d2 d3
   * y1 y2 y3 = b1 b2 b3 x e1 e2 e3
   * z1 z2 z3   c1 c2 c3   f1 f2 f3
   *
   * use boost libs ... and make it general ...
  */
  extern bool mat_mul_3x3(vector3_t a, vector3_t b, vector3_t c, vector3_t d, vector3_t e, vector3_t f, vector3_t& x, vector3_t& y, vector3_t& z);
  extern bool mat_mul_3x3(const real_vec_t a, const real_vec_t d, real_vec_t& x);
  extern bool mat_mul_3x3(const real_vec_t a, const real_t* d, real_t*& x);
  extern bool mat_mul_3x3(const real_t* a, const real_t* d, real_t*& x);

  /** matrix vector product for matrix of size 3x3 and vector of size 1x3
   * operation is:
   * x1   a1 a2 a3   d1
   * x2 = b1 b2 b3 x d2
   * x3   c1 c2 c3   d3
   * note: transpose of d is used
   *
   * use boost libs ...
   */
  extern bool mat_mul_3x1(vector3_t a, vector3_t b, vector3_t c, vector3_t d, vector3_t& x);

} // namespace hig

#endif /* __MATMUL_HPP__ */
