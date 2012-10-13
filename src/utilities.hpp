/***
  *  $Id: utilities.hpp 33 2012-08-06 16:22:01Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: utilities.hpp
  *  Created: Jun 25, 2012
  *  Modified: Fri 12 Oct 2012 12:34:34 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _UTILITIES_HPP_
#define _UTILITIES_HPP_


#include "globals.hpp"
#include "typedefs.hpp"

namespace hig {

	/**
	 * various min and max functions
	 */

	/** compute the minimum of a and b
	 * requires operator '>' to be defined for type_t
	 */
	template<typename type_t>
	type_t min(type_t a, type_t b) {
		return (a > b) ? b : a;
	} // min()

	/** compute the maximum of a and b
	 * requires operator '<' to be defined for type_t
	 */
	template<typename type_t>
	type_t max(type_t a, type_t b) {
		return (a < b) ? b : a;
	} // max()

	/** specialization for vector3_t
	 */
	template <>
	vector3_t min <vector3_t> (vector3_t a, vector3_t b);

	/** specialization for vector3_t
	 */
	template <>
	vector3_t max <vector3_t> (vector3_t a, vector3_t b);

	/** compute the minimum of a, b and c
	 * requires operator '>' to be defined for type_t
	 */
	template <typename type_t>
	type_t min(type_t a, type_t b, type_t c) {
		type_t d = min(a, b);
		return (c > d) ? d : c;
	} // min()

	/** compute the maximum of a, b and c
	 * requires operator '<' to be defined for type_t
	 */
	template<typename type_t>
	type_t max(type_t a, type_t b, type_t c) {
		type_t d = max(a, b);
		return (c < d) ? d : c;
	} // max()


	/**
	 * operators
	 */

	/** multiply each element of given matrix or vector by a scalar
	 * requires iterator to be defined
	 */
	template <typename scalar_t>		// how to restrict scalar_t to just scalars? ...
	std::vector<float_t>& operator*(scalar_t scalar, std::vector<float_t>& vec) {
		for(std::vector<float_t>::iterator i = vec.begin(); i != vec.end(); ++ i) {
			(*i) = (*i) * scalar;
		} // for
		return vec;
	} // operator*()

	/** comparison operator for two complex numbers
	 */
	extern bool operator<(complex_t a, complex_t b);

	/**
	 * complex operators
	 */

	extern complex_t operator*(complex_t c, float_t s);
	extern complex_t operator*(float_t s, complex_t c);
	extern complex_t operator*(float_t s, cucomplex_t c);
	extern complex_t operator*(cucomplex_t c, float_t s);
	extern complex_t operator*(complex_t s, cucomplex_t c);
	extern complex_t operator*(cucomplex_t c, complex_t s);

	extern std::complex<long double> operator*(std::complex<long double> c, long double s);

	extern complex_t operator+(float_t s, cucomplex_t c);
	extern complex_t operator+(cucomplex_t c, float_t s);
	extern complex_t operator+(complex_t s, cucomplex_t c);
	extern complex_t operator+(cucomplex_t c, complex_t s);

	/**
	 * matrix and vector operation functions
	 * use boost libs ...
	 */

	extern bool mat_log10_2d(unsigned int x_size, unsigned int y_size, float_t* &data);
	extern vector3_t floor(vector3_t a);

	extern complex_vec_t& mat_sqr(complex_vec_t&);
	extern bool mat_sqr(const complex_vec_t&, complex_vec_t&);
	extern bool mat_sqr_in(complex_vec_t&);
	extern complex_vec_t& mat_sqrt(complex_vec_t&);
	extern bool mat_sqrt(const complex_vec_t&, complex_vec_t&);
	extern bool mat_sqrt_in(complex_vec_t&);
	extern bool mat_exp(complex_vec_t& matrix, complex_vec_t& result);
	extern bool mat_exp_in(complex_vec_t& matrix);
	extern complex_vec_t& mat_besselj(int, unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern bool mat_besselj(int, unsigned int, unsigned int, unsigned int, const complex_vec_t&,
						complex_vec_t&);
	extern bool mat_besselj_in(int, unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern complex_t besselj(int, complex_t);

	extern complex_vec_t& mat_add(unsigned int, unsigned int, unsigned int, complex_vec_t&,
						unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern bool mat_add(unsigned int, unsigned int, unsigned int, const complex_vec_t&,
						unsigned int, unsigned int, unsigned int, const complex_vec_t&,
						complex_vec_t&);
	extern bool mat_add_in(unsigned int, unsigned int, unsigned int, complex_vec_t&,
						unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern complex_vec_t& mat_mul(float_t scalar, std::vector<complex_t>& matrix);
	extern complex_vec_t& mat_mul(complex_t scalar, std::vector<complex_t>& matrix);
	extern complex_vec_t& mat_mul(std::vector<complex_t>& matrix, float_t scalar);
	extern complex_vec_t& mat_mul(std::vector<complex_t>& matrix, complex_t scalar);
	extern bool mat_mul(float_t, const std::vector<complex_t>&, complex_vec_t&);
	extern bool mat_mul(complex_t, const std::vector<complex_t>&, complex_vec_t&);
	extern bool mat_mul(const std::vector<complex_t>&, float_t, complex_vec_t&);
	extern bool mat_mul(const std::vector<complex_t>&, complex_t, complex_vec_t&);
	extern bool mat_mul_in(float_t scalar, std::vector<complex_t>& matrix);
	extern bool mat_mul_in(complex_t scalar, std::vector<complex_t>& matrix);
	extern bool mat_mul_in(std::vector<complex_t>& matrix, float_t scalar);
	extern bool mat_mul_in(std::vector<complex_t>& matrix, complex_t scalar);
	extern complex_vec_t& mat_dot_prod(unsigned int, unsigned int, unsigned int, complex_vec_t&,
						unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern bool mat_dot_prod(unsigned int, unsigned int, unsigned int, const complex_vec_t&,
						unsigned int, unsigned int, unsigned int, const complex_vec_t&, complex_vec_t&);
	extern bool mat_dot_prod_in(unsigned int, unsigned int, unsigned int, complex_vec_t&,
						unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern complex_vec_t& mat_dot_div(unsigned int, unsigned int, unsigned int, complex_vec_t&,
						unsigned int, unsigned int, unsigned int, complex_vec_t&);
	extern bool mat_dot_div(unsigned int, unsigned int, unsigned int, const complex_vec_t&,
						unsigned int, unsigned int, unsigned int, const complex_vec_t&, complex_vec_t&);
	extern bool mat_dot_div_in(unsigned int, unsigned int, unsigned int, complex_vec_t&,
						unsigned int, unsigned int, unsigned int, complex_vec_t&);
	//extern std::vector<complex_t>& mat_sinc(unsigned int, unsigned int, unsigned int, complex_vec_t&);

	/** compute the transpose of a matrix
	 */
	extern bool transpose(unsigned int x_size, unsigned int y_size, const float_t *matrix, float_t* &transp);

	/** matrix multiplication for two 3x3 matrices
	 * operation is:
	 * x1 x2 x3   a1 a2 a3   d1 d2 d3
	 * y1 y2 y3 = b1 b2 b3 x e1 e2 e3
	 * z1 z2 z3   c1 c2 c3   f1 f2 f3
	 *
	 * use boost libs ... and make it general ...
	*/
	extern bool mat_mul_3x3(vector3_t a, vector3_t b, vector3_t c, vector3_t d, vector3_t e, vector3_t f,
					vector3_t& x, vector3_t& y, vector3_t& z);

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

	extern int count_naninfs(int, int, int, const complex_t*);
	extern int count_naninfs(int, int, int, const cucomplex_t*);
	extern int count_naninfs(int, int, int, const std::vector<complex_t>&);

} // namespace hig

#endif /* _UTILITIES_HPP_ */
