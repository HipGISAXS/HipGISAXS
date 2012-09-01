/***
  *  $Id: utilities.hpp 33 2012-08-06 16:22:01Z asarje $
  *
  *  Project: HipGISAXS
  *
  *  File: utilities.hpp
  *  Created: Jun 25, 2012
  *  Modified: Tue 28 Aug 2012 05:39:03 PM PDT
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

	extern complex_t operator+(float_t s, cucomplex_t c);
	extern complex_t operator+(cucomplex_t c, float_t s);
	extern complex_t operator+(complex_t s, cucomplex_t c);
	extern complex_t operator+(cucomplex_t c, complex_t s);

	/**
	 * matrix and vector operation functions
	 * use boost libs ...
	 */

	/** apply log10 to all elements of the 2D matrix
	 */
	extern bool mat_log10_2d(unsigned int x_size, unsigned int y_size, float_t* &data);

	/** compute floor for each element in a vector
	 */ 
	extern vector3_t floor(vector3_t a);


	extern std::vector<complex_t>& mat_sqr(//unsigned int, unsigned int, unsigned int,
											std::vector<complex_t>&);

	extern std::vector<complex_t>& mat_sqrt(//unsigned int, unsigned int, unsigned int,
											std::vector<complex_t>&);

	extern std::vector<complex_t>& mat_besselj(int, unsigned int, unsigned int, unsigned int,
											std::vector<complex_t>&);

	extern complex_t besselj(int, complex_t);

	extern std::vector<complex_t>& mat_add(unsigned int, unsigned int, unsigned int, std::vector<complex_t>&,
										unsigned int, unsigned int, unsigned int, std::vector<complex_t>&);

	/** mat_mul for a scalar and matrix
	 */
	extern std::vector<complex_t>& mat_mul(float_t scalar,
										//unsigned int x_size, unsigned int y_size, unsigned int z_size,
										std::vector<complex_t>& matrix);

	extern std::vector<complex_t>& mat_mul(complex_t scalar,
										//unsigned int x_size, unsigned int y_size, unsigned int z_size,
										std::vector<complex_t>& matrix);

	/** mat_mul for a scalar and matrix
	 */
	extern std::vector<complex_t>& mat_mul(//unsigned int x_size, unsigned int y_size, unsigned int z_size,
										std::vector<complex_t>& matrix, float_t scalar);

	extern std::vector<complex_t>& mat_mul(//unsigned int x_size, unsigned int y_size, unsigned int z_size,
										std::vector<complex_t>& matrix, complex_t scalar);

	/** matrix-matrix element-by-element product
	 */
	extern std::vector<complex_t>& mat_dot_prod(unsigned int, unsigned int, unsigned int,
												std::vector<complex_t>&,
												unsigned int, unsigned int, unsigned int,
												std::vector<complex_t>&);

	extern std::vector<complex_t>& mat_dot_div(unsigned int, unsigned int, unsigned int,
												std::vector<complex_t>&,
												unsigned int, unsigned int, unsigned int,
												std::vector<complex_t>&);

	/** a special compute function on a matrix
	 */
	extern std::vector<complex_t>& mat_sinc(unsigned int, unsigned int, unsigned int,
											std::vector<complex_t>& matrix);

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

	//extern std::vector<float_t>& mat_mul(float_t,
	//									unsigned int, unsigned int, unsigned int, std::vector<float_t>&);

	//extern std::vector<float_t>& mat_mul(unsigned int, unsigned int, unsigned int, std::vector<float_t>&,
	//									float_t);

} // namespace hig

#endif /* _UTILITIES_HPP_ */
