/***
  *  $Id: utilities.cpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS
  *
  *  File: utilities.cpp
  *  Created: Jun 25, 2012
  *  Modified: Mon 27 Aug 2012 10:44:20 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <cmath>

#include "utilities.hpp"

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
	bool mat_log10_2d(unsigned int x_size, unsigned int y_size, float_t* &data) {
		for(unsigned int i = 0; i < x_size * y_size; ++ i) {
			if(data[i] <= 0) {
				if(data[i] == 0) {
					std::cerr << "warning: matrix has a zero value. cannot calculate logarithm. keeping zero."
							<< std::endl;
					data[i] = 0;
					continue;
				} else {
					std::cerr << "error: matrix has a negative value. cannot calculate logarithm"
							<< std::endl;
					return false;
				} // if-else
			} // if
			data[i] = log10(data[i]);
		} // for
		return true;
	} // mat_log10()


	/**
	 * compute the transpose of a matrix
	 * use boost libs ...
	 */
	bool transpose(unsigned int x_size, unsigned int y_size, const float_t *matrix, float_t* &transp) {
		if(matrix == NULL) {
			return false;
		} // if

		transp = new (std::nothrow) float_t[x_size * y_size];

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
	bool mat_mul_3x3(vector3_t a, vector3_t b, vector3_t c,
					vector3_t d, vector3_t e, vector3_t f,
					vector3_t& x, vector3_t& y, vector3_t& z) {
		float_t *A = new (std::nothrow) float_t[9];
		float_t *B = new (std::nothrow) float_t[9];
		float_t *C = new (std::nothrow) float_t[9];

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

		delete[] C;
		delete[] B;
		delete[] A;
		return true;
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


	/*cucomplex_t operator*(float_t s, cucomplex_t c) {
		cucomplex_t temp;
		temp.x = s * c.x;
		temp.y = s * c.y;
		return temp;
	} // operator*()

	cucomplex_t operator*(cucomplex_t c, float_t s) {
		cucomplex_t temp;
		temp.x = s * c.x;
		temp.y = s * c.y;
		return temp;
	} // operator*()

	cucomplex_t operator*(complex_t s, cucomplex_t c) {
		cucomplex_t temp;
		temp.x = s.real() * c.x - s.imag() * c.y;
		temp.y = s.real() * c.y + s.imag() * c.x;
		return temp;
	} // operator*()

	cucomplex_t operator*(cucomplex_t c, complex_t s) {
		cucomplex_t temp;
		temp.x = s.real() * c.x - s.imag() * c.y;
		temp.y = s.real() * c.y + s.imag() * c.x;
		return temp;
	} // operator*() */

	complex_t operator*(complex_t c, float_t s) {
		return complex_t(c.real() * s, c.imag() * s);
	} // operator*()

	complex_t operator*(float_t s, complex_t c) {
		return complex_t(c.real() * s, c.imag() * s);
	} // operator*()

	complex_t operator*(float_t s, cucomplex_t c) {
		return complex_t(s * c.x, s * c.y);
	} // operator*()

	complex_t operator*(cucomplex_t c, float_t s) {
		return complex_t(s * c.x, s * c.y);
	} // operator*()

	complex_t operator*(complex_t s, cucomplex_t c) {
		return complex_t(s.real() * c.x - s.imag() * c.y, s.real() * c.y + s.imag() * c.x);
	} // operator*()

	complex_t operator*(cucomplex_t c, complex_t s) {
		return complex_t(s.real() * c.x - s.imag() * c.y, s.real() * c.y + s.imag() * c.x);
	} // operator*()

	complex_t operator+(float_t s, cucomplex_t c) {
		return complex_t(s + c.x, c.y);
	} // operator+()

	complex_t operator+(cucomplex_t c, float_t s) {
		return complex_t(s + c.x, c.y);
	} // operator+()

	complex_t operator+(complex_t s, cucomplex_t c) {
		return complex_t(s.real() + c.x, s.imag() + c.y);
	} // operator+()

	complex_t operator+(cucomplex_t c, complex_t s) {
		return complex_t(s.real() + c.x, s.imag() + c.y);
	} // operator+()

	/**
	 * returns element-by-element sum of two matrices
	 */
	std::vector<complex_t>& mat_add(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
									std::vector<complex_t>& matrix1,
									unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
									std::vector<complex_t>& matrix2) {
		if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
				|| matrix1.size() != matrix2.size()) {
			std::cerr << "error: matrix sizes are not the same for addition operation" << std::endl;
			return matrix1;
		} // if
		std::vector<complex_t>::iterator i1 = matrix1.begin();
		std::vector<complex_t>::iterator i2 = matrix2.begin();
		for(; i1 != matrix1.end(); ++ i1, ++ i2) {
			(*i1) = (*i1) + (*i2);
			//(*i1).x += (*i2).x;
			//(*i1).y += (*i2).y;
		} // for

		return matrix1;
	} // mat_add()


	/**
	 * scalar-matrix multiplication
	 */
	std::vector<complex_t>& mat_mul(float_t scalar,
									//unsigned int x_size, unsigned int y_size, unsigned int z_size,
									std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			(*i) = (*i) * scalar;
			//(*i).x = (*i).x * scalar;
			//(*i).y = (*i).y * scalar;
		} // for
		return matrix;
	} // mat_mul()


	std::vector<complex_t>& mat_mul(complex_t scalar,
									//unsigned int x_size, unsigned int y_size, unsigned int z_size,
									std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			(*i) = (*i) * scalar;
			//complex_t temp;
			//temp.x = (*i).x * (scalar).x - (*i).y * (scalar).y;
			//temp.y = (*i).x * (scalar).y + (*i).y * (scalar).x;
			//(*i).x = temp.x;
			//(*i).y = temp.y;
		} // for
		return matrix;
	} // mat_mul()


	std::vector<complex_t>& mat_mul(//unsigned int x_size, unsigned int y_size, unsigned int z_size,
									std::vector<complex_t>& matrix, float_t scalar) {
		return mat_mul(scalar, /*x_size, y_size, z_size,*/ matrix);
	} // mat_mul()


	std::vector<complex_t>& mat_mul(//unsigned int x_size, unsigned int y_size, unsigned int z_size,
									std::vector<complex_t>& matrix, complex_t scalar) {
		return mat_mul(scalar, /*x_size, y_size, z_size,*/ matrix);
	} // mat_mul()


	/**
	 * returns element-by-element product of two matrices
	 * it modifies the first matrix to store the result
	 */
	std::vector<complex_t>& mat_dot_prod(unsigned int x1_size, unsigned int y1_size, unsigned int z1_size,
										std::vector<complex_t>& matrix1,
										unsigned int x2_size, unsigned int y2_size, unsigned int z2_size,
										std::vector<complex_t>& matrix2) {
		if(x1_size != x2_size || y1_size != y2_size || z1_size != z2_size
				|| matrix1.size() != matrix2.size()) {
			std::cerr << "error: matrix sizes are not the same for dot product operation" << std::endl;
			return matrix1;
		} // if

		std::vector<complex_t>::iterator i1 = matrix1.begin();
		std::vector<complex_t>::iterator i2 = matrix2.begin();
		for(; i1 != matrix1.end(); ++ i1, ++ i2) {
			(*i1) = (*i1) * (*i2);
			//complex_t temp;
			//temp.x = (*i1).x * (*i2).x - (*i1).y * (*i2).y;
			//temp.y = (*i1).x * (*i2).y + (*i1).y * (*i2).x;
			//(*i1).x = temp.x;
			//(*i1).y = temp.y;
		} // for

		return matrix1;
	} // mat_dot_prod()


	std::vector<complex_t>& mat_dot_div(unsigned int nx1, unsigned int ny1, unsigned int nz1,
										std::vector<complex_t>& matrix1,
										unsigned int nx2, unsigned int ny2, unsigned int nz2,
										std::vector<complex_t>& matrix2) {
		if(nx1 != nx2 || ny1 != ny2 || nz1 != nz2 || matrix1.size() != matrix2.size()) {
			std::cerr << "error: matrix sizes are not the same for dot division operation"
						<< std::endl;
			return matrix1;
		} // if

		std::vector<complex_t>::iterator i1 = matrix1.begin();
		std::vector<complex_t>::iterator i2 = matrix2.begin();
		for(; i1 != matrix1.end(); ++ i1, ++ i2) {
			(*i1) = (*i1) / (*i2);
			//complex_t temp;
			//temp.x = ((*i1).x * (*i2).x + (*i1).y * (*i2).y) / ((*i2).x * (*i2).x + (*i2).y * (*i2).y);
			//temp.y = ((*i1).y * (*i2).x - (*i1).x * (*i2).y) / ((*i2).x * (*i2).x + (*i2).y * (*i2).y);
			//(*i1).x = temp.x;
			//(*i1).y = temp.y;
		} // for

		return matrix1;
	} // mat_dot_div()


	std::vector<complex_t>& mat_sqr(//unsigned int nx, unsigned int ny, unsigned int nz,
									std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			(*i) = (*i) * (*i);
			//complex_t temp;
			//temp.x = (*i).x * (*i).x - (*i).y * (*i).y;
			//temp.y = (*i).x * (*i).y + (*i).y * (*i).x;
			//(*i).x = temp.x;
			//(*i).y = temp.y;
		} // for

		return matrix;
	} // mat_sqr()


	std::vector<complex_t>& mat_sqrt(//unsigned int nx, unsigned int ny, unsigned int nz,
									std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			(*i) = sqrt(*i);
			//std::complex<float_t> temp((*i).x, (*i).y);
			//temp = sqrt(temp);
			//(*i).x = temp.real();
			//(*i).y = temp.imag();
		} // for

		return matrix;
	} // mat_sqrt()


	std::vector<complex_t>& mat_besselj(int j, unsigned int nx, unsigned int ny, unsigned int nz,
										std::vector<complex_t>& matrix) {
		// ...
	} // mat_besselj()


	float_t besselj(int j, complex_t m) {
		// ...
	} // besselj()


} // namespace hig
