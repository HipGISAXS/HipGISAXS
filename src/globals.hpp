/***
  *  $Id: globals.hpp 33 2012-08-06 16:22:01Z asarje $
  *
  *  Project: HipGISAXS - High Performance GISAXS
  *
  *  File: globals.hpp
  *  Created: Jun 05, 2012
  *  Modified: Mon 27 Aug 2012 10:38:54 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _GLOBALS_HPP_
#define _GLOBALS_HPP_

#include <boost/array.hpp>
#include <unordered_map>
#include <vector>
#include <cmath>

#include "tokens.hpp"
#include "typedefs.hpp"


namespace hig {

	const unsigned int LIGHT_SPEED_ = 3e+8;		/* speed of light in m/s */
	const float_t PI_ = 3.141592653589793;		/* PI correct upto 15 decimal places */

	const unsigned int MAX_DEPTH_ = 500;		/* maximum depth allowed */


	typedef struct vector2_t {
		boost::array <float_t, 2> vec_;

		/* constructors */

		vector2_t() {
			vec_[0] = 0; vec_[1] = 0;
		} // vector2_t()

		vector2_t(float_t a, float_t b) {
			vec_[0] = a; vec_[1] = b;
		} // vector2_t()

		vector2_t(vector2_t& a) {
			vec_[0] = a[0]; vec_[1] = a[1];
		} // vector2_t()

		vector2_t(const vector2_t& a) {
			vec_ = a.vec_;
		} // vector2_t()

		/* operators */

		vector2_t& operator=(const vector2_t& a) {
			vec_ = a.vec_;
			return *this;
		} // operator=

		vector2_t& operator=(vector2_t& a) {
			vec_ = a.vec_;
			return *this;
		} // operator=

		float_t& operator[](int i) {
			return vec_[i];
		} // operator[]
	} vector2_t;


	typedef struct vector3_t {
		boost::array <float_t, 3> vec_;

		/* constructors */

		vector3_t() {
			vec_[0] = 0; vec_[1] = 0; vec_[2] = 0;
		} // vector3_t()

		vector3_t(float_t a, float_t b, float_t c) {
			vec_[0] = a; vec_[1] = b; vec_[2] = c;
		} // vector3_t()

		vector3_t(vector3_t& a) {
			vec_[0] = a[0]; vec_[1] = a[1]; vec_[2] = a[2];
		} // vector3_t()

		vector3_t(const vector3_t& a) {
			vec_ = a.vec_;
		} // vector3_t()

		/* operators */

		vector3_t& operator=(const vector3_t& a) {
			vec_ = a.vec_;
			return *this;
		} // operator=

		vector3_t& operator=(vector3_t& a) {
			vec_ = a.vec_;
			return *this;
		} // operator=

		float_t& operator[](int i) {
			return vec_[i];
		} // operator[]

		vector3_t operator+(int toadd) {
			return vector3_t(vec_[0] + toadd, vec_[1] + toadd, vec_[2] + toadd);
		} // operator+()

		vector3_t operator+(vector3_t toadd) {
			return vector3_t(vec_[0] + toadd[0], vec_[1] + toadd[1], vec_[2] + toadd[2]);
		} // operator+()

		vector3_t operator-(int tosub) {
			return vector3_t(vec_[0] - tosub, vec_[1] - tosub, vec_[2] - tosub);
		} // operator-()

		vector3_t operator-(vector3_t tosub) {
			return vector3_t(vec_[0] - tosub[0], vec_[1] - tosub[1], vec_[2] - tosub[2]);
		} // operator-()

		vector3_t operator/(vector3_t todiv) {
			return vector3_t(vec_[0] / todiv[0], vec_[1] / todiv[1], vec_[2] / todiv[2]);
		} // operator/()
	} vector3_t;


	typedef struct matrix3x3_t {
		typedef boost::array<float_t, 3> mat3_t;
		mat3_t mat_[3];

		/* constructors */

		matrix3x3_t() {
			for(int i = 0; i < 3; ++ i)
				for(int j = 0; j < 3; ++ j)
					mat_[i][j] = 0.0;
		} // matrix3x3_t()

		matrix3x3_t(const matrix3x3_t& a) {
			mat_[0] = a.mat_[0];
			mat_[1] = a.mat_[1];
			mat_[2] = a.mat_[2];
		} // matrix3x3_t

		/* operators */

		mat3_t& operator[](unsigned int index) {
			return mat_[index];
		} // operator[]()

		mat3_t& operator[](int index) {
			return mat_[index];
		} // operator[]()

		matrix3x3_t& operator=(matrix3x3_t& a) {
			for(int i = 0; i < 3; ++ i)
				for(int j = 0; j < 3; ++ j)
					mat_[i][j] = a[i][j];
			return *this;
		} // operstor=()

		matrix3x3_t operator+(int toadd) {
			matrix3x3_t sum;
			for(int i = 0; i < 3; ++ i)
				for(int j = 0; j < 3; ++ j)
					sum.mat_[i][j] = mat_[i][j] + toadd;
			return sum;
		} // operator+()

		matrix3x3_t operator+(matrix3x3_t& toadd) {
			matrix3x3_t sum;
			for(int i = 0; i < 3; ++ i)
				for(int j = 0; j < 3; ++ j)
					sum.mat_[i][j] = mat_[i][j] + toadd[i][j];
			return sum;
		} // operator+()

		matrix3x3_t operator*(int tomul) {
			matrix3x3_t prod;
			for(int i = 0; i < 3; ++ i)
				for(int j = 0; j < 3; ++ j)
					prod.mat_[i][j] = mat_[i][j] * tomul;
			return prod;
		} // operator*()

		matrix3x3_t operator*(matrix3x3_t& tomul) {
			matrix3x3_t prod;
			for(int i = 0; i < 3; ++ i)
				for(int j = 0; j < 3; ++ j)
					prod.mat_[i][j] = mat_[i][j] * tomul[i][j];
			return prod;
		} // operator*()

	} matrix3x3_t;

} // namespace hig

#endif /* _GLOBALS_HPP_ */
