/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: globals.hpp
 *  Created: Jun 05, 2012
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

#ifndef _GLOBALS_HPP_
#define _GLOBALS_HPP_

#include <boost/array.hpp>
#include <vector>
#include <cmath>

#include <common/typedefs.hpp>
#include <common/cudafy.hpp>

namespace hig {

  typedef struct vector2_t {
    boost::array <real_t, 2> vec_;

    /* constructors */

    vector2_t() {
      vec_[0] = 0; vec_[1] = 0;
    } // vector2_t()

    vector2_t(real_t a, real_t b) {
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

    real_t& operator[](int i) {
      return vec_[i];
    } // operator[]
  } vector2_t;


  typedef struct vector3_t {
    real_t vec_[3];

    /* constructors */

    CUDAFY vector3_t() {
      vec_[0] = 0; vec_[1] = 0; vec_[2] = 0;
    } // vector3_t()

    CUDAFY vector3_t(real_t a, real_t b, real_t c) {
      vec_[0] = a; vec_[1] = b; vec_[2] = c;
    } // vector3_t()

    INLINE vector3_t(vector3_t& a) {
      vec_[0] = a.vec_[0]; 
      vec_[1] = a.vec_[1]; 
      vec_[2] = a.vec_[2];
    } // vector3_t()

    INLINE vector3_t(const std::vector<real_t>& a) {
      if(a.size() != 3) return;
      vec_[0] = a[0]; 
      vec_[1] = a[1]; 
      vec_[2] = a[2];
    } // vector3_t()

    CUDAFY vector3_t(const vector3_t& a) {
      vec_[0] = a.vec_[0]; 
      vec_[1] = a.vec_[1]; 
      vec_[2] = a.vec_[2]; 
    } // vector3_t()

    INLINE real_t norm(){
      real_t norm = 0.;
      for (int i = 0; i < 3; i++) norm += vec_[i] * vec_[i];
      return norm;
    }

    INLINE real_t abs(){
      real_t v = this->norm();
#ifdef DOUBLEP
      return sqrt(v);
#else
      return sqrtf(v);
#endif
    }

    /* operators */
    INLINE vector3_t& operator=(const vector3_t& a) {
      vec_[0] = a.vec_[0]; 
      vec_[1] = a.vec_[1]; 
      vec_[2] = a.vec_[2];
      return *this;
    } // operator=

    INLINE vector3_t& operator=(vector3_t& a) {
      vec_[0] = a.vec_[0]; 
      vec_[1] = a.vec_[1]; 
      vec_[2] = a.vec_[2]; 
      return *this;
    } // operator=

    INLINE vector3_t& operator+=(vector3_t& a) {
      vec_[0] += a.vec_[0]; 
      vec_[1] += a.vec_[1]; 
      vec_[2] += a.vec_[2]; 
      return *this;
    } // operator=

    INLINE real_t & operator[](int i) {
      return vec_[i];
    } // operator[]

    INLINE real_t operator[](int i) const {
      return vec_[i];
    } // read_only

    INLINE vector3_t operator+(int toadd) {
      return vector3_t(vec_[0] + toadd, vec_[1] + toadd, vec_[2] + toadd);
    } // operator+()

    INLINE vector3_t operator+(vector3_t toadd) {
      return vector3_t(vec_[0] + toadd[0], vec_[1] + toadd[1], vec_[2] + toadd[2]);
    } // operator+()

    INLINE vector3_t operator-(int tosub) {
      return vector3_t(vec_[0] - tosub, vec_[1] - tosub, vec_[2] - tosub);
    } // operator-()

    INLINE vector3_t operator-(vector3_t tosub) {
      return vector3_t(vec_[0] - tosub[0], vec_[1] - tosub[1], vec_[2] - tosub[2]);
    } // operator-()

    INLINE vector3_t operator*(real_t tomul) {
      return vector3_t(vec_[0] * tomul, vec_[1] * tomul, vec_[2] * tomul);
    } // operator*()

    INLINE vector3_t operator*(vector3_t rhs) {
      return vector3_t(vec_[0]*rhs[0], vec_[1]*rhs[1], vec_[2]*rhs[2]);
    } // operator*()

    INLINE vector3_t operator/(real_t todiv) {
      return vector3_t(vec_[0] / todiv, vec_[1] / todiv, vec_[2] / todiv);
    }

    INLINE vector3_t operator/(vector3_t todiv) {
      return vector3_t(vec_[0] / todiv[0], vec_[1] / todiv[1], vec_[2] / todiv[2]);
    } // operator/()
  } vector3_t;


  typedef struct matrix3x3_t {
    typedef boost::array <real_t, 3> mat3_t;
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
