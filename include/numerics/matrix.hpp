/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: globals.hpp
 *  Created: Jun 05, 2012
 *  Modified: Thu 13 Mar 2014 04:40:43 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef ROT_MATRIX__HPP
#define ROT_MATRIX__HPP

#include <vector>
#include <cmath>
#include <stdexcept>

#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <common/cudafy.hpp>

#ifdef FF_CPU_OPT
#include <mkl_cblas.h>
#endif // FF_CPU_OPT


namespace hig {
  class RotMatrix_t {
    private:
      INLINE void swapit(real_t & a, real_t & b){
        real_t t = a; a = b; b = t;
      }
      real_t data_[9];
      INLINE void set(real_t v){ for (int i=0; i<9; i++) data_[i] = v; }

    public:
      // constructor 1
      CUDAFY RotMatrix_t() {
        data_[0] = 1.; data_[1] = 0.; data_[2] = 0.;
        data_[3] = 0.; data_[4] = 1.; data_[5] = 0.;
        data_[6] = 0.; data_[7] = 0.; data_[8] = 1.;
      }

      //constructor 2
      CUDAFY RotMatrix_t(int axis, real_t a){
        real_t c = std::cos(a);
        real_t s = std::sin(a);
        switch (axis){
          case 0:
            data_[0] = 1.; data_[1] = 0.; data_[2] = 0.;
            data_[3] = 0.; data_[4] = c;  data_[5] = -s;
            data_[6] = 0.; data_[7] = s;  data_[8] =  c;
            break;
          case 1:
            data_[0] = c;  data_[1] = 0.; data_[2] =-s;
            data_[3] = 0.; data_[4] = 1.; data_[5] = 0.;
            data_[6] = s;  data_[7] = 0.; data_[8] = c;
            break;
          case 2:
            data_[0] = c;  data_[1] =-s;  data_[2] = 0.;
            data_[3] = s;  data_[4] = c;  data_[5] = 0.;
            data_[6] = 0.; data_[7] = 0.; data_[8] = 1.;
            break;
        }
      }

      // constructor 3
      CUDAFY RotMatrix_t(vector3_t & a, vector3_t & b, vector3_t & c){
        data_[0] = a[0]; data_[1] = a[1]; data_[2] = a[2];
        data_[3] = b[0]; data_[4] = b[1]; data_[5] = b[2];
        data_[6] = c[0]; data_[7] = c[1]; data_[8] = c[2];
      }

      // copy constructor
      CUDAFY RotMatrix_t(const RotMatrix_t & obj){
        for(int i=0; i<9; i++) data_[i] = obj.data_[i];
      } 

      // transpose
      CUDAFY void transpose(){
        swapit(data_[1],data_[3]);
        swapit(data_[2],data_[6]);
        swapit(data_[5],data_[7]);
      }

      // assignment operator
      CUDAFY RotMatrix_t operator= (const RotMatrix_t & rhs){
        for(int i=0; i<9; i++) data_[i] = rhs.data_[i];
        return *this;
      }

      // DEBUG -- print 
      void print() const {
        for (int i = 0; i < 9; i++){
          std::cout << "  " << data_[i];
          if (i % 3 == 2) std::cout << std::endl;
        }
      }

      // Matrix-Matrix multiplication
      INLINE RotMatrix_t operator*(const RotMatrix_t & rhs){
        RotMatrix_t res;
        res.set(0.);
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                res.data_[i * 3 + j] += data_[i * 3 + k] * rhs.data_[j + 3 * k];
        return res;
      }

      // rotate real vector
      CUDAFY vector3_t operator*(const vector3_t & rhs){
        vector3_t lhs;
        lhs[0] = data_[0] * rhs.vec_[0] + data_[1] * rhs.vec_[1] + data_[2] * rhs.vec_[2];
        lhs[1] = data_[3] * rhs.vec_[0] + data_[4] * rhs.vec_[1] + data_[5] * rhs.vec_[2];
        lhs[2] = data_[6] * rhs.vec_[0] + data_[7] * rhs.vec_[1] + data_[8] * rhs.vec_[2];
        return lhs;
      }

      // Matrix-vector multiplication (Rotation)
      inline std::vector<complex_t> rotate(real_t x, real_t y, complex_t z) {
        std::vector<complex_t> res;
        res.resize(3);
        res[0] = data_[0] * x + data_[1] * y + data_[2] * z;
        res[1] = data_[3] * x + data_[4] * y + data_[5] * z;
        res[2] = data_[6] * x + data_[7] * y + data_[8] * z;
        return res;
      } // rotate()

      inline void rotate_new(real_t x, real_t y, complex_t z, complex_t* res) {
        res[0] = data_[0] * x + data_[1] * y + data_[2] * z;
        res[1] = data_[3] * x + data_[4] * y + data_[5] * z;
        res[2] = data_[6] * x + data_[7] * y + data_[8] * z;
      } // rotate()

#ifdef FF_CPU_OPT
  #ifdef FF_CPU_OPT_MKL
      inline void rotate_vec(const int VEC_LEN, real_t* x, real_t* y, complex_t* z, complex_t* res) {
        complex_t *res0 = res, *res1 = res + VEC_LEN, *res2 = res + 2 * VEC_LEN;
        real_t temp[VEC_LEN];

        cblas_dcopy(VEC_LEN, y, 1, temp, 1);
        cblas_daxpby(VEC_LEN, data_[0], x, 1, data_[1], temp, 1);
        for(int i = 0; i < VEC_LEN; ++ i) res0[i] = complex_t(temp[i], (real_t) 0.0);
        cblas_zaxpy(VEC_LEN, &data_[2], z, 1, res0, 1);

        cblas_dcopy(VEC_LEN, y, 1, temp, 1);
        cblas_daxpby(VEC_LEN, data_[3], x, 1, data_[4], temp, 1);
        for(int i = 0; i < VEC_LEN; ++ i) res1[i] = complex_t(temp[i], (real_t) 0.0);
        cblas_zaxpy(VEC_LEN, &data_[5], z, 1, res1, 1);

        cblas_dcopy(VEC_LEN, y, 1, temp, 1);
        cblas_daxpby(VEC_LEN, data_[6], x, 1, data_[7], temp, 1);
        for(int i = 0; i < VEC_LEN; ++ i) res2[i] = complex_t(temp[i], (real_t) 0.0);
        cblas_zaxpy(VEC_LEN, &data_[8], z, 1, res2, 1);
      } // rotate()
  #endif // FF_CPU_OPT_MKL
#endif // FF_CPU_OPT

#ifdef USE_GPU
      __device__ void rotate(real_t x, real_t y, cucomplex_t z, 
              cucomplex_t &mx, cucomplex_t &my, cucomplex_t &mz){
          mx.x = data_[0] * x + data_[1] * y + data_[2] * z.x; mx.y = data_[2] * z.y;
          my.x = data_[3] * x + data_[4] * y + data_[5] * z.x; my.y = data_[5] * z.y;
          mz.x = data_[6] * x + data_[7] * y + data_[8] * z.x; mz.y = data_[8] * z.y;
      }
#endif
  };
} // namespace hig
#endif // _MATRIX_HPP_
