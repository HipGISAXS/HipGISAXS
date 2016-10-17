/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: cu_utilities.cuh
 *  Created: Feb 19, 2013
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

#ifndef __CU_UTILITIES_CUH__
#define __CU_UTILITIES_CUH__

#include <common/typedefs.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>

namespace hig {

  __device__ static __inline__ cucomplex_t fq_inv(cucomplex_t value, real_t y) {
    cucomplex_t temp1 = value * y / (real_t) 2.0;
    cucomplex_t temp;
    if(cuC_abs(value) < 1e-14) temp = make_cuC(y, (real_t) 0.0);
    else temp = 2.0 * cuCexpi(temp1) * cuCsin(temp1) / value;
    return temp;
  } // fq_inv()


  __device__ static __inline__ void compute_meshpoints(const real_t qx, const real_t qy,
                    const cucomplex_t qz, const real_t* rot,
                    cucomplex_t& mqx, cucomplex_t& mqy, cucomplex_t& mqz) {
    mqx = make_cuC(qx * rot[0] + qy * rot[3] + qz.x * rot[6], qz.y * rot[6]);
    mqy = make_cuC(qx * rot[1] + qy * rot[4] + qz.x * rot[7], qz.y * rot[7]);
    mqz = make_cuC(qx * rot[2] + qy * rot[5] + qz.x * rot[8], qz.y * rot[8]);
    //mqx = make_cuC(qx * rot[0] + qy * rot[1] + cu_real(qz) * rot[2], cu_imag(qz) * rot[2]);
    //mqy = make_cuC(qx * rot[3] + qy * rot[4] + cu_real(qz) * rot[5], cu_imag(qz) * rot[5]);
    //mqz = make_cuC(qx * rot[6] + qy * rot[7] + cu_real(qz) * rot[8], cu_imag(qz) * rot[8]);
  } // compute_meshpoints()

  // vector norm
  __device__ __inline__ real_t norm2 (real_t * v) {
      return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  }

  // dot product
  __device__ __inline__ real_t dot_prod (real_t * v1, real_t * v2) {
      return (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]);
  }

  // cross product
  __device__ __inline__ void cross_prod (real_t * v1, real_t * v2, real_t * res) {
      res[0] = v1[1] * v2[2] - v1[2] * v2[1];
      res[1] = v1[2] * v2[0] - v1[0] * v2[2];
      res[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }


} // namespace hig

#endif // __CU_UTILITIES_CUH__
