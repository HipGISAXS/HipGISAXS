/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: cu_utilities.cuh
 *  Created: Feb 19, 2013
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

#ifndef __CU_UTILITIES_CUH__
#define __CU_UTILITIES_CUH__

#include <common/typedefs.hpp>
#include <numerics/gpu/cu_complex_numeric.cuh>

namespace hig {

  __device__ static __inline__ cucomplex_t fq_inv(cucomplex_t value, float_t y) {
    cucomplex_t temp1 = value * y / (float_t) 2.0;
    cucomplex_t temp;
    if(cuCabsf(value) < 1e-14) temp = make_cuC(y, (float_t) 0.0);
    else temp = 2.0 * cuCexpi(temp1) * cuCsin(temp1) / value;
    //if(cuCabsf(temp) < 1e-38) temp = make_cuC(y, (float_t) 0.0);
    return temp;
  } // fq_inv()


  __device__ static __inline__ void compute_meshpoints(const float_t qx, const float_t qy,
                    const cucomplex_t qz, const float_t* rot,
                    cucomplex_t& mqx, cucomplex_t& mqy, cucomplex_t& mqz) {
    //mqx = make_cuC(qx * rot[0] + qy * rot[1] + qz.x * rot[2], qz.y * rot[2]);
    //mqy = make_cuC(qx * rot[3] + qy * rot[4] + qz.x * rot[5], qz.y * rot[5]);
    //mqz = make_cuC(qx * rot[6] + qy * rot[7] + qz.x * rot[8], qz.y * rot[8]);
    mqx = make_cuC(qx * rot[0] + qy * rot[1] + cu_real(qz) * rot[2], cu_imag(qz) * rot[2]);
    mqy = make_cuC(qx * rot[3] + qy * rot[4] + cu_real(qz) * rot[5], cu_imag(qz) * rot[5]);
    mqz = make_cuC(qx * rot[6] + qy * rot[7] + cu_real(qz) * rot[8], cu_imag(qz) * rot[8]);
  } // compute_meshpoints()

} // namespace hig

#endif // __CU_UTILITIES_CUH__
