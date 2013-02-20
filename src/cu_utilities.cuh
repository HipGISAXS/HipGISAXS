/***
  *  Project:
  *
  *  File: cu_utilities.cuh
  *  Created: Feb 19, 2013
  *  Modified: Tue 19 Feb 2013 08:05:09 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __CU_UTILITIES_CUH__
#define __CU_UTILITIES_CUH__

#include "typedefs.hpp"
#include "cu_complex_numeric.cuh"

namespace hig {

	__device__ static __inline__ cucomplex_t fq_inv(cucomplex_t value, float_t y) {
		cucomplex_t temp1 = value * y / (float_t) 2.0;
		cucomplex_t temp = 2.0 * cuCexpi(temp1) * cuCsin(temp1) / value;
		if(cuCabsf(temp) < 1e-14) temp = make_cuC(y, (float_t) 0.0);
		return temp;
	} // fq_inv()

} // namespace hig

#endif // __CU_UTILITIES_CUH__
