/**
 * $Id: reduction.cuh 38 2012-08-09 23:01:20Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 */

#ifndef _REDUCTION_CUH_
#define _REDUCTION_CUH_

namespace hig {

	#ifdef GPUR		// GPU reductions

	__global__ void reduction_kernel(cuFloatComplex*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuFloatComplex*);

	__global__ void reduction_kernel_parallel(cuFloatComplex*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuFloatComplex*);

	__global__ void reduction_kernel(cuDoubleComplex*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuDoubleComplex*);

	__global__ void reduction_kernel_parallel(cuDoubleComplex*,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuDoubleComplex*);

	#else			// CPU reductions (OpenMP)

	void reduction_kernel(unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned long int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuFloatComplex*, cuFloatComplex*);

	void reduction_kernel(unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned long int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuDoubleComplex*, cuDoubleComplex*);

	#endif // GPUR

} // namespace hig

#endif // _REDUCTION_CUH_
