/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: reduction.cuh
 *  Created: Aug 25, 2012
 *  Modified: Tue 16 Jul 2013 12:18:09 PM PDT
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
