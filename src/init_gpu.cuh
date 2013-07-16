/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: init_gpu.cuh
 *  Created: Feb 22, 2013
 *  Modified: Tue 16 Jul 2013 11:51:03 AM PDT
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

#ifndef __INIT_GPU_CUH__
#define __INIT_GPU_CUH__

#include <cuda.h>
#include <cuda_runtime.h>

namespace hig {

	void init_gpu() {
		std::cout << "-- Waking up GPU(s) ..." << std::flush << std::endl;
		cudaSetDevice(3);
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
		#ifdef FF_NUM_GPU_DYNAMICP
			cudaDeviceSetLimit(cudaLimitMallocHeapSize, 4 * sizeof(cucomplex_t));
		#endif
		cudaFree(0);
	} // init_gpu()

} // namespace hig

#endif // __INIT_GPU_CUH__
