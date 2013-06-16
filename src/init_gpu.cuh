/***
  *  Project:
  *
  *  File: init_gpu.cuh
  *  Created: Feb 22, 2013
  *  Modified: Mon 06 May 2013 10:36:10 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __INIT_GPU_CUH__
#define __INIT_GPU_CUH__

#include <cuda.h>
#include <cuda_runtime.h>

namespace hig {

	void init_gpu() {
		std::cout << "-- Waking up GPU(s) ..." << std::flush << std::endl;
		//cudaSetDevice(3);
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
		#ifdef FF_NUM_GPU_DYNAMICP
			cudaDeviceSetLimit(cudaLimitMallocHeapSize, 4 * sizeof(cucomplex_t));
		#endif
		cudaFree(0);
	} // init_gpu()

} // namespace hig

#endif // __INIT_GPU_CUH__
