/***
  *  Project:
  *
  *  File: init_gpu.cuh
  *  Created: Feb 22, 2013
  *  Modified: Wed 01 May 2013 06:02:04 PM PDT
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
		cudaSetDevice(2);
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
		#ifdef FF_NUM_GPU_DYNAMICP
			cudaDeviceSetLimit(cudaLimitMallocHeapSize, 4 * sizeof(cucomplex_t));
		#endif
		cudaFree(0);
	} // init_gpu()

} // namespace hig

#endif // __INIT_GPU_CUH__
