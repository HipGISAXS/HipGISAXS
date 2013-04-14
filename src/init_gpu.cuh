/***
  *  Project:
  *
  *  File: init_gpu.cuh
  *  Created: Feb 22, 2013
  *  Modified: Sat 13 Apr 2013 08:46:19 PM PDT
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
		//cudaSetDevice(2);
		cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
		cudaFree(0);
	} // init_gpu()

} // namespace hig

#endif // __INIT_GPU_CUH__
