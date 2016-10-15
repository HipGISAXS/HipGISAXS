/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: init_gpu.cuh
 *  Created: Feb 22, 2013
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

#ifndef __INIT_GPU_CUH__
#define __INIT_GPU_CUH__

#include <cuda.h>
#include <cuda_runtime.h>

namespace hig {

  void init_gpu() {
    //cudaSetDevice(1);
    cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
    #ifdef FF_NUM_GPU_DYNAMICP
      cudaDeviceSetLimit(cudaLimitMallocHeapSize, 4 * sizeof(cucomplex_t));
    #endif
    cudaFree(0);
  } // init_gpu()

} // namespace hig

#endif // __INIT_GPU_CUH__
