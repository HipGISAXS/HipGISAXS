/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: cudafy.hpp
 *  Created: Aug 25, 2012
 *  Modified: Wed 08 Oct 2014 12:13:01 PM PDT
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

#ifndef __CUDAFY_HPP__
#define __CUDAFY_HPP__

#ifdef USE_GPU
#define CUDAFY __host__ __device__
#define INLINE __host__ __device__ __inline__
#else
#define CUDAFY
#define INLINE inline
#endif



#endif // __CUDAFY_HPP__

