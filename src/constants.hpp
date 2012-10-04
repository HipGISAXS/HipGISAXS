/**
 * $Id: constants.hpp 38 2012-08-09 23:01:20Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 */

#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

// some extra buffer when estimating device memory need
const int DEVICE_RESERVE_MEM_ = 100 * 1024 * 1024;

// 4 decomposition hyperblock parameters to auto-tune later.
// Y and Z are most important for tuning.
const int BLOCK_X_ = 100;
const int BLOCK_Y_ = 20;
const int BLOCK_Z_ = 15;
const int BLOCK_T_ = 2000;

// CUDA block sizes (kernel 1)
//const int BLOCK_CUDA_ = 64;

// CUDA block sizes (kernel 2)
//const int BLOCK_FF_T_ = 4;
//const int BLOCK_FF_Y_ = 2;
//const int BLOCK_FF_Z_ = 2;

// CUDA block size (reduction kernel)
const int BLOCK_REDUCTION_ = 8;

#endif // _CONSTANTS_H_
