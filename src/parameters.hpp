/**
 * $Id: constants.hpp 38 2012-08-09 23:01:20Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 *
 * File: parameters.hpp
 * Created: June 5, 2012
 * Modified: Mon 26 Nov 2012 11:53:42 AM PST
 *
 * Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

// put these into namespace too ...
namespace hig {

	// some extra buffer when estimating device memory need
	const int DEVICE_RESERVE_MEM_ = 100 * 1024 * 1024;

	// padding triangle inputs to 8 entries (8 * 4 = 32 bytes)
	// instead of 7 entries (7 * 4 = 28 bytes) in the hope to
	// have better efficiency with cache line alignments
	const int T_PROP_SIZE_ = 8;
	const float T_PROP_SIZE_F_ = T_PROP_SIZE_ * 1.0f;

	// 4 decomposition hyperblock parameters to auto-tune later.
	// Y and Z are most important for tuning.
	const int BLOCK_X_ = 100;
	const int BLOCK_Y_ = 32;
	const int BLOCK_Z_ = 16;
	const int BLOCK_T_ = 2048;

	// CUDA block sizes (kernel 1)
	//const int BLOCK_CUDA_ = 64;

	// CUDA block sizes (kernel 2)
	//const int BLOCK_FF_T_ = 4;
	//const int BLOCK_FF_Y_ = 2;
	//const int BLOCK_FF_Z_ = 2;

	// CUDA block size (reduction kernel)
	const int BLOCK_REDUCTION_ = 8;

	// at most these many entries are possible for each of t, y, z in a bock
	const int MAX_SHARED_SIZE_ = 64;   // not used in shared2
	const int MAX_NUM_THREADS_ = 32;

	// the number of entries to copy from shared to device at a time
	const int FQ_COPY_SIZE_ = 32;
	const float FQ_COPY_SIZE_F_ = FQ_COPY_SIZE_ * 1.0f;
	// offset to avoid bank conflicts when FQ_COPY_SIZE_ is multiple of # banks (32)
	const int BANK_OFF_ = 1;

	// constants used only by kernels:
	/*
	// at most these many entries are possible for each of t, y, z in a bock
	__constant__ const int MAX_SHARED_SIZE_ = 64;   // not used in shared2
	__constant__ const int MAX_NUM_THREADS_ = 32;

	// the number of entries to copy from shared to device at a time
	__constant__ const int FQ_COPY_SIZE_ = 32;
	__constant__ const float FQ_COPY_SIZE_F_ = FQ_COPY_SIZE_ * 1.0f;
	// offset to avoid bank conflicts when FQ_COPY_SIZE_ is multiple of # banks (32)
	__constant__ const int BANK_OFF_ = 1;
	*/

} // namespace hig

#endif // _PARAMETERS_H_
