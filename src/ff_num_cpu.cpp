/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
 *  Modified: Sat 20 Apr 2013 12:47:15 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */
	
#include <iostream>
#include <cmath>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef PROFILE_PAPI
#include <papi.h>
#endif

#ifdef __SSE3__
#include "sse3_numerics.hpp"
#endif

#include "woo/timer/woo_boostchronotimers.hpp"

#include "constants.hpp"
#include "parameters.hpp"
#include "parameters_cpu.hpp"

#include "ff_num_cpu.hpp"
	
namespace hig {
	
	NumericFormFactorC::NumericFormFactorC() { }

	NumericFormFactorC::~NumericFormFactorC() { }

	// TODO: ...
	bool NumericFormFactorC::init() { return true; }

	/**
	 * The main host function called from outside, as part of the API for a single node.
	 */
	unsigned int NumericFormFactorC::compute_form_factor(int rank,
						#ifndef __SSE3__
							float_vec_t &shape_def,
						#else
							float_t* shape_def, unsigned int num_triangles,
						#endif
						complex_t* &ff,
						float_t* &qx, int nqx, float_t* &qy, int nqy, complex_t* &qz, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
						#ifdef FINDBLOCK
							, const int block_x, const int block_y, const int block_z, const int block_t
						#endif
						) {
		double kernel_time = 0.0, reduce_time = 0.0, total_kernel_time = 0.0, total_reduce_time = 0.0,
				temp_mem_time = 0.0, total_mem_time = 0.0;
		#ifdef _OPENMP
			if(rank == 0)
				std::cout << "++      Number of OpenMP threads: " << omp_get_max_threads() << std::endl;
		#endif
	
		#ifndef __SSE3__
			unsigned int num_triangles = shape_def.size() / CPU_T_PROP_SIZE_;
		#endif
		if(num_triangles < 1) return 0;

		#ifdef __SSE3__
			unsigned int shape_padding = (16 - (num_triangles & 15)) & 15;
		#endif
	
		//#ifndef FF_NUM_CPU_PADDING
			unsigned long int total_qpoints = nqx * nqy * nqz;
			unsigned long int host_mem_usage = ((unsigned long int) nqx + nqy) * sizeof(float_t) +
												nqz * sizeof(complex_t);
		//#else
			// padding to 16 bytes
			//const unsigned int PAD_LINE_ = 16;
			//unsigned int pad_x = 0;
			//if(nqx != 1) pad_x = (PAD_LINE_ - (nqx % PAD_LINE_)) % PAD_LINE_;
			//unsigned int pad_y = (PAD_LINE_ - (nqy % PAD_LINE_)) % PAD_LINE_;
			//unsigned int pad_z = (PAD_LINE_ - (nqz % PAD_LINE_)) % PAD_LINE_;
			//unsigned int pnqx = nqx + pad_x, pnqy = nqy + pad_y, pnqz = nqz + pad_z;
			//unsigned long int total_qpoints = pnqx * pnqy * pnqz;
			//unsigned long int host_mem_usage = ((unsigned long int) pnqx + pnqy) * sizeof(float_t) +
			//									pnqz * sizeof(complex_t);
		//#endif
	
		// allocate memory for the final FF 3D matrix
		ff = new (std::nothrow) complex_t[total_qpoints];	// allocate and initialize to 0
		memset(ff, 0, total_qpoints * sizeof(complex_t));
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(complex_t) << " b" << std::endl;
			return 0;
		} // if
		host_mem_usage += total_qpoints * sizeof(complex_t);
	
		//unsigned long int matrix_size = (unsigned long int) nqx * nqy * nqz * num_triangles;
		
		// do hyperblocking to use less memory
		unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
		//#ifndef FF_NUM_CPU_PADDING
			compute_block_size(nqx, nqy, nqz, num_triangles,
								b_nqx, b_nqy, b_nqz, b_num_triangles
								#ifdef FINDBLOCK
									, block_x, block_y, block_z, block_t
								#endif
								);
		//#else
		//	// lets have hyperblock dimensions a multiple of 16 (TODO: autotune ...)
		//	compute_block_size(pnqx, pnqy, pnqz, num_triangles,
		//						b_nqx, b_nqy, b_nqz, b_num_triangles
		//						#ifdef FINDBLOCK
		//							, block_x, block_y, block_z, block_t
		//						#endif
		//						);
		//#endif
	
		unsigned long int blocked_3d_matrix_size = (unsigned long int) b_nqx * b_nqy * b_nqz;
		
		//size_t estimated_host_mem_need = host_mem_usage + blocked_matrix_size * sizeof(complex_t);
		//if(rank == 0) {
		//	std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
		//				<< " MB" << std::endl;
		//} // if
		#ifndef FF_NUM_CPU_FUSED
			unsigned long int blocked_matrix_size =
										(unsigned long int) blocked_3d_matrix_size * b_num_triangles;
			host_mem_usage += blocked_matrix_size * sizeof(complex_t);
			complex_t *fq_buffer = new (std::nothrow) complex_t[blocked_matrix_size]();
			if(fq_buffer == NULL) {
				std::cerr << "Memory allocation failed for fq_buffer. blocked_matrix_size = "
							<< blocked_matrix_size << std::endl
							<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
							<< std::endl;
				delete[] ff;
				return 0;
			} // if
		#endif
		if(rank == 0) {
			std::cout << "++             Host memory usage: " << (float) host_mem_usage / 1024 / 1024
						<< " MB" << std::endl << std::flush;
		} // if

		// compute the number of sub-blocks, along each of the 4 dimensions
		// formulate loops over each dimension, to go over each sub block
		unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
		unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
		unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
		unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);

		unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
		unsigned int curr_b_num_triangles = b_num_triangles;
		unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;

		if(rank == 0) {
			std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
						<< " x " << b_nqz << " x " << b_num_triangles << std::endl;
			std::cout << "++  Number of decomposed Hblocks: " << num_blocks
						<< " [" << nb_x << " x " << nb_y << " x " << nb_z << " x " << nb_t << "]"
						<< std::endl;
		} // if

		unsigned int block_num = 0;

		#ifdef PROFILE_PAPI
			long long int papi_total_cycles = 0, papi_total_inst = 0, papi_total_flop = 0;
			double overall_ipc = 0.0;
		#endif

		if(rank == 0) std::cout << "-- Computing form factor on CPU ... " << std::flush;

		woo::BoostChronoTimer kernel_timer;
		kernel_timer.start();

		// compute for each hyperblock
		curr_b_nqx = b_nqx;
		for(unsigned int ib_x = 0; ib_x < nb_x; ++ ib_x) {
			if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
			curr_b_nqy = b_nqy;
			for(unsigned int ib_y = 0; ib_y < nb_y; ++ ib_y) {
				if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
				curr_b_nqz = b_nqz;
				for(unsigned int ib_z = 0; ib_z < nb_z; ++ ib_z) {
					if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
					curr_b_num_triangles = b_num_triangles;
					for(unsigned int ib_t = 0; ib_t < nb_t; ++ ib_t) {
						if(ib_t == nb_t - 1)
							curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

						#ifdef VERBOSITY_3
							if(rank == 0) {
								std::cout << "-- Processing hyperblock " << ++ block_num << " of "
										<< num_blocks << ": Kernel... " << std::flush;
							} // if
						#endif

						#ifdef PROFILE_PAPI
							// PAPI_L1_DCM  0x80000000  No   Level 1 data cache misses
							// PAPI_L1_ICM  0x80000001  No   Level 1 instruction cache misses
							// PAPI_L2_DCM  0x80000002  No   Level 2 data cache misses
							// PAPI_L2_ICM  0x80000003  No   Level 2 instruction cache misses
							// PAPI_L1_TCM  0x80000006  Yes  Level 1 cache misses
							// PAPI_L2_TCM  0x80000007  No   Level 2 cache misses
							// PAPI_FPU_IDL 0x80000012  No   Cycles floating point units are idle
							// PAPI_TLB_DM  0x80000014  No   Data translation lookaside buffer misses
							// PAPI_TLB_IM  0x80000015  No   Instruction translation lookaside buffer misses
							// PAPI_TLB_TL  0x80000016  Yes  Total translation lookaside buffer misses
							// PAPI_STL_ICY 0x80000025  No   Cycles with no instruction issue
							// PAPI_HW_INT  0x80000029  No   Hardware interrupts
							// PAPI_BR_TKN  0x8000002c  No   Conditional branch instructions taken
							// PAPI_BR_MSP  0x8000002e  No   Conditional branch instructions mispredicted
							// PAPI_TOT_INS 0x80000032  No   Instructions completed
							// PAPI_FP_INS  0x80000034  No   Floating point instructions
							// PAPI_BR_INS  0x80000037  No   Branch instructions
							// PAPI_VEC_INS 0x80000038  No   Vector/SIMD instructions (could include integer)
							// PAPI_RES_STL 0x80000039  No   Cycles stalled on any resource
							// PAPI_TOT_CYC 0x8000003b  No   Total cycles
							// PAPI_L1_DCH  0x8000003e  Yes  Level 1 data cache hits
							// PAPI_L2_DCH  0x8000003f  Yes  Level 2 data cache hits
							// PAPI_L1_DCA  0x80000040  No   Level 1 data cache accesses
							// PAPI_L2_DCA  0x80000041  No   Level 2 data cache accesses
							// PAPI_L1_ICH  0x80000049  Yes  Level 1 instruction cache hits
							// PAPI_L2_ICH  0x8000004a  No   Level 2 instruction cache hits
							// PAPI_L1_ICA  0x8000004c  No   Level 1 instruction cache accesses
							// PAPI_L2_ICA  0x8000004d  No   Level 2 instruction cache accesses
							// PAPI_L1_ICR  0x8000004f  No   Level 1 instruction cache reads
							// PAPI_L1_TCH  0x80000055  Yes  Level 1 total cache hits
							// PAPI_L2_TCH  0x80000056  Yes  Level 2 total cache hits
							// PAPI_L1_TCA  0x80000058  Yes  Level 1 total cache accesses
							// PAPI_L2_TCA  0x80000059  No   Level 2 total cache accesses
							// PAPI_FML_INS 0x80000061  No   Floating point multiply instructions
							// PAPI_FAD_INS 0x80000062  No   Floating point add instructions
							//                               (Also includes subtract instructions)
							// PAPI_FDV_INS 0x80000063  No   Floating point divide instructions
							//                               (Counts both divide and square root instructions)
							// PAPI_FSQ_INS 0x80000064  No   Floating point square root instructions
							//                               (Counts both divide and square root instructions)
							// PAPI_FP_OPS  0x80000066  No   Floating point operations
							// PAPI_SP_OPS  0x80000067  No   Floating point operations; optimized to count
							//                               scaled single precision vector operations
							// PAPI_DP_OPS  0x80000068  No   Floating point operations; optimized to count
							//                               scaled double precision vector operations

							int papi_events[3] = { PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_OPS };
							//int papi_events[3] = { PAPI_FML_INS, PAPI_FAD_INS, PAPI_FDV_INS };
							//int papi_events[3] = { PAPI_FP_OPS, PAPI_SP_OPS, PAPI_DP_OPS };
							long long  papi_counter_values[3];
							PAPI_start_counters(papi_events, 3);
						#endif

						// call the main kernel
						#ifndef FF_NUM_CPU_FUSED
							form_factor_kernel(qx, qy, qz, shape_def,
									curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
									b_nqx, b_nqy, b_nqz, b_num_triangles,
									ib_x, ib_y, ib_z, ib_t,
									fq_buffer);
						#else
							if(nqx == 1) {
								form_factor_kernel_fused_nqx1(qx, qy, qz, shape_def,
								//form_factor_kernel_fused_nqx1_unroll4(qx, qy, qz, shape_def,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										nqx, nqy, nqz, num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										ff);
							} else {
								#ifdef __SSE3__
									if(rank == 0)
										std::cout << "uh-oh: no SSE3 version!" << std::endl;
								#else
									form_factor_kernel_fused_unroll4(qx, qy, qz, shape_def,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										nqx, nqy, nqz, num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										ff);
								#endif // __SSE3__
							} // if-else
						#endif

						#ifdef VERBOSITY_3
							if(rank == 0) {
								std::cout << "done [" << temp_kernel_time << "s]. Reduction... "
										<< std::flush;
							} // if
						#endif

						#ifndef FF_NUM_CPU_FUSED
							// call the reduction kernel
							reduction_kernel(curr_b_nqx, curr_b_nqy, curr_b_nqz,
									curr_b_num_triangles, blocked_matrix_size,
									b_nqx, b_nqy, b_nqz, num_triangles,
									nqx, nqy, nqz,
									ib_x, ib_y, ib_z, ib_t,
									fq_buffer, ff);
						#endif

						#ifdef VERBOSITY_3
							if(rank == 0) {
								std::cout << "done [" << temp_reduce_time << "s]." << std::endl;
							} // if
						#endif

						#ifdef PROFILE_PAPI
							PAPI_stop_counters(papi_counter_values, 3);
							papi_total_cycles += papi_counter_values[0];
							papi_total_inst += papi_counter_values[1];
							papi_total_flop += papi_counter_values[2];
						#endif
					} // for ib_t
				} // for ib_z
			} // for ib_y
		} // for ib_x

		kernel_timer.stop();
		float_t ktime = kernel_timer.elapsed_msec();

		#ifndef FF_NUM_CPU_FUSED
			delete[] fq_buffer;
		#endif

		if(rank == 0) std::cout << "done." << std::endl;

		#ifdef PROFILE_PAPI
			if(rank == 0) {
				std::cout << "++                  PAPI_TOT_CYC: " << papi_total_cycles << std::endl;
				std::cout << "++                  PAPI_TOT_INS: " << papi_total_inst << std::endl;
				std::cout << "++                   PAPI_FP_OPS: " << papi_total_flop << std::endl;
				std::cout << "++                           IPC: "
							<< (double) papi_total_inst / papi_total_cycles << std::endl;
				std::cout << "++                        Flop/s: "
							<< (double) papi_total_flop / (1000000 * ktime) << " GFlop/s" << std::endl;
			} // if
		#endif

		if(rank == 0) {
			std::cout << "**                FF kernel time: " << ktime << " ms."
						<< std::endl;
		} // if
	
		pass_kernel_time = total_kernel_time;
		red_time = total_reduce_time;
		mem_time = total_mem_time;
	
		return num_triangles;
	} // NumericFormFactorC::compute_form_factor()
		

	#ifndef FF_NUM_CPU_FUSED

	/**
	 * the main Form Factor kernel function - for one hyperblock.
	 * needs reduction kernel after this.
	 */
	void NumericFormFactorC::form_factor_kernel(float_t* qx, float_t* qy, complex_t* qz,
					float_vec_t& shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					complex_t* fq_buffer) {
		if(fq_buffer == NULL || qx == NULL || qy == NULL || qz == NULL) return;
	
		#pragma omp parallel
		{
			#pragma omp for nowait //schedule(auto)
			for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
				unsigned int shape_off = (ib_t * b_num_triangles + i_t) * CPU_T_PROP_SIZE_;
				float_t s = shape_def[shape_off];
				float_t nx = shape_def[shape_off + 1];
				float_t ny = shape_def[shape_off + 2];
				float_t nz = shape_def[shape_off + 3];
				float_t x = shape_def[shape_off + 4];
				float_t y = shape_def[shape_off + 5];
				float_t z = shape_def[shape_off + 6];
	
				unsigned long int xy_size = (unsigned long int) curr_nqx * curr_nqy;
				unsigned long int matrix_off = xy_size * curr_nqz * i_t;
				unsigned int start_z = b_nqz * ib_z;
				unsigned int start_y = b_nqy * ib_y;
				unsigned int start_x = b_nqx * ib_x;
	
				for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
						++ i_z, ++ global_i_z) {
					unsigned long int off_start = matrix_off + xy_size * i_z;
					complex_t temp_z = qz[global_i_z];
					complex_t qz2 = temp_z * temp_z;
					complex_t qzn = temp_z * nz;
					complex_t qzt = temp_z * z;
	
					for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
							++ i_y, ++ global_i_y) {
						unsigned long int xy_off_start = (unsigned long int) curr_nqx * i_y;
						float_t temp_y = qy[global_i_y];
						float_t qy2 = temp_y * temp_y;
						float_t qyn = temp_y * ny;
						float_t qyt = temp_y * y;
	
						for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
								++ i_x, ++ global_i_x) {
							unsigned long int off = off_start + xy_off_start + i_x;
							float_t temp_x = qx[global_i_x];
							complex_t q2 = temp_x * temp_x + qy2 + qz2;
							complex_t qt = temp_x * x + qyt + qzt;
							complex_t qn = (temp_x * nx + qyn + qzn) / q2;
	
							fq_buffer[off] = compute_fq(s, qt, qn);
						} // for z
					} // for y
				} // for x
			} // for t
		} // pragma omp parallel
	} // NumericFormFactorC::form_factor_kernel()
	
	
	/**
	 * Reduction kernel
	 */
	void NumericFormFactorC::reduction_kernel(
			unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
			unsigned int curr_b_num_triangles, unsigned long int blocked_matrix_size,
			unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
			unsigned int num_triangles, unsigned int nqx, unsigned int nqy, unsigned int nqz,
			unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
			complex_t* fq_buffer, complex_t* ff) {
		if(fq_buffer == NULL || ff == NULL) return;
		
		unsigned long int curr_b_xyz = (unsigned long int) curr_b_nqx * curr_b_nqy * curr_b_nqz;
		unsigned long int curr_b_nqxy = curr_b_nqx * curr_b_nqy;
		complex_t temp_complex(0.0, -1.0);
	
		unsigned long int super_nqxy = nqx * nqy;
		unsigned long int temp_1 = super_nqxy * ib_z * b_nqz + nqx * ib_y * b_nqy;
	
		// reduction over all triangles
		for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
			unsigned long int temp_3 = ib_x * b_nqx + i_x;
			#pragma omp parallel
			{
				#pragma omp for nowait //schedule(auto)
				for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
					unsigned long int temp_nqxiyix = curr_b_nqx * i_y + i_x;
					unsigned long int temp_2 = nqx * i_y;
	
					for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
						unsigned long int temp_i = curr_b_nqxy * i_z + temp_nqxiyix;
						complex_t total(0.0, 0.0);
	
						for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
							total += fq_buffer[curr_b_xyz * i_t + temp_i];
						} // for i_t
	
						unsigned long int super_i = temp_1 + temp_2 + temp_3 + super_nqxy * i_z;
						//unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
						//							nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
						ff[super_i] += total * temp_complex;
					} // for i_z
				} // for i_y
			} // omp parallel
		} // for i_x
	} // NumericFormFactorC::reduction_kernel()
	
	#else
	
	/**
	 * the main Form Factor kernel with fused reduction function - for one hyperblock.
	 * using this as baseline.
	 */
	void NumericFormFactorC::form_factor_kernel_fused(float_t* qx, float_t* qy, complex_t* qz,
					float_vec_t& shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					complex_t* ff) {
		if(ff == NULL || qx == NULL || qy == NULL || qz == NULL) return;
	
		unsigned long int xy_size = (unsigned long int) curr_nqx * curr_nqy;
		unsigned int start_z = b_nqz * ib_z;
		unsigned int start_y = b_nqy * ib_y;
		unsigned int start_x = b_nqx * ib_x;

		#pragma omp parallel
		{
			#pragma omp for collapse(3)
			for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
				for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
					for(int i_x = 0; i_x < curr_nqx; ++ i_x) {
						for(int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
							unsigned int shape_off = (ib_t * b_num_triangles + i_t) * CPU_T_PROP_SIZE_;
							float_t s = shape_def[shape_off];
							float_t nx = shape_def[shape_off + 1];
							float_t ny = shape_def[shape_off + 2];
							float_t nz = shape_def[shape_off + 3];
							float_t x = shape_def[shape_off + 4];
							float_t y = shape_def[shape_off + 5];
							float_t z = shape_def[shape_off + 6];
	
							complex_t temp_z = qz[start_z + i_z];
							complex_t qz2 = temp_z * temp_z;
							complex_t qzn = temp_z * nz;
							complex_t qzt = temp_z * z;
	
							float_t temp_y = qy[start_y + i_y];
							float_t qy2 = temp_y * temp_y;
							float_t qyn = temp_y * ny;
							float_t qyt = temp_y * y;
	
							float_t temp_x = qx[start_x + i_x];
							complex_t q2 = temp_x * temp_x + qy2 + qz2;
							complex_t qt = temp_x * x + qyt + qzt;
							complex_t qn = (temp_x * nx + qyn + qzn) / q2;

							unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
														nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
							ff[super_i] += compute_fq(s, qt, qn);
						} // for t
					} // for x
				} // for y
			} // for z
		} // pragma omp parallel
	} // NumericFormFactorC::form_factor_kernel_fused()
	

	/**
	 * Include kernels with various optimizations.
	 */
	#include "ff_num_cpu_kernels.hpp"
	
	#endif // FF_NUM_CPU_FUSED
	

	/**
	 * Computational kernel function.
	 */
	inline complex_t NumericFormFactorC::compute_fq(float_t s, complex_t qt, complex_t qn) {
		complex_t v1 = qn * complex_t(cos(qt.real()), sin(qt.real()));
		float_t v2 = s * exp(qt.imag());
		return v1 * v2;
	} // NumericFormFactorC::compute_fq()
	
	
	/**
	 * Function to compute the decomposition block size
	 * TODO: Improve it later ...
	 */
	void NumericFormFactorC::compute_block_size(int nqx, int nqy, int nqz, int num_triangles,
			unsigned int& b_nqx, unsigned int& b_nqy, unsigned int& b_nqz, unsigned int& b_num_triangles
			#ifdef FINDBLOCK
				, const int block_x, const int block_y, const int block_z, const int block_t
			#endif
			) {
		b_nqx = (unsigned int) nqx; b_nqy = (unsigned int) nqy; b_nqz = (unsigned int) nqz;
		b_num_triangles = (unsigned int) num_triangles;
	
		#ifdef FINDBLOCK
			b_nqx = (b_nqx > block_x) ? block_x : b_nqx;
			b_nqy = (b_nqy > block_y) ? block_y : b_nqy;
			b_nqz = (b_nqz > block_z) ? block_z : b_nqz;
			b_num_triangles = (b_num_triangles > block_t) ? block_t : b_num_triangles;
		#else
			b_nqx = (b_nqx > CPU_BLOCK_X_) ? CPU_BLOCK_X_ : b_nqx;
			b_nqy = (b_nqy > CPU_BLOCK_Y_) ? CPU_BLOCK_Y_ : b_nqy;
			b_nqz = (b_nqz > CPU_BLOCK_Z_) ? CPU_BLOCK_Z_ : b_nqz;
			b_num_triangles = (b_num_triangles > CPU_BLOCK_T_) ? CPU_BLOCK_T_ : b_num_triangles;
		#endif
	} // NumericFormFactorC::compute_block_size()
	
	
} // namespace hig
