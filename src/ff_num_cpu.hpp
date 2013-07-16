/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
 *  Modified: Tue 16 Jul 2013 11:50:27 AM PDT
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

#ifndef __NUMERIC_FF_CPU_HPP__
#define __NUMERIC_FF_CPU_HPP__
	
//#include <iostream>
//#include <vector>
//#include <cmath>
//#include <complex>
//#ifdef _OPENMP
//#include <omp.h>
//#endif

//#include "constants.hpp"
//#include "parameters_cpu.hpp"

#include "typedefs.hpp"
	
namespace hig {
	
	/**
	 * Class for computing Form Factor in either single or double precision on a single CPU node
	 */
	class NumericFormFactorC {
		public:

			NumericFormFactorC();
			~NumericFormFactorC();

			bool init();	// TODO ...
	
		private:
			#ifndef FF_NUM_CPU_FUSED
			
			void form_factor_kernel(float_t*, float_t*, complex_t*, float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			void reduction_kernel(unsigned int, unsigned int, unsigned int,	unsigned int,
									unsigned long int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*, complex_t*);
			#else

			void form_factor_kernel_fused(float_t*, float_t*, complex_t*, float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			#ifndef __SSE3__
			void form_factor_kernel_fused_unroll4(float_t*, float_t*, complex_t*,
									float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									#ifdef FF_NUM_CPU_PADDING
										unsigned int, unsigned int, unsigned int,
									#endif
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);
			#endif // __SSE3__

			void form_factor_kernel_fused_nqx1(const float_t*, const float_t*, const complex_t*,
									#ifndef __SSE3__
										float_vec_t&,
									#else
										float_t*,
									#endif // __SSE3__
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									#ifdef FF_NUM_CPU_PADDING
										unsigned int, unsigned int, unsigned int,
									#endif
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			void form_factor_kernel_fused_nqx1_unroll4(float_t*, float_t*, complex_t*, float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									#ifdef FF_NUM_CPU_PADDING
										unsigned int, unsigned int, unsigned int,
									#endif
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			#endif // FF_NUM_CPU_FUSED


			complex_t compute_fq(float_t, complex_t, complex_t);

			#ifdef INTEL_SB_AVX

			avx_m256c_t avx_compute_fq(avx_m256_t, avx_m256c_t, avx_m256c_t);

			#elif defined __SSE3__

			sse_m128c_t sse_compute_fq(sse_m128_t, sse_m128c_t, sse_m128c_t);

			#endif // __SSE3__

			void compute_block_size(int, int, int, int, unsigned int&, unsigned int&, unsigned int&,
									unsigned int&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);

		public:

			unsigned int compute_form_factor(int,
									#ifndef __SSE3__
										float_vec_t&,
									#else
										float_t*, unsigned int,
									#endif
									complex_t*&,
									float_t*&, int, float_t*&, int, complex_t* &qz, int,
									float_t&, float_t&, float_t&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);
			
	}; // class NumericFormFactorC
		
} // namespace hig
		
#endif // __NUMERIC_FF_CPU_HPP__
