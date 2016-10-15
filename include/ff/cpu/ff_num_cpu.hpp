/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
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

#ifndef __NUMERIC_FF_CPU_HPP__
#define __NUMERIC_FF_CPU_HPP__
	
#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <numerics/matrix.hpp>	

namespace hig {
	
	/**
	 * Class for computing Form Factor in either single or double precision on a single CPU node
	 */
	class NumericFormFactorC {
		public:

			NumericFormFactorC();
			~NumericFormFactorC();

			bool init();	// TODO ...
	
            unsigned int compute_exact_triangle(triangle_t *, int,
                    complex_t *&, 
                    int, real_t *, real_t *, int, complex_t *,
                    RotMatrix_t &, real_t &);

            unsigned int compute_approx_triangle(real_vec_t &,
                    complex_t *&,
                    int, real_t *, real_t *,
                    int, complex_t *, RotMatrix_t &, real_t &); 
		private:
			#ifndef FF_NUM_CPU_FUSED			
				void form_factor_kernel(real_t*, real_t*, complex_t*, real_vec_t&,
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
				void form_factor_kernel_fused(real_t*, real_t*, complex_t*, real_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									complex_t*);
//				#ifndef __SSE3__
					void form_factor_kernel_fused_unroll4(real_t*, real_t*, complex_t*,
									real_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									#ifdef FF_NUM_CPU_PADDING
										unsigned int, unsigned int, unsigned int,
									#endif
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									complex_t*);
//				#endif // __SSE3__
				void form_factor_kernel_fused_nqx1(const real_t*, const real_t*, const complex_t*,
//									#ifndef __SSE3__
										real_vec_t&,
//									#else
//										real_t*,
//									#endif // __SSE3__
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									#ifdef FF_NUM_CPU_PADDING
										unsigned int, unsigned int, unsigned int,
									#endif
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									complex_t*);

				void form_factor_kernel_fused_nqx1_unroll4(real_t*, real_t*, complex_t*, real_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									#ifdef FF_NUM_CPU_PADDING
										unsigned int, unsigned int, unsigned int,
									#endif
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									complex_t*);

			#endif // FF_NUM_CPU_FUSED


			complex_t compute_fq(real_t, complex_t, complex_t);

/*			#ifdef INTEL_SB_AVX
				avx_m256c_t avx_compute_fq(avx_m256_t, avx_m256c_t, avx_m256c_t);
			#elif defined __SSE3__
				sse_m128c_t sse_compute_fq(sse_m128_t, sse_m128c_t, sse_m128c_t);
			#endif // __SSE3__
*/
			void compute_block_size(int, int, int, int, unsigned int&, unsigned int&, unsigned int&,
									unsigned int&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);

		public:

			unsigned int compute_form_factor(int,
//									#ifndef __SSE3__
										real_vec_t&,
//									#else
//										real_t*, unsigned int,
//									#endif
									complex_t*&,
									real_t*&, int, real_t*&, int, complex_t* &qz, int,
									real_t*&,
									real_t&, real_t&, real_t&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);
			
	}; // class NumericFormFactorC
		
} // namespace hig
		
#endif // __NUMERIC_FF_CPU_HPP__
