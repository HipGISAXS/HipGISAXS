/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_mic.hpp
 *  Created: Apr 02, 2013
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

#ifndef __NUMERIC_FF_MIC_HPP__
#define __NUMERIC_FF_MIC_HPP__
	
#include <common/typedefs.hpp>
	
namespace hig {
	
	/**
	 * Class for computing Form Factor in either single or double precision on a CPU + MIC node
	 * TODO: this should be merged with the pure CPU version
	 */
	class NumericFormFactorM {
		public:

			NumericFormFactorM();
			~NumericFormFactorM();

			bool init();	// TODO ...
	
		private:
			void form_factor_kernel(real_t*, real_t*, scomplex_t*, real_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									scomplex_t*);
			__attribute__((target(mic:0)))
			void form_factor_kernel_db(real_t*, real_t*, scomplex_t*, real_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									scomplex_t*);
			#ifndef FF_NUM_MIC_SWAP
			__attribute__((target(mic:0)))
			void form_factor_kernel_opt(real_t*, real_t*, scomplex_t*, real_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									scomplex_t*, scomplex_t*);
			#else
			__attribute__((target(mic:0)))
			void form_factor_kernel_loopswap(real_t*, real_t*, scomplex_t*, real_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									scomplex_t*);
			__attribute__((target(mic:0)))
			void form_factor_kernel_loopswap_nqx1(real_t*, real_t*, scomplex_t*, real_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									scomplex_t*);
			__attribute__((target(mic:0)))
			void form_factor_kernel_loopswap_vec_nqx1(real_t*, real_t*, scomplex_t*, real_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									real_t*,
									scomplex_t*);

			__attribute__((target(mic:0))) float2_t compute_fq_nqx1(float, float2_t, float2_t);
			__attribute__((target(mic:0))) double2_t compute_fq_nqx1(double, double2_t, double2_t);
			__attribute__((target(mic:0))) mic_m512c_t compute_fq_vec_nqx1(mic_m512_t, mic_m512c_t, mic_m512c_t);
			//__attribute__((target(mic:0))) double2_t compute_fq_vec_nqx1(double, double2_t, double2_t);
			#endif

			__attribute__((target(mic:0))) float2_t compute_fq(float, float2_t, float2_t);
			__attribute__((target(mic:0))) double2_t compute_fq(double, double2_t, double2_t);

			void reduction_kernel(unsigned int, unsigned int, unsigned int,	unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									scomplex_t*, scomplex_t*);
			__attribute__((target(mic:0)))
			void reduction_kernel_db(unsigned int, unsigned int, unsigned int,	unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									scomplex_t*, scomplex_t*);

			void compute_hyperblock_size(int, int, int, int, unsigned int&, unsigned int&, unsigned int&,
									unsigned int&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);

			void move_to_main_ff(scomplex_t*,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

		public:

			unsigned int compute_form_factor(int, real_vec_t&, complex_t*&,
									real_t*, int, real_t*, int, complex_t* qz, int,
									real_t*,
									real_t&, real_t&, real_t&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);
			unsigned int compute_form_factor_db(int, real_vec_t&, complex_t*&,
									real_t*, int, real_t*, int, complex_t* qz, int,
									real_t*,
									real_t&, real_t&, real_t&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);
			unsigned int compute_form_factor_kb(int, real_t*, unsigned int, complex_t*&,
									real_t*, int, real_t*, int, complex_t* qz, int, int,
									real_t*,
									real_t&, real_t&, real_t&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);
			
	}; // class NumericFormFactorM
		
} // namespace hig
		
#endif // __NUMERIC_FF_MIC_HPP__
