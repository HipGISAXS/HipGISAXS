/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
 *  Modified: Sat 02 Mar 2013 09:30:13 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
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
			void form_factor_kernel(float_t*, float_t*, complex_t*, float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			void form_factor_kernel_fused(float_t*, float_t*, complex_t*, float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			void form_factor_kernel_fused_unroll4(float_t*, float_t*, complex_t*, float_vec_t&,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*);

			complex_t compute_fq(float_t, complex_t, complex_t);

			void reduction_kernel(unsigned int, unsigned int, unsigned int,	unsigned int,
									unsigned long int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									complex_t*, complex_t*);

			void compute_block_size(int, int, int, int, unsigned int&, unsigned int&, unsigned int&,
									unsigned int&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);

		public:

			unsigned int compute_form_factor(int, float_vec_t&, complex_t*&,
									float_t*&, int, float_t*&, int, complex_t* &qz, int,
									float_t&, float_t&, float_t&
									#ifdef FINDBLOCK
										, const int, const int, const int, const int
									#endif
									);
			
	}; // class NumericFormFactorC
		
} // namespace hig
		
#endif // __NUMERIC_FF_CPU_HPP__
