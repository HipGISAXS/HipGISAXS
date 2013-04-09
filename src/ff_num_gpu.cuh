/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_gpu.cuh
 *  Created: Nov 05, 2011
 *  Modified: Tue 09 Apr 2013 11:51:53 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __FF_NUM_GPU_CUH__
#define __FF_NUM_GPU_CUH__

//#include <iostream>
//#include <fstream>
#include <vector>
//#include <complex>
//#include <cuComplex.h>

#include "typedefs.hpp"
//#include "constants.hpp"
//#include "parameters.hpp"

namespace hig {

	/**
	 * Class for computing Form Factor in either single or double precision on a single GPU.
	 */
	//template<typename float_t, typename complex_t>
	class NumericFormFactorG {
		public:
			NumericFormFactorG(int block_cuda):
				block_cuda_(block_cuda), block_cuda_t_(0), block_cuda_y_(0), block_cuda_z_(0) { }

			NumericFormFactorG(int block_cuda_t, int block_cuda_y, int block_cuda_z): block_cuda_(0),
				block_cuda_t_(block_cuda_t), block_cuda_y_(block_cuda_y), block_cuda_z_(block_cuda_z) { }

			#ifdef FF_NUM_GPU_FUSED
				NumericFormFactorG(int block_cuda_y, int block_cuda_z): block_cuda_(0),
					block_cuda_t_(0), block_cuda_y_(block_cuda_y), block_cuda_z_(block_cuda_z) { }
			#endif

			NumericFormFactorG():	// called when not using this GPU version
				block_cuda_(0), block_cuda_t_(0), block_cuda_y_(0), block_cuda_z_(0) { }

			~NumericFormFactorG() { }

			/* original */
			unsigned int compute_form_factor(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					cucomplex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					cucomplex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
					#ifdef FINDBLOCK
						, const int, const int, const int, const int
					#endif
					);

			/* with double buffering - default */
			unsigned int compute_form_factor_db(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					cucomplex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					cucomplex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
					#ifdef FINDBLOCK
						, const int, const int, const int, const int
					#endif
					);

			/* with fused kernels and double buffering */
			unsigned int compute_form_factor_db_fused(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					cucomplex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					cucomplex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
					#ifdef FINDBLOCK
						, const int, const int, const int, const int
					#endif
					);

			/* with double buffering and optimized memory (incomplete ... TODO)*/
			unsigned int compute_form_factor_db_mem(int,
					std::vector<float_t> &shape_def, std::vector<short int> &axes,
					cucomplex_t* &ff,
					float_t* &qx_h, int nqx,
					float_t* &qy_h, int nqy,
					cucomplex_t* &qz_h, int nqz,
					float_t&, float_t&, float_t&
					#ifdef FINDBLOCK
						, const int, const int, const int, const int
					#endif
					);

		private:
			// TODO make these conditional ...

			// for kernel 1
			int block_cuda_;

			// for kernel 2 and up
			int block_cuda_t_;
			int block_cuda_y_;
			int block_cuda_z_;

			void compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
					unsigned long int estimated_device_mem_need, unsigned long int device_mem_avail,
					unsigned int& b_nqx, unsigned int& b_nqy, unsigned int& b_nqz,
					unsigned int& b_num_triangles
					#ifdef FINDBLOCK
						, const int, const int, const int, const int
					#endif
					);

			void move_to_main_ff(cucomplex_t* fq_buffer,
					unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
					unsigned int nqx, unsigned int nqy, unsigned int nqz,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z,
					cucomplex_t* ff);
	}; // class NumericFormFactorG

} // namespace hig

#endif // __FF_NUM_GPU_CUH__
