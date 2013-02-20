/***
  *  $Id$
  *
  *  Project:
  *
  *  File: ff_ana_gpu.cuh
  *  Created: Oct 16, 2012
  *  Modified: Wed 20 Feb 2013 02:22:48 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _FF_ANA_GPU_CUH_
#define _FF_ANA_GPU_CUH_

#include <vector>

#include "typedefs.hpp"
#include "constants.hpp"

namespace hig {

	class AnalyticFormFactorG {
		private:
			unsigned int nqx_;
			unsigned int nqy_;
			unsigned int nqz_;

			// device buffers

			float_t* qx_;
			float_t* qy_;
			cucomplex_t* qz_;
			cucomplex_t* ff_;

			float_t* transvec_;
			float_t* rot_;

			bool construct_output_ff(std::vector<complex_t>&);

		public:
			AnalyticFormFactorG(unsigned int, unsigned int, unsigned int);
			AnalyticFormFactorG();
			~AnalyticFormFactorG();

			bool init(unsigned int, unsigned int, unsigned int);
			bool run_init(const float_t*, const std::vector<float_t>&);
			bool destroy();

			void grid_size(unsigned int, unsigned int, unsigned int);

			bool compute_box(const float_t, const float_t,
							const std::vector<float_t>&, const std::vector<float_t>&,
							const std::vector<float_t>&, const std::vector<float_t>&,
							const std::vector<float_t>&, const std::vector<float_t>&,
							//const float_t*, const float_t*, const cucomplex_t*,
							const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_cylinder(const float_t, const float_t, const std::vector<float_t>&,
									const std::vector<float_t>&, const std::vector<float_t>&,
									const std::vector<float_t>&,
									//const float_t*, const float_t*, const cucomplex_t*,
									const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_horizontal_cylinder(const float_t, const float_t,
									const std::vector<float_t>&, const std::vector<float_t>&,
									const std::vector<float_t>&, const std::vector<float_t>&,
									const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_random_cylinders();
            bool compute_sphere(const std::vector<float_t>& r,
								const std::vector<float_t>& distr_r,
								//const float_t* qx,
								//const float_t* qy,
								//const cucomplex_t* qz,
								const float_t* rot,
								const std::vector<float_t>& transvec,
								std::vector<complex_t>& ff);
            bool compute_prism(/*shape_param_list_t&, std::vector<complex_t>&,
                                float_t, float_t, vector3_t*/);
            bool compute_prism6(/*shape_param_list_t&, std::vector<complex_t>&,
                                float_t, float_t, vector3_t*/);
            bool compute_prism3x();
            bool compute_sawtooth_up();
            bool compute_sawtooth_down();
            bool compute_pyramid();
            bool compute_truncated_pyramid();
            bool compute_truncated_cone();

	}; // class AnalyticFormFactorG

} // namespace


#endif // _FF_ANA_GPU_CUH_

