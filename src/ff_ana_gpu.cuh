/***
  *  $Id$
  *
  *  Project:
  *
  *  File: ff_ana_gpu.cuh
  *  Created: Oct 16, 2012
  *  Modified: Thu 18 Oct 2012 04:28:46 PM PDT
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

		public:
			AnalyticFormFactorG(unsigned int, unsigned int, unsigned int);
			AnalyticFormFactorG();
			~AnalyticFormFactorG();

			void grid_size(unsigned int, unsigned int, unsigned int);

			bool compute_box(/*unsigned int nqx, unsigned int nqy, unsigned int nqz,
                            std::vector<complex_t>& ff,
                            ShapeName shape, shape_param_list_t& params,
                            float_t tau, float_t eta, vector3_t &transvec,
                            vector3_t &rot1, vector3_t &rot2, vector3_t &rot3*/);
            bool compute_cylinder(/*shape_param_list_t&, float_t, float_t,
                                    std::vector<complex_t>&, vector3_t*/);
            bool compute_horizontal_cylinder(/*shape_param_list_t&, vector3_t, std::vector<complex_t>&*/);
            bool compute_random_cylinders();
            /*bool compute_sphere(const std::vector<float_t>& r,
								const std::vector<float_t>& distr_r,
								const std::vector<complex_t>& mesh_qx,
								const std::vector<complex_t>& mesh_qy,
								const std::vector<complex_t>& mesh_qz,
								const std::vector<float_t>& transvec,
								std::vector<complex_t>& ff);*/
            bool compute_sphere(const std::vector<float_t>& r,
								const std::vector<float_t>& distr_r,
								const cucomplex_t* mesh_qx,
								const cucomplex_t* mesh_qy,
								const cucomplex_t* mesh_qz,
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

