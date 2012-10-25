/***
  *  $Id: ff_ana.hpp 38 2012-08-09 23:01:20Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana.hpp
  *  Created: Jul 12, 2012
  *  Modified: Sat 20 Oct 2012 05:02:03 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _FF_ANA_HPP_
#define _FF_ANA_HPP_

#include <vector>
#include <mpi.h>

#include "typedefs.hpp"
#include "globals.hpp"
#include "enums.hpp"
#include "shape.hpp"
#include "ff_ana_gpu.cuh"

namespace hig {

	class AnalyticFormFactor {	// make this and numerical ff inherited from class FormFactor ...
		private:
			//std::vector<complex_t> ff_;

			unsigned int nqx_;
			unsigned int nqy_;
			unsigned int nqz_;

			std::vector<complex_t> mesh_qx_;
			std::vector<complex_t> mesh_qy_;
			std::vector<complex_t> mesh_qz_;

			float_t *rot_;

			AnalyticFormFactorG ff_gpu_;

		public:
			AnalyticFormFactor()/*: mqx_(NULL), mqy_(NULL), mqz_(NULL)*/ { }
			~AnalyticFormFactor() {
				/*if(mqx_ != NULL) delete[] mqx_;
				if(mqy_ != NULL) delete[] mqy_;
				if(mqz_ != NULL) delete[] mqz_;*/
			} // ~AnalyticFormFactor()

			bool init(vector3_t&, vector3_t&, vector3_t&, std::vector<complex_t> &ff);
			void clear();
			bool compute(ShapeName shape, float_t tau, float_t eta, vector3_t transvec,
						std::vector<complex_t>&,
						shape_param_list_t& params, float_t single_layer_thickness_,
						vector3_t rot1, vector3_t rot2, vector3_t rot3,
						MPI::Intracomm& world_comm);

		private:
			/* compute ff for various shapes */
			bool compute_box(unsigned int nqx, unsigned int nqy, unsigned int nqz,
							std::vector<complex_t>& ff,
							ShapeName shape, shape_param_list_t& params,
							float_t tau, float_t eta, vector3_t &transvec,
							vector3_t &rot1, vector3_t &rot2, vector3_t &rot3);
			bool compute_cylinder(shape_param_list_t&, float_t, float_t, std::vector<complex_t>&, vector3_t);
			bool compute_horizontal_cylinder(shape_param_list_t&, vector3_t, std::vector<complex_t>&);
			bool compute_random_cylinders();
			bool compute_sphere(shape_param_list_t&, std::vector<complex_t>&, vector3_t);
			bool compute_prism(shape_param_list_t&, std::vector<complex_t>&, float_t, float_t, vector3_t);
			bool compute_prism6(shape_param_list_t&, std::vector<complex_t>&, float_t, float_t, vector3_t);
			bool compute_prism3x();
			bool compute_sawtooth_up();
			bool compute_sawtooth_down();
			bool compute_pyramid();
			bool compute_truncated_pyramid();
			bool compute_truncated_cone();

			/* compute ff on gpu */
			/*bool compute_box_gpu(unsigned int nqx, unsigned int nqy, unsigned int nqz,
							std::vector<complex_t>& ff,
							ShapeName shape, shape_param_list_t& params,
							float_t tau, float_t eta, vector3_t &transvec,
							vector3_t &rot1, vector3_t &rot2, vector3_t &rot3);
			bool compute_cylinder_gpu(shape_param_list_t&, float_t, float_t,
									std::vector<complex_t>&, vector3_t);
			bool compute_horizontal_cylinder_gpu(shape_param_list_t&, vector3_t, std::vector<complex_t>&);
			bool compute_random_cylinders_gpu();
			bool compute_sphere_gpu(shape_param_list_t&, std::vector<complex_t>&, vector3_t);
			bool compute_prism_gpu(shape_param_list_t&, std::vector<complex_t>&,
								float_t, float_t, vector3_t);
			bool compute_prism6_gpu(shape_param_list_t&, std::vector<complex_t>&,
								float_t, float_t, vector3_t);
			bool compute_prism3x_gpu();
			bool compute_sawtooth_up_gpu();
			bool compute_sawtooth_down_gpu();
			bool compute_pyramid_gpu();
			bool compute_truncated_pyramid_gpu();
			bool compute_truncated_cone_gpu(); */

			/* other helpers */ // check if they should be private ...
			bool param_distribution(ShapeParam&, std::vector<float_t>&, std::vector<float_t>&);
			bool mat_fq_inv_in(unsigned int, unsigned int, unsigned int, complex_vec_t&, float_t);
			bool mat_fq_inv(unsigned int, unsigned int, unsigned int, const complex_vec_t&,
							float_t, complex_vec_t&);
			complex_t fq_inv(complex_t, float_t);
			bool mat_sinc(unsigned int, unsigned int, unsigned int,	const complex_vec_t&, complex_vec_t&);
			bool mat_sinc_in(unsigned int, unsigned int, unsigned int, complex_vec_t&);
			float_t sinc(complex_t value);


	}; // class AnalyticFormFactor

} // namespace hig

#endif /* _FF_ANA_HPP */
