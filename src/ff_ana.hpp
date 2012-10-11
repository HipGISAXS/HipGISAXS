/***
  *  $Id: ff_ana.hpp 38 2012-08-09 23:01:20Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana.hpp
  *  Created: Jul 12, 2012
  *  Modified: Tue 09 Oct 2012 11:54:15 AM PDT
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
			/*float_t* mesh_qx_;
			float_t* mesh_qy_;
			float_t* mesh_qz_;*/

		public:
			AnalyticFormFactor() { }
			~AnalyticFormFactor() { }

			bool init(vector3_t&, vector3_t&, vector3_t&, std::vector<complex_t> &ff);
			void clear();
			bool compute(ShapeName shape, float_t tau, float_t eta, vector3_t transvec,
						std::vector<complex_t>&,
						shape_param_list_t& params, float_t single_layer_thickness_,
						vector3_t rot1, vector3_t rot2, vector3_t rot3,
						MPI::Intracomm& world_comm);

		private:
			/* computing ff for various shapes */
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
			bool compute_prism6();
			bool compute_prism3x();
			bool compute_sawtooth_up();
			bool compute_sawtooth_down();
			bool compute_pyramid();
			bool compute_truncated_pyramid();
			bool compute_truncated_cone();

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
