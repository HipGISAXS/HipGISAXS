/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana.hpp
 *  Created: Jul 12, 2012
 *  Modified: Wed 08 Oct 2014 12:13:01 PM PDT
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

#ifndef _FF_ANA_HPP_
#define _FF_ANA_HPP_

#include <vector>

#ifdef USE_MPI
  #include <mpi.h>
  #include <woo/comm/multi_node_comm.hpp>
#endif

#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <common/enums.hpp>
#include <model/shape.hpp>

#ifdef USE_GPU
  #include <ff/gpu/ff_ana_gpu.cuh>
#endif

namespace hig {

  class AnalyticFormFactor {  // make this and numerical ff inherited from class FormFactor ...
    private:
      unsigned int nqx_;
      unsigned int nqy_;
      unsigned int nqz_;

      float_t *rot_;

      #ifdef FF_ANA_GPU
        AnalyticFormFactorG gff_;
      #endif // FF_ANA_GPU

    public:
      AnalyticFormFactor() { }
      ~AnalyticFormFactor() { }

      bool init(vector3_t&, vector3_t&, vector3_t&, std::vector<complex_t> &ff);
      void clear();

      bool compute(ShapeName shape, float_t tau, float_t eta, vector3_t transvec,
            std::vector<complex_t>&,
            shape_param_list_t& params, float_t single_layer_thickness_,
            vector3_t rot1, vector3_t rot2, vector3_t rot3
            #ifdef USE_MPI
              , woo::MultiNode& multi_node, std::string comm_key
            #endif
            );

    private:
      /* compute ff for various shapes */
      bool compute_box(unsigned int nqx, unsigned int nqy, unsigned int nqz,
              std::vector<complex_t>& ff,
              ShapeName shape, shape_param_list_t& params,
              float_t tau, float_t eta, vector3_t &transvec,
              vector3_t &rot1, vector3_t &rot2, vector3_t &rot3);
      bool compute_cylinder(shape_param_list_t&, float_t, float_t,
              std::vector<complex_t>&, vector3_t);
      bool compute_horizontal_cylinder(float_t, float_t, shape_param_list_t&, vector3_t,
              std::vector<complex_t>&);
      bool compute_random_cylinders(shape_param_list_t&, std::vector<complex_t>&,
              float_t, float_t, vector3_t);
      bool compute_sphere(shape_param_list_t&, std::vector<complex_t>&, vector3_t);
      bool compute_prism(shape_param_list_t&, std::vector<complex_t>&,
              float_t, float_t, vector3_t);
      bool compute_prism6(shape_param_list_t&, std::vector<complex_t>&,
              float_t, float_t, vector3_t);
      bool compute_prism3x(shape_param_list_t&, std::vector<complex_t>&,
              float_t, float_t, vector3_t);
      bool compute_truncated_pyramid(shape_param_list_t&, std::vector<complex_t>&, vector3_t);
      complex_t truncated_pyramid_core(int, float_t, complex_t, complex_t, complex_t,
              float_t, float_t, float_t, float_t, vector3_t);
      bool compute_rotation_matrix(int, float_t, vector3_t&, vector3_t&, vector3_t&);
      bool compute_truncated_cone(shape_param_list_t&, float_t, float_t, std::vector<complex_t>&,
              vector3_t);
      bool compute_sawtooth_up();
      bool compute_sawtooth_down();
      bool compute_pyramid();

      /* other helpers */ // check if they should be private ...
      bool param_distribution(ShapeParam&, std::vector<float_t>&, std::vector<float_t>&);
      bool mat_fq_inv_in(unsigned int, unsigned int, unsigned int, complex_vec_t&, float_t);
      bool mat_fq_inv(unsigned int, unsigned int, unsigned int, const complex_vec_t&,
              float_t, complex_vec_t&);
      complex_t fq_inv(complex_t, float_t);
      bool mat_sinc(unsigned int, unsigned int, unsigned int,  const complex_vec_t&, complex_vec_t&);
      bool mat_sinc_in(unsigned int, unsigned int, unsigned int, complex_vec_t&);
      complex_t sinc(complex_t value);
      void compute_meshpoints(const float_t, const float_t, const complex_t, const float_t*,
              complex_t&, complex_t&, complex_t&);

  }; // class AnalyticFormFactor

} // namespace hig

#endif /* _FF_ANA_HPP */
