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

#include <numerics/matrix.hpp>

#ifdef USE_GPU
  #include <ff/gpu/ff_ana_gpu.cuh>
#endif

namespace hig {

  class AnalyticFormFactor {  // make this and numerical ff inherited from class FormFactor ...
    private:
      unsigned int nqx_;
      unsigned int nqy_;
      unsigned int nqz_;
      RotMatrix_t rot_;
      

      #ifdef FF_ANA_GPU
        AnalyticFormFactorG gff_;
      #endif // FF_ANA_GPU

    public:
      AnalyticFormFactor() { }
      ~AnalyticFormFactor() { }

      bool init(RotMatrix_t &, std::vector<complex_t> &);
      void clear();

      bool compute(ShapeName , real_t , real_t , vector3_t ,
            std::vector<complex_t>&,
            shape_param_list_t& , real_t , RotMatrix_t &
            #ifdef USE_MPI
              , woo::MultiNode& multi_node, std::string comm_key
            #endif
            );

      // optimized for CPU
      #ifdef FF_CPU_OPT
        bool cylinder_opt(std::vector<real_t>&, std::vector<real_t>&, vector3_t, std::vector<complex_t>&);
        bool cylinder_kernel_opt(std::vector<real_t>&, std::vector<real_t>&, vector3_t, std::vector<complex_t>&);
        complex_t ff_cylinder_kernel_opt(complex_t qpar, complex_t qz, real_t radius, real_t height);
        void ff_cylinder_kernel_opt_vec(complex_t*, const complex_t*, const complex_t*, real_t, real_t, complex_t*);
      #endif

      // helper kernels
      complex_t ff_cylinder_kernel(complex_t qpar, complex_t qz, real_t radius, real_t height);

    private:
      /* compute ff for various shapes */
      bool compute_box(unsigned int nqx, unsigned int nqy, unsigned int nqz,
              std::vector<complex_t>& ff,
              ShapeName shape, shape_param_list_t& params,
              real_t tau, real_t eta, vector3_t &transvec);
      bool compute_cylinder(shape_param_list_t&, real_t, real_t,
              std::vector<complex_t>&, vector3_t);
      bool compute_horizontal_cylinder(real_t, real_t, shape_param_list_t&, vector3_t,
              std::vector<complex_t>&);
      bool compute_random_cylinders(shape_param_list_t&, std::vector<complex_t>&,
              real_t, real_t, vector3_t);
      bool compute_sphere(shape_param_list_t&, std::vector<complex_t>&, vector3_t);
      bool compute_prism(shape_param_list_t&, std::vector<complex_t>&,
              real_t, real_t, vector3_t);
      bool compute_prism6(shape_param_list_t&, std::vector<complex_t>&,
              real_t, real_t, vector3_t);
      bool compute_prism3x(shape_param_list_t&, std::vector<complex_t>&,
              real_t, real_t, vector3_t);
      bool compute_pyramid(shape_param_list_t&, std::vector<complex_t>&,
              real_t, real_t, vector3_t);
      bool compute_rotation_matrix(int, real_t, vector3_t&, vector3_t&, vector3_t&);
      bool compute_truncated_cone(shape_param_list_t&, real_t, real_t, std::vector<complex_t>&, 
              vector3_t);
      bool compute_truncated_sphere(shape_param_list_t&, real_t, real_t, std::vector<complex_t>&, 
              vector3_t);
      bool compute_sawtooth_up();
      bool compute_sawtooth_down();

      /* other helpers */ // check if they should be private ...
      bool param_distribution(ShapeParam&, std::vector<real_t>&, std::vector<real_t>&);
      bool mat_fq_inv_in(unsigned int, unsigned int, unsigned int, complex_vec_t&, real_t);
      bool mat_fq_inv(unsigned int, unsigned int, unsigned int, const complex_vec_t&,
              real_t, complex_vec_t&);
      complex_t fq_inv(complex_t, real_t);
      bool mat_sinc(unsigned int, unsigned int, unsigned int,  const complex_vec_t&, complex_vec_t&);
      bool mat_sinc_in(unsigned int, unsigned int, unsigned int, complex_vec_t&);
      complex_t sinc(complex_t value);
      void compute_meshpoints(const real_t, const real_t, const complex_t, const real_t*,
              complex_t&, complex_t&, complex_t&);

  }; // class AnalyticFormFactor

} // namespace hig

#endif /* _FF_ANA_HPP */
