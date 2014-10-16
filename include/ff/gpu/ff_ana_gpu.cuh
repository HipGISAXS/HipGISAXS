/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_gpu.cuh
 *  Created: Oct 16, 2012
 *  Modified: Wed 08 Oct 2014 12:15:06 PM PDT
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

#ifndef _FF_ANA_GPU_CUH_
#define _FF_ANA_GPU_CUH_

#include <vector>

#include <common/typedefs.hpp>
#include <common/constants.hpp>

namespace hig {

  class AnalyticFormFactorG {
    private:
      unsigned int nqx_;
      unsigned int nqy_;
      unsigned int nqz_;

      // device buffers

      float_t* transvec_;
      float_t* rot_;

      float_t* qx_;
      float_t* qy_;
      cucomplex_t* qz_;
      cucomplex_t* ff_;

      unsigned int b_nqx_;    // hyperblock size
      unsigned int b_nqy_;
      unsigned int b_nqz_;
      unsigned int nb_x_;
      unsigned int nb_y_;
      unsigned int nb_z_;
      cucomplex_t* ff_buff_d_;  // device buffer
      cucomplex_t* ff_buff_h_;  // host buffer

      bool construct_output_ff(std::vector<complex_t>&);
      bool compute_hyperblock_size(unsigned int, unsigned int);
      bool move_ff_buff_to_host_ff(std::vector<complex_t>&,
                unsigned int, unsigned int, unsigned int,
                unsigned int, unsigned int, unsigned int,
                int, int, int);

    public:
      AnalyticFormFactorG(unsigned int, unsigned int, unsigned int);
      AnalyticFormFactorG();
      ~AnalyticFormFactorG();

      bool init(unsigned int, unsigned int, unsigned int);
      bool run_init(const float_t*, const std::vector<float_t>&);
      bool run_init(const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
      bool clear();
      bool destroy();

      void grid_size(unsigned int, unsigned int, unsigned int);

      bool compute_box(const float_t, const float_t,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_cylinder(const float_t, const float_t,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_horizontal_cylinder(const float_t, const float_t,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_random_cylinders(const float_t, const float_t,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_sphere(const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_prism3(const float_t, const float_t,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_prism6(const float_t, const float_t,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const std::vector<float_t>&, const std::vector<float_t>&,
                const float_t*, const std::vector<float_t>&, std::vector<complex_t>&);
            bool compute_prism3x();
            bool compute_sawtooth_up();
            bool compute_sawtooth_down();
            bool compute_pyramid();
            bool compute_truncated_pyramid();
            bool compute_truncated_cone();

  }; // class AnalyticFormFactorG

} // namespace


#endif // _FF_ANA_GPU_CUH_

