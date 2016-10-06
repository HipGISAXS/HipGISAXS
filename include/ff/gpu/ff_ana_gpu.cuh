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
#include <numerics/matrix.hpp>

namespace hig {

  class AnalyticFormFactorG {
    private:
      unsigned int nqy_;
      unsigned int nqz_;

      // device buffers

      real_t* qx_;
      real_t* qy_;
      cucomplex_t* qz_;
      cucomplex_t* ff_;
      bool construct_output_ff(std::vector<complex_t>&);

    public:
      AnalyticFormFactorG(unsigned int, unsigned int);
      AnalyticFormFactorG();
      ~AnalyticFormFactorG();

      bool init(unsigned int, unsigned int);
      bool run_init(const real_t*, const std::vector<real_t>&);
      bool run_init(const real_t*, const std::vector<real_t>&, std::vector<complex_t>&);
      bool clear();
      bool destroy();

      void grid_size(unsigned int, unsigned int);

      bool compute_cube(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);

      bool compute_box(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_cube(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_cylinder(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_horizontal_cylinder(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_random_cylinders(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_sphere(const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_prism3(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_prism6(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_prism3x();
      bool compute_sawtooth_up();
      bool compute_sawtooth_down();
      bool compute_pyramid(const real_t, const real_t,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const std::vector<real_t>&, const std::vector<real_t>&,
                const RotMatrix_t &, const std::vector<real_t>&, std::vector<complex_t>&);
      bool compute_truncated_cone();

  }; // class AnalyticFormFactorG
} // namespace
#endif // _FF_ANA_GPU_CUH_
