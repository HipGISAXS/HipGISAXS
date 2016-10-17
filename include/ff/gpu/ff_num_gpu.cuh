/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_gpu.cuh
 *  Created: Nov 05, 2011
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __FF_NUM_GPU_CUH__
#define __FF_NUM_GPU_CUH__

#include <vector>

#include <common/typedefs.hpp>
#include <numerics/matrix.hpp>

namespace hig {

  /**
   * Class for computing Form Factor in either single or double precision on a single GPU.
   */
  //template<typename real_t, typename complex_t>
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

      NumericFormFactorG():  // called when not using this GPU version
        block_cuda_(0), block_cuda_t_(0), block_cuda_y_(0), block_cuda_z_(0) { }

      ~NumericFormFactorG() {}

      /* Spherical Q_grid */
      unsigned int compute_exact_triangle(triangle_t *, int, 
              cucomplex_t * &, 
              int, real_t *, real_t *, int,
              cucomplex_t *, RotMatrix_t &, real_t &);

      unsigned int compute_approx_triangle(std::vector<real_t> &, 
              cucomplex_t * &,
              int, real_t *, real_t *, int, 
              cucomplex_t *, RotMatrix_t &, real_t &);

      
      /* original */
      unsigned int compute_form_factor(int,
          std::vector<real_t> &shape_def, std::vector<short int> &axes,
          cucomplex_t* &ff,
          real_t* &qx_h, int nqx,
          real_t* &qy_h, int nqy,
          cucomplex_t* &qz_h, int nqz,
          real_t* &rot,
          real_t&, real_t&, real_t&
          #ifdef FINDBLOCK
            , const int, const int, const int, const int
          #endif
          );

      /* with double buffering - default */
      unsigned int compute_form_factor_db(int,
          std::vector<real_t> &shape_def, std::vector<short int> &axes,
          cucomplex_t* &ff,
          real_t* &qx_h, int nqx,
          real_t* &qy_h, int nqy,
          cucomplex_t* &qz_h, int nqz,
          real_t* &rot,
          real_t&, real_t&, real_t&
          #ifdef FINDBLOCK
            , const int, const int, const int, const int
          #endif
          );

      /* with fused kernels and double buffering */
      unsigned int compute_form_factor_db_fused(int,
          std::vector<real_t> &shape_def, std::vector<short int> &axes,
          cucomplex_t* &ff,
          real_t* &qx_h, int nqx,
          real_t* &qy_h, int nqy,
          cucomplex_t* &qz_h, int nqz,
          real_t* &rot,
          real_t&, real_t&, real_t&
          #ifdef FINDBLOCK
            , const int, const int, const int, const int
          #endif
          );

      /* with fused kernels and triple buffering */
      unsigned int compute_form_factor_kb_fused(int,
          std::vector<real_t> &shape_def, std::vector<short int> &axes,
          cucomplex_t* &ff,
          real_t* &qx_h, int nqx,
          real_t* &qy_h, int nqy,
          cucomplex_t* &qz_h, int nqz, int k,
          real_t* &rot,
          real_t&, real_t&, real_t&
          #ifdef FINDBLOCK
            , const int, const int, const int, const int
          #endif
          );

      /* with double buffering and optimized memory (incomplete ... TODO)*/
      unsigned int compute_form_factor_db_mem(int,
          std::vector<real_t> &shape_def, std::vector<short int> &axes,
          cucomplex_t* &ff,
          real_t* &qx_h, int nqx,
          real_t* &qy_h, int nqy,
          cucomplex_t* &qz_h, int nqz,
          real_t* &rot,
          real_t&, real_t&, real_t&
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
