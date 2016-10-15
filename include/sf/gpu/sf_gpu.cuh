/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf_gpu.cuh
 *  Created: Oct 16, 2012
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

#ifndef __SF_GPU_CUH__
#define __SF_GPU_CUH__

#include <vector>

#include <common/typedefs.hpp>
#include <common/constants.hpp>
#include <model/structure.hpp>
#include <numerics/matrix.hpp>

namespace hig {

#ifdef SF_GPU
  class StructureFactorG {
    private:
      bool inited_;
      unsigned int nqx_;
      unsigned int nqy_;
      unsigned int nqz_;
      unsigned int nrow_;
      unsigned int ncol_;

      // device buffers

      cucomplex_t* sf_;
      real_t* qx_;
      real_t* qy_;
      cucomplex_t* qz_;
      real_t* rot_;
      real_t* repet_;
      real_t* center_;
      real_t* transvec_;

    public:
      StructureFactorG(unsigned int, unsigned int, unsigned int);
      StructureFactorG();
      ~StructureFactorG();

      bool init(unsigned int, unsigned int, unsigned int);
      //bool run_init(const real_t*, const std::vector<real_t>&);
      bool clear();
      bool destroy();

      void grid_size(unsigned int, unsigned int, unsigned int);
      void get_sf(complex_t*&);

      bool compute(std::string, vector3_t, Lattice*, vector3_t, vector3_t,
                   RotMatrix_t
                   //#ifdef USE_MPI
                   // , woo::MultiNode&, std::string
                   //#endif
                   );

  }; // class StructureFactorG

  __global__ void structure_factor_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                                          real_t* qx, real_t* qy, cucomplex_t* qz,
                                          real_t* rot, real_t* repet,
                                          real_t* center, real_t* transvec,
                                          cucomplex_t* sf);
  __global__ void structure_factor_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                                          real_t* qx, real_t* qy, cucomplex_t* qz,
                                          RotMatrix_t r, vector3_t a, vector3_t b, vector3_t c,
                                          real_t* repet,
                                          real_t* center, real_t* transvec,
                                          cucomplex_t* sf);

#endif // SF_GPU

} // namespace


#endif // __SF_GPU_CUH__

