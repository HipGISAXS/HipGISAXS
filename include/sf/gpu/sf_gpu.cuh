/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf_gpu.cuh
 *  Created: Oct 16, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Abhinav Sarje <asarje@lbl.gov>
 *              Dinesh Kumar <dkumar@lbl.gov>
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

#ifndef __SF_GPU_CUH__
#define __SF_GPU_CUH__

#include <vector>

#include <common/typedefs.hpp>
#include <common/constants.hpp>

namespace hig {

#ifdef SF_GPU
  class StructureFactorG {
    private:
      unsigned int nqx_;
      unsigned int nqy_;
      unsigned int nqz_;

      // device buffers

      cucomplex_t* sf_;
      float_t* qx_;
      float_t* qy_;
      cucomplex_t* qz_;
      float_t* rot_;
      float_t* repet_;
      float_t* center_;
      float_t* transvec_;

    public:
      StructureFactorG(unsigned int, unsigned int, unsigned int);
      StructureFactorG();
      ~StructureFactorG();

      bool init(unsigned int, unsigned int, unsigned int);
      //bool run_init(const float_t*, const std::vector<float_t>&);
      bool clear();
      bool destroy();

      void grid_size(unsigned int, unsigned int, unsigned int);
      void get_sf(complex_t*&);

      bool compute(std::string, vector3_t, Lattice*, vector3_t, vector3_t,
                   vector3_t, vector3_t, vector3_t
                   #ifdef USE_MPI
                    , woo::MultiNode&, std::string
                   #endif
                   );

  }; // class StructureFactorG

  __global__ void structure_factor_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                                          float_t* qx, float_t* qy, cucomplex_t* qz,
                                          float_t* rot, float_t* repet,
                                          float_t* center, float_t* transvec,
                                          cucomplex_t* sf);

#endif // SF_GPU

} // namespace


#endif // __SF_GPU_CUH__

