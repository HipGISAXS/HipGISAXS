/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf.hpp
 *  Created: Jun 18, 2012
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

#ifndef __HIG_SF_HPP__
#define __HIG_SF_HPP__

#ifdef USE_MPI
#include <woo/comm/multi_node_comm.hpp>
#endif // USE_MPI

#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <model/structure.hpp>
#include <numerics/matrix.hpp>

#ifdef USE_GPU
  #include <sf/gpu/sf_gpu.cuh>
#endif


namespace hig {
  
  class StructureFactor {
    private:
      complex_t *sf_;
      unsigned int ny_;
      unsigned int nz_;
      StructureType type_;

      #ifdef SF_GPU
        StructureFactorG gsf_;
      #endif

    public:
      StructureFactor();
      ~StructureFactor();
      void putStructureType(StructureType d){ type_ = d; }
      void clear(void);

      bool compute_structure_factor(std::string, vector3_t, Lattice*, vector3_t, vector3_t,
                      RotMatrix_t &, std::shared_ptr<Paracrystal> , std::shared_ptr<PercusYevick> 
                      #ifdef USE_MPI
                        , woo::MultiNode&, std::string
                      #endif
                      );
      #ifdef SF_GPU
      bool compute_structure_factor_gpu(std::string, vector3_t, Lattice*, vector3_t, vector3_t,
                      RotMatrix_t &
                      #ifdef USE_MPI
                        , woo::MultiNode&, std::string
                      #endif
                      );
      #endif
      complex_t & operator[](unsigned int i) const { return sf_[i]; }

      void save_sf(const std::string& filename);
      void save(const char* filename);

      StructureFactor & operator=(const StructureFactor & rhs);
      StructureFactor & operator+=(const StructureFactor & rhs);
      StructureFactor & operator*(const real_t val) {
        for (int i = 0; i < nz_; i++) sf_[i] *= val;
        return *this;
      }

      //! 1-D Paracrystal
      bool paracrystal1d(real_t, real_t, real_t);

      //! 2-D Paracrystal Cubic
      bool paracrystal2d_cubic(real_t, real_t, real_t, real_t, real_t);

      //! 2-D Paracrystal Hex
      bool paracrystal2d_hex(real_t, real_t, real_t, real_t, real_t);

      //! Percus-Yevik
      real_t PYFcn(real_t, real_t, real_t, real_t);
      bool percus_yevik(real_t, real_t, int);

      // only for testing - remove it ...
      // complex_t* sf() { return sf_; }
 
  }; // class StructureFactor

} // namespace hig


#endif /* __HIG_SF_HPP__ */
