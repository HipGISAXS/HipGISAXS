/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: qgrid.hpp
 *  Created: Jun 17, 2012
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

#ifndef __QGRID_HPP__
#define __QGRID_HPP__

#include <vector>

#include <common/globals.hpp>
#include <model/compute_params.hpp>

namespace hig {

  class QGrid {       
    /* data types for qgrid data */
    typedef std::vector<real_t> qvec_t;
    typedef std::vector<complex_t> cqvec_t;
    /* iterators for qgrid data types */
    typedef qvec_t::const_iterator qvec_iter_t;
    typedef cqvec_t::const_iterator cqvec_iter_t;

    private:
      int nrow_;
      int ncol_;

      vector2_t qmin_;
      vector2_t qmax_;
      qvec_t qx_;
      qvec_t qy_;
      qvec_t qz_;
      qvec_t alpha_;
      cqvec_t qz_extended_;

      /* singleton */
      QGrid() { }
      QGrid(const QGrid&);
      QGrid& operator=(const QGrid&);

      bool pixel_to_kspace(vector2_t, real_t, real_t, real_t, real_t, vector2_t, vector3_t&);
      vector3_t pixel_to_kspace(vector2_t, real_t, real_t, real_t, real_t, vector2_t);
      bool kspace_to_pixel();    // not implemented yet ...

    public:
      static QGrid& instance() {
        static QGrid qgrid;
        return qgrid;
      } // instance()

      /* create the Q-grid */
      bool create(const ComputeParams &, real_t, real_t, int);
      bool create_qz_extended(real_t, real_t, complex_t); 
      bool create_test();

      /* for fitting */
      bool update(unsigned int, unsigned int, real_t, real_t, real_t, real_t,
            real_t, real_t, real_t, int);

      /* temporary for steepest descent fitting */
      bool create_z_cut(real_t, real_t, real_t, real_t);

      /* sizes */
      int nqx() const          { return qx_.size();          }
      int nqy() const          { return qy_.size();          }
      int nqz() const          { return qz_.size();          }
      int nqz_extended() const { return qz_extended_.size(); }
      int nrows() const        { return nrow_;               }
      int ncols() const        { return ncol_;               }
      int nalpha() const       { return alpha_.size();       }
      vector2_t qmin() const   { return qmin_;               }
      vector2_t qmax() const   { return qmax_;               }


      /* NOTE: delta are not constants in q-space */
      /* deltas */ 
      //real_t delta_x() const { return (qx_[qx_.size() - 1] - qx_[0]) / (qx_.size() - 1); }
      //real_t delta_y() const { return (qy_[qy_.size() - 1] - qy_[0]) / (qy_.size() - 1); }
      //complex_t delta_z() const { return (qz_[qz_.size() - 1] - qz_[0]) / (qz_.size() - 1); }

      /* value accessors */
      real_t qx(int i) const { return qx_[i]; }  // do some error checking also ...
      real_t qy(int i) const { return qy_[i]; }
      real_t qz(int i) const { return qz_[i]; }
      complex_t qz_extended(int i) const { return qz_extended_[i]; }
      real_t qx(unsigned int i) const { return qx_[i]; }  // do some error checking also ...
      real_t qy(unsigned int i) const { return qy_[i]; }
      real_t qz(unsigned int i) const { return qz_[i]; }
      complex_t qz_extended(unsigned int i) const { return qz_extended_[i]; }
      real_t alpha(unsigned int i) const { return alpha_[i]; }

      /* iterator helpers */
      qvec_iter_t qx_begin() const { return qx_.begin(); }
      qvec_iter_t qy_begin() const { return qy_.begin(); }
      qvec_iter_t qz_begin() const { return qz_.begin(); }
      cqvec_iter_t qz_extended_begin() const { return qz_extended_.begin(); }
      qvec_iter_t qx_end() const { return qx_.end(); }
      qvec_iter_t qy_end() const { return qy_.end(); }
      qvec_iter_t qz_end() const { return qz_.end(); }
      cqvec_iter_t qz_extended_end() const { return qz_extended_.end(); }

      /* temporary, just for performance testing purposes */
      real_t* qx(void) { return &qx_[0]; }
      real_t* qy(void) { return &qy_[0]; }
      real_t* qz(void) { return &qz_[0]; }
      complex_t* qz_extended(void) { return &qz_extended_[0]; }

      /* debug */
      void save (const char *);

  }; // class QGrid

} // namespace hig

#endif // __QGRID_HPP__ */
