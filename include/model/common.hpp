/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: common.hpp
 *  Created: Jun 13, 2012
 *  Modified: Wed 08 Oct 2014 12:13:02 PM PDT
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

#ifndef __MODEL_COMMON_HPP__
#define __MODEL_COMMON_HPP__

#include <common/typedefs.hpp>

namespace hig {

  class RefractiveIndex {
    private:
      real_t delta_;
      real_t beta_;

    public:
      RefractiveIndex(): delta_(0.0), beta_(0.0) { }
      RefractiveIndex(real_t delta, real_t beta): delta_(delta), beta_(beta) { }
      ~RefractiveIndex() { }

      real_t delta() { return delta_; }
      real_t beta() { return beta_; }
      complex_t one_minus_n2(){ return (complex_t(2*delta_, 2*beta_)); }

      void delta(real_t a) { delta_ = a; }
      void beta(real_t a) { beta_ = a; }

      void init();
      void clear();
  }; // class RefractiveIndex

} // namespace hig

#endif /* __MODEL_COMMON_HPP__ */
