/**
 *  Project: HipGISAXS
 *
 *  File: hipgisaxs_fit_lmvm.hpp
 *  Created: Dec 26, 2013
 *
 *  Author: Dinesh Kumar <dkumar@lbl.gov>
 *          Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __HIPGISAXS_FIT_LMVM_HPP__
#define __HIPGISAXS_FIT_LMVM_HPP__

#include <analyzer/analysis_algorithm.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
*/

namespace hig {

  class FitLMVMAlgo : public AnalysisAlgorithm {

    private:
      unsigned int num_obs_;
      std::vector<std::pair<hig::real_t, hig::real_t> > plimits_;   // parameter limits/bounds
      std::vector<hig::real_t> psteps_;                             // parameter step sizes

    public:
      FitLMVMAlgo() { name_= algo_lmvm; max_iter_ = 200; max_hist_ = 200; tol_ = 1e-6; }
      FitLMVMAlgo(int narg, char** args, ObjectiveFunction* obj, unsigned int algo_num) {
        name_= algo_lmvm; obj_func_ = obj; max_iter_ = 200; max_hist_ = 200;
        tol_ = (*obj_func_).analysis_tolerance(algo_num);
        num_obs_ = (*obj_func_).data_size();
        num_params_ = (*obj_func_).num_fit_params();
        x0_ = (*obj_func_).fit_param_init_values();
      } // FitLMVMAlgo()

      ~FitLMVMAlgo() { }

      bool run(int argc,char **argv, int, int);
      void print();

  }; // class FitLMVMAlgo

} // namespace hig

#endif /* __HIPGISAXS_FIT_LMVM_HPP__ */
