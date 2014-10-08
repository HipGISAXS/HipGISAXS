/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: FitLMVMAlgo.hpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 08 Oct 2014 12:12:05 PM PDT
 *  Description: Defines the TAO LMVM Algorithm
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _FITLMVMALGO_HPP_
#define _FITLMVMALGO_HPP_

#include <analyzer/analysis_algorithm.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
*/

namespace hig{

  class FitLMVMAlgo : public AnalysisAlgorithm {
  private:
    unsigned int num_obs_;


  public:
    FitLMVMAlgo() { name_= algo_lmvm; max_iter_ = 50; max_hist_ = 100; tol_ = 1e-10; }
    FitLMVMAlgo(ObjectiveFunction* obj) {
    name_= algo_lmvm; obj_func_ = obj; max_iter_ = 50; max_hist_ = 100; tol_ = 1e-10;
    num_obs_ = (*obj_func_).data_size();
    num_params_ = (*obj_func_).num_fit_params();
    x0_ = (*obj_func_).fit_param_init_values();
  } // FitLMVMAlgo()

    ~FitLMVMAlgo(){}

    bool run(int argc,char **argv, int);
    void print();

  }; /* class FitLMVMAlgo  */

} /* namespace hig */

#endif /* FITLMVMALGO_HPP_ */
