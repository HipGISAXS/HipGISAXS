/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: FitPOUNDERSAlgo.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 03 Feb 2014 12:35:02 PM PST
 *  Description: Defines the TAO POUNDERS Algorithm
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

#ifndef _FITPOUNDERSALGO_HPP_
#define _FITPOUNDERSALGO_HPP_

#include <analyzer/analysis_algorithm.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
*/

namespace hig{

  class FitPOUNDERSAlgo : public AnalysisAlgorithm {
  private:
	  unsigned int num_obs_;


  public:
    FitPOUNDERSAlgo() { name_= algo_pounders; max_iter_ = 200; max_hist_ = 100; tol_ = 1e-10; }
    FitPOUNDERSAlgo(ObjectiveFunction* obj) {
		name_= algo_pounders; obj_func_ = obj; max_iter_ = 200; max_hist_ = 100; tol_ = 1e-10;
		num_obs_ = (*obj_func_).data_size();
		num_params_ = (*obj_func_).num_fit_params();
		x0_ = (*obj_func_).fit_param_init_values();
	} // FitPOUNDERSAlgo()

    ~FitPOUNDERSAlgo(){}

    bool run(int argc,char **argv, int);
    void print();

  }; /* class FitPOUNDERSAlgo  */

} /* namespace hig */

#endif /* FITPOUNDERSALGO_HPP_ */
