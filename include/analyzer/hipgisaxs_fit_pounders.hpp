/**
 *  Project: HipGISAXS
 *
 *  File: hipgisaxs_fit_pounders.hpp
 *  Created: Dec 26, 2013
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *          Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __HIPGISAXS_FIT_POUNDERS_HPP__
#define __HIPGISAXS_FIT_POUNDERS_HPP__

#include <analyzer/analysis_algorithm.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
*/

namespace hig {

  class FitPOUNDERSAlgo : public AnalysisAlgorithm {

    private:
      unsigned int num_obs_;

    public:
      FitPOUNDERSAlgo();
      FitPOUNDERSAlgo(int narg, char** args, ObjectiveFunction* obj, unsigned int algo_num);
      ~FitPOUNDERSAlgo();

      bool run(int argc,char **argv, int, int);
      void print();

  }; /* class FitPOUNDERSAlgo  */

} /* namespace hig */

#endif /* __HIPGISAXS_FIT_POUNDERS_HPP_ */
