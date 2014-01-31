/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: FitPOUNDERSAlgo.hpp
 *  Created: Dec 26, 2013
 *  Modified: Thu 30 Jan 2014 08:43:39 PM PST
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

#include <analyzer/AnaAlgorithm.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
*/

namespace hig{

  class FitPOUNDERSAlgo : public AnaAlgorithm{
  private:


  public:
    FitPOUNDERSAlgo(){type_= fit_pounders; max_iter_=200; max_hist_=100; tol_=1e-4; }
    FitPOUNDERSAlgo(float_t* params){type_= fit_pounders  ;params_= params; max_iter_=200; max_hist_=100; tol_=1e-4;}
    ~FitPOUNDERSAlgo(){}
    virtual bool init(){return false;}
    virtual bool parse_params(){return false;}
    //virtual Analysis run();
    virtual Analysis run(int argc,char **argv);
    virtual void print();

  }; /* class FitPOUNDERSAlgo  */

} /* namespace hig */

#endif /* FITPOUNDERSALGO_HPP_ */
