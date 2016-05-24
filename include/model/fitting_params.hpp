/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File:
 *
 *  Author: Dinesh Kumar
 *  Email:  dkumar@lbl.gov
 *  Date create: Apr-19, 2016
 *  Date modified: 
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


#ifndef FITTING_PARAMS__H
#define FITTING_PARAMS__H

#include <config/temp_helpers.hpp>
#include <common/typedefs.hpp>

namespace hig {
  struct ParamSpace {
    real_t min_;
    real_t max_;
    real_t step_;
    ParamSpace(): min_(0), max_(0), step_(-1) {}
    ParamSpace(real_t a, real_t b): min_(a), max_(b), step_(-1) {}
    ParamSpace(real_t a, real_t b, real_t c): min_(a), max_(b), step_(c) {}
    void clear(){ min_ = 0; max_ = 0; step_ = -1; }
  };

  struct FitParam {
    std::string key_;
    std::string variable_;
    ParamSpace range_;
    real_t init_;

    void clear() { key_ = ""; variable_ = ""; range_.clear(); init_ = 0; }
    void init() { clear(); }
  };

  typedef std::map<std::string, FitParam> fit_params_map_t;

  class FittingParams {
    private:
      fit_params_map_t fit_params_;
      AnalysisAlgorithmData fit_algo_;
      AnalysisAlgorithmParamData fit_algo_param_;
      FitReferenceData ref_data_;
    public:
      real_t reference_region_min_x(int i);
      real_t reference_region_min_y(int i);
      real_t reference_region_max_x(int i);
      real_t reference_region_max_y(int i);
  };
}

  

#endif // FITTING_PARAMS__H
