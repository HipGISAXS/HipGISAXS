/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: fitting_params.hpp
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


#ifndef __FITTING_PARAMS_H__
#define __FITTING_PARAMS_H__

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
  }; // struct ParamSpace

  struct FitParam {
    std::string key_;
    std::string variable_;
    ParamSpace range_;
    real_t init_;

    void clear() { key_ = ""; variable_ = ""; range_.clear(); init_ = 0; }
    void init() { clear(); }
  }; // struct FitParam

  typedef std::map<std::string, FitParam> fit_params_map_t;

  // TODO ... 
  class FittingParams {
    private:
      fit_params_map_t fit_params_;
      AnalysisAlgorithmData fit_algo_;
      AnalysisAlgorithmParamData fit_algo_param_;
      FitReferenceData ref_data_;
    public:
      real_t reference_region_min_x(int i) const { return 0; }
      real_t reference_region_min_y(int i) const { return 0; }
      real_t reference_region_max_x(int i) const { return 0; }
      real_t reference_region_max_y(int i) const { return 0; }

      std::vector<std::string> fit_param_keys() const {
        std::vector<std::string> keys;
        return keys;
      }
      std::vector<std::pair<real_t, real_t> > fit_param_limits() const {
        std::vector<std::pair<real_t, real_t> > lims;
        return lims;
      }
      std::vector<real_t> fit_param_step_values() const {
        std::vector<real_t> steps;
        return steps;
      }
      std::string reference_data_path() const {
        return std::string();
      }
      std::string reference_data_mask() const {
        return std::string();
      }
      int num_fit_params() const {
        return 0;
      }
      std::vector<real_t> fit_param_init_values() const {
        std::vector<real_t> init_vals;
        return init_vals;
      }
      real_t param_space_mean(const std::string& key) { return 0; }
  }; // class FittingParams

}

#endif // __FITTING_PARAMS_H__
