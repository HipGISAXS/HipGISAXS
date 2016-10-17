/**
 *  Project: HipGISAXS
 *
 *  File: hipgisaxs_compute_objective.cpp
 *  Created: May 28, 2016
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>

#include <analyzer/hipgisaxs_compute_objective.hpp>


namespace hig {

  /* constructor to set objective function */
  ComputeObjectiveFunction::ComputeObjectiveFunction(int narg, char** args,
                                   ObjectiveFunction* obj, unsigned int algo_num) {
    name_ = algo_none_pounders;
    obj_func_ = obj;
    num_obs_ = (*obj_func_).data_size();
    num_params_ = (*obj_func_).num_fit_params();
    x0_ = (*obj_func_).fit_param_init_values();
  } // FitPOUNDERSAlgo::FitPOUNDERSAlgo()


  bool ComputeObjectiveFunction::run(int argc, char **argv, int algo_num, int img_num) {
    if(!(*obj_func_).set_reference_data(img_num)) return false;

    static char help[] = "** computing objective function...";
    std::cout << help << " [ " << img_num << " ]" << std::endl;

    // initial vector = x0_
    std::vector<std::pair<hig::real_t, hig::real_t> > plimits = obj_func_->fit_param_limits();

    std::cout << "** [none] parameter vector: [ ";
    for(real_vec_t::iterator i = x0_.begin(); i != x0_.end(); ++ i) std::cout << *i << " ";
    std::cout << "]" << std::endl;

    // run objective function
    std::cout << "++ [none] evaluating objective function..." << std::endl;
    real_vec_t temp = (*obj_func_)(x0_);

    // compute the distance and set residual vector
    real_t dist = 0.0;      // square of gradient norm
    for(unsigned int i = 0; i < num_obs_; ++ i) {
      dist += temp[i] * temp[i];
    } // for
    std::cout << "** [none] distance (objective function value) = " << dist << std::endl;

    return true;
  } // ComputeObjectiveFunction::run()

} // namespace hig

