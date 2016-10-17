/**
 *  Project:
 *
 *  File: hipgisaxs_fit_bruteforce.cpp
 *  Created: Feb 07, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */


#include <analyzer/hipgisaxs_fit_bruteforce.hpp>

namespace hig {

  BruteForceOptimization::BruteForceOptimization(int narg, char** args, ObjectiveFunction* obj, unsigned int algo_num) {
    name_ = algo_bruteforce;
    obj_func_ = obj;
    num_params_ = (*obj_func_).num_fit_params();
    params_ = (*obj_func_).fit_param_keys();
    real_vec_t temp_step = (*obj_func_).fit_param_step_values();
    std::vector<std::pair<real_t, real_t> > temp_minmax = (*obj_func_).fit_param_limits();

    for(int i = 0; i < num_params_; ++ i) {
      x_min_.push_back(temp_minmax[i].first);
      x_max_.push_back(temp_minmax[i].second);
      real_t step = (temp_step[i] < 0) ? 1 : temp_step[i];
      x_step_.push_back(step);
    } // for

    tol_ = 0;        // not used
    max_iter_ = 0;   // not used
    max_hist_ = 0;   // not used

    x0_ = x_min_;    // initial vector is min of all
    xn_.clear();
    error_list_.clear();
  } // BruteForceOptimization::BruteForceOptimization()

  BruteForceOptimization::~BruteForceOptimization() {

  } // BruteForceOptimization::~BruteForceOptimization()


  bool BruteForceOptimization::run(int narg, char** args, int algo_num, int img_num) {
    std::cout << "Running Brute Force Optimization ..." << std::endl;

    (*obj_func_).set_reference_data(img_num);

    real_vec_t curr_vec;
    loop_over_param(0, curr_vec);

    if(narg == 3) save_history(args[2]);
    else save_history("brute_force.dat");

    return true;
  } // BruteForceOptimization::run()


  void BruteForceOptimization::loop_over_param(int param_i, real_vec_t& curr_vec) {
    if(param_i == num_params_) {
      real_vec_t curr_err = (*obj_func_)(curr_vec);
      real_t tot_err = 0;      // if a vector, then sum over all values
      for(int i = 0; i < curr_err.size(); ++ i) tot_err += curr_err[i];
      error_list_.push_back(std::pair<real_vec_t, real_t>(curr_vec, tot_err));
    } else {
      for(real_t p = x_min_[param_i]; p <= x_max_[param_i]; p += x_step_[param_i]) {
        curr_vec.push_back(p);
        loop_over_param(param_i + 1, curr_vec);
        curr_vec.pop_back();
      } // for
    } // if-else
  } // BruteForceOptimization::loop_over_param()


  bool BruteForceOptimization::save_history(char* filename) {
    std::ofstream f(filename);
    if(!f.is_open()) {
      std::cerr << "error: could not open history file " << filename << " for writing" << std::endl;
      return false;
    } // if
    for(int i = 0; i < error_list_.size(); ++ i) {
      for(int j = 0; j < num_params_; ++ j) {
        f << error_list_[i].first[j] << "\t";
      } // for
      f << error_list_[i].second << std::endl;
    } // for
    f.close();
  } // BruteForceOptimization::save_history()

} // namespace hig
