/**
 *  Project:
 *
 *  File: objective_func_hipgisaxs.hpp
 *  Created: Feb 02, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __OBJECTIVE_FUNC_HIPGISAXS_HPP__
#define __OBJECTIVE_FUNC_HIPGISAXS_HPP__

#include <analyzer/objective_func.hpp>
#include <hipgisaxs.hpp>

namespace hig {

  class HipGISAXSObjectiveFunction : public ObjectiveFunction {
    private:
      HipGISAXS hipgisaxs_;   // hipgisaxs object
      unsigned int n_par_;    // nqy
      unsigned int n_ver_;    // nqz

      real_t* mean_data_;     // buffer to store simulated data with mean parameter vector
      real_t reg_alpha_;      // alpha for regularization

    public:
      HipGISAXSObjectiveFunction(int, char**, DistanceMeasure*);
      HipGISAXSObjectiveFunction(int, char**, std::string);
      ~HipGISAXSObjectiveFunction();

      bool set_distance_measure(DistanceMeasure*);
      bool set_reference_data(int);
      bool set_reference_data(char*) { }
      bool set_mean_data(void);
      bool set_regularization(real_t r) { reg_alpha_ = r; }
      bool read_mask_data(string_t);
      bool read_edf_mask_data(string_t);

      real_vec_t operator()(const real_vec_t&);

      int num_fit_params() const { return hipgisaxs_.num_fit_params(); }
      unsigned int n_par() const { return n_par_; }
      unsigned int n_ver() const { return n_ver_; }
      unsigned int data_size() const { return n_par_ * n_ver_; }
      real_t analysis_tolerance(int n) const { return hipgisaxs_.analysis_tolerance(n); }
      real_t analysis_regularization(int n) const { return hipgisaxs_.analysis_regularization(n); }
      std::vector <std::string> fit_param_keys() const { return hipgisaxs_.fit_param_keys(); }
      std::vector <real_pair_t> fit_param_limits() const { return hipgisaxs_.fit_param_limits(); }
      real_vec_t fit_param_step_values() const { return hipgisaxs_.fit_param_step_values(); }
      real_vec_t fit_param_init_values() const { return hipgisaxs_.fit_param_init_values(); }
      int num_analysis_algos() const { return hipgisaxs_.num_analysis_algos(); }
      FittingAlgorithmName analysis_algo(int n) const { return hipgisaxs_.analysis_algo(n); }
      bool analysis_algo_param(int n, const std::string name, real_t& value) const {
        return hipgisaxs_.analysis_algo_param(n, name, value); } 
      FittingDistanceMetric analysis_distance_metric(int n) const {
        return hipgisaxs_.analysis_distance_metric(n); }

      std::string param_pathprefix() const { return hipgisaxs_.path(); }
      std::string runname() const { return hipgisaxs_.runname(); }

      #ifdef USE_MPI
        woo::MultiNode* multi_node_comm() { return hipgisaxs_.multi_node_comm(); }
        bool update_sim_comm(woo::comm_t comm) { return hipgisaxs_.update_sim_comm(comm); }
      #endif

      // for testing
      bool simulate_and_set_ref(const real_vec_t&);
  }; // class HipGISAXSObjectiveFunction

} // namespace hig

#endif // __OBJECTIVE_FUNC_HIPGISAXS_HPP__
