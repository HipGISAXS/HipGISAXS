/**
 *  Project:
 *
 *  File: analysis_algorithm.hpp
 *  Created: Feb 02, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 */

#ifndef __ANALYSIS_ALGORITHM_HPP__
#define __ANALYSIS_ALGORITHM_HPP__

#include <common/typedefs.hpp>
//#include <analyzer/enums.hpp>
#include <analyzer/objective_func.hpp>


namespace hig {

  class AnalysisAlgorithm {

    protected:
      bool is_valid_;
      FittingAlgorithmName name_;   // algorithm name
      ObjectiveFunction* obj_func_; // the objective function
      real_t tol_;                  // error tolerance
      int max_iter_;                // max num of iterations
      int max_hist_;                // max history
      int num_params_;              // number of parameters
      real_vec_t x0_;               // initial param values
      real_vec_t xn_;               // final param values

    public:
      AnalysisAlgorithm(): max_iter_(200), max_hist_(200), tol_(1e-6), num_params_(0) { }
      ~AnalysisAlgorithm() { }

      bool init_params(const real_vec_t& X0);
      void set_objective_function(ObjectiveFunction* func) { obj_func_= func; }

      real_vec_t get_param_values() const { return xn_; }
      real_t tolerance() const { return tol_; }

      virtual bool run(int, char**, int, int) = 0;

    }; // class AnalysisAlgorithm

} // namespace hig

#endif // __ANALYSIS_ALGORITHM_HPP__
