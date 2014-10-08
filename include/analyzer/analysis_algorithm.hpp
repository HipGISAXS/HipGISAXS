/**
 *  Project:
 *
 *  File: analysis_algorithm.hpp
 *  Created: Feb 02, 2014
 *  Modified: Wed 08 Oct 2014 12:10:40 PM PDT
 */

#ifndef __ANALYSIS_ALGORITHM_HPP__
#define __ANALYSIS_ALGORITHM_HPP__

#include <common/typedefs.hpp>
#include <analyzer/enums.hpp>
#include <analyzer/objective_func.hpp>


namespace hig {

  class AnalysisAlgorithm {

    protected:
      bool is_valid_;
      FittingAlgorithmName name_;    // algorithm name
      ObjectiveFunction* obj_func_;  // the objective function
      float tol_;            // error tolerance
      int max_iter_;          // max num of iterations
      int max_hist_;          // max history
      int num_params_;        // number of parameters
      float_vec_t x0_;        // initial param values
      float_vec_t xn_;        // final param values

    public:
      AnalysisAlgorithm(): max_iter_(200), max_hist_(100), tol_(1e-4), num_params_(0) { }
      /*AnalysisAlgorithm(ObjectiveFunction* func):
          max_iter_(200), max_hist_(100), tol_(1e-4) {
        obj_func_ = func;
        num_params_ = (*obj_func_).num_fit_params();
        x0_ = (*obj_func_).fit_param_init_values();
      } // AnalysisAlgorithm()*/

      ~AnalysisAlgorithm() { }

      bool init_params(const float_vec_t& X0);
      void set_objective_function(ObjectiveFunction* func) { obj_func_= func;  }

      float_vec_t get_param_values() { return xn_; }

      virtual bool run(int, char**, int) = 0;

    }; // class AnalysisAlgorithm

} // namespace hig

#endif // __ANALYSIS_ALGORITHM_HPP__
