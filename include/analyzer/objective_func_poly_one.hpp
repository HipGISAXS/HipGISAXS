/**
 *  Project:
 *
 *  File: objective_func.hpp
 *  Created: Feb 02, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __OBJECTIVE_FUNC_POLY_ONE_HPP__
#define __OBJECTIVE_FUNC_POLY_ONE_HPP__

#include <tao.h>
#include <common/typedefs.hpp>
#include <analyzer/ImageData.hpp>
#include <analyzer/distance_functions.hpp>
#include <analyzer/objective_func.hpp>


namespace hig{

  class PolyOneObjectiveFunction : public ObjectiveFunction {
    private:
      real_t q_min_, q_max_, q_step_;    // the q-range
      real_t x1_, x2_;            // function parameters
      real_t x1_init_, x2_init_;        // init parameter values
      real_t x1_min_, x1_max_, x2_min_, x2_max_;    // parameter value ranges

      real_t func(real_t, real_t, real_t);    // the main function

    public:
      PolyOneObjectiveFunction(int, char**, DistanceMeasure*);
      ~PolyOneObjectiveFunction();

      real_vec_t operator()(const real_vec_t&);
      bool set_reference_data(char*);
      bool set_reference_data(int) { }

      int num_fit_params() const { return 2; }
      std::vector <std::string> fit_param_keys() const;
      std::vector <real_pair_t> fit_param_limits() const;
      real_vec_t fit_param_init_values() const;
      unsigned int data_size() const { return ref_data_->size(); }
  }; // class PolyOneObjectiveFunction


} // namespace hig

#endif // __OBJECTIVE_FUNC_HPP__
