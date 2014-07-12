/**
 *  Project:
 *
 *  File: objective_func.hpp
 *  Created: Feb 02, 2014
 *  Modified: Fri 11 Jul 2014 09:05:18 AM PDT
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
			float_t q_min_, q_max_, q_step_;		// the q-range
			float_t x1_, x2_;						// function parameters
			float_t x1_init_, x2_init_;				// init parameter values
			float_t x1_min_, x1_max_, x2_min_, x2_max_;		// parameter value ranges

			float_t func(float_t, float_t, float_t);		// the main function

		public:
			PolyOneObjectiveFunction(int, char**, DistanceMeasure*);
			~PolyOneObjectiveFunction();

			float_vec_t operator()(const float_vec_t&);
			bool set_reference_data(char*);
			bool set_reference_data(int) { }

			int num_fit_params() const { return 2; }
			std::vector <std::string> fit_param_keys() const;
			std::vector <float_pair_t> fit_param_limits() const;
			float_vec_t fit_param_init_values() const;
			unsigned int data_size() const { return ref_data_->size(); }
	}; // class PolyOneObjectiveFunction


} // namespace hig

#endif // __OBJECTIVE_FUNC_HPP__