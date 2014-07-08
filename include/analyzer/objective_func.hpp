/**
 *  Project:
 *
 *  File: objective_func.hpp
 *  Created: Feb 02, 2014
 *  Modified: Wed 26 Feb 2014 10:02:54 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __OBJECTIVE_FUNC_HPP__
#define __OBJECTIVE_FUNC_HPP__

#include <analyzer/ImageData.hpp>
#include <analyzer/typedefs.hpp>
#include <analyzer/distance_functions.hpp>
#include <hipgisaxs.hpp>

#include <tao.h>

namespace hig{

	class ObjectiveFunction {
		protected:
			DistanceMeasure* pdist_;	// distance function
			ImageData* ref_data_;		// reference data
			float_vec_t curr_dist_;		// current computed distance output

		public:
			virtual float_vec_t operator()(const float_vec_t&) = 0;
			virtual int num_fit_params() const = 0;
			virtual std::vector <std::string> fit_param_keys() const = 0;
			virtual std::vector <float_pair_t> fit_param_limits() const = 0;
			virtual std::vector <float_t> fit_param_step_values() const { }
			virtual float_vec_t fit_param_init_values() const = 0;
			virtual bool set_distance_measure(DistanceMeasure*) = 0;
			virtual bool set_reference_data(int) = 0;
			virtual bool set_reference_data(char*) = 0;
			virtual unsigned int data_size() const = 0;
			//virtual unsigned int n_par() const { }
			//virtual unsigned int n_ver() const { }

			// for testing
			//virtual bool update_params(const float_vec_t&);
			virtual bool simulate_and_set_ref(const float_vec_t&) = 0;
	}; // class ObjectiveFunction


	class HipGISAXSObjectiveFunction : public ObjectiveFunction {
		private:
			HipGISAXS hipgisaxs_;		// the hipgisaxs object
			unsigned int n_par_;		// nqy
			unsigned int n_ver_;		// nqz

		public:
			HipGISAXSObjectiveFunction(int, char**, DistanceMeasure*);
			HipGISAXSObjectiveFunction(int, char**, std::string);
			~HipGISAXSObjectiveFunction();

			bool set_distance_measure(DistanceMeasure*);
			bool set_reference_data(int);
			bool set_reference_data(char*) { }
			float_vec_t operator()(const float_vec_t&);

			int num_fit_params() const { return hipgisaxs_.num_fit_params(); }
			unsigned int n_par() const { return n_par_; }
			unsigned int n_ver() const { return n_ver_; }
			unsigned int data_size() const { return n_par_ * n_ver_; }
			std::vector <std::string> fit_param_keys() const { return hipgisaxs_.fit_param_keys(); }
			std::vector <float_pair_t> fit_param_limits() const { return hipgisaxs_.fit_param_limits(); }
			float_vec_t fit_param_step_values() const { return hipgisaxs_.fit_param_step_values(); }
			float_vec_t fit_param_init_values() const { return hipgisaxs_.fit_param_init_values(); }

			// for testing
			bool simulate_and_set_ref(const float_vec_t&);
	}; // class HipGISAXSObjectiveFunction


	PetscReal EvaluateFunction(TaoSolver , float_vec_t , void *);
	PetscErrorCode EvaluateFunction(TaoSolver , Vec , Vec , void *);
	PetscErrorCode EvaluateJacobian(TaoSolver , Vec , Mat *, Mat *, MatStructure*,void *);

} // namespace hig

#endif // __OBJECTIVE_FUNC_HPP__
