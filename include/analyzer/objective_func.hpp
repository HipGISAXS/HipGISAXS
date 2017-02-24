/**
 *  Project: HipGISAXS
 *
 *  File: objective_func.hpp
 *  Created: Feb 02, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __OBJECTIVE_FUNC_HPP__
#define __OBJECTIVE_FUNC_HPP__

#if defined PETSC_35 || defined PETSC_36 || defined PETSC_37
# include <petsctao.h>
# define TaoSolver Tao
# define TaoInitialize(...)
# define TaoFinalize(...)
# define TaoSolverTerminationReason TaoConvergedReason
# define TaoGetTerminationReason TaoGetConvergedReason
# if defined PETSC_36 || defined PETSC_37
#   define TaoSetHistory TaoSetConvergenceHistory
#   define TaoGetHistory TaoGetConvergenceHistory
# endif // PETSC_36
#else   // PETSC_34 or lower
# include <tao.h>
#endif // PETSC_35 || PETSC_36 || PETSC_37

#include <hipgisaxs.hpp>
#include <analyzer/ImageData.hpp>
#include <analyzer/distance_functions.hpp>
#include <woo/comm/multi_node_comm.hpp>

namespace hig {

  /**
   * The abstract objective function class
   */
  class ObjectiveFunction {

    protected:

      DistanceMeasure* pdist_;  // distance function
      ImageData* ref_data_;     // reference data
      bool mask_set_;           // whether mask data is set or not
      uint_vec_t mask_data_;    // mask with 0s and 1s
      //real_vec_t curr_dist_;  // current computed distance output

    public:

      virtual real_vec_t operator()(const real_vec_t&) = 0;
      virtual int num_fit_params() const = 0;
      virtual std::vector <std::string> fit_param_keys() const = 0;
      virtual std::vector <real_pair_t> fit_param_limits() const = 0;
      virtual std::vector <real_t> fit_param_step_values() const { }
      virtual real_vec_t fit_param_init_values() const = 0;
      virtual real_t analysis_tolerance(int) const = 0;
      virtual real_t analysis_regularization(int) const = 0;
      virtual bool set_distance_measure(DistanceMeasure*) = 0;
      virtual bool set_reference_data(int) = 0;
      virtual bool set_reference_data(char*) = 0;
      virtual bool set_regularization(real_t) = 0;
      virtual unsigned int data_size() const = 0;
      real_t* get_reference_data() { return ref_data_->data(); }
      unsigned int* get_mask_data() { return &(mask_data_[0]); }
      //virtual unsigned int n_par() const { }
      //virtual unsigned int n_ver() const { }

      virtual bool analysis_algo_param(int, std::string, real_t&) const = 0;
      virtual std::string param_pathprefix() const = 0;
      virtual std::string runname() const = 0;

      #ifdef USE_MPI
        virtual woo::MultiNode* multi_node_comm() = 0;
        virtual bool update_sim_comm(std::string) { }
      #endif

      // for testing
      //virtual bool update_params(const real_vec_t&);
      virtual bool simulate_and_set_ref(const real_vec_t&) = 0;

  }; // class ObjectiveFunction


  // lmvm
  PetscReal EvaluateFunction(TaoSolver , real_vec_t , void *);
  // pounders
  PetscErrorCode EvaluateFunction(TaoSolver, Vec, Vec, void *);
  PetscErrorCode EvaluateJacobian(TaoSolver, Vec, Mat *, Mat *, MatStructure *, void *);
  PetscErrorCode convergence_test(Tao, void *);

} // namespace hig

#endif // __OBJECTIVE_FUNC_HPP__
