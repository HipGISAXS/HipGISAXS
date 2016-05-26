/**
 *  Project:
 *
 *  File: objective_func.cpp
 *  Created: Feb 02, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <map>

#include <analyzer/objective_func.hpp>

namespace hig {

  /* evaluate function used in pounders */
  //PetscErrorCode EvaluateFunction(TaoSolver tao, Vec X, Vec F, void *ptr) {
  PetscErrorCode EvaluateFunction(Tao tao, Vec X, Vec F, void* ptr) {
    PetscErrorCode ierr;
    PetscReal *x, *f;
    ObjectiveFunction* obj_func = (ObjectiveFunction *) ptr;
    int data_size = obj_func->data_size();
    int num_params = obj_func->num_fit_params();

    ierr = VecGetArray(X, &x); CHKERRQ(ierr);
    ierr = VecGetArray(F, &f); CHKERRQ(ierr);

    // construct the parameter vector
    real_vec_t params;
    for(int i = 0; i < num_params; ++ i) params.push_back(x[i]);

    // run objective function
    std::cout << "++ [pounders] evaluating objective function..." << std::endl;
    real_vec_t temp = (*obj_func)(params);

    // compute the distance and set residual vector
    real_t dist = 0.0;      // square of gradient norm
    for(unsigned int i = 0; i < data_size; ++ i) {
      f[i] = temp[i];
      dist += f[i] * f[i];
    } // for
    std::cout << "** [pounders] distance (gradient norm square) = " << dist << std::endl;
    ierr = VecRestoreArray(F, &f); CHKERRQ(ierr);

    return 0;
  } // EvaluateFunction()


  /* evaluate function used in lmvm */
  PetscReal EvaluateFunction(TaoSolver tao, real_vec_t params, void *ptr) {
    ObjectiveFunction* obj_func = (ObjectiveFunction *) ptr;
    std::cout << "++ [lmvm] evaluating objective function..." << std::endl;
    real_vec_t temp = (*obj_func)(params);
    real_t dist = temp[0];
    std::cout << "** [lmvm] distance = " << dist << std::endl;

    return (PetscReal) dist;
  } // EvaluateFunction()


} // namespace hig
