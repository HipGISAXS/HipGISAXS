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

  /* evaluate function used in pounders algo */
  PetscErrorCode EvaluateFunction(TaoSolver tao, Vec X, Vec F, void *ptr) {
    PetscErrorCode ierr;
    PetscReal *x, *ff;
    int data_size = ((ObjectiveFunction*)ptr)->data_size();
    int num_params = ((ObjectiveFunction*)ptr)->num_fit_params();

    ierr = VecGetArray(X, &x); CHKERRQ(ierr);
    ierr = VecGetArray(F, &ff); CHKERRQ(ierr);

    real_vec_t params;
    for(int i = 0; i < num_params; ++ i) params.push_back(x[i]);
    // run objective function
    std::cout << "++ [pounders] evaluating objective function..." << std::endl;
    real_vec_t temp = (*(ObjectiveFunction*)ptr)(params);
    real_t err = 0.0;
    for(int i = 0; i < data_size; ++ i) {
      ff[i] = temp[i];
      err += ff[i] * ff[i];
    } // for

    ierr = VecRestoreArray(F, &ff); CHKERRQ(ierr);
    std::cout << "** [pounders] distance = " << err << std::endl;

    return 0;
  } // EvaluateFunction()


  /* evaluate function used in lmvm */
  PetscReal EvaluateFunction(TaoSolver tao, real_vec_t X, void *ptr) {
    std::cout << "++ [lmvm] evaluating objective function..." << std::endl;
    real_vec_t temp = (*(ObjectiveFunction*) ptr)(X);
    real_t dist = temp[0];
    std::cout << "** Distance = " << dist << std::endl;

    return (PetscReal) dist;
  } // EvaluateFunction()


} // namespace hig
