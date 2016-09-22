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


  /* convergence test */
  PetscErrorCode convergence_test(Tao tao, void * ctx) {
    PetscErrorCode ierr;
    PetscInt iter;    // iteration number
    PetscReal f,      // function value
              gnorm,  // square of gradient norm (distance)
              cnorm,  // infeasibility
              xdiff;  // trust region step length
    TaoConvergedReason reason;
    PetscReal gatol, grtol, gttol, fatol, frtol;
    PetscInt maxiter;

    TaoView(tao, PETSC_VIEWER_STDOUT_SELF);

    ierr = TaoGetSolutionStatus(tao, &iter, &f, &gnorm, &cnorm, &xdiff, &reason); CHKERRQ(ierr);
    std::cout << "** [status] iteration number                  : " << iter << std::endl
              << "**          function value (f)                : " << f << std::endl
              << "**          gradient norm square (gnorm)      : " << gnorm << std::endl
              << "**          infeasibility (cnorm)             : " << cnorm << std::endl
              << "**          trust region step length (xdiff)  : " << xdiff << std::endl;

    Vec xi; PetscInt size; PetscReal y;
    //VecCreateSeq(PETSC_COMM_SELF, num_params_, &xi);
    TaoGetSolutionVector(tao, &xi);
    VecGetSize(xi, &size);
    std::cout << "**          current parameter vector          : [ ";
    for(PetscInt j = 0; j < size; ++ j) {
      VecGetValues(xi, 1, &j, &y);
      std::cout << y << " ";
    } // for
    std::cout << "]" << std::endl;

    #ifdef PETSC_37
      ierr = TaoGetTolerances(tao, &gatol, &grtol, &gttol); CHKERRQ(ierr);
    #else
      ierr = TaoGetTolerances(tao, &fatol, &frtol, &gatol, &grtol, &gttol); CHKERRQ(ierr);
    #endif
    ierr = TaoGetMaximumIterations(tao, &maxiter); CHKERRQ(ierr);
    std::cout << "++ optimization iteration " << iter << " / " << maxiter << " done." << std::endl;
    //if(gnorm <= gatol) {
    if(f <= gatol) {
      #ifdef PETSC_37
        TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL);
      #else
        TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL);
      #endif
    } else if(iter >= maxiter) {
      TaoSetConvergedReason(tao, TAO_DIVERGED_MAXITS);
    } else {
      TaoSetConvergedReason(tao, TAO_CONTINUE_ITERATING);
    } // if-else
    
    return 0;
  } // convergence_test()


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
