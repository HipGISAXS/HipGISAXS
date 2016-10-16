/**
 *  Project: HipGISAXS
 *
 *  File: hipgisaxs_fit_lmvm.cpp
 *  Created: Dec 26, 2013
 *
 *  Author: Dinesh Kumar <dkumar@lbl.gov>
 *          Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>

#include <analyzer/hipgisaxs_fit_lmvm.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
 */

namespace hig {

  // Declare user function for tao
  PetscErrorCode HipGISAXSFormFunctionGradient(TaoSolver, Vec, PetscReal *, Vec, void *);

  // context for lmvm
  typedef struct {
    ObjectiveFunction* obj_func_;
    std::vector<hig::real_t> psteps_;
  } lmvm_ctx_t;


  bool FitLMVMAlgo::run(int argc, char **argv, int algo_num, int img_num) {

    if(!(*obj_func_).set_reference_data(img_num)) return false;

    static char help[] = "** Attempting fitting using LMVM algorithm...";
    std::cout << help << " [ " << img_num << " ]" << std::endl;

    int size, rank;
    PetscErrorCode ierr;
    PetscInitialize(&argc, &argv, (char*) 0, help);

    //std::vector<std::pair<hig::real_t, hig::real_t> > plimits = HiGInput::instance().fit_param_limits();
    //std::vector<hig::real_t> psteps = HiGInput::instance().fit_param_steps();
    plimits_ = (*obj_func_).fit_param_limits();
    psteps_ = (*obj_func_).fit_param_step_values();

    Vec x0, xmin, xmax;
    double y;
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &x0);
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &xmin);
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &xmax);
    for(PetscInt i = 0; i < num_params_; ++ i) {
      y = (double) x0_[i];
      VecSetValues(x0, 1, &i, &y, INSERT_VALUES);
      std::cout << "** " << y << " [ ";
      y = (double) plimits_[i].first;
      VecSetValues(xmin, 1, &i, &y, INSERT_VALUES);
      std::cout << y << " ";
      y = (double) plimits_[i].second;
      VecSetValues(xmax, 1, &i, &y, INSERT_VALUES);
      std::cout << y << " ] " << std::endl;
    } // for

    real_t reg_init = (*obj_func_).analysis_regularization(algo_num);
    real_t reg_factor = reg_init;

    for(int reg_iter = 0; reg_iter < MAX_ITER_REG_; ++ reg_iter) {

    std::cout << "++ Regularization optimization iteration " << reg_iter + 1 << std::endl;
    (*obj_func_).set_regularization(reg_factor);

    Vec f;            // gradient vector
    PetscReal hist[max_hist_], resid[max_hist_];
    PetscInt nhist = max_hist_;

    Tao tao;    // TaoSolver solver context
    TaoConvergedReason reason;

    ierr = VecCreateSeq(PETSC_COMM_SELF, num_obs_, &f);

    ierr = TaoCreate(PETSC_COMM_SELF, &tao);
    //ierr = TaoSetType(tao, TAOLMVM);
    ierr = TaoSetType(tao, TAOBLMVM);

    ierr = TaoSetFromOptions(tao);

    lmvm_ctx_t ctx;
    ctx.obj_func_ = obj_func_;
    ctx.psteps_ = psteps_;

    //ierr = TaoSetObjectiveAndGradientRoutine(tao, HipGISAXSFormFunctionGradient, obj_func_);
    ierr = TaoSetObjectiveAndGradientRoutine(tao, HipGISAXSFormFunctionGradient, (void*) &ctx);
    ierr = TaoSetConvergenceTest(tao, &convergence_test, NULL);

    TaoSetMaximumIterations(tao, max_iter_);
    TaoDefaultMonitor(tao, PETSC_NULL);
    #ifdef PETSC_37
      ierr = TaoSetConvergenceHistory(tao, hist, resid, NULL, NULL, nhist, PETSC_TRUE);
    #elif defined PETSC_36
      ierr = TaoSetHistory(tao, hist, resid, NULL, NULL, nhist, PETSC_TRUE);
    #else
      ierr = TaoSetHistory(tao, hist, resid, NULL, nhist, PETSC_TRUE);
    #endif // PETSC_36

    std::cout << "++ [lmvm] Setting tolerances = " << tol_ << std::endl;
    #ifdef PETSC_37
      ierr = TaoSetTolerances(tao, tol_, tol_, tol_); CHKERRQ(ierr);
    #else
      ierr = TaoSetTolerances(tao, tol_, tol_, tol_, tol_, tol_); CHKERRQ(ierr);
    #endif

    ierr = TaoSetVariableBounds(tao, xmin, xmax); CHKERRQ(ierr);
    ierr = TaoSetInitialVector(tao, x0);

    ierr = TaoSolve(tao); CHKERRQ(ierr);

    TaoView(tao, PETSC_VIEWER_STDOUT_SELF);

    ierr = TaoGetConvergedReason(tao, &reason);
    std::cout << "** [lmvm] converged reason: " << reason << std::endl;

    #if defined PETSC_36 || defined PETSC_37
      TaoGetConvergenceHistory(tao, NULL, NULL, NULL, NULL, &nhist);
    #else
      TaoGetHistory(tao, 0, 0, 0, &nhist);
    #endif // PETSC_36

    PetscPrintf(PETSC_COMM_WORLD, "** [lmvm] history: [ iter\tobj_val\tresidual ]\n");
    for(int j = 0; j < nhist; ++ j)
      PetscPrintf(PETSC_COMM_WORLD, ">> \t%d\t%g\t%g\n", j, hist[j], resid[j]);

    PetscInt iterate; // current iterate number
    PetscReal f_cv,   // current function value
              gnorm,  // square of gradient norm (distance)
              cnorm,  // infeasibility of current solution w.r.t. constraints
              xdiff;  // current trust region step length
    ierr = TaoGetSolutionStatus(tao, &iterate, &f_cv, &gnorm, &cnorm, &xdiff, &reason);
    std::cout << "** [lmvm] converged reason    : " << reason   << std::endl
              << "** [lmvm] iterate number      : " << iterate  << std::endl
              << "** [lmvm] function value      : " << f_cv     << std::endl
              << "** [lmvm] distance            : " << gnorm    << std::endl
              << "** [lmvm] infeasibility       : " << cnorm    << std::endl
              << "** [lmvm] trust region length : " << xdiff    << std::endl;

    TaoGetSolutionVector(tao, &x0); xn_.clear();
    for(PetscInt j = 0; j < num_params_; ++ j) {
      VecGetValues(x0, 1, &j, &y);
      xn_.push_back(y);
    } // for

    std::cout << "** [lmvm] final parameter vector: [ ";
    for(real_vec_t::iterator i = xn_.begin(); i != xn_.end(); ++ i) std::cout << *i << " ";
    std::cout << "]" << std::endl;

    ierr = TaoDestroy(&tao);
    ierr = VecDestroy(&f);

    reg_factor /= 5.0;
    if(reg_factor <= TINY_) break;

    } // for reg_iter

    ierr = VecDestroy(&x0);
    ierr = VecDestroy(&xmin);
    ierr = VecDestroy(&xmax);

    PetscFinalize();

    return true;
  } // FitLMVMAlgo::run()


  void FitLMVMAlgo::print() {
    // ...
  } // FitLMVMAlgo::print()


  PetscErrorCode HipGISAXSFormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ptr) {
    lmvm_ctx_t* ctx = (lmvm_ctx_t*) ptr;
    PetscInt i, j, size;
    PetscReal fxp, fxm, dx = 0.04;
    PetscReal *x, *g;
    real_vec_t xvec, xpvec, xmvec;
    PetscErrorCode ierr;

    /* Get pointers to vector data */
    ierr = VecGetArray(X, &x);
    ierr = VecGetArray(G, &g);
    ierr = VecGetSize(X, &size);

    /* Compute Objective Funtion at X */
    for(i = 0; i < size; ++ i) xvec.push_back(x[i]);
    *f = EvaluateFunction(tao, xvec, ctx->obj_func_);

    /* Evaluate gradients */
    for(i = 0; i < size; ++ i) {
      xpvec = xmvec = xvec;
      real_t pstep = ctx->psteps_[i];
      pstep = (pstep > TINY_) ? pstep : dx;
      xpvec[i] = xvec[i] + 0.5 * pstep;
      xmvec[i] = xvec[i] - 0.5 * pstep;
      fxp = EvaluateFunction(tao, xpvec, ctx->obj_func_);
      fxm = EvaluateFunction(tao, xmvec, ctx->obj_func_);
      g[i] = (fxp - fxm) / dx;
    } // for

    /* Restore vectors */
    ierr = VecRestoreArray(G, &g);

    return 0;
  } // HipGISAXSFormFunctionGradient()

} // namespace hig
