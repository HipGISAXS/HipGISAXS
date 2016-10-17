/**
 *  Project: HipGISAXS
 *
 *  File: hipgisaxs_fit_pounders.cpp
 *  Created: Dec 26, 2013
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *          Slim Chourou <stchourou@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>

#include <analyzer/hipgisaxs_fit_pounders.hpp>


namespace hig {

  /* pounders: Finds the nonlinear least-squares solution to the model
   *           y = exp[-b1 * x] / (b2 + b3 * x) + e
   */

  /*
   * Concepts: TAO^Solving a system of nonlinear equations, nonlinear least squares
   * Routines: TaoCreate();
   * Routines: TaoSetType();
   * Routines: TaoSetSeparableObjectiveRoutine();
   * Routines: TaoSetJacobianRoutine();
   * Routines: TaoSetInitialVector();
   * Routines: TaoSetFromOptions();
   * Routines: TaoSetConvergenceHistory(); TaoGetConvergenceHistory();
   * Routines: TaoSolve();
   * Routines: TaoView(); TaoDestroy();
   * Processors: 1
   */

  /* default constructor */
  FitPOUNDERSAlgo::FitPOUNDERSAlgo() {
        name_ = algo_pounders;
        max_iter_ = 200;
        max_hist_ = 200;
        tol_ = 1e-6;
  } // FitPOUNDERSAlgo::FitPOUNDERSAlgo()


  /* constructor to set objective function */
  FitPOUNDERSAlgo::FitPOUNDERSAlgo(int narg, char** args,
                                   ObjectiveFunction* obj, unsigned int algo_num) {
    name_ = algo_pounders;
    obj_func_ = obj;
    max_iter_ = 200;
    max_hist_ = 200;
    tol_ = (*obj_func_).analysis_tolerance(algo_num);
    num_obs_ = (*obj_func_).data_size();
    num_params_ = (*obj_func_).num_fit_params();
    x0_ = (*obj_func_).fit_param_init_values();
  } // FitPOUNDERSAlgo::FitPOUNDERSAlgo()


  FitPOUNDERSAlgo::~FitPOUNDERSAlgo() { }


  bool FitPOUNDERSAlgo::run(int argc, char **argv, int algo_num, int img_num) {
    if(!(*obj_func_).set_reference_data(img_num)) return false;

    static char help[] = "** Attempting fitting using Pounders algorithm...";
    std::cout << help << " [ " << img_num << " ]" << std::endl;

    const int str_max = 100;

    real_t pdelta, pnpmax, pgqt;
    bool isdelta = false;
    int newnarg = argc;
    //char* newargs[argc + 3];     // possibly add arguments for the tao routines
    //char newargs[argc + 3][50];
    char** newargs = new char*[argc + 6];     // possibly add arguments for the tao routines
    for(int i = 0; i < argc; ++ i) {
      newargs[i] = new char[str_max];
      strncpy(newargs[i], argv[i], str_max);
    } // for

    if(img_num >= 0) {

      if((*obj_func_).analysis_algo_param(algo_num, "pounders_delta", pdelta)) {
        std::stringstream arg1; arg1 << "-tao_pounders_delta";
        //newargs[newnarg] = new char[arg1.str().size() + 1];
        newargs[newnarg] = new char[str_max];
        strncpy(newargs[newnarg], arg1.str().c_str(), str_max);
        ++ newnarg;
        std::stringstream arg2; arg2 << pdelta;
        //newargs[newnarg] = new char[arg2.str().size() + 1];
        newargs[newnarg] = new char[str_max];
        strncpy(newargs[newnarg], arg2.str().c_str(), str_max);
        ++ newnarg;
        isdelta = true;
      } else {
        std::cerr << "warning: default pounders_delta being used" << std::endl;
      } // if-else

      if((*obj_func_).analysis_algo_param(algo_num, "pounders_npmax", pnpmax)) {
        std::stringstream arg1; arg1 << "-tao_pounders_npmax";
        //newargs[newnarg] = new char[arg1.str().size() + 1];
        newargs[newnarg] = new char[str_max];
        strncpy(newargs[newnarg], arg1.str().c_str(), str_max);
        ++ newnarg;
        std::stringstream arg2; arg2 << pnpmax;
        //newargs[newnarg] = new char[arg2.str().size() + 1];
        newargs[newnarg] = new char[str_max];
        strncpy(newargs[newnarg], arg2.str().c_str(), str_max);
        ++ newnarg;
      } else {
        std::cerr << "warning: default pounders_npmax being used" << std::endl;
      } // if-else

      if((*obj_func_).analysis_algo_param(algo_num, "pounders_gqt", pgqt)) {
        std::stringstream arg1; arg1 << "-tao_pounders_gqt";
        //newargs[newnarg] = new char[arg1.str().size() + 1];
        newargs[newnarg] = new char[str_max];
        strncpy(newargs[newnarg], arg1.str().c_str(), str_max);
        ++ newnarg;
        std::stringstream arg2; arg2 << pgqt;
        //newargs[newnarg] = new char[arg2.str().size() + 1];
        newargs[newnarg] = new char[str_max];
        strncpy(newargs[newnarg], arg2.str().c_str(), str_max);
        ++ newnarg;
      } else {
        std::cerr << "warning: default pounders_gqt being used" << std::endl;
      } // if-else

    } // if

    // temp ...
    for(int i = 0; i < newnarg; ++ i) std::cout << newargs[i] << std::endl;

    int size, rank;
    PetscErrorCode ierr;
    PetscInitialize(&newnarg, &newargs, (char*) 0, help);

    // need to free the newargs memory ... TODO ...

    std::vector<std::pair<hig::real_t, hig::real_t> > plimits = (*obj_func_).fit_param_limits();

    Vec x0,         // initial parameter vector
        xmin, xmax; // parameter min and max limits (bounds)
    double y;
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &x0);
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &xmin);
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &xmax);
    for(PetscInt i = 0; i < num_params_; ++ i) {
      y = (double) x0_[i];
      VecSetValues(x0, 1, &i, &y, INSERT_VALUES);
      std::cout << "** " << y << "\t[ ";
      y = (double) plimits[i].first;
      VecSetValues(xmin, 1, &i, &y, INSERT_VALUES);
      std::cout << y << "\t";
      y = (double) plimits[i].second;
      if(isdelta) y += pdelta;
      VecSetValues(xmax, 1, &i, &y, INSERT_VALUES);
      std::cout << y << "\t]" << std::endl;
    } // for

    real_t reg_init = (*obj_func_).analysis_regularization(algo_num);
    real_t reg_factor = reg_init;

    for(int reg_iter = 0; reg_iter < MAX_ITER_REG_; ++ reg_iter) {

      std::cout << "++ Regularization optimization iteration " << reg_iter + 1 << std::endl;

      (*obj_func_).set_regularization(reg_factor);

      Vec f;
      PetscReal hist[max_hist_], resid[max_hist_];
      PetscInt nhist = max_hist_;

      Tao tao;
      TaoConvergedReason reason;

      // allocate vectors
      ierr = VecCreateSeq(PETSC_COMM_SELF, num_obs_, &f);

      // create TAO solver with pounders
      ierr = TaoCreate(PETSC_COMM_SELF, &tao);
      ierr = TaoSetType(tao, TAOPOUNDERS); CHKERRQ(ierr);

      // check for command line options
      ierr = TaoSetFromOptions(tao);

      // set objective function
      ierr = TaoSetSeparableObjectiveRoutine(tao, f, EvaluateFunction, (void*) obj_func_);
      // set jacobian function
      // ierr = TaoSetJacobianRoutine(tao, J, J, EvaluateJacobian, (void*) &user); CHKERRQ(ierr);

      // set the convergence test function
      ierr = TaoSetConvergenceTest(tao, &convergence_test, NULL);

      TaoSetMaximumIterations(tao, max_iter_);
      #ifdef PETSC_37
        ierr = TaoSetConvergenceHistory(tao, hist, resid, NULL, NULL, max_hist_, PETSC_TRUE); CHKERRQ(ierr);
      #elif defined PETSC_36
        ierr = TaoSetHistory(tao, hist, resid, NULL, NULL, max_hist_, PETSC_TRUE); CHKERRQ(ierr);
      #else
        ierr = TaoSetHistory(tao, hist, resid, NULL, max_hist_, PETSC_TRUE); CHKERRQ(ierr);
      #endif // PETSC_36

      std::cout << "++ [pounders] Setting tolerances = " << tol_ << std::endl;
      #ifdef PETSC_37
        ierr = TaoSetTolerances(tao, tol_, tol_, tol_); CHKERRQ(ierr);
      #else
        ierr = TaoSetTolerances(tao, tol_, tol_, tol_, tol_, tol_); CHKERRQ(ierr);
      #endif

      ierr = TaoSetVariableBounds(tao, xmin, xmax); CHKERRQ(ierr);

      // set the initial parameter vector (initial guess)
      ierr = TaoSetInitialVector(tao, x0); CHKERRQ(ierr);

      // perform the solve
      ierr = TaoSolve(tao); CHKERRQ(ierr);

      // temporary, to check the details
      TaoView(tao, PETSC_VIEWER_STDOUT_SELF);

      ierr = TaoGetConvergedReason(tao, &reason);
      std::cout << "** [pounders] converged reason: " << reason << std::endl;

      #if defined PETSC_36 || defined PETSC_37
        ierr = TaoGetConvergenceHistory(tao, NULL, NULL, NULL, NULL, &nhist); CHKERRQ(ierr);
      #else
        TaoGetHistory(tao, 0, 0, 0, &nhist);
      #endif

      PetscPrintf(PETSC_COMM_WORLD, "** [pounders] history: [ iter\tobj_val\tresidual ]\n");
      for(int i = 0; i < nhist; ++ i)
        PetscPrintf(PETSC_COMM_WORLD, ">> \t%d\t%g\t%g\n", i, hist[i], resid[i]);

      PetscInt iterate; // current iterate number
      PetscReal f_cv,   // current function value
                gnorm,  // square of gradient norm (distance)
                cnorm,  // infeasibility of current solution w.r.t. constraints
                xdiff;  // current trust region step length
      ierr = TaoGetSolutionStatus(tao, &iterate, &f_cv, &gnorm, &cnorm, &xdiff, &reason);
      std::cout << "** [pounders] converged reason    : " << reason   << std::endl
                << "** [pounders] iterate number      : " << iterate  << std::endl
                << "** [pounders] function value      : " << f_cv     << std::endl
                << "** [pounders] distance            : " << gnorm    << std::endl
                << "** [pounders] infeasibility       : " << cnorm    << std::endl
                << "** [pounders] trust region length : " << xdiff    << std::endl;

      // obtain the parameter solution vector
      TaoGetSolutionVector(tao, &x0); xn_.clear();
      for(PetscInt j = 0; j < num_params_; ++ j) {
        VecGetValues(x0, 1, &j, &y);
        xn_.push_back(y);
      } // for

      std::cout << "** [pounders] final parameter vector: [ ";
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
  } // FitPOUNDERSAlgo::run()


  /*PetscErrorCode convergence_test(Tao tao, void * ctx) {
    PetscErrorCode ierr;
    PetscInt iter;    // iteration number
    PetscReal f,      // function value
              gnorm,  // square of gradient norm (distance)
              cnorm,  // infeasibility
              xdiff;  // trust region step length
    TaoConvergedReason reason;

    ierr = TaoGetSolutionStatus(tao, &iter, &f, &gnorm, &cnorm, &xdiff, &reason); CHKERRQ(ierr);
    if(gnorm <= tol_) {
      #ifdef PETSC_37
        TaoSetConvergedReason(tao, TAO_CONVERGED_ATOL);
      #else
        TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL);
      #endif
    } else if(iter >= max_iter_) {
      TaoSetConvergedReason(tao, TAO_DIVERGED_MAXITS);
    } else {
      TaoSetConvergedReason(tao, TAO_CONTINUE_ITERATING);
    } // if-else
    
    return 0;
  } // convergence_test()*/


  void FitPOUNDERSAlgo::print() {
    // ...
  } // FitPOUNDERSAlgo::print()

} // namespace hig

