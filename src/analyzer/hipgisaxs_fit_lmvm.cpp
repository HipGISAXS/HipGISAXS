/**
 *  Project: HipGISAXS
 *
 *  File: AnaAlgorithm.cpp
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


  bool FitLMVMAlgo::run(int argc, char **argv, int img_num) {

    if(!(*obj_func_).set_reference_data(img_num)) return false;

    static char help[] = "Running LMVM fitting...";
    std::cout << help << " [" << img_num << "]" << std::endl;

    PetscErrorCode ierr;
    int size, rank;
    PetscInitialize(&argc, &argv, (char*) 0, help);
    TaoInitialize(&argc, &argv, (char*) 0, help);

    Vec x0;
    double y[1] = { 0 };
    VecCreateSeq(PETSC_COMM_SELF, num_params_, &x0);
    for(PetscInt i = 0; i < num_params_; ++ i) {
        y[0] = (double) x0_[i];
        VecSetValues(x0, 1, &i, y, INSERT_VALUES);
    } // for

    PetscReal zero = 0.0;
    PetscReal hist[max_hist_], resid[max_hist_];
    PetscInt nhist = max_hist_;
    Vec G;            // gradient vector
    TaoSolver tao;    // TaoSolver solver context
    TaoSolverTerminationReason reason;

    ierr = VecCreateSeq(PETSC_COMM_SELF, num_obs_, &G);
    ierr = TaoCreate(PETSC_COMM_SELF, &tao);
    //ierr = TaoSetType(tao, "tao_lmvm");
    ierr = TaoSetType(tao, TAOLMVM);

    ierr = TaoSetObjectiveAndGradientRoutine(tao, HipGISAXSFormFunctionGradient, obj_func_);
    ierr = TaoSetFromOptions(tao);
    TaoSetMaximumIterations(tao, max_iter_);

    TaoSetTolerances(tao, tol_, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    TaoDefaultMonitor(tao, PETSC_NULL);

    ierr = TaoSetInitialVector(tao, x0);
    #ifdef PETSC_36
      ierr = TaoSetHistory(tao, hist, resid, NULL, NULL, nhist, PETSC_TRUE);
    #else
      ierr = TaoSetHistory(tao, hist, resid, 0, nhist, PETSC_TRUE);
    #endif // PETSC_36
    ierr = TaoSolve(tao);
    ierr = TaoGetTerminationReason(tao, &reason);
    #ifdef PETSC_36
      TaoGetHistory(tao, NULL, NULL, NULL, NULL, &nhist);
    #else
      TaoGetHistory(tao, 0, 0, 0, &nhist);
    #endif // PETSC_36

    // print history and converged values to file
    char filename[30];
    sprintf(filename, "output_%.2f_%.2f.txt", x0_[0], x0_[1]);
    std::fstream out(filename, std::ios::out);
    if(!out.is_open()) {
      std::cerr << "Error: unable to create new file" << std::endl;
      exit(1);
    } // if
    for(int j = 0; j < nhist; ++ j) {
      //PetscPrintf(MPI_COMM_SELF, "History: %G\t%G\n", hist[j], resid[j]);
      PetscPrintf(MPI_COMM_SELF, "History: %g\t%g\n", hist[j], resid[j]);
      out << "History: " << j << "\t" << hist[j] << "\t" << resid[j] << std::endl;
    } // for

    TaoGetSolutionVector(tao, &x0);
    xn_.clear();
    for(PetscInt j = 0; j < num_params_; ++ j) {
      VecGetValues(x0, 1, &j, y);
      xn_.push_back(y[0]);
    } // for

    std::cout << "Converged vector: ";
    out << "Converged vector: ";
    for(real_vec_t::iterator i = xn_.begin(); i != xn_.end(); ++i) {
      std::cout << *i << " ";
      out << *i << " ";
    } // for
    std::cout << std::endl;
    out << std::endl;
    out.close();

    ierr = VecDestroy(&x0);
    ierr = VecDestroy(&G);
    ierr = TaoDestroy(&tao);
    TaoFinalize();

    return true;
  } // FitLMVMAlgo::run()


  void FitLMVMAlgo::print() {
    //std::cout << get_type_string() << " - Parameters: default." <<std::endl;
  } // FitLMVMAlgo::print()


  PetscErrorCode HipGISAXSFormFunctionGradient(TaoSolver tao, Vec X, PetscReal *f, Vec G, void *ptr) {
    PetscInt i, j, size;
    PetscReal fxp, fxm, dx = 0.05;
    PetscReal *x, *g;
    real_vec_t xvec, xpvec, xmvec;
    PetscErrorCode ierr;

    /* Get pointers to vector data */
    ierr = VecGetArray(X, &x);
    ierr = VecGetArray(G, &g);
    ierr = VecGetSize(X, &size);

    /* Compute Objective Funtion at X */
    for(i = 0; i < size; ++ i) xvec.push_back(x[i]);
    *f = EvaluateFunction(tao, xvec, ptr);

    /* Evaluate gradients */
    for(i = 0; i < size; ++ i) {
      xpvec = xmvec = xvec;
      xpvec[i] = xvec[i] + 0.5 * dx;
      xmvec[i] = xvec[i] - 0.5 * dx;
      fxp = EvaluateFunction(tao, xpvec, ptr);
      fxm = EvaluateFunction(tao, xmvec, ptr);
      g[i] = (fxp - fxm) / dx;
    } // for

    /* Restore vectors */
    ierr = VecRestoreArray(G, &g);

    return 0;
  } // HipGISAXSFormFunctionGradient()

} // namespace hig
