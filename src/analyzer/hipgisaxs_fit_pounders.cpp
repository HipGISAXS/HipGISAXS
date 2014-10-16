/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaAlgorithm.cpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 08 Oct 2014 12:17:42 PM PDT
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
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

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
 */

namespace hig{

  bool FitPOUNDERSAlgo::run(int argc,char **argv, int img_num){

    //Analysis ana_out;
    std::cout << "Running POUNDERS fitting..." <<std::endl;

  if(!(*obj_func_).set_reference_data(img_num)) return false;

    static char help[]= "Running POUNDERS fitting...";

    /* PETSC & TAO INITIALIZATION  */
    PetscErrorCode  ierr;
    int        size,rank;     /* number of processes running */
    PetscInitialize(&argc,&argv,(char *)0,help);
    TaoInitialize(&argc,&argv,(char*)0,help);
//    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);// CHKERRQ(ierr);
//    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);// CHKERRQ(ierr);

    /* Variables to read from context   */
    Vec x0;

    double y[1]= {0};
    VecCreateSeq(PETSC_COMM_SELF, num_params_ , &x0);
    for(int i=0 ; i<num_params_ ;i++){
      y[0] =(double) x0_[i];
      VecSetValues(x0, 1, &i, y, INSERT_VALUES);
    }

    PetscReal tr_rad =10;

    /* Variables needed for TAO run   */
    Vec f   ;//, zeros;
    PetscReal zero = 0.0;
    PetscReal  hist[max_hist_],resid[max_hist_];
    PetscInt   nhist;
    //Mat        H;        /* Hessian matrix */
    //Mat        J;        /* Jacobian matrix */
    TaoSolver  tao;      /* TaoSolver solver context */
    TaoSolverTerminationReason reason;

    /* Allocate vectors for the solution, gradient, Hessian, etc. */
    ierr = VecCreateSeq(PETSC_COMM_SELF, num_obs_ ,&f);

    /* The TAO code begins here */

    /* Create TAO solver with desired solution method */
    ierr = TaoCreate(PETSC_COMM_SELF,&tao);
    ierr = TaoSetType(tao, "tao_pounders");
    /* Set routines for function, gradient, hessian evaluation */
    ierr = TaoSetSeparableObjectiveRoutine(tao, f, EvaluateFunction, obj_func_);

    /* Check for TAO command line options */
    ierr = TaoSetFromOptions(tao);
    TaoSetMaximumIterations( tao, max_iter_);
    //TaoSetInitialTrustRegionRadius(tao,  tr_rad);
    TaoSetHistory(tao, hist, resid, 0, max_hist_, PETSC_TRUE);
    TaoSetTolerances(tao, tol_ , PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

    /* Set initial vector  */
    ierr = TaoSetInitialVector(tao, x0);

    /* Run solver */
    ierr = TaoSolve(tao);

    /* Get termination information */
    ierr = TaoGetTerminationReason(tao, &reason);

    /*  get solution  */
    TaoGetHistory(tao, 0, 0, 0, &nhist);
    for(int it = 0; it < nhist; ++ it) PetscPrintf(PETSC_COMM_WORLD, "%G\t%G\n", hist[it], resid[it]);

    PetscInt iterate;
    PetscReal f_cv;
    PetscReal gnorm;
    PetscReal cnorm;
    PetscReal xdiff;
    PetscReal *x_cv;
    TaoGetSolutionStatus(tao, &iterate, &f_cv, &gnorm, &cnorm, &xdiff, &reason);
    TaoGetSolutionVector(tao, &x0);
    VecGetArray(x0, &x_cv);

    /*  Write to output vector    */
    xn_.clear();
    for(int j = 0; j < num_params_; ++ j){
      VecGetValues(x0, 1 , &j , y);
      xn_.push_back((float) y[0]);
    }

    std::cout << "Converged vector: ";
  for(float_vec_t::iterator i = xn_.begin(); i != xn_.end(); ++ i) std::cout << *i << " ";
  std::cout << std::endl;

    ierr = TaoDestroy(&tao);
    ierr = VecDestroy(&x0);
    TaoFinalize();

    return true;
  }


  void FitPOUNDERSAlgo::print(){
    //std::cout << get_type_string() << " - Parameters: default." <<std::endl;
  }

}
