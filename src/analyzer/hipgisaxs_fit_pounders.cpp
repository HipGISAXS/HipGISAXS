/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaAlgorithm.cpp
 *  Created: Dec 26, 2013
 *  Modified: Fri 11 Jul 2014 10:20:12 AM PDT
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

    static  char help[]= "Running POUNDERS fitting...";

    /* PETSC & TAO INITIALIZATION  */
    PetscErrorCode  ierr;
    int        size,rank;     /* number of processes running */
    PetscInitialize(&argc,&argv,(char *)0,help);
    TaoInitialize(&argc,&argv,(char*)0,help);
//    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);// CHKERRQ(ierr);
//    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);// CHKERRQ(ierr);
//    if (size >1) {
//      if (rank == 0) {
//	PetscPrintf(PETSC_COMM_SELF,"This example is intended for single processor use!\n");
	//SETERRQ(PETSC_COMM_SELF,1,"Incorrect number of processors");
//      }
//    }

    /* Variables to read from context   */
    Vec x0;

    double y[1]= {0};
    VecCreateSeq(PETSC_COMM_SELF, num_params_ , &x0);
    for(int i=0 ; i<num_params_ ;i++){
      y[0] =(double) x0_[i];
      VecSetValues(x0, 1, &i, y, INSERT_VALUES);
    }
    //VecView(x0, PETSC_VIEWER_STDOUT_WORLD);
    /*******************************/

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

    /*      int argc =2;

	    char opt[50]="-tao_pounders_delta";
	    char val[10]="5";
	    char* opt_ch= new char[strlen(opt)+1] ;//= opt.c_str();
	    char* val_ch= new char[strlen(val)+1] ;//= val.c_str();
	    strcpy(opt_ch,  opt);
	    strcpy(val_ch,  val);


	    char **argv= new char*[2];
	    argv[0]= opt_ch;
	    argv[1]= val_ch;

	    std::cout << argv[0] << " = " << argv[1] << std::endl;

	    // argv[0]=meth;
	    //argv[1]=val;
	    */

    /* Allocate vectors for the solution, gradient, Hessian, etc. */
    ierr = VecCreateSeq(PETSC_COMM_SELF, num_obs_ ,&f);// CHKERRQ(ierr);

    /* The TAO code begins here */

    /* Create TAO solver with desired solution method */
    ierr = TaoCreate(PETSC_COMM_SELF,&tao);// CHKERRQ(ierr);
    ierr = TaoSetType(tao, "tao_pounders");// CHKERRQ(ierr);  // tao_nls
    /* Set routines for function, gradient, hessian evaluation */
    ierr = TaoSetSeparableObjectiveRoutine(tao, f, EvaluateFunction, obj_func_);// CHKERRQ(ierr);

    /* Check for TAO command line options */
    ierr = TaoSetFromOptions(tao); // CHKERRQ(ierr);
    TaoSetMaximumIterations( tao, max_iter_);
    //      TaoSetInitialTrustRegionRadius(tao,  tr_rad);
    TaoSetHistory(tao, hist, resid, 0, max_hist_, PETSC_TRUE);
    //TaoSetTolerances(tao, tol_ , tol_, tol_, PETSC_DEFAULT, PETSC_DEFAULT);
    TaoSetTolerances(tao, tol_ , PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    //TaoDefaultMonitor(tao, PETSC_NULL);

    /* Set initial vector  */
    ierr = TaoSetInitialVector(tao, x0);// CHKERRQ(ierr);

    /* Run solver */
    ierr = TaoSolve(tao); //CHKERRQ(ierr);
    //f=      pobj_fct_->tao_compute(x0);
    //VecView(f, PETSC_VIEWER_STDOUT_WORLD);

    /* Get termination information */
    ierr = TaoGetTerminationReason(tao,&reason); // CHKERRQ(ierr);

    /*  get solution  */
    TaoGetHistory(tao,0,0,0,&nhist);
    for (int it=0;it<nhist;it++) {
      PetscPrintf(PETSC_COMM_WORLD,"%G\t%G\n",hist[it],resid[it]);
    }

    PetscInt iterate;
    PetscReal f_cv;
    PetscReal gnorm;
    PetscReal cnorm;
    PetscReal xdiff;
    PetscReal *x_cv;
    TaoGetSolutionStatus(tao, &iterate, &f_cv, &gnorm, &cnorm, &xdiff, &reason);
    TaoGetSolutionVector(tao, &x0);
    VecGetArray(x0,&x_cv);

    /*  Write to output vector    */
    xn_.clear();
    for(int j=0; j<num_params_; j++){
      VecGetValues(x0, 1 , &j , y);
      xn_.push_back( (float) y[0]  );
    }

    std::cout << "Converged vector: ";
  //  VecView(x0, PETSC_VIEWER_STDOUT_WORLD);
	for(float_vec_t::iterator i = xn_.begin(); i != xn_.end(); ++ i)
		std::cout << *i << " ";
	std::cout << std::endl;

    /* Compose output analysis */

    /* Destroy & finalize  */
    ierr = TaoDestroy(&tao); //CHKERRQ(ierr);
    /* Free PETSc data structures */
    ierr = VecDestroy(&x0); //CHKERRQ(ierr);
    TaoFinalize();

    return true;
  }


  void FitPOUNDERSAlgo::print(){
    //std::cout << get_type_string() << " - Parameters: default." <<std::endl;
  }

}
