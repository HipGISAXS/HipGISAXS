/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaAlgorithm.cpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 28 Jul 2014 01:23:24 PM PDT
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

		//Analysis ana_out;
		std::cout << "Running LMVM fitting ... [" << img_num << "]" << std::endl;

		if(!(*obj_func_).set_reference_data(img_num)) return false;

		static char help[] = "Running LMVM fitting...";

		/* PETSC & TAO INITIALIZATION  */
		PetscErrorCode ierr;
		int size, rank;		/* number of processes running */
		PetscInitialize(&argc, &argv, (char *) 0, help);
		TaoInitialize(&argc, &argv, (char *) 0, help);

		/* Variables to read from context   */
		Vec x0;

		double y[1] = { 0 };
		VecCreateSeq(PETSC_COMM_SELF, num_params_, &x0);
		for(int i = 0; i < num_params_; i++) {
		    y[0] = (double) x0_[i];
		    VecSetValues(x0, 1, &i, y, INSERT_VALUES);
		}

		/*******************************/
		/* Variables needed for TAO run   */
		PetscReal zero = 0.0;
		PetscReal hist[max_hist_], resid[max_hist_];
		PetscInt nhist = max_hist_;
		Vec G;			/* Gradient vector */
		TaoSolver tao;		/* TaoSolver solver context */
		TaoSolverTerminationReason reason;

		/* Allocate vectors for the solution, gradient, Hessian, etc. */
		ierr = VecCreateSeq(PETSC_COMM_SELF, num_obs_, &G);	//CHKERRQ(ierr);

		/****************************/
		/* The TAO code begins here */
		/****************************/

		/* Create TAO solver with desired solution method */
		ierr = TaoCreate(PETSC_COMM_SELF, &tao);	//CHKERRQ(ierr);
		ierr = TaoSetType(tao, "tao_lmvm");	//CHKERRQ(ierr);  // tao_nls

		/* Set routines for function, gradient, hessian evaluation */
		ierr = TaoSetObjectiveAndGradientRoutine(tao,
						      HipGISAXSFormFunctionGradient,
						      obj_func_);

		/* Check for TAO command line options */
		ierr = TaoSetFromOptions(tao);
		TaoSetMaximumIterations(tao, 100);

		//TaoSetTolerances(tao, tol_ , tol_, tol_, PETSC_DEFAULT, PETSC_DEFAULT);
		TaoSetTolerances(tao, tol_, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
		TaoDefaultMonitor(tao, PETSC_NULL);

		/* Set initial vector  */
		ierr = TaoSetInitialVector(tao, x0);	// CHKERRQ(ierr);

		/* Record History */
		ierr = TaoSetHistory (tao, hist, resid, 0, nhist, PETSC_TRUE);

		/* Run solver */
		ierr = TaoSolve(tao);	//CHKERRQ(ierr);

		std::cout << "---------------- DONE TAO SOLVE" << std::endl;

		/* Get termination information */
		ierr = TaoGetTerminationReason(tao, &reason);	// CHKERRQ(ierr);

		/*  get solution  */
		TaoGetHistory(tao, 0, 0, 0, &nhist);

		/* print history and converged values to file */
		char filename[30];
		sprintf(filename, "output%4.2f_%4.2f.txt", x0_[0], x0_[1]);
		std::fstream out(filename, std::ios::out);
		if (!out.is_open()){
		    std::cerr << "Error: unable to create new file" << std::endl;
	    	exit(1);
		}
		for (int j = 0; j < nhist; j++)
		{
		    PetscPrintf(MPI_COMM_SELF, "History : %G\t%G\n", hist[j], resid[j]);
	    	out << "History: " << j << "\t" << hist[j] << "\t" << resid[j] << std::endl;
		}

		TaoGetSolutionVector(tao, &x0);

		/*  Write to output vector    */
		xn_.clear();
		for (int j = 0; j < num_params_; j++) {
		    VecGetValues(x0, 1, &j, y);
		    xn_.push_back((float) y[0]);
		}

		/* Compose output analysis */
		std::cout << "Converged vector: ";
		out << "Converged vector: ";
		for (float_vec_t::iterator i = xn_.begin(); i != xn_.end(); ++i)
		{
		    std::cout << *i << " ";
	    out << *i << " ";
		}
		std::cout << std::endl;
		out << std::endl;
		out.close();

		/* Free PETSc data structures */
		ierr = VecDestroy(&x0);
		ierr = VecDestroy(&G);

		/* Destroy & finalize  */
		ierr = TaoDestroy(&tao);
		TaoFinalize();

		return true;
	}

	void FitLMVMAlgo::print() {
		//std::cout << get_type_string() << " - Parameters: default." <<std::endl;
	}

	PetscErrorCode HipGISAXSFormFunctionGradient(TaoSolver tao, Vec X,
							 PetscReal * f, Vec G, void *ptr) {
		PetscInt i, j, size;
		PetscReal fxp, fxm, dx = 0.05;
		PetscReal *x, *g;
		float_vec_t xvec, xpvec, xmvec;
		PetscErrorCode ierr;


		/* Get pointers to vector data */
		ierr = VecGetArray(X, &x);
		ierr = VecGetArray(G, &g);
		ierr = VecGetSize(X, &size);

		/* Compute Objective Funtion at X */
		for (i = 0; i < size; i++) 
		    xvec.push_back(x[i]);
		*f = EvaluateFunction(tao, xvec, ptr);

		/* Evaluate gradients */
		for (i = 0; i < size; i++) {
		    xpvec = xmvec = xvec;
		    xpvec[i] = xvec[i] + 0.5 * dx;
	    	xmvec[i] = xvec[i] - 0.5 * dx;
		    fxp = EvaluateFunction(tao, xpvec, ptr);
		    fxm = EvaluateFunction(tao, xmvec, ptr);
		    g[i] = (fxp - fxm) / dx;
		}

		/* Restore vectors */
		ierr = VecRestoreArray(G, &g);
		return 0;
	}
}
