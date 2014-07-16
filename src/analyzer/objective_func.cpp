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

	PetscErrorCode EvaluateFunction(TaoSolver tao, Vec X, Vec F, void *ptr) {
		PetscFunctionBegin;
		PetscErrorCode ierr;
		PetscReal *x, *ff;

		ierr = VecGetArray(X, &x); CHKERRQ(ierr);
        ierr = VecGetArray(F, &ff); CHKERRQ(ierr);

		int data_size = ((ObjectiveFunction*) ptr)->data_size();

		float_vec_t params;
		int num_params = ((ObjectiveFunction*) ptr)->num_fit_params();
		for(int i = 0; i < num_params; ++ i) params.push_back(x[i]);
		// run objective function
		float_vec_t temp = (*(ObjectiveFunction*) ptr)(params);
		for(int i = 0; i < data_size; ++ i) ff[i] = temp[i];

		ierr = VecRestoreArray(F, &ff); CHKERRQ(ierr);

		std::cout << "Eval X =\n" ;
		VecView(X, PETSC_VIEWER_STDOUT_WORLD);

		PetscFunctionReturn(0);
		return 0;
	} // EvaluateFunction()


	PetscReal EvaluateFunction(TaoSolver tao, float_vec_t X, void *ptr) {
		// Compute F(X)
		float_vec_t temp  =  (*(ObjectiveFunction*) ptr)(X);
        return  (PetscReal) temp[0];
	} // EvaluateFunction()
} // namespace hig
