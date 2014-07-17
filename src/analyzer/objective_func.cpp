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
		int data_size = ((ObjectiveFunction*) ptr)->data_size();
		int num_params = ((ObjectiveFunction*) ptr)->num_fit_params();

		ierr = VecGetArray(X, &x); CHKERRQ(ierr);
        ierr = VecGetArray(F, &ff); CHKERRQ(ierr);

		float_vec_t params;
		for(int i = 0; i < num_params; ++ i) params.push_back(x[i]);
		// run objective function
		float_vec_t temp = (*(ObjectiveFunction*) ptr)(params);
		float_t err = 0.0;
		float_t* ref_data = (*(ObjectiveFunction*) ptr).get_reference_data();
		unsigned int* mask_data = (*(ObjectiveFunction*) ptr).get_mask_data();
		for(int i = 0; i < data_size; ++ i) {
			ff[i] = mask_data[i] * temp[i];
			//err += mask_data[i] * (temp[i] * temp[i] / ref_data[i]);
			err += mask_data[i] * (temp[i] * temp[i]);
		} // for

		ierr = VecRestoreArray(F, &ff); CHKERRQ(ierr);

		std::cout << "Distance = " << err << std::endl;
		std::cout << "Eval X =\n" ;
		VecView(X, PETSC_VIEWER_STDOUT_WORLD);

		PetscFunctionReturn(0);
		return 0;
	} // EvaluateFunction()


	PetscReal EvaluateFunction(TaoSolver tao, float_vec_t X, void *ptr) {
		std::cout << "evaluate function ..." << std::endl;
		// Compute F(X)
		float_vec_t temp = (*(ObjectiveFunction*) ptr)(X);
		std::cout << "*********************** " << temp[0] << std::endl;
        return  (PetscReal) temp[0];
	} // EvaluateFunction()
} // namespace hig
