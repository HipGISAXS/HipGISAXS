/**
 *  Project:
 *
 *  File: objective_func.cpp
 *  Created: Feb 02, 2014
 *  Modified: Thu 27 Feb 2014 11:01:33 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <map>

#include <analyzer/objective_func.hpp>

namespace hig{

	HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction(int narg, char** args, DistanceMeasure* d) :
			hipgisaxs_(narg, args) {
		if(!hipgisaxs_.construct_input(args[1])) {
			std::cerr << "error: failed to construct HipGISAXS input containers" << std::endl;
			exit(1);
		} // if

		if(!hipgisaxs_.fit_init()) {
			std::cerr << "error: failed to initialize HipGISAXS for fitting" << std::endl;
			exit(1);
		} // if

		n_par_ = hipgisaxs_.nqy();
		n_ver_ = hipgisaxs_.nqz();

		ref_data_ = NULL;
		pdist_ = d;
		curr_dist_.clear();

	} // HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction()


	HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction(int narg, char** args, std::string config) :
			hipgisaxs_(narg, args) {
		//if(!hipgisaxs_.construct_input(args[1])) {
		if(!hipgisaxs_.construct_input(config.c_str())) {
			std::cerr << "error: failed to construct HipGISAXS input containers" << std::endl;
			exit(1);
		} // if

		if(!hipgisaxs_.fit_init()) {
			std::cerr << "error: failed to initialize HipGISAXS for fitting" << std::endl;
			exit(1);
		} // if

		n_par_ = hipgisaxs_.nqy();
		n_ver_ = hipgisaxs_.nqz();

		ref_data_ = NULL;
		pdist_ = NULL;
		curr_dist_.clear();

	} // HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction()


	HipGISAXSObjectiveFunction::~HipGISAXSObjectiveFunction() {
		delete ref_data_;
	} // HipGISAXSObjectiveFunction::~HipGISAXSObjectiveFunction()


	bool HipGISAXSObjectiveFunction::set_distance_measure(DistanceMeasure* dist) {
		pdist_ = dist;
		return true;
	} // HipGISAXSObjectiveFunction::set_distance_measure()


	bool HipGISAXSObjectiveFunction::set_reference_data(int i) {
		if(i > 0) {
			if(ref_data_ != NULL) delete ref_data_;
			ref_data_ = new ImageData(hipgisaxs_.reference_data_path(i));
		} // if
		return true;
	} // HipGISAXSObjectiveFunction::set_reference_data()


	float_vec_t HipGISAXSObjectiveFunction::operator()(const float_vec_t& x) {
		float_t *gisaxs_data = NULL;
		// construct param_vals
		std::vector <std::string> params = hipgisaxs_.fit_param_keys();
		// TODO check if param values are within range ...
		std::map <std::string, float_t> param_vals;
		for(int i = 0; i < x.size(); ++ i) param_vals[params[i]] = x[i];

		for(std::map<std::string, float_t>::iterator i = param_vals.begin(); i != param_vals.end(); ++ i)
			std::cout << (*i).first << ": " << (*i).second << "  ";
		std::cout << std::endl;

		// update and compute gisaxs
		hipgisaxs_.update_params(param_vals);
		hipgisaxs_.compute_gisaxs(gisaxs_data);

		// compute error/distance
		float_t* ref_data = (*ref_data_).data();
		(*pdist_)(gisaxs_data, ref_data, n_par_ * n_ver_, curr_dist_);

		// write to output file
		std::string prefix(HiGInput::instance().param_pathprefix()+"/"+HiGInput::instance().runname());
		std::ofstream out(prefix + "/convergance.dat", std::ios::app);
		for(float_vec_t::const_iterator i = curr_dist_.begin(); i != curr_dist_.end(); ++ i)
			out << (*i) << " ";
		out << std::endl;
		out.close();

		return curr_dist_;
	} // ObjectiveFunction::operator()()


	bool HipGISAXSObjectiveFunction::simulate_and_set_ref(const float_vec_t& x) {
		float_t *gisaxs_data = NULL;
		if(x.size() > 0) {
			// construct param_vals
			std::vector <std::string> params = hipgisaxs_.fit_param_keys();
			std::map <std::string, float_t> param_vals;
			for(int i = 0; i < x.size(); ++ i) param_vals[params[i]] = x[i];

			for(std::map<std::string, float_t>::iterator i = param_vals.begin();
					i != param_vals.end(); ++ i)
				std::cout << (*i).first << ": " << (*i).second << "  ";
			std::cout << std::endl;

		// update and compute gisaxs
			hipgisaxs_.update_params(param_vals);
		} // if

		hipgisaxs_.compute_gisaxs(gisaxs_data);
		if(ref_data_ == NULL) ref_data_ = new ImageData(n_par_, n_ver_);
		(*ref_data_).set_data(gisaxs_data);

		return true;
	} // ObjectiveFunction::operator()()


	// from Slim's original code
	PetscErrorCode EvaluateFunction(TaoSolver tao, Vec X, Vec F, void *ptr) {
		// Compute F(X)
		PetscFunctionBegin;
		VecView(X, PETSC_VIEWER_STDOUT_WORLD);

		PetscErrorCode ierr;
		PetscReal *x, *f;

		ierr = VecGetArray(X,&x);
		ierr = VecGetArray(F,&f);

		// either ff or f can be eliminated ...
		int data_size = ((ObjectiveFunction*) ptr)->data_size();
		PetscReal* ff = new PetscReal[data_size];
		float_vec_t params;
		int num_params = ((ObjectiveFunction*) ptr)->num_fit_params();
		for(int i = 0; i < num_params; ++ i) params.push_back(x[i]);
		float_vec_t temp = (*(ObjectiveFunction*) ptr)(params);
		for(int i = 0; i < data_size; ++ i) { ff[i] = temp[i]; f[i] = temp[i]; }

		ierr = VecRestoreArray(X, &x); CHKERRQ(ierr);
		ierr = VecRestoreArray(F, &ff); CHKERRQ(ierr);

		std::cout << "Eval X =\n" ;
		VecView(X, PETSC_VIEWER_STDOUT_WORLD);
		//std::cout << "F = \n" ;
		//VecView(F, PETSC_VIEWER_STDOUT_WORLD);

		PetscFunctionReturn(0);
		return 0;
	} // EvaluateFunction()


} // namespace hig
