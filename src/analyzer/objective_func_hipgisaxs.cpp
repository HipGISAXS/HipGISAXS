/**
 *  Project:
 *
 *  File: objective_func.cpp
 *  Created: Feb 02, 2014
 *  Modified: Wed 16 Jul 2014 10:32:30 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <map>

#include <analyzer/objective_func_hipgisaxs.hpp>

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
		mask_data_.clear();
		mask_set_ = false;
		pdist_ = d;
		curr_dist_.clear();

	} // HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction()


	HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction(int narg, char** args, std::string config) :
			hipgisaxs_(narg, args) {
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
		mask_data_.clear();
		mask_set_ = false;
		pdist_ = NULL;
		curr_dist_.clear();

	} // HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction()


	HipGISAXSObjectiveFunction::~HipGISAXSObjectiveFunction() {
		if(ref_data_ != NULL) delete ref_data_;
	} // HipGISAXSObjectiveFunction::~HipGISAXSObjectiveFunction()


	bool HipGISAXSObjectiveFunction::set_distance_measure(DistanceMeasure* dist) {
		pdist_ = dist;
		return true;
	} // HipGISAXSObjectiveFunction::set_distance_measure()


	bool HipGISAXSObjectiveFunction::set_reference_data(int i) {
        if(i >= 0) {
            if(ref_data_ != NULL) delete ref_data_;
            ref_data_ = new ImageData(hipgisaxs_.reference_data_path(i));
            if(n_par_ != ref_data_->n_par() || n_ver_ != ref_data_->n_ver()) {
                std::cerr << "warning: reference and simulation data dimension sizes do not match [ "
                            << ref_data_->n_par() << " * " << ref_data_->n_ver() << " ] != [ "
                            << n_par_ << " * " << n_ver_ << " ]"
                            << std::endl;
				std::cerr << "warning: overriding simulation region with reference data region"
							<< std::endl;
                n_par_ = ref_data_->n_par();
                n_ver_ = ref_data_->n_ver();
                hipgisaxs_.override_qregion(n_par_, n_ver_, i);
            } // if
            if(!read_mask_data(hipgisaxs_.reference_data_mask(i))) return false;
            if(mask_data_.size() != n_par_ * n_ver_) {
                std::cerr << "error: mask and reference data dimension sizes do not match [ "
                            << n_par_ << " * " << n_ver_ << " ] != " << mask_data_.size()
                            << std::endl;
                return false;
            } // if
            mask_set_ = true;
        } // if
		if(!mask_set_) {
			mask_data_.clear();
			mask_data_.resize(n_par_ * n_ver_, 1);
		} // if
        return true;
	} // HipGISAXSObjectiveFunction::set_reference_data()


	bool HipGISAXSObjectiveFunction::read_mask_data(string_t filename) {
		if(filename.empty()) {
			mask_data_.resize(n_par_ * n_ver_, 1);
			return true;
		} // if
		std::ifstream maskf(filename);
		if(!maskf.is_open()) {
			std::cerr << "error: could not open mask data file " << filename << std::endl;
			return false;
		} // if
		while(true) {
			unsigned int val = 0;
			maskf >> val;
			if(maskf.eof()) break;
			mask_data_.push_back(val);
		} // while
		maskf.close();
		return true;
	} // HipGISAXSObjectiveFunction::read_mask_data()


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
		if(gisaxs_data == NULL) {
			std::cerr << "error: something went wrong in compute_gisaxs. gisaxs_data == NULL."
						<< std::endl;
			exit(-1);
		} // if

		// compute error/distance
		std::cout << "+++++ computing distance ..." << std::endl;
		float_t* ref_data = (*ref_data_).data();
		if(ref_data == NULL) std::cerr << "woops: ref_data is NULL" << std::endl;
		unsigned int* mask_data = NULL;
		//if(mask_set_) mask_data = &mask_data_[0];
		mask_data = &mask_data_[0];
		if(mask_data == NULL) std::cerr << "warning: mask_data is NULL" << std::endl;
		(*pdist_)(gisaxs_data, ref_data, mask_data, n_par_ * n_ver_, curr_dist_);

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
		std::cout << "++ Reference data set after simulation" << std::endl;

		return true;
	} // ObjectiveFunction::operator()()

} // namespace hig
