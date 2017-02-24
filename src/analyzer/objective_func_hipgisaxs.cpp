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
#include <boost/math/special_functions/fpclassify.hpp>

#include <analyzer/objective_func_hipgisaxs.hpp>
#include <file/edf_reader.hpp>

namespace hig {

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

    n_par_ = hipgisaxs_.ncol();
    n_ver_ = hipgisaxs_.nrow();

    ref_data_ = NULL;
    mean_data_ = NULL;
    mask_data_.clear();
    mask_set_ = false;
    pdist_ = d;

    reg_alpha_ = 0.0;

    set_mean_data();

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

    n_par_ = hipgisaxs_.ncol();
    n_ver_ = hipgisaxs_.nrow();

    ref_data_ = NULL;
    mean_data_ = NULL;
    mask_data_.clear();
    mask_set_ = false;
    pdist_ = NULL;

    reg_alpha_ = 0.0;

    set_mean_data();

  } // HipGISAXSObjectiveFunction::HipGISAXSObjectiveFunction()


  HipGISAXSObjectiveFunction::~HipGISAXSObjectiveFunction() {
    if(ref_data_ != NULL) delete ref_data_;
    if(mean_data_ != NULL) delete[] mean_data_;
  } // HipGISAXSObjectiveFunction::~HipGISAXSObjectiveFunction()


  bool HipGISAXSObjectiveFunction::set_distance_measure(DistanceMeasure* dist) {
    pdist_ = dist;
    return true;
  } // HipGISAXSObjectiveFunction::set_distance_measure()


  ReferenceFileType get_reference_file_type(const std::string& fname) {
    //std::cout << "input file = " << fname << std::endl;
    size_t len = fname.size();
    char temp_ext[16];
    int extlen = 0;
    while(len > 0) {
      char temp = fname[len - 1];
      if(temp == '.') break;
      temp_ext[extlen ++] = temp;
      if(extlen == 15) break;  // failsafe
      -- len;
    } // while()
    temp_ext[extlen] = '\0';
    std::string temp_ext2(temp_ext);
    std::string ext(temp_ext2.rbegin(), temp_ext2.rend());
    if(ext.compare(std::string("edf")) == 0) return reference_file_edf;
    if(ext.compare(std::string("EDF")) == 0) return reference_file_edf;
    if(ext.compare(std::string("txt")) == 0) return reference_file_ascii;
    if(ext.compare(std::string("TXT")) == 0) return reference_file_ascii;
    if(ext.compare(std::string("dat")) == 0) return reference_file_ascii;
    if(ext.compare(std::string("DAT")) == 0) return reference_file_ascii;
    if(ext.compare(std::string("out")) == 0) return reference_file_ascii;
    if(ext.compare(std::string("OUT")) == 0) return reference_file_ascii;
    return reference_file_error;
  } // get_reference_file_type()


  bool HipGISAXSObjectiveFunction::set_reference_data(int i) {
    if(i >= 0) {
      if(ref_data_ != NULL) delete ref_data_;
      std::string ref_filename = hipgisaxs_.reference_data_path(i);
      ReferenceFileType ref_type = get_reference_file_type(ref_filename);
      real_t* temp_data = NULL;
      unsigned int temp_n_par = 0, temp_n_ver = 0;
      EDFReader* edfreader = NULL;
      switch(ref_type) {
        case reference_file_ascii:
          ref_data_ = new ImageData(ref_filename);
          break;

        case reference_file_edf:
          ref_data_ = new ImageData();
          edfreader = new EDFReader(ref_filename.c_str());
          edfreader->get_data(temp_data, temp_n_par, temp_n_ver);
          ref_data_->set_data(temp_data, temp_n_par, temp_n_ver);
          delete edfreader;
          break;

        case reference_file_null:
        case reference_file_error:
        default:
          std::cerr << "error: error in reference data file format (i.e. extension!!!)"
                << std::endl;
          return false;
      } // switch
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

      std::string mask_filename = hipgisaxs_.reference_data_mask(i);
      ref_type = get_reference_file_type(mask_filename);
      if(ref_type == reference_file_edf) {
        if(!read_edf_mask_data(hipgisaxs_.reference_data_mask(i))) return false;
      } else {
        if(!read_mask_data(hipgisaxs_.reference_data_mask(i))) return false;
      } // if-else
      if(mask_data_.size() != n_par_ * n_ver_) {
        std::cerr << "error: mask and reference data dimension sizes do not match [ "
                  << n_par_ << " * " << n_ver_ << " ] != " << mask_data_.size()
                  << std::endl;
        return false;
      } // if
    } // if
    if(!mask_set_) {
      mask_data_.clear();
      mask_data_.resize(n_par_ * n_ver_, 1);
    } // if
    return true;
  } // HipGISAXSObjectiveFunction::set_reference_data()


  bool HipGISAXSObjectiveFunction::read_mask_data(string_t filename) {
    mask_data_.clear();
    if(filename.empty()) {
      mask_data_.resize(n_par_ * n_ver_, 1);
      return true;
    } // if
    //std::cout << "+++++++++++++ reading mask: " << filename << std::endl;
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
    mask_set_ = true;
    return true;
  } // HipGISAXSObjectiveFunction::read_mask_data()


  bool HipGISAXSObjectiveFunction::read_edf_mask_data(string_t filename) {
    mask_data_.clear();
    if(filename.empty()) {
      mask_data_.resize(n_par_ * n_ver_, 1);
      return true;
    } // if
    //std::cout << "-- Reading mask data from " << filename << "..." << std::endl;
    EDFReader* edfreader = new EDFReader(filename.c_str());
    real_t* temp_data = NULL;
    unsigned int temp_n_par = 0, temp_n_ver = 0;
    edfreader->get_data(temp_data, temp_n_par, temp_n_ver);
    if(temp_data == NULL) {
      std::cerr << "error: failed to get edf mask data" << std::endl;
      return false;
    } // if
    for(unsigned int i = 0; i < temp_n_par * temp_n_ver; ++ i) mask_data_.push_back(temp_data[i]);
    mask_set_ = true;
    delete edfreader;
    return true;
  } // HipGISAXSObjectiveFunction::read_edf_mask_data()


  // computes and sets the simulated data with mean parameter vector values
  bool HipGISAXSObjectiveFunction::set_mean_data(void) {
    mean_data_ = NULL;

    // construct mean param_vals
    std::vector <std::string> params = hipgisaxs_.fit_param_keys();
    std::map <std::string, real_t> param_vals;
    for(int i = 0; i < params.size(); ++ i) {
      param_vals[params[i]] = hipgisaxs_.param_space_mean(params[i]);
    } // for

    std::cout << "-- Mean/estimated parameter values: ";
    for(std::map<std::string, real_t>::iterator i = param_vals.begin(); i != param_vals.end(); ++ i)
      std::cout << (*i).first << ": " << (*i).second << "  ";
    std::cout << std::endl;

    hipgisaxs_.update_params(param_vals);
    hipgisaxs_.compute_gisaxs(mean_data_);

    return true;
  } // HipGISAXSObjectiveFunction::set_mean_data()


  real_vec_t HipGISAXSObjectiveFunction::operator()(const real_vec_t& x) {

    //bool use_mean = true;    // temporary
    bool use_mean = false;    // temporary

    real_t *gisaxs_data = NULL;
    // construct param_vals
    std::vector <std::string> params = hipgisaxs_.fit_param_keys();
    std::map <std::string, real_t> param_vals;
    for(int i = 0; i < x.size(); ++ i) param_vals[params[i]] = x[i];

    std::cout << "-- Parameter values input to objective function: ";
    for(std::map<std::string, real_t>::iterator i = param_vals.begin(); i != param_vals.end(); ++ i)
      std::cout << (*i).first << ": " << (*i).second << "  ";
    std::cout << std::endl;

    real_vec_t curr_dist;

    /**
     * update and compute gisaxs
     */
    hipgisaxs_.update_params(param_vals);
    hipgisaxs_.compute_gisaxs(gisaxs_data);

    // only the master process does the following
    if(hipgisaxs_.is_master()) {
      if(gisaxs_data == NULL) {
        std::cerr << "error: something went wrong in compute_gisaxs. gisaxs_data == NULL."
                  << std::endl;
        exit(-1);
      } // if

      #ifdef MEMDEBUG
        std::cerr << "** MEMDEBUG enabled" << std::endl;
        if(!hipgisaxs_.check_finite((*ref_data_).data(), n_par_ * n_ver_))
          std::cerr << "** ARRAY CHECK ** ref_data_ failed check" << std::endl;
        if(!hipgisaxs_.check_finite(gisaxs_data, n_par_ * n_ver_))
          std::cerr << "** ARRAY CHECK ** gisaxs_data failed check" << std::endl;
      #endif // MEMDEBUG

      /**
       * compute error/distance
       */
      std::cout << "-- Computing distance..." << std::endl;
      real_t* ref_data = (*ref_data_).data();
      if(ref_data == NULL) std::cerr << "error: ref_data is NULL" << std::endl;
      unsigned int* mask_data = &(mask_data_[0]);

      // distance function
      if(use_mean) (*pdist_)(ref_data, gisaxs_data, mask_data, n_par_ * n_ver_, curr_dist, mean_data_);
      else (*pdist_)(ref_data, gisaxs_data, mask_data, n_par_ * n_ver_, curr_dist);

      /** regularization **/
      // compute regularization = (alpha / 2) * || x - x_mean || ^ 2
      double pmean = 0.0;
      for(std::map<std::string, real_t>::const_iterator i = param_vals.begin();
          i != param_vals.end(); ++ i) {
        double x_mean = hipgisaxs_.param_space_mean((*i).first);
        double ptemp = ((*i).second - x_mean);
        ptemp *= ptemp;
        pmean += ptemp;
      } // for
      // for better alpha selection, plot pmean vs. curr_dist
      real_t reg = (reg_alpha_ / 2) * pmean;
      std::cout << "## reg_const: " << reg_alpha_ << ", param_norm: " << pmean
                << ", reg_value: " << reg << ", dist[0]: " << curr_dist[0] << std::endl;
      real_t reg_dist = reg / curr_dist.size();   // distribute it across the distance vector
      for(auto i = 0; i < curr_dist.size(); ++ i) {
        if((boost::math::isfinite)(curr_dist[i])) curr_dist[i] += reg_dist;
        else curr_dist[i] = 1e8;    // some large number
      } // for

      // write to output file
      // do something better ...
      // as here all procs for different particles will have same ranks within their hipgisaxs object
      int myrank = hipgisaxs_.rank();
      std::stringstream cfilename;
      cfilename << "distance." << myrank << ".dat";
      std::string prefix(hipgisaxs_.path() + "/" + hipgisaxs_.runname());
      std::ofstream out(prefix + "/" + cfilename.str(), std::ios::app);
      out.precision(10);
      for(real_vec_t::const_iterator i = curr_dist.begin(); i != curr_dist.end(); ++ i) out << *i << " ";
      out << std::endl;
      out.close();
    } // if

    return curr_dist;
  } // ObjectiveFunction::operator()()


  // this is used only in the test mode
  bool HipGISAXSObjectiveFunction::simulate_and_set_ref(const real_vec_t& x) {
    real_t *gisaxs_data = NULL;
    if(x.size() > 0) {
      // construct param_vals
      std::vector <std::string> params = hipgisaxs_.fit_param_keys();
      std::map <std::string, real_t> param_vals;
      for(int i = 0; i < x.size(); ++ i) param_vals[params[i]] = x[i];

      std::cout << "-- Reference parameter values: ";
      for(std::map<std::string, real_t>::iterator i = param_vals.begin(); i != param_vals.end(); ++ i)
        std::cout << (*i).first << " = " << (*i).second << "\t";
      std::cout << std::endl;

      // update and compute gisaxs
      hipgisaxs_.update_params(param_vals);
    } // if

    hipgisaxs_.compute_gisaxs(gisaxs_data);
    #ifdef USE_MPI
      woo::MultiNode* temp_multi = hipgisaxs_.multi_node_comm();
      woo::comm_t root = (*temp_multi).universe_key();
      if((*temp_multi).size(root) > 1) {
        // the root processor has the computed data. it needs to send it to all other procs.
        if(!(*temp_multi).is_master(root)) {
          if(gisaxs_data != NULL) delete[] gisaxs_data;
          gisaxs_data = new (std::nothrow) real_t[n_par_ * n_ver_];
        } // if
        if(!(*temp_multi).broadcast(root, gisaxs_data, n_par_ * n_ver_)) {
          std::cerr << "error: failed to broadcast gisaxs_data to all" << std::endl;
          return false;
        } // if
      } // if
    #endif // USE_MPI
    if(ref_data_ == NULL) ref_data_ = new ImageData(n_par_, n_ver_);
    (*ref_data_).set_data(gisaxs_data);
    std::cout << "** Reference data set after simulation" << std::endl;

    return true;
  } // ObjectiveFunction::operator()()

} // namespace hig
