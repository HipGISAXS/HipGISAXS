/**
 *  Project:
 *
 *  File:
 *  Created: Feb 03, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <cmath>

#include <analyzer/objective_func_poly_one.hpp>

namespace hig {

  PolyOneObjectiveFunction::PolyOneObjectiveFunction(int narg, char** args, DistanceMeasure* d) {

    pdist_ = d;

    ref_data_ = new (std::nothrow) ImageData;
    set_reference_data(args[1]);

    q_min_ = atof(args[2]);
    q_max_ = atof(args[3]);
    q_step_ = atof(args[4]);
    x1_min_ = atof(args[5]);
    x1_max_ = atof(args[6]);
    x2_min_ = atof(args[7]);
    x2_max_ = atof(args[8]);
    x1_init_ = atof(args[9]);
    x2_init_ = atof(args[10]);

    std::cout << " Q: [ " << q_min_ << " " << q_max_ << " " << q_step_ << " ]" << std::endl;
    std::cout << "x1: [ " << x1_min_ << " " << x1_max_ << " ]" << std::endl;
    std::cout << "x2: [ " << x2_min_ << " " << x2_max_ << " ]" << std::endl;

  } // PolyOneObjectiveFunction::PolyOneObjectiveFunction()


  PolyOneObjectiveFunction::~PolyOneObjectiveFunction() { delete ref_data_; }


  real_vec_t PolyOneObjectiveFunction::operator()(const real_vec_t& x) {
    x1_ = x[0];
    x2_ = x[1];

    std::cout << " ** params: " << x1_ << " " << x2_ << std::endl;

    real_vec_t data;
    
    std::cout << " ** computed:\t";
    for(real_t i = q_min_; i < q_max_; i += q_step_) {
      real_t val = func(i, x1_, x2_);
      data.push_back(val);
      std::cout << val << " ";
    } // for
    std::cout << std::endl;
  
    //curr_dist_.clear();
    real_vec_t curr_dist;
    real_t* ref_data = ref_data_->data();
    real_t* d = &data[0];
    (*pdist_)(ref_data, d, data.size(), curr_dist);

    std::cout << " ** reference:\t";
    for(int i = 0; i < ref_data_->size(); ++ i) std::cout << ref_data[i] << " ";
    std::cout << std::endl;
    std::cout << " ** residual:\t";
    for(int i = 0; i < curr_dist.size(); ++ i) std::cout << curr_dist[i] << " ";
    std::cout << std::endl;

    return curr_dist;
  } // PolyOneObjectiveFunction::operator()


  bool PolyOneObjectiveFunction::set_reference_data(char* filename) {
         std::ifstream fin(filename);
         if(!fin.is_open()) {
             std::cerr << "error: failed to open reference data file " << filename << std::endl;
             return false;
         } // if
     real_vec_t data;
         while(true) {
             real_t temp = -1.0;
             fin >> temp;
             if(fin.eof() || !fin.good()) break;
             data.push_back(temp);
         } // while
         fin.close();
     ref_data_->set_data(data);
         return true;
  } // PolyOneObjectiveFunction::set_reference_data()


  std::vector <std::string> PolyOneObjectiveFunction::fit_param_keys() const {
    std::vector <std::string> keys;
    keys.push_back("x1"); keys.push_back("x2");
    return keys;
  } // PolyOneObjectiveFunction::fit_param_keys()


  std::vector <float_pair_t> PolyOneObjectiveFunction::fit_param_limits() const {
    std::vector <float_pair_t> limits;
    limits.push_back(float_pair_t(x1_min_, x1_max_));
    limits.push_back(float_pair_t(x2_min_, x2_max_));
    return limits;
  } // PolyOneObjectiveFunction::fit_param_limits()


  real_vec_t PolyOneObjectiveFunction::fit_param_init_values() const {
    std::vector <real_t> init_vals;
    init_vals.push_back(x1_init_);
    init_vals.push_back(x2_init_);
    return init_vals;
  } // PolyOneObjectiveFunction::fit_param_init_values()


  real_t PolyOneObjectiveFunction::func(real_t q, real_t x1, real_t x2) {
         real_t q3 = q * q * q;
         real_t q4 = q3 * q;
         float val = fabs(x1 * q4 + x2 * q3);
         return val;
     } // PolyOneFitness::func()

} // namespace hig
mespace hig
e hig
