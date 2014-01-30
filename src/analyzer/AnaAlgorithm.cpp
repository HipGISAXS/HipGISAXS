/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaAlgorithm.cpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 29 Jan 2014 03:52:22 PM PST
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
#include <analyzer/AnaAlgorithm.hpp>

namespace hig{

  void AnaAlgorithm::set_initial_vec(AnaInput X0){
    int nvars = X0.size();
    X0_.clear();
    for(int i=0; i<nvars; i++){
      X0_.push_back( X0.get_var_init_val(i) );
    }

    if(pobj_fct_)
      pobj_fct_->set_dim(nvars);
  }


  /*
  Vec  AnaAlgorithm::tao_get_vec( float_vec_t x , int dim ){
    PetscFunctionBegin;
    Vec X;
    double y[1]= {0};
    int n = x.size() <= dim ?  x.size() : dim ;
    VecCreateSeq(PETSC_COMM_SELF, n , &X);
    for(int i=0 ; i<n ;i++){
      y[0] =(double) x[i];
      VecSetValues(X, 1, &i, y, INSERT_VALUES);
    }
    return X;
  }
  */

  AlgoType AnaAlgorithm::get_type(string_t algo_str){

    if(algo_str.compare("pounders")==0)
      {
	return fit_pounders;
      }
    else
      {
	return fit_brute_force;
      }
  }

  string_t AnaAlgorithm::get_type_string(){
    switch ( type_ )
      {
      case fit_pounders : return "TAO POUNDerS Fit";
      case fit_cg : return "TAO Conjugate Gradient Fit";
      default: return "Brute Force Fit";
      }
  }

  Analysis AnaAlgorithm::run(int argc,char **argv){
    Analysis ana_out;
    std::cout << "Running brute force analysis..." <<std::endl;
    return ana_out;
  }

  void AnaAlgorithm::print(){
    std::cout << get_type_string() << " - Parameters: default." <<std::endl;
  }

}
