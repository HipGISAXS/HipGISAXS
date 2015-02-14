/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: SimModel.cpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 29 Jan 2014 03:53:07 PM PST
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

#include <analyzer/SimModel.hpp>

namespace hig{

  int SimModel::compute(){

    return 0;
  }

  int SimModel::update_vars(real_vec_t X){
    return 0;
  }

  void SimModel::generate_rand_vars(float min, float max, int dim){
    fit_vars_.clear();
    //srand(time(NULL));
    for(int p=0; p<dim; p++){
      fit_vars_.push_back(min+ rand()% abs(max-min));
    }
  }

  void SimModel::print(){
    //data_sim_.print();
    /*  print vars   */
    std::cout << "Fit vars = \n";
    for(int p =0; p< fit_vars_.size();p++){
      std::cout << fit_vars_[p] << "      " ;
    }
    std::cout <<std::endl;
  }

  void SimModel::save(string_t filename){
    data_sim_.save(filename);
  }

}

