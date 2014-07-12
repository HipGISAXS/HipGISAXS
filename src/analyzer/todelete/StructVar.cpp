/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: StructVar.cpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 29 Jan 2014 03:53:18 PM PST
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
#include <analyzer/StructVar.hpp>

namespace hig{

  StructVar::StructVar(float init_val, string_t uid, var_attr_t attr , var_fit_params_t params){
    init_val_=init_val;
    uid_ = uid;
    set_attr(attr);
    set_fit_params(params);
  }

  StructVar::StructVar(float init_val, string_t uid, var_fit_params_t params){
    init_val_=init_val;
    uid_ = uid;
    set_fit_params(params);
  }

  void StructVar::print(){
    std::cout << "[\n";
    std::cout << "Variable UID: "  << uid_  << std::endl;
    std::cout << "Current Value: "  << val_  << std::endl;
    std::cout << "Initial Value: "  << init_val_  << std::endl;
    std::cout << "Final Value: "  << fin_val_  << std::endl;
    std::cout << "Value Min: "  << fit_params_.min_  << std::endl;
    std::cout << "Value Max: "  << fit_params_.max_  << std::endl;
    std::cout << "Nbr Values: "  << fit_params_.n_vals_  << std::endl;
    std::cout << "]\n";
  }

}
