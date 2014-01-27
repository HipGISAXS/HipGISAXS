/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaInput.cpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:20:20 AM PST
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

#include <time.h>
#include <iostream>
#include <analyzer/AnaInput.hpp>

namespace ana_hig{

  StructVar AnaInput::get_var(string_t uid){
    /* Linear access, may improve to cst time with hash map  */
    ana_struct_vec_it it;
    for(it=struct_vec_.begin(); it!= struct_vec_.end(); it++){
      if( uid.compare( it->get_uid() ) == 0 ){
	break;
      }
    }
    if( uid.compare( it->get_uid() ) == 0 )
      return *it;
    else
      return StructVar::dummy();
  }

  bool AnaInput::set_var_val(int i, float val){
    if ( i >= size() ) return false;
    struct_vec_[i].set_val(val);
    return true;
  }

  bool AnaInput::set_var_final_val(int i, float val){
    if ( i >= size() ) return false;
    struct_vec_[i].set_final_val(val);
    return true;
  }

  bool AnaInput::set_var_val(string_t uid, float val){
    /* Linear access, improve to cst time with hash map  */
    ana_struct_vec_it it;
    for(it=struct_vec_.begin(); it!= struct_vec_.end(); it++){
      if( uid.compare( it->get_uid() ) == 0 ){
	it->set_val(val);
        return true;
      }
    }
    return false;
  }

  void AnaInput::generate_random_vars(float min, float max, int dim){
    struct_vec_.clear();
    srand (time(NULL));
    for(int p=0; p<dim; p++){
      StructVar sv( min + rand()% abs(max-min) , "" );
      sv.set_fit_params(true,  min, max, -1, -1);
      struct_vec_.push_back(sv);
    }
  }

  void AnaInput::print(){

    std::cout << "Nbr Fit Variables : " << size() << "\n";
    ana_struct_vec_it it;
    int cnt =1;
    for(it=struct_vec_.begin(); it!= struct_vec_.end(); it++){
      std::cout << "Fit Var " << cnt << " : \n";
      it->print();
      cnt++;
    }
  }

}
