/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaInput.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:21:47 AM PST
 *  Description: Abstract class for gisaxs data
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

#ifndef _ANAINPUT_HPP_
#define _ANAINPUT_HPP_

#include <analyzer/StructVar.hpp>

namespace hig{
  typedef std::vector<StructVar> ana_struct_vec_t;
  typedef std::vector<StructVar>::iterator  ana_struct_vec_it;

  class AnaInput{
  private :
    ana_struct_vec_t struct_vec_;

  public:
    AnaInput(){}
    ~AnaInput(){}

    /*modifiers*/
    void push(StructVar sv){ struct_vec_.push_back(sv);}

    /*getters*/
    StructVar get_var(int i){return struct_vec_[i];}
    StructVar get_var(string_t uid);
    float get_var_val(int i){return struct_vec_[i].get_val();}
    float get_var_init_val(int i){return struct_vec_[i].get_init_val();}

    /*setters*/
    bool set_var_val(int i, float val);
    bool set_var_val(string_t uid, float val);
    bool set_var_final_val(int i, float val);

    void generate_random_vars(float min, float max, int dim);

    void clear(){ struct_vec_.clear();}
    int size(){return struct_vec_.size();}
    void print();

  }; /* class AnaInput */

} /* namespace hig */

#endif /* ANAINPUT_HPP_ */
