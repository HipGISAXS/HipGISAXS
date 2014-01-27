/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: SimHipGISAXS.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:13 AM PST
 *  Description: Class that computes a forward HipGISAXS simulation model.
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

#ifndef _SIMHIPGISAXS_HPP_
#define _SIMHIPGISAXS_HPP_

#include <iostream>
#include <analyzer/SimModel.hpp>
//#include <hipgisaxs_main.hpp>

namespace hig{

  /* Idea: A variable can be accessed via a path down a syntax tree  */
  struct VarPath{
    string_t key_;
    int ind_;
    VarPath* next_;
  };

  class SimHipGISAXS : public SimModel{

  protected :
    //    hig::HipGISAXS *phig_; /* hipgisaxs obj used in fwd sim.   */

  public:
    SimHipGISAXS(){}
    SimHipGISAXS(string_t higfile); /* generates the hipgisaxs obj. for fwd sim.   */

    ~SimHipGISAXS(){}
    virtual  bool init(){is_valid_=false; return is_valid_;}

    /* computer  */
    virtual int compute();

    /* accessors  */
    virtual int update_vars(float_vec_t X){}
    //    hig::HipGISAXS* get_hipgisaxs() {  return phig_; }

    /* MOST IMPORTANT PARAM ACCESSORS  */
    int set_struct_param(string_t uid, float val); // parses uid to get path and sets value
    int set_value(VarPath vp , float val); //uses path to set value
    float get_value(VarPath vp); // uses path to read value
    /********************************/

  }; /* class SimHipGISAXS */

} /* namespace hig */

#endif /* SIMHIPGISAXS_HPP_ */
