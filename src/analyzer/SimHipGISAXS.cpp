/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: SimHipGISAXS.cpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:20:53 AM PST
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

#include <analyzer/SimHipGISAXS.hpp>

namespace ana_hig{

  SimHipGISAXS::SimHipGISAXS(string_t higfile){


  }

  /* MAPPER METHODS BEWTEEN STRUCTVAR & HIPGISAXS PARAMS */
  /*
  int SimHipGISAXS::set_struct_param(string_t uid, float val){
    VarPath vp = parse_param_uid(uid);
    return set_value(vp , val);
  }

  int SimHipGISAXS::set_value(VarPath vp, float val){

    return 0;
  }

  float SimHipGISAXS::get_value(VarPath vp){

    return 0;
  }

  VarPath SimHipGISAXS::parse_param_uid(string_t uid){

  }
  */
  /******************************************************/

  int SimHipGISAXS::compute(){

    return 0;
  }

}
