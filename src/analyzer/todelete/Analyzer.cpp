/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Analyzer.cpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 08 Oct 2014 12:17:47 PM PDT
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
#include <analyzer/Analyzer.hpp>

namespace hig{

  bool Analyzer::setup(){
    return false;
  }

  bool Analyzer::init(){
    return false;
  }

  bool Analyzer::set_obj_fct(ObjFct* pof){
    pobj_fct_ = pof;
    return true;
  }

  bool Analyzer::set_workflow(Workflow wf){
    wf_.set(wf.get_wfq());
    return true;
  }

  bool Analyzer::set_workflow(string_t wf_str){
    //  return    wf_.init(wf_str);
  }

  bool Analyzer::set_data(ImageDataProc data){
    data_ = data;
  }

  bool Analyzer::set_data(string_t filename){
    //data_.init(filename);
  }

  bool Analyzer::set_input(AnaInput inp){
    inputs_ = inp;
  }

  bool Analyzer::set_ref_data( ImageData*  pimg){
    pobj_fct_->set_ref_data(pimg);
  }

  void Analyzer::write_final_vec(float_vec_t XN){

    for(int i=0; i<inputs_.size();i++){
      inputs_.set_var_final_val(i, XN[i]);
      inputs_.set_var_val(i, XN[i]);
    }
  }

  int Analyzer::analyze(int argc,char **argv){
    img_data_vec_t img_data_set = data_.get_data();
    float_vec_t XN;
    for(img_data_vec_t::iterator it=img_data_set.begin(); it!=img_data_set.end(); it++){
      std::cout << (*it).get_name()   << std::endl;
      set_ref_data(&(*it));
      Workflow wf_c = wf_;
      while (!wf_c.empty()){
  AnaAlgorithm* palgo = wf_c.deq();
  palgo->run(argc,argv);
  XN = palgo->get_final_vec();
      } // run the workflow queue
    }// loop over image data set

    /* write the final vector into input_  */
    write_final_vec(XN);
  }

  void Analyzer::print_output(){

    inputs_.print();

  }


}
