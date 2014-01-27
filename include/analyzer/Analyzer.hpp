/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Analyzer.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:21:52 AM PST
 *  Description: The main analysis class that executes the workflows defined therein wrt
 *  the inputs (a priori structural info)  and datasets (expt. data) provided.
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

#ifndef _ANALYZER_HPP_
#define _ANALYZER_HPP_

#include <analyzer/Workflow.hpp>
#include <analyzer/ImageDataProc.hpp>
#include <analyzer/AnaInput.hpp>
#include <analyzer/AnaOutput.hpp>
#include <analyzer/ObjFct.hpp>

namespace hig{

  class Analyzer{

  private :
    string_t title_;
    Workflow wf_;
    ImageDataProc data_;
    AnaInput inputs_;
    ObjFct* pobj_fct_;
    AnaOutput outputs_;
    bool is_valid;

  public:
    Analyzer(){}
    ~Analyzer(){}

    bool setup();
    bool init();

    void set_title(string_t title) { title_ = title;}

    bool set_obj_fct(ObjFct* pof);

    bool set_workflow(Workflow wf);
    bool set_workflow(string_t wf_str);

    bool set_data(ImageDataProc data);
    bool set_data(string_t filename);
    bool set_ref_data(ImageData* pimg);

    bool set_input(AnaInput inp);

    AnaOutput get_output(){return outputs_;}
    void write_final_vec(float_vec_t XN);

    int analyze(int argc,char **argv);

    void print_output();

  }; /* class Analyzer */

} /* namespace hig */

#endif /* ANALYZER_HPP_ */
