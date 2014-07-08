/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Analyzer.hpp
 *  Created: Dec 26, 2013
 *  Modified: Thu 27 Feb 2014 11:02:21 AM PST
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

#include <vector>

#include <analyzer/typedefs.hpp>
#include <analyzer/analysis_algorithm.hpp>

//#include <analyzer/Workflow.hpp>
//#include <analyzer/ImageDataProc.hpp>
//#include <analyzer/objective_functions.hpp>

namespace hig {

  class HipGISAXSAnalyzer {

  private :
    string_t title_;
    //Workflow wf_;
    //ImageDataProc data_;
    //ObjFct* pobj_fct_;
    //AnaOutput outputs_;
	std::vector <AnalysisAlgorithm*> wf_;		// the workflow
    bool is_valid;

  public:
    HipGISAXSAnalyzer() { }
    ~HipGISAXSAnalyzer() { }

    bool init();
    void set_title(string_t title) { title_ = title; }

    //bool set_obj_fct(ObjFct* pof);

    //bool set_workflow(Workflow wf);
    //bool set_workflow(string_t wf_str);

    //bool set_data(ImageDataProc data);
    //bool set_data(string_t filename);
    //bool set_ref_data(ImageData* pimg);

    //bool set_input(AnaInput inp);
	bool add_analysis_algo(AnalysisAlgorithm* algo) { wf_.push_back(algo); return true; }

    int analyze(int argc,char **argv, int);
    //AnaOutput get_output(){return outputs_;}
    //void write_final_vec(float_vec_t XN);

    void print_output();

  }; /* class HipGISAXSAnalyzer */

} /* namespace hig */

#endif /* ANALYZER_HPP_ */
