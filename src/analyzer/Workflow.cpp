/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Workflow.cpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:21:04 AM PST
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
#include <analyzer/Workflow.hpp>
//#include "FitPOUNDERSAlgo.hpp"

namespace ana_hig{

  Workflow::Workflow(string_t wf_str){
    clear();
    AlgoType at =  AnaAlgorithm::get_type(wf_str);
    enq(create_algo(at));
  }

  Workflow::Workflow(wfq_t wfq){
    set(wfq);
  }

  Workflow::Workflow(const Workflow& wf){
   this->set(wf.get_wfq());
  }

  Workflow& Workflow::operator=(const Workflow& wf){
    this->set(wf.get_wfq());
  }

  AnaAlgorithm* Workflow::create_algo(AlgoType at){
    //AnaAlgorithm ana;
    switch(at){
    case fit_pounders:{
      //std::cout << "POUNDERS\n" ;
      return new  AnaAlgorithm() ;//FitPOUNDERSAlgo();
    }
    default:{
      //std::cout << "BASIC\n" ;
      return new AnaAlgorithm();
    }
    }
  }

  bool Workflow::set(wfq_t wfq){
    wfq_ = wfq;
    return true;
  }

  AnaAlgorithm* Workflow::deq(){
    AnaAlgorithm* fr_algo = wfq_.front();
    //std::cout << fr_algo->get_type_string() << std::endl ;
    wfq_.pop();
    return fr_algo;
  }

  void Workflow::enq(string_t algo_str){
    AlgoType at =  AnaAlgorithm::get_type(algo_str);
    enq(create_algo(at));
  }

  void Workflow::enq(AnaAlgorithm* algo){
    wfq_.push(algo);
  }

  void Workflow::clear(){
    while(!wfq_.empty())
      {
	deq();
      }
  }


  void Workflow::print(){
    std::cout << "Workflow:\n";
    wfq_t wfq_c = wfq_;
    int cnt=1;
    while(!wfq_c.empty())
      {
	AnaAlgorithm* fr_algo = wfq_c.front();
	std::cout << cnt << ": ";
	fr_algo->print();
	wfq_c.pop();
	cnt++;
      }
  }
}
