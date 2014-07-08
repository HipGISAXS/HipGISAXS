/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Workflow.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:26 AM PST
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

#ifndef _WORKFLOW_HPP_
#define _WORKFLOW_HPP_

//#include "FitPOUNDERSAlgo.hpp"
#include <analyzer/AnaAlgorithm.hpp>

#include <queue>

namespace hig{

  typedef std::queue<AnaAlgorithm*> wfq_t;

  class Workflow{
  private :
    wfq_t wfq_;
    int is_void_;
  public:
    Workflow(){}
    Workflow(string_t qf_str);
    Workflow(wfq_t wfq);
    Workflow(const Workflow& wf);
    ~Workflow(){}

    /* copy */
    Workflow& operator=(const Workflow& wf);

    AnaAlgorithm* create_algo(AlgoType at);
    bool init(){clear();return true;}
    bool set(wfq_t wfq);
    void enq(string_t algo_str);
    void enq(AnaAlgorithm* algo);
    AnaAlgorithm* deq();
    void clear();
    bool empty(){return wfq_.empty();}

    /* accessor */
    wfq_t get_wfq() const {return wfq_;}

    void print();

  }; /* class Workflow */

} /* namespace hig */

#endif /* WORKFLOW_HPP_ */
