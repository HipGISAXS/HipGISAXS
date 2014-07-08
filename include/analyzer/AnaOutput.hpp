/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaOutput.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:21:54 AM PST
 *  Description: Aggregates all analysis output from all workflows
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

#ifndef _ANAOUTPUT_HPP_
#define _ANAOUTPUT_HPP_

#include <analyzer/Analysis.hpp>

namespace hig{

  typedef std::vector<Analysis> ana_output_vec_t;
  //typedef std::vector<Analysis>::const_iterator  ana_output_vec_cit;

  class AnaOutput{
  private :
    ana_output_vec_t outputs_;
  public:
    AnaOutput(){}
    ~AnaOutput(){}
    bool aggregate();
    bool write(string_t filename);
    bool plot();
    bool print();
  }; /* class AnaOutput */

} /* namespace hig */

#endif /* ANAOUTPUT_HPP_ */
