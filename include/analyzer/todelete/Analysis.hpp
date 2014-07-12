/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Analysis.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:21:49 AM PST
 *  Description: Stores and manipulates analysis output
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

#ifndef _ANALYSIS_HPP_
#define _ANALYSIS_HPP_

#include <sstream>
#include <analyzer/typedefs.hpp>

namespace hig{

  class Analysis{

  private :

    //    std::ostringstream  log_;
    //std::ostringstream  out_;

  public:
    Analysis(){}
    ~Analysis(){}
    bool write(string_t filename);
    bool plot();
    bool print();
  }; /* class Analysis */

} /* namespace hig */

#endif /* ANALYSIS_HPP_ */
