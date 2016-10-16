/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_ana.hpp
 *  Created: Dec 26, 2013
 *  Description: The main analysis class that executes the workflows defined therein wrt
 *  the inputs (a priori structural info)  and datasets (expt. data) provided.
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __HIPGISAXS_ANA_HPP__
#define __HIPGISAXS_ANA_HPP__

#include <vector>

#include <common/typedefs.hpp>
#include <analyzer/analysis_algorithm.hpp>
#include <config/input.hpp>


namespace hig {

  class HipGISAXSAnalyzer {

    private :
      std::vector <AnalysisAlgorithm*> wf_;   // the workflow
      int num_algo_;                          // the input container

    public:
      HipGISAXSAnalyzer() { }
      ~HipGISAXSAnalyzer() { }

      bool add_analysis_algo(AnalysisAlgorithm* algo) {
        wf_.push_back(algo);
        num_algo_ = wf_.size();
        return true;
      } // add_analysis_algo()

      bool analyze(int argc,char **argv, int);

      real_t tolerance(int n) const { return wf_[n]->tolerance(); }

  }; // class HipGISAXSAnalyzer

} // namespace hig

#endif // __HIPGISAXS_ANA_HPP__
