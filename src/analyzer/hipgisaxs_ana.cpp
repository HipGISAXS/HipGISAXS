/**
 *  Project: HipGISAXS
 *
 *  File: Analyzer.cpp
 *  Created: Dec 26, 2013
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
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

#include <analyzer/hipgisaxs_ana.hpp>
#include <hipgisaxs.hpp>

namespace hig {

  bool HipGISAXSAnalyzer::analyze(int argc, char **argv, int flag) {

    std::vector <real_vec_t> all_results;
    for(int i = 0; i < HiGInput::instance().num_analysis_data(); ++ i) {
      for(int j = 0; j < wf_.size(); ++ j) {
        if(flag < 0) wf_[j]->run(argc, argv, -1);    // ref data to be computed
        else wf_[j]->run(argc, argv, i);        // ref data to be read
        all_results.push_back(wf_[j]->get_param_values());
      } // for
    } // for

    return true;
  } // HipGISAXSAnalyzer::analyze()

} // namespace hig

