/**
 *  Project: HipGISAXS
 *
 *  File: hipgisaxs_compute_objective.hpp
 *  Created: May 28, 2016
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __HIPGISAXS_COMPUTE_OBJECTIVE_HPP__
#define __HIPGISAXS_COMPUTE_OBJECTIVE_HPP__

#include <analyzer/analysis_algorithm.hpp>

namespace hig {

  class ComputeObjectiveFunction: public AnalysisAlgorithm {

    private:

      unsigned int num_obs_;

    public:
      ComputeObjectiveFunction(int, char**, ObjectiveFunction*, unsigned int);
      ~ComputeObjectiveFunction() { }

      bool run(int argc, char **argv, int, int);

  }; /* class ComputeObjectiveFunction  */

} /* namespace hig */

#endif /* __HIPGISAXS_COMPUTE_OBJECTIVE_HPP__ */
