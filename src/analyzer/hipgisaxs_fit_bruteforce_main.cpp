/**
 *  Project:
 *
 *  File: hipgisaxs_fit_bruteforce_main.cpp
 *  Created: Jan 13, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <analyzer/hipgisaxs_fit_bruteforce.hpp>
#include <analyzer/objective_func.hpp>
#include <analyzer/hipgisaxs_ana.hpp>

int main(int narg, char** args) {
  if(narg < 2 || narg > 3) {
    std::cout << "usage: analyze <input_config> [<history_output_filename>]" << std::endl;
    return 1;
  } // if

  //AbsoluteDifferenceError err;
  //AbsoluteDifferenceNorm err;
  AbsoluteDifferenceSquareNorm err;
  hig::HipGISAXSObjectiveFunction hip_func(narg, args, &err);
  hig::BruteForceOptimization my_bfo(narg, args, &hip_func);

  hig::HipGISAXSAnalyzer ana;
  ana.add_analysis_algo(&my_bfo);

  ana.analyze(narg, args);

  return 0;
} // main()
