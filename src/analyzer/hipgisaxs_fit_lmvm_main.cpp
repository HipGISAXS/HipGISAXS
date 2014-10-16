#include <iostream>

#include <analyzer/distance_functions.hpp>
#include <analyzer/hipgisaxs_fit_lmvm.hpp>
#include <analyzer/objective_func_hipgisaxs.hpp>
#include <analyzer/hipgisaxs_ana.hpp>


/* The main for LMVM HipGISAXS Analyzer
 */
int main(int narg, char** args) {

  if(narg < 2 || narg > 4) {
    std::cout << "usage: hipgisaxs <input_config>" << std::endl;
    return 1;
  } // if

  AbsoluteDifferenceSquare err;                  // define the distance measure
  hig::HipGISAXSObjectiveFunction hip_func(narg, args, &err);    // define the objective function
  hig::FitLMVMAlgo lmvm (&hip_func);                    // define analysis algo

  hig::HipGISAXSAnalyzer ana;                      // define the analyzer
  ana.add_analysis_algo(&lmvm);                    // add analysis algo
                                    // similarly add more if needed
  ana.analyze(narg, args, 1);                    // perform the analysis
                                // takes flag >= 0 for reading ref data
                                // and flag < 0 for computing ref data
  return 0;
} // main()
