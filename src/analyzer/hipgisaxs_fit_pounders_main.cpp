#include <iostream>

#include <analyzer/distance_functions.hpp>
#include <analyzer/hipgisaxs_fit_pounders.hpp>
#include <analyzer/objective_func.hpp>
#include <analyzer/hipgisaxs_ana.hpp>


/* The main for POUNDers HipGISAXS Analyzer
 */
int main(int narg, char** args) {

	if(narg < 2 || narg > 4) {
		std::cout << "usage: hipgisaxs <input_config>" << std::endl;
		return 1;
	} // if

	ResidualVector err;												// define the distance measure
	hig::HipGISAXSObjectiveFunction hip_func(narg, args, &err);		// define the objective function
	hig::FitPOUNDERSAlgo pounders(&hip_func);						// define analysis algo

	hig::HipGISAXSAnalyzer ana;										// define the analyzer
	ana.add_analysis_algo(&pounders);								// add analysis algo
																	// similarly add more if needed

	ana.analyze(narg, args);										// perform the analysis

	return 0;
} // main()
