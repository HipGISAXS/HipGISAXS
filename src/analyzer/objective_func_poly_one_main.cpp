#include <iostream>

#include <analyzer/distance_functions.hpp>
#include <analyzer/FitPOUNDERSAlgo.hpp>
#include <analyzer/objective_func_poly_one.hpp>
#include <analyzer/hipgisaxs_ana.hpp>


/* The main
 */
int main(int narg, char** args) {

	if(narg < 11) {
		std::cout << "usage: poly_one <reference file> <q_min> <q_max> <a_step> <x1_min> <x1_max> <x2_min> <x2_max> <x1_init> <x2_init>" << std::endl;
		return 1;
	} // if

	ResidualVector err;												// define the distance measure
	hig::PolyOneObjectiveFunction func(narg, args, &err);		// define the objective function
	hig::FitPOUNDERSAlgo pounders(&func);						// define analysis algo

	hig::HipGISAXSAnalyzer ana;										// define the analyzer
	ana.add_analysis_algo(&pounders);								// add analysis algo
																	// similarly add more if needed

	ana.analyze(narg, args);										// perform the analysis

	return 0;
} // main()
