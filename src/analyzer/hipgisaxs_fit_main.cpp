/**
 *  Project:
 *
 *  File: hipgisaxs_fit_main.cpp
 *  Created: Feb 25, 2014
 *  Modified: Wed 09 Jul 2014 11:56:33 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */
#include <iostream>

#include <analyzer/hipgisaxs_fit_pounders.hpp>
#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <analyzer/hipgisaxs_fit_bruteforce.hpp>

#include <analyzer/distance_functions.hpp>
#include <analyzer/objective_func.hpp>
#include <analyzer/hipgisaxs_ana.hpp>

//#include <config/hig_input.hpp>


/**
 * The common main for all fitting. Uses algorithm data from the input file
 */
int main(int narg, char** args) {
	if(narg < 2) {
		std::cout << "usage: hipgisaxs <input_config>" << std::endl;
		return 1;
	} // if

	hig::HipGISAXSObjectiveFunction hip_func(narg, args, args[1]);
	hig::HipGISAXSAnalyzer ana;

	for(int i = 0; i < hig::HiGInput::instance().num_analysis_algos(); ++ i) {
		hig::FittingAlgorithmName algo = hig::HiGInput::instance().analysis_algo(i);

		if(i > 0) {
			std::cout << "warning: currently hipgisaxs supports only one analysis algorithm at a time"
						<< std::endl;
			break;
		} // if

		if(algo == hig::algo_pounders) {
			hip_func.set_distance_measure(new ResidualVector());
			ana.add_analysis_algo(new hig::FitPOUNDERSAlgo(&hip_func));
		} else if(algo == hig::algo_pso) {
			hip_func.set_distance_measure(new AbsoluteDifferenceSquareNorm());
			ana.add_analysis_algo(new hig::ParticleSwarmOptimization(narg, args, &hip_func, i));
		} else if(algo == hig::algo_bruteforce) {
			hip_func.set_distance_measure(new AbsoluteDifferenceSquareNorm());
			ana.add_analysis_algo(new hig::BruteForceOptimization(narg, args, &hip_func));
		} else if(algo == hig::algo_error) {
			std::cerr << "error: unknown optimization algorithm encountered" << std::endl;
			return -1;
		} else {
			std::cerr << "error: NULL optimization algorithm encountered" << std::endl;
			return -1;
		} // if-else

	} // for

	ana.analyze(narg, args, 1);		// perform the analysis

	return 0;
} // main()
