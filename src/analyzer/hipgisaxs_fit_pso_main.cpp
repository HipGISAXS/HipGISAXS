/**
 *  Project:
 *
 *  File: pso.cpp
 *  Created: Jan 13, 2014
 *  Modified: Wed 05 Feb 2014 10:29:25 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <analyzer/objective_func.hpp>
#include <analyzer/hipgisaxs_ana.hpp>

int main(int narg, char** args) {
	if(narg != 7) {
		std::cout << "usage: hipgisaxs_pso <input_config> <num_particles> <num_generations> "
			<< "<omega> <phi1> <phi2>"
			<< std::endl;
		return 1;
	} // if

	//AbsoluteDifferenceError err;
	AbsoluteDifferenceNorm err;
	hig::HipGISAXSObjectiveFunction hip_func(narg, args, &err);
	hig::ParticleSwarmOptimization my_pso(narg, args, &hip_func);

	hig::HipGISAXSAnalyzer ana;
	ana.add_analysis_algo(&my_pso);

	ana.analyze(narg, args);

	return 0;
} // main()
