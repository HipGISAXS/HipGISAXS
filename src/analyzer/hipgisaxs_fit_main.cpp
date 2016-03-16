/**
 *  Project:
 *
 *  File: hipgisaxs_fit_main.cpp
 *  Created: Feb 25, 2014
 *  Modified: Wed 08 Oct 2014 12:17:42 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */
#include <iostream>

#include <analyzer/hipgisaxs_fit_pounders.hpp>
#include <analyzer/hipgisaxs_fit_lmvm.hpp>
#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <analyzer/hipgisaxs_fit_bruteforce.hpp>

#include <analyzer/distance_functions.hpp>
#include <analyzer/objective_func.hpp>
#include <analyzer/objective_func_hipgisaxs.hpp>
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

  // for PSO
  /*bool tune_omega = false;
  if(narg > 7) tune_omega = (atoi(args[7]) == 1);
  int type = 0;
  if(narg > 8) {
    std::string type_str(args[8]);
    if(type_str.compare("base") == 0) {
      type = 0;
    } else if(type_str.compare("fips") == 0) {
      type = 1;
    } else if(type_str.compare("foresee") == 0) {
      type = 2;
    } else if(type_str.compare("fdr") == 0) {
      type = 3;
    } else if(type_str.compare("bb") == 0) {
      std::cerr << "WARNING: barebones is not complete yet!" << std::endl;
      type = 4;
    } else if(type_str.compare("lbest") == 0) {
      type = 5;
    } else if(type_str.compare("von") == 0) {
      type = 6;
    } else if(type_str.compare("random") == 0) {
      type = 7;
    } else {
      type = -1;
      std::cerr << "error: invalid type given. valid types are: base fips foresee" << std::endl;
      return -1;
    } // if-else
  } // if*/

  for(int i = 0; i < hig::HiGInput::instance().num_analysis_algos(); ++ i) {
    hig::FittingAlgorithmName algo = hig::HiGInput::instance().analysis_algo(i);

    if(i > 0) {
      std::cout << "warning: currently hipgisaxs supports only one analysis algorithm at a time"
            << std::endl;
      break;
    } // if

    if(algo == hig::algo_pounders) {
      hip_func.set_distance_measure(new ResidualVector());
      ana.add_analysis_algo(new hig::FitPOUNDERSAlgo(narg, args, &hip_func, i));
    } else if(algo == hig::algo_lmvm) {
      hip_func.set_distance_measure(new AbsoluteDifferenceSquare());
      ana.add_analysis_algo(new hig::FitLMVMAlgo(narg, args, &hip_func, i));
    } else if(algo == hig::algo_pso) {
      hip_func.set_distance_measure(new AbsoluteDifferenceSquareNorm());
      ana.add_analysis_algo(new hig::ParticleSwarmOptimization(narg, args, &hip_func, i, false, 0));
    } else if(algo == hig::algo_bruteforce) {
      hip_func.set_distance_measure(new AbsoluteDifferenceSquareNorm());
      ana.add_analysis_algo(new hig::BruteForceOptimization(narg, args, &hip_func, i));
    } else if(algo == hig::algo_error) {
      std::cerr << "error: unknown optimization algorithm encountered" << std::endl;
      return -1;
    } else {
      std::cerr << "error: NULL optimization algorithm encountered" << std::endl;
      return -1;
    } // if-else

  } // for

  ana.analyze(narg, args, 1);    // perform the analysis

  return 0;
} // main()
