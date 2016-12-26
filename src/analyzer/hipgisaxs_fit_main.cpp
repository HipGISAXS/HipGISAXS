/**
 *  Project:
 *
 *  File: hipgisaxs_fit_main.cpp
 *  Created: Feb 25, 2014
 *  Modified: Tue 24 May 2016 08:04:50 PM EDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */
#include <iostream>

#include <analyzer/hipgisaxs_fit_pounders.hpp>
#include <analyzer/hipgisaxs_fit_lmvm.hpp>
#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <analyzer/hipgisaxs_fit_bruteforce.hpp>
#include <analyzer/hipgisaxs_compute_objective.hpp>

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

  for(int i = 0; i < hip_func.num_analysis_algos(); ++ i) {
    if(i > 0) {
      std::cout << "warning: currently hipgisaxs supports only one analysis algorithm at a time"
            << std::endl;
      break;
    } // if

    hig::FittingAlgorithmName algo = hip_func.analysis_algo(i);
    hig::FittingDistanceMetric dist_metric = hip_func.analysis_distance_metric(i);

    // set the distance metric
    switch(dist_metric) {
      case hig::metric_sqrt_unit_norm_l2:
        hip_func.set_distance_measure(new SqrtUnitVectorNormL2DistanceSquare());
        break;
      case hig::metric_sqrt_c_norm_l2:
        hip_func.set_distance_measure(new SqrtCNormL2DistanceSquare());
        break;
      case hig::metric_cbrt_unit_norm_l2:
        hip_func.set_distance_measure(new CbrtUnitVectorNormL2DistanceSquare());
        break;
      case hig::metric_cbrt_c_norm_l2:
        hip_func.set_distance_measure(new CbrtCNormL2DistanceSquare());
        break;
      case hig::metric_sqrt_unit_norm_l2_residual:
        hip_func.set_distance_measure(new SqrtUnitVectorNormL2DistanceSquareResidualVector());
        break;
      case hig::metric_sqrt_c_norm_l2_residual:
        hip_func.set_distance_measure(new SqrtCNormL2DistanceSquareResidualVector());
        break;
      default:
        std::cerr << "error: unsupported distance metric encountered" << std::endl;
        return -1;
    } // switch

    // set the analysis algorithm
    switch(algo) {
      case hig::algo_pounders:      // pounders fitting (uses petsc/tao)
        ana.add_analysis_algo(new hig::FitPOUNDERSAlgo(narg, args, &hip_func, i));
        break;
      case hig::algo_lmvm:          // lmvm fitting (uses petsc/tao)
        ana.add_analysis_algo(new hig::FitLMVMAlgo(narg, args, &hip_func, i));
        break;
      case hig::algo_pso:           // particle swarm optimization
        ana.add_analysis_algo(new hig::ParticleSwarmOptimization(narg, args, &hip_func, i, false, 0));
        break;
      case hig::algo_bruteforce:    // brute force: try all possibilities
        ana.add_analysis_algo(new hig::BruteForceOptimization(narg, args, &hip_func, i));
        break;
      case hig::algo_none_pounders: // compute the objective function, do not fit
        ana.add_analysis_algo(new hig::ComputeObjectiveFunction(narg, args, &hip_func, i));
        break;
      case hig::algo_error:
        std::cerr << "error: unknown optimization algorithm encountered" << std::endl;
        return -1;
      default:
        std::cerr << "error: unsupported optimization algorithm encountered" << std::endl;
        return -1;
    } // switch

  } // for

  ana.analyze(narg, args, 1);    // perform the analysis

  return 0;
} // main()
