/**
 *  Project:
 *
 *  File: pso.cpp
 *  Created: Jan 13, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <analyzer/objective_func_hipgisaxs.hpp>
#include <analyzer/hipgisaxs_ana.hpp>

int main(int narg, char** args) {
  if(narg < 7 || narg > 9) {
    std::cout << "usage: hipgisaxs_pso <input_config> <num_particles> <num_generations> "
      << "<omega> <phi1> <phi2> [<tune omega>] [<type>]"
      << std::endl;
    return 1;
  } // if

  bool tune_omega = false;
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
  } // if

  //AbsoluteDifferenceError err;
  //AbsoluteDifferenceNorm err;
  AbsoluteDifferenceSquareNorm err;
  hig::HipGISAXSObjectiveFunction hip_func(narg, args, &err);
  hig::ParticleSwarmOptimization my_pso(narg, args, &hip_func,
                      atof(args[4]), atof(args[5]), atof(args[6]),
                      atoi(args[2]), atoi(args[3]), tune_omega, type);
  hig::HipGISAXSAnalyzer ana;
  ana.add_analysis_algo(&my_pso);

  woo::BoostChronoTimer maintimer;

  maintimer.start();
  ana.analyze(narg, args, 1);
  maintimer.stop();

  hig::parameter_map_t result = my_pso.get_best_values();
  if(my_pso.is_master()) {
    std::cout << "** ** Final parameter values: " << std::endl;
    for(hig::parameter_map_t::const_iterator i = result.begin(); i != result.end(); ++ i)
      std::cout << "      ++ " << (*i).first << " = " << (*i).second << std::endl;
    std::cout << "** ** TOTAL ANALYSIS TIME: " << maintimer.elapsed_msec() << " ms. ** **"
          << std::endl;
  } // if

  return 0;
} // main()
