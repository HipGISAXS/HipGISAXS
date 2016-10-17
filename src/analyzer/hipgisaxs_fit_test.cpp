/**
 *  Project:
 *
 *  File: hipgisaxs_fit_test.cpp
 *  Created: Feb 23, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <string>
#include <map>
#include <vector>

typedef std::map<std::string, std::vector<std::string> > arg_map_t;

#include <hipgisaxs_fit.hpp>


bool parse_arguments(int narg, char** args, arg_map_t& arg_map) {
  arg_map.clear();

  std::string key;
  for(int i = 1; i < narg; ++ i) {
    if(args[i][0] == '-') {
      key = std::string(args[i]);
      //arg_map[key].clear();
    } else {
      arg_map[key].push_back(std::string(args[i]));
    } // if-else
  } // for

  std::cout << "** Input arguments:" << std::endl;
  for(arg_map_t::const_iterator i = arg_map.begin(); i != arg_map.end(); ++ i) {
    std::cout << "    " << (*i).first << " =";
    for(std::vector<std::string>::const_iterator j = (*i).second.begin();
        j != (*i).second.end(); ++ j) {
      std::cout << " " << (*j);
    } // for
    std::cout << std::endl;
  } // for
} // parse_arguments()


int main(int narg, char** args) {
  if(narg < 2) {
    std::cout << "usage: analyze -model <model_name> -algo <ana_algo> "
          << "[-i <input_config>] "
          << "[-init <val1> <val2> ...] "
          << "[-ref <val1> <val2> ...] " << std::endl;
    std::cout << "  valid values for <model_name>: hipgisaxs, f1, f2, f3" << std::endl;
    std::cout << "  valid values for <ana_algo>  : pounders, pso, bruteforce" << std::endl;
    return 0;
  } // if

  arg_map_t arg_map;
  parse_arguments(narg, args, arg_map);

  DistanceMeasure* dist = NULL;
  hig::HipGISAXSAnalyzer hig_ana;

  // choose model
  hig::ObjectiveFunction* obj_func = NULL;

  std::string model = arg_map["-model"][0];
  if(model.compare(std::string("hipgisaxs")) == 0) {
    obj_func = new hig::HipGISAXSObjectiveFunction(narg, args, arg_map["-i"][0]);
  } else if(model.compare(std::string("f1")) == 0) {
  } else if(model.compare(std::string("f2")) == 0) {
  } else if(model.compare(std::string("f3")) == 0) {
  } else {
    std::cerr << "error: unknown model specified" << std::endl;
    std::cerr << "valid values are: hipgisaxs, f1, f2, f3" << std::endl;
    return -1;
  } // if-else
  
  // set reference parameters
  std::vector<hig::real_t> ref_param_vals;
  if(arg_map.count("-ref") == 1) {
    for(std::vector<std::string>::const_iterator i = arg_map.at("-ref").begin();
        i != arg_map.at("-ref").end(); ++ i) {
      ref_param_vals.push_back(atof((*i).c_str()));
    } // for
    //(*obj_func).update_params(ref_param_vals);
  } // if

  // simulate reference and set ref data
  (*obj_func).simulate_and_set_ref(ref_param_vals);

  // choose algo
  hig::AnalysisAlgorithm* ana_algo = NULL;

  if(arg_map.count("-algo") < 1) {
    std::cerr << "error: analysis algorithm not specified" << std::endl;
    std::cerr << "valid values are: pounders, pso, bruteforce" << std::endl;
    return -1;
  } // if
  std::string algo(arg_map.at("-algo")[0]);
  if(algo.compare("pounders") == 0) {
    //(*obj_func).set_distance_measure(new ResidualVector());
    (*obj_func).set_distance_measure(new RelativeResidualVector());
    ana_algo = new hig::FitPOUNDERSAlgo(narg, args, obj_func, 0);
  } else if(algo.compare("pso") == 0) {
    if(arg_map.count("-pso_omega") < 1 ||
        arg_map.count("-pso_phi1") < 1 ||
        arg_map.count("-pso_phi2") < 1 ||
        arg_map.count("-pso_npart") < 1 ||
        arg_map.count("-pso_ngen") < 1) {
      std::cerr << "error: missing parameters for PSO. " << std::endl
            << "  parameters are: -pso_omega, -pso_phi1, -pso_phi2, -pso_npart, -pso_ngen"
            << std::endl;
      return -1;
    } // if
    //(*obj_func).set_distance_measure(new AbsoluteDifferenceSquareNorm());
    (*obj_func).set_distance_measure(new ScaledRelativeAbsoluteDifferenceSquare());
    ana_algo = new hig::ParticleSwarmOptimization(narg, args, obj_func,
                            atof(arg_map.at("-pso_omega")[0].c_str()),
                            atof(arg_map.at("-pso_phi1")[0].c_str()),
                            atof(arg_map.at("-pso_phi2")[0].c_str()),
                            atoi(arg_map.at("-pso_npart")[0].c_str()),
                            atoi(arg_map.at("-pso_ngen")[0].c_str()),
                            false, 0);
  } else if(algo.compare("bruteforce") == 0) {
    (*obj_func).set_distance_measure(new AbsoluteDifferenceSquareNorm());
    ana_algo = new hig::BruteForceOptimization(narg, args, obj_func, 0);
  } else if(algo.compare ("lmvm") == 0) {
    //(*obj_func).set_distance_measure(new AbsoluteDifferenceSquareNorm());
    (*obj_func).set_distance_measure(new RelativeAbsoluteDifferenceSquare());
    ana_algo = new hig::FitLMVMAlgo(narg, args, obj_func, 0);
  } else {
    std::cerr << "error: unknown analysis algorithm specified" << std::endl;
    std::cerr << "valid values are: pounders, pso, bruteforce" << std::endl;
    return -1;
  } // if-else

  // set initial parameter values
  std::vector<hig::real_t> param_vals;
  if(arg_map.count(std::string("-init")) == 1) {
    for(std::vector<std::string>::const_iterator i = arg_map.at(std::string("-init")).begin();
        i != arg_map.at(std::string("-init")).end(); ++ i) {
      param_vals.push_back(atof((*i).c_str()));
    } // for
    (*ana_algo).init_params(param_vals);
  } // if


  // set analysis algo
  hig_ana.add_analysis_algo(ana_algo);

  // perform the analysis
  hig_ana.analyze(narg, args, -1);

  delete ana_algo;
  delete obj_func;

  return 0;
} // main()
