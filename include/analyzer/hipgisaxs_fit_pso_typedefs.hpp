/**
 *  Project:
 *
 *  File: typedefs.hpp
 *  Created: Jan 14, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __HIPGISAXS_FIT_PSO_TYPEDEFS_HPP__
#define __HIPGISAXS_FIT_PSO_TYPEDEFS_HPP__

#include <vector>
#include <map>
#include <string>
#include <common/typedefs.hpp>  // from hipgisaxs library

namespace hig {

  typedef std::map <std::string, real_t>    parameter_map_t;
  typedef std::vector <real_t>        parameter_data_list_t;
  typedef std::vector <std::string>      parameter_name_list_t;

  enum PSOParameterDistribution {
    PSO_UNIFORM,                  // randomized with uniform distribution
    PSO_GAUSSIAN,                  // randomized with gaussian distribution
    PSO_SINGLE,                    // values set to given values
    PSO_DEFAULT                    // currently same as PSO_UNIFORM
  }; // enum PSOParameterDistribution
  typedef PSOParameterDistribution pso_parameter_dist_t;

  // PSO algorithm types
  enum PSOAlgorithmType {
    pso_algo_base = 0,      // base algorithm (inertia weight)
    pso_algo_fips,        // FIPS algorithm
    pso_algo_soothsayer,    // Look-ahead/soothsayer algorithm
    pso_algo_fdr,        // FDR (Fitness Distance Ratio) algorithm
    pso_algo_barebones,      // barebones algorithm
    pso_algo_lbest,        // LBest topology (K = 2)
    pso_algo_vonnewmann,    // Von Newmann topology (K = 4)
    pso_algo_random,      // Random topology with 20% connectivity
    pso_algo_null,        // no algorithm
    pso_algo_error        // error case
  };

} // namespace hig

#endif // __HIPGISAXS_FIT_PSO_TYPEDEFS_HPP__
