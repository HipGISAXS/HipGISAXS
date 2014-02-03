/**
 *  Project:
 *
 *  File: typedefs.hpp
 *  Created: Jan 14, 2014
 *  Modified: Tue 14 Jan 2014 10:23:37 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __HIPGISAXS_FIT_PSO_TYPEDEFS_HPP__
#define __HIPGISAXS_FIT_PSO_TYPEDEFS_HPP__

#include <vector>
#include <map>
#include <string>
#include <common/typedefs.hpp>	// from hipgisaxs library

namespace hig {

	typedef std::map <std::string, float_t>		parameter_map_t;
	typedef std::vector <float_t>				parameter_data_list_t;
	typedef std::vector <std::string>			parameter_name_list_t;

	enum PSOParameterDistribution {
		PSO_UNIFORM,									// randomized with uniform distribution
		PSO_GAUSSIAN,									// randomized with gaussian distribution
		PSO_SINGLE,										// values set to given values
		PSO_DEFAULT										// currently same as PSO_UNIFORM
	}; // enum PSOParameterDistribution
	typedef PSOParameterDistribution pso_parameter_dist_t;

} // namespace hig

#endif // __HIPGISAXS_FIT_PSO_TYPEDEFS_HPP__
