/**
 *  Project:
 *
 *  File: analysis_algorithm.cpp
 *  Created: Feb 02, 2014
 *  Modified: Wed 09 Jul 2014 11:40:24 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <analyzer/analysis_algorithm.hpp>

namespace hig {

	bool AnalysisAlgorithm::init_params(const float_vec_t& x0) {
		x0_.clear();
		x0_.insert(x0_.end(), x0.begin(), x0.end());
		num_params_ = x0_.size();

		return true;
	} // AnalysisAlgorithm::init_params()

} // namespace hig
