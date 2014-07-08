/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_fit_bruteforce.hpp
 *  Created: Feb 07, 2014
 *  Modified: Fri 07 Feb 2014 10:03:27 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __HIPGISAXS_FIT_BRUTEFORCE_HPP__
#define __HIPGISAXS_FIT_BRUTEFORCE_HPP__

#include <analyzer/analysis_algorithm.hpp>

namespace hig{

	class BruteForceOptimization : public AnalysisAlgorithm {
		private:
			std::vector<std::string> params_;	// list of parameter keys
			float_vec_t x_min_;		// range min for all parameters
			float_vec_t x_max_;		// range max for all parameters
			float_vec_t x_step_;	// stepping delta for all parameters. defaults to 1

			std::vector <std::pair <float_vec_t, float_t> > error_list_;
									// list of errors for corresponding parameter values

			void loop_over_param(int, float_vec_t&);
			bool save_history(char*);

		public:
			BruteForceOptimization(int, char**, ObjectiveFunction*);
			~BruteForceOptimization();

			bool run(int argc,char **argv, int);
	}; // class BruteForceOptimization

} // namespace hig

#endif // __HIPGISAXS_FIT_BRUTEFORCE__
