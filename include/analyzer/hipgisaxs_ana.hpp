/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Analyzer.hpp
 *  Created: Dec 26, 2013
 *  Modified: Fri 11 Jul 2014 09:03:25 AM PDT
 *  Description: The main analysis class that executes the workflows defined therein wrt
 *  the inputs (a priori structural info)  and datasets (expt. data) provided.
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
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

#ifndef __HIPGISAXS_ANA_HPP__
#define __HIPGISAXS_ANA_HPP__

#include <vector>

#include <common/typedefs.hpp>
#include <analyzer/analysis_algorithm.hpp>


namespace hig {

	class HipGISAXSAnalyzer {

		private :
			std::vector <AnalysisAlgorithm*> wf_;		// the workflow

		public:
			HipGISAXSAnalyzer() { }
			~HipGISAXSAnalyzer() { }

			bool add_analysis_algo(AnalysisAlgorithm* algo) { wf_.push_back(algo); return true; }
			bool analyze(int argc,char **argv, int);

	}; // class HipGISAXSAnalyzer

} // namespace hig

#endif // __HIPGISAXS_ANA_HPP__
