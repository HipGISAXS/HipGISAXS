/**
 *  Project:
 *
 *  File: temp_helpers.hpp
 *  Created: Jan 29, 2014
 *  Modified: Wed 29 Jan 2014 03:28:59 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <string>
#include <vector>

namespace hig {

	class FitReferenceData {
		private:
			std::string image_path_;
			unsigned int np_parallel_;
			unsigned int np_perpendicular_;
			OutputRegionType region_type_;
			float_t x_min_, y_min_;
			float_t x_max_, y_max_;

		public:
			FitReferenceData() { }
			FitReferenceData(std::string p) : image_path_(p) { }
			~FitReferenceData() { }

			bool path(std::string p) { image_path_ = p; return true; }
			bool region_min(float_t x, float_t y) { x_min_ = x; y_min_ = y; return true; }
			bool region_max(float_t x, float_t y) { x_max_ = x; y_max_ = y; return true; }
			bool npoints_parallel(unsigned int n) { np_parallel_ = n; return true; }
			bool npoints_perpendicular(unsigned int n) { np_perpendicular_ = n; return true; }
			bool region_type(OutputRegionType r) { region_type_ = r; return true; }

			std::string image_path() const { return image_path_; }
			OutputRegionType get_region_type() const { return region_type_; }

	}; // class FitReferenceData


	class AnalysisAlgorithmParam {
		private:
			float_t value_;
			FitAlgorithmParamType type_;
			std::string type_name_;

		public:
			AnalysisAlgorithmParam() { }
			~AnalysisAlgorithmParam() { }

			bool init() { return clear(); }
			bool clear() { return true; }
			bool value(float_t v) { value_ = v; return true; }
			bool type(FitAlgorithmParamType t) { type_ = t; return true; }
			bool type_name(std::string t) { type_name_ = t; return true; }
	}; // class AnalysisAlgorithmParam


	class AnalysisAlgorithm {
		private:
			std::vector <AnalysisAlgorithmParam> params_;
			int order_;
			float_t tolerance_;
			FittingAlgorithmName name_;
			std::string name_str_;
			bool restart_;

		public:
			AnalysisAlgorithm() { }
			~AnalysisAlgorithm() { }

			bool init() { return clear(); }
			bool clear() { params_.clear(); return true; }
			bool add_param(const AnalysisAlgorithmParam& p) { params_.push_back(p); return true; }
			bool order(int o) { order_ = o; return true; }
			bool tolerance(float_t t) { tolerance_ = t; return true; }
			bool name(FittingAlgorithmName n) { name_ = n; return true; }
			bool name_str(std::string n) { name_str_ = n; return true; }
			bool restart(bool r) { restart_ = r; return true; }
	}; // class AnalysisAlgorithm

	typedef std::vector <AnalysisAlgorithm> analysis_algo_list_t;


} // namespace hig
