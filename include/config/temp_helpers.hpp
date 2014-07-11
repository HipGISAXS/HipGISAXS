/**
 *  Project:
 *
 *  File: temp_helpers.hpp
 *  Created: Jan 29, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <string>
#include <vector>
#include <map>

namespace hig {

	class FitReferenceData {
		private:
			std::string image_path_;
			std::string image_mask_;
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
			bool mask(std::string m) { image_mask_ = m; return true; }
			bool region_min(float_t x, float_t y) { x_min_ = x; y_min_ = y; return true; }
			bool region_max(float_t x, float_t y) { x_max_ = x; y_max_ = y; return true; }
			bool npoints_parallel(unsigned int n) { np_parallel_ = n; return true; }
			bool npoints_perpendicular(unsigned int n) { np_perpendicular_ = n; return true; }
			bool region_type(OutputRegionType r) { region_type_ = r; return true; }

			std::string image_path() const { return image_path_; }
			std::string image_mask() const { return image_mask_; }
			OutputRegionType get_region_type() const { return region_type_; }

			void print() {
				std::cout << "Reference Data:" << std::endl;
				std::cout << "  Image: " << image_path_ << std::endl;
				std::cout << "  Mask: " << image_mask_ << std::endl;
				std::cout << "  min: [" << x_min_ << " " << y_min_ << "]" << std::endl;
				std::cout << "  max: [" << x_max_ << " " << y_max_ << "]" << std::endl;
				std::cout << "  n: " << np_parallel_ << " x " << np_perpendicular_ << std::endl;
			} // print()

	}; // class FitReferenceData


	class AnalysisAlgorithmParamData {
		private:
			float_t value_;
			FitAlgorithmParamType type_;
			std::string type_name_;

		public:
			AnalysisAlgorithmParamData() { }
			~AnalysisAlgorithmParamData() { }

			bool init() { return clear(); }
			bool clear() { value_ = 0.0; type_ = algo_param_null; type_name_ = ""; return true; }
			bool value(float_t v) { value_ = v; return true; }
			bool type(FitAlgorithmParamType t) { type_ = t; return true; }
			bool type_name(std::string t) { type_name_ = t; return true; }

			FitAlgorithmParamType type() const { return type_; }
			float_t value() const { return value_; }

			void print() const {
				std::cout << "  " << type_name_ << " [" << type_ << "] = " << value_ << std::endl;
			} // print()
	}; // class AnalysisAlgorithmParam

	typedef std::map <FitAlgorithmParamType, AnalysisAlgorithmParamData> analysis_algo_param_map_t;

	class AnalysisAlgorithmData {
		private:
			//std::vector <AnalysisAlgorithmParamData> params_;
			analysis_algo_param_map_t params_map_;
			int order_;
			float_t tolerance_;
			FittingAlgorithmName name_;
			std::string name_str_;
			bool restart_;

		public:
			AnalysisAlgorithmData() { }
			~AnalysisAlgorithmData() { }

			bool init() { return clear(); }

			bool clear() {
				//params_.clear();
				params_map_.clear();
				return true;
			} // clear()

			bool add_param(const AnalysisAlgorithmParamData& p) {
				//params_.push_back(p);
				params_map_[p.type()] = p;
				return true;
			} // add_param()

			bool order(int o) { order_ = o; return true; }
			bool tolerance(float_t t) { tolerance_ = t; return true; }
			bool name(FittingAlgorithmName n) { name_ = n; return true; }
			bool name_str(std::string n) { name_str_ = n; return true; }
			bool restart(bool r) { restart_ = r; return true; }

			FittingAlgorithmName name() const { return name_; }

			float_t param(const std::string pstr) const {
				FitAlgorithmParamType type = TokenMapper::instance().get_fit_algorithm_param_token(pstr);
				if(type == algo_param_error) {
					std::cerr << "error: invalid analysis algorithm parameter type encountered"
								<< std::endl;
					exit(2);
				} // if
				return params_map_.at(type).value();
			} // param()

			void print() const {
				std::cout << order_ << ": " << name_str_ << " [" << name_ << "]" << std::endl;
				std::cout << "  Tolerance: " << tolerance_ << std::endl;
				std::cout << "  Algorithm Parameters: " << std::endl;
				for(analysis_algo_param_map_t::const_iterator i = params_map_.begin();
						i != params_map_.end(); ++ i) {
					std::cout << "  "; (*i).second.print();
				} // for
				//for(std::vector <AnalysisAlgorithmParamData>::const_iterator i = params_.begin();
				//		i != params_.end(); ++ i) {
				//	std::cout << "  ";
				//	(*i).print();
				//} // for
			} // print()
	}; // class AnalysisAlgorithm

	typedef std::vector <AnalysisAlgorithmData> analysis_algo_list_t;


} // namespace hig
