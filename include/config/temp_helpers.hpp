/**
 *  Project:
 *
 *  File: temp_helpers.hpp
 *  Created: Jan 29, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef TEMP_HELPERS__H
#define TEMP_HELPERS__H

#include <string>
#include <vector>
#include <map>

#include <config/tokens.hpp>
#include <config/token_mapper.hpp>
#include <common/typedefs.hpp>
namespace hig {

  class FitReferenceData {
    private:
      std::string image_path_;
      std::string image_mask_;
      unsigned int np_parallel_;
      unsigned int np_perpendicular_;
      OutputRegionType region_type_;
      real_t x_min_, y_min_;
      real_t x_max_, y_max_;

    public:
      FitReferenceData() { }
      FitReferenceData(std::string p) : image_path_(p) { }
      ~FitReferenceData() { }

      bool path(std::string p) { image_path_ = p; return true; }
      bool mask(std::string m) { image_mask_ = m; return true; }
      bool region_min(real_t x, real_t y) { x_min_ = x; y_min_ = y; return true; }
      bool region_max(real_t x, real_t y) { x_max_ = x; y_max_ = y; return true; }
      bool npoints_parallel(unsigned int n) { np_parallel_ = n; return true; }
      bool npoints_perpendicular(unsigned int n) { np_perpendicular_ = n; return true; }
      bool region_type(OutputRegionType r) { region_type_ = r; return true; }

      std::string image_path() const { return image_path_; }
      std::string image_mask() const { return image_mask_; }
      OutputRegionType get_region_type() const { return region_type_; }
      real_t region_min_x() const { return x_min_; }
      real_t region_min_y() const { return y_min_; }
      real_t region_max_x() const { return x_max_; }
      real_t region_max_y() const { return y_max_; }

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
      real_t value_;
      FitAlgorithmParamType type_;
      std::string type_name_;

    public:
      AnalysisAlgorithmParamData() { }
      ~AnalysisAlgorithmParamData() { }

      bool init() { return clear(); }
      bool clear() { value_ = 0.0; type_ = algo_param_null; type_name_ = ""; return true; }
      bool value(real_t v) { value_ = v; return true; }
      bool type(FitAlgorithmParamType t) { type_ = t; return true; }
      bool type_name(std::string t) { type_name_ = t; return true; }

      FitAlgorithmParamType type() const { return type_; }
      real_t value() const { return value_; }

      void print() const {
        std::cout << "  " << type_name_ << " [" << type_ << "] = " << value_ << std::endl;
      } // print()
  }; // class AnalysisAlgorithmParam

  typedef std::map <FitAlgorithmParamType, AnalysisAlgorithmParamData> analysis_algo_param_map_t;


  class AnalysisAlgorithmData {

    private:

      analysis_algo_param_map_t params_map_;
      int order_;
      real_t tolerance_;
      real_t regularization_;
      FittingAlgorithmName name_;
      std::string name_str_;
      FittingDistanceMetric distance_metric_;
      bool restart_;

    public:

      AnalysisAlgorithmData() { regularization_ = 0.0; }
      ~AnalysisAlgorithmData() { }

      bool init() { return clear(); }

      bool clear() {
        params_map_.clear();
        return true;
      } // clear()

      bool add_param(const AnalysisAlgorithmParamData& p) {
        params_map_[p.type()] = p;
        return true;
      } // add_param()

      bool order(int o) { order_ = o; return true; }
      bool tolerance(real_t t) { tolerance_ = t; return true; }
      bool regularization(real_t r) { regularization_ = r; return true; }
      bool name(FittingAlgorithmName n) { name_ = n; return true; }
      bool name_str(std::string n) { name_str_ = n; return true; }
      bool distance_metric(FittingDistanceMetric m) { distance_metric_ = m; return true; }
      bool restart(bool r) { restart_ = r; return true; }

      FittingAlgorithmName name() const { return name_; }
      real_t tolerance() const { return tolerance_; }
      real_t regularization() const { return regularization_; }
      FittingDistanceMetric distance_metric() const { return distance_metric_; }

      bool param(const std::string pstr, real_t& val) const {
        FitAlgorithmParamType type = TokenMapper::instance().get_fit_algorithm_param_token(pstr);
        if(type == algo_param_error) {
          std::cerr << "warning: non-existent analysis algorithm parameter type encountered ["
                << pstr << "]" << std::endl;
          return false;
        } // if
        val = params_map_.at(type).value();
        return true;
      } // param()

      void print() const {
        std::cout << order_ << ": " << name_str_ << " [" << name_ << "]" << std::endl;
        std::cout << "  Tolerance: " << tolerance_ << std::endl;
        std::cout << "  Regularization: " << regularization_ << std::endl;
        std::cout << "  Algorithm Parameters: " << std::endl;
        for(analysis_algo_param_map_t::const_iterator i = params_map_.begin();
            i != params_map_.end(); ++ i) {
          std::cout << "  "; (*i).second.print();
        } // for
        std::cout << "  Distance Metric: " << distance_metric_ << std::endl;
      } // print()
  }; // class AnalysisAlgorithmData

  typedef std::vector <AnalysisAlgorithmData> analysis_algo_list_t;


} // namespace hig

#endif // TEMP_HELPERS__H
