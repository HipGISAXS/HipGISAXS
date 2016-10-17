/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: compute_params.hpp
 *  Created: Jun 05, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef __COMPUTE_PARAMS_HPP__
#define __COMPUTE_PARAMS_HPP__

#include <iostream>
#include <string>

#include <common/globals.hpp>
#include <common/enums.hpp>

namespace hig {

  class ComputeParams {  // just one instance of this ... TODO ...
    friend class HiGInput;

    private:
      std::string pathprefix_;
      std::string runname_;
      std::string method_;  // TODO: ... change to enum - "dwba" ?
      StructCorrelationType correlation_;    /* grain/ensemble correlation type */
      struct OutputRegion {
        OutputRegionType type_;
        vector2_t minpoint_;
        vector2_t maxpoint_;
        OutputRegion() { }
      } output_region_;
      std::vector<int> resolution_;
      std::string palette_;
      unsigned int nslices_;
      std::string timestamp();
      bool saveff_;
      bool savesf_;

    public:
      ComputeParams();

      void init();
      void clear();

      const std::string& pathprefix() const { return pathprefix_; }
      const std::string& runname() const { return runname_; }
      bool saveff() const { return saveff_; }
      bool savesf() const { return savesf_; }
      StructCorrelationType param_structcorrelation() const { return correlation_; }

      /* setters */

      void pathprefix(std::string s) { pathprefix_ = s; }
      void runname(std::string s) { runname_ = s + "_" + timestamp(); }
      void method(std::string s) { method_ = s; }
      void saveff(bool b) { saveff_ = b; }
      void savesf(bool b) { savesf_ = b; }

      void output_region_type(OutputRegionType o) { output_region_.type_ = o; }
      void output_region_minpoint(vector2_t v) { output_region_.minpoint_ = v; }
      void output_region_maxpoint(vector2_t v) { output_region_.maxpoint_ = v; }
      void output_region_minpoint(real_t v, real_t w) {
        output_region_.minpoint_[0] = v; output_region_.minpoint_[1] = w; }
      void output_region_maxpoint(real_t v, real_t w) {
        output_region_.maxpoint_[0] = v; output_region_.maxpoint_[1] = w; }
      void resolution(int a, int b) { resolution_[0] = a; resolution_[1] = b; }

      void palette(std::string p) { palette_ = p; }
      void nslices(real_t d) { nslices_ = (unsigned int) d; }
      void structcorrelation(StructCorrelationType c) { correlation_ = c; }

      /* getters */
      OutputRegion output_region() const { return output_region_; }
      OutputRegionType output_region_type() const { return output_region_.type_; }
      vector2_t output_minpoint() const { return output_region_.minpoint_; }
      vector2_t output_maxpoint() const { return output_region_.maxpoint_; }
      std::vector<int> resolution() const { return resolution_; }
      std::string palette() const { return palette_; }
      int nslices() const { return nslices_; }
      bool save_ff() const { return saveff_; }
      bool save_sf() const { return savesf_; }

      /* modifiers (updates) */
      bool update_param(const std::string&, real_t);

      void print() {
        std::cout << " pathprefix_ = " << pathprefix_ << std::endl
              << " runname_ = " << runname_ << std::endl
              << " method_ = " << method_ << std::endl
              << " output_region_: type = " << output_region_.type_
              << ", minpoint = [" << output_region_.minpoint_[0] << ", "
              << output_region_.minpoint_[1] << "], maxpoint = ["
              << output_region_.maxpoint_[0] << ", "
              << output_region_.maxpoint_[1] << "]" << std::endl
              << " resolution_ = [" << resolution_[0] << ", "
              << resolution_[1] << "]" << std::endl
              << " nslices_ = " << nslices_ << std::endl
              << " palette_ = " << palette_ << std::endl
              << std::endl;
      } // print()

  }; // class ComputeParams

} // namespace hig

#endif /* __COMPUTE_PARAMS_HPP__ */
