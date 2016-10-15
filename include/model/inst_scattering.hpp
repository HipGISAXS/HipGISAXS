/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: inst_scattering.hpp
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

#ifndef __INST_SCATTERING_HPP__
#define __INST_SCATTERING_HPP__

#include <iostream>
#include <string>

#include <common/globals.hpp>

namespace hig {

  class ScatteringParams {
    private:
      typedef struct StepRange_ {
        real_t min_;
        real_t max_;
        real_t step_;
      } StepRange;

      std::string expt_;    /* the kind of experiment: saxs, gisaxs, waxs, etc. */
      StepRange alpha_i_;
      StepRange inplane_rot_;
      StepRange tilt_;
      struct Photon {
        real_t value_;
        std::string unit_;
      } photon_;
      std::string polarization_;  // takes values "s", "p", and "sp"
      real_t coherence_;
      real_t spot_area_;
      //vector3_t smearing_;
      real_t smearing_;

    public:
      ScatteringParams();
      ~ScatteringParams();

      void init();
      void clear();

      void expt(std::string s) { expt_ = s; }
      void experiment(std::string s) { expt_ = s; }

      void coherence(real_t d) { coherence_ = d; }
      void spot_area(real_t d) { spot_area_ = d; }

      //void smearing(vector3_t v) { smearing_ = v; }
      //void smearing(real_t v, real_t w, real_t x) {
      //  smearing_[0] = v; smearing_[1] = w; smearing_[2] = x; }
      void smearing(real_t s) { smearing_ = s; }

      void alphai_min(real_t d) { alpha_i_.min_ = d; }
      void alphai_max(real_t d) { alpha_i_.max_ = d; }
      void alphai_step(real_t d) { alpha_i_.step_ = d; }

      void inplane_rot_min(real_t d) { inplane_rot_.min_ = d; }
      void inplane_rot_max(real_t d) { inplane_rot_.max_ = d; }
      void inplane_rot_step(real_t d) { inplane_rot_.step_ = d; }

      void photon_value(real_t d) { photon_.value_ = d; }
      void photon_unit(std::string s) { photon_.unit_ = s; }
      void polarization(std::string s) { polarization_ = s; }

      void tilt_min(real_t d) { tilt_.min_ = d; }
      void tilt_max(real_t d) { tilt_.max_ = d; }
      void tilt_step(real_t d) { tilt_.step_ = d; }

      // getters
      real_t smearing() const { return smearing_; }

      void tilt(real_t & vmin, real_t & vmax, real_t & vstep) const {
        vmin = tilt_.min_; vmax = tilt_.max_; vstep= tilt_.step_;
      }

      real_t spot_area() const { return spot_area_; }
      Photon photon_energy() const { return photon_; }
      real_t energy() const { return photon_.value_; }
      std::string unit() const { return photon_.unit_; }

      real_t alphai_min() const { return alpha_i_.min_; }
      void alphai(real_t & ai_min, real_t & ai_max, real_t & ai_step) const {
        ai_min = alpha_i_.min_;
        ai_max = alpha_i_.max_;
        ai_step = alpha_i_.step_;
      } 
      void inplanerot(real_t & phi_min, real_t & phi_max, real_t & phi_step) const {
        phi_min = inplane_rot_.min_;
        phi_max = inplane_rot_.max_;
        phi_step= inplane_rot_.step_;
      }

      // return experiment string
      std::string experiment() const { return expt_; }
      bool update_param(const std::string&, real_t);


      void print() {
        std::cout << " expt_ = " << expt_ << std::endl
              << " alpha_i_: min_ = " << alpha_i_.min_ << ", max_ = "
              << alpha_i_.max_ << ", step = " << alpha_i_.step_ << std::endl
              << " inplane_rot_: min_ = " << inplane_rot_.min_ << ", max_ = "
              << inplane_rot_.max_ << ", step_ = " << inplane_rot_.step_ << std::endl
              << " tilt_: min_ = " << tilt_.min_ << ", max_ = "
              << tilt_.max_ << ", step_ = " << tilt_.step_ << std::endl
              << " photon_: value_ = " << photon_.value_
              << ", unit = " << photon_.unit_ << std::endl
              << " polarization_ = " << polarization_ << std::endl
              << " coherence_ = " << coherence_ << std::endl
              << " spot_area_ = " << spot_area_ << std::endl
              << " smearing_ = " << smearing_ << std::endl
              //<< " smearing_ = [" << smearing_[0] << ", " << smearing_[1] << ", "
              //<< smearing_[2] << "]" << std::endl
              << std::endl;
      } // print()

      friend class HiGInput;

  }; // class Scattering

} // namespace hig

#endif // __INST_SCATTERING_HPP__
