/***
  *  $Id: inst_scattering.hpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: inst_scattering.hpp
  *  Created: Jun 05, 2012
  *  Modified: Mon 01 Oct 2012 11:14:44 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _INST_SCATTERING_HPP_
#define _INST_SCATTERING_HPP_

#include <iostream>
#include <string>

#include "globals.hpp"

namespace hig {

	class ScatteringParams {
		private:
			typedef struct StepRange_ {
				float_t min_;
				float_t max_;
				float_t step_;
			} StepRange;

			std::string expt_;		/* the kind of experiment: saxs, gisaxs, waxs, etc. */
			StepRange alpha_i_;
			StepRange inplane_rot_;
			StepRange tilt_;
			struct Photon {
				float_t value_;
				std::string unit_;
			} photon_;
			std::string polarization_;	// takes values "s", "p", and "sp"
			float_t coherence_;
			float_t spot_area_;
			vector3_t smearing_;

		public:
			ScatteringParams();
			~ScatteringParams();

			void init();
			void clear();

			void expt(std::string s) { expt_ = s; }

			void coherence(float_t d) { coherence_ = d; }
			void spot_area(float_t d) { spot_area_ = d; }

			void smearing(vector3_t v) { smearing_ = v; }
			void smearing(float_t v, float_t w, float_t x) {
				smearing_[0] = v; smearing_[1] = w; smearing_[2] = x; }

			void alphai_min(float_t d) { alpha_i_.min_ = d; }
			void alphai_max(float_t d) { alpha_i_.max_ = d; }
			void alphai_step(float_t d) { alpha_i_.step_ = d; }

			void inplane_rot_min(float_t d) { inplane_rot_.min_ = d; }
			void inplane_rot_max(float_t d) { inplane_rot_.max_ = d; }
			void inplane_rot_step(float_t d) { inplane_rot_.step_ = d; }

			void photon_value(float_t d) { photon_.value_ = d; }
			void photon_unit(std::string s) { photon_.unit_ = s; }
			void polarization(std::string s) { polarization_ = s; }

			void tilt_min(float_t d) { tilt_.min_ = d; }
			void tilt_max(float_t d) { tilt_.max_ = d; }
			void tilt_step(float_t d) { tilt_.step_ = d; }

			float_t spot_area() const { return spot_area_; }
			Photon photon_energy() const { return photon_; }


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
							<< " smearing_ = [" << smearing_[0] << ", " << smearing_[1] << ", "
							<< smearing_[2] << "]" << std::endl
							<< std::endl;
			} // print()

			friend class HiGInput;

	}; // class Scattering

} // namespace hig

#endif /* _INST_SCATTERING_HPP_ */
