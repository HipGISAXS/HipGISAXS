/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_hcylinder.cpp
  *  Created: Jul 12, 2012
  *  Modified: Tue 19 Feb 2013 11:41:51 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/math/special_functions/fpclassify.hpp>
//#include <boost/timer/timer.hpp>

#include "woo/timer/woo_boostchronotimers.hpp"

#include "ff_ana.hpp"
#include "shape.hpp"
#include "enums.hpp"
#include "qgrid.hpp"
#include "utilities.hpp"
#include "numeric_utils.hpp"

namespace hig {

	/**
	 * horizontal cylinder
	 */
	bool AnalyticFormFactor::compute_horizontal_cylinder(shape_param_list_t& params,
														vector3_t transvec,
														std::vector<complex_t>& ff) {
		std::vector<float_t> r, distr_r;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in horizontal cylinder"
								<< std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: unknown or invalid parameter given for horizontal cylinder"
								<< std::endl;
					return false;
			} // switch
		} // for

		if(r.size() < 1 || h.size() < 1) {
			std::cerr << "error: both radius and height parameters are required for horizontal cylinder"
						<< std::endl;
			return false;
		} // if

		// in slims code, why not doing range of r and h ???

		complex_t unitc(0, 1);

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqy_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					complex_t mqx = QGrid::instance().qy(y) * rot_[0] +
									QGrid::instance().qx(x) * rot_[1] +
									QGrid::instance().qz_extended(z) * rot_[2];
					complex_t mqy = QGrid::instance().qy(y) * rot_[3] +
									QGrid::instance().qx(x) * rot_[4] +
									QGrid::instance().qz_extended(z) * rot_[5];
					complex_t mqz = QGrid::instance().qy(y) * rot_[6] +
									QGrid::instance().qx(x) * rot_[7] +
									QGrid::instance().qz_extended(z) * rot_[8];
					complex_t temp_qpar = sqrt(mqz * mqz + mqy * mqy);
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
								temp_ff += 2 * PI_ * h[i_h] * r[i_r] * r[i_r] *
											(cbessj(temp_qpar * r[i_r], 1) / (temp_qpar * r[i_r])) *
											exp(complex_t(-(mqz * r[i_r]).imag(), (mqz * r[i_r]).real())) *
											sinc(mqx * h[i_h] / (float_t)2.0);
						} // for h
					} // for r
					complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					ff[index] = temp_ff * exp(complex_t(-temp1.imag(), temp1.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_horizontal_cylinder()

} // namespace hig
