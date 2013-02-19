/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_prism.cpp
  *  Created: Jul 12, 2012
  *  Modified: Tue 19 Feb 2013 11:42:33 AM PST
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
	 * prism - 3 face
	 */
	bool AnalyticFormFactor::compute_prism(shape_param_list_t& params, std::vector<complex_t>& ff,
											float_t tau, float_t eta, vector3_t transvec) {
		std::vector<float_t> r, distr_r;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: invalid parameter given for prism shape" << std::endl;
					return false;
			} // switch
		} // for

		if(h.size() < 1 || r.size() < 1) {
			std::cerr << "error: need radius and height for prism shape" << std::endl;
			return false;
		} // if

		// why not doing for a range of r, h for this? ... ???

		float_t sqrt3 = sqrt(3.0);
		complex_t unit(0, 1.0);

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqz_);
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
					complex_t temp1 = mqx * (mqx * mqx - 3.0 * mqy * mqy);
					complex_t temp2 = mqz + tan(tau) * (mqx * sin(eta) + mqy * cos(eta));
					complex_t temp_ff(0, 0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							complex_t temp3 = mqy * r[i_r] / sqrt3;
							complex_t temp4 = mqx * exp(unit * mqy * r[i_r] * sqrt3);
							complex_t temp5 = mqx * cos(mqx * r[i_r]);
							complex_t temp6 = unit * sqrt3 * mqy * sin(mqx * r[i_r]);
							complex_t temp7 = fq_inv(temp2, h[i_h]);
							complex_t temp8 = (temp4 - temp5 - temp6) * temp7;
							complex_t temp9 = 2.0 * sqrt3 * exp(complex_t(-temp3.imag(), temp3.real()));
							temp_ff += temp9 / temp1 * temp8;
						} // for h
					} // for r
					complex_t temp10 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					ff[index] = temp_ff * exp(complex_t(-temp10.imag(), temp10.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_prism()

} // namespace hig
