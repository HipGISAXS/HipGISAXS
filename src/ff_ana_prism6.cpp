/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_prism6.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 05:11:14 PM PST
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
	 * six faceted prism
	 */
	bool AnalyticFormFactor::compute_prism6(shape_param_list_t& params, std::vector<complex_t>& ff,
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

/*            R = dims(1);
            H = dims(2);
            sq3 = sqrt(3);
            qm = tan(tau) *( qx *sin(eta) + qy *cos(eta));
            FF = 4*sq3./( 3*qy.^2 - qx.^2 ) .* (R^2 * qy.^2 .*  SINC_Matrix(qx*R/sq3) .* SINC_Matrix(qy*R) + cos(2*qx*R/sq3) - cos(qy*R) .* cos(qx*R/sq3) )  .* Fq_Inv_Matrix(qz+qm, H) .* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3)))  ; */

		float_t sqrt3 = sqrt(3.0);

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
					complex_t qm = tan(tau) * (mqx * sin(eta) + mqy * cos(eta));
					complex_t temp1 = ((float_t)4.0 * sqrt3) / (3.0 * mqy * mqy - mqx * mqx);
					complex_t temp_ff(0, 0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							complex_t rmqx = r[i_r] * mqx / sqrt3;
							complex_t rmqy = r[i_r] * mqy;
							complex_t temp2 = rmqy * rmqy *	sinc(rmqx) * sinc(rmqy);
							complex_t temp3 = cos(2.0 * rmqx);
							complex_t temp4 = cos(rmqy) * cos(rmqx);
							complex_t temp5 = temp1 * (temp2 + temp3 - temp4);
							complex_t temp6 = fq_inv(mqz + qm, h[i_h]);
							temp_ff += temp5 * temp6;
						} // for h
					} // for r
					complex_t temp7 = (mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
					ff[index] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_prism6()

} // namespace hig
