/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_prism3x.cpp
 *  Created: Jul 12, 2012
 *  Modified: Sat 28 Sep 2013 11:07:09 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <boost/math/special_functions/fpclassify.hpp>

#include "../woo/timer/woo_boostchronotimers.hpp"

#include "ff_ana.hpp"
#include "../common/enums.hpp"
#include "../model/shape.hpp"
#include "../model/qgrid.hpp"
#include "../utils/utilities.hpp"
#include "../numerics/numeric_utils.hpp"

namespace hig {

	/**
	 * triangular grating in the x-direction
	 */
	bool AnalyticFormFactor::compute_prism3x(shape_param_list_t& params, std::vector<complex_t>& ff,
												float_t tau, float_t eta, vector3_t transvec) {
		std::vector <float_t> lx, distr_lx;
		std::vector <float_t> ly, distr_ly;
		std::vector <float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_xsize:
					param_distribution((*i).second, lx, distr_lx);
					break;
				case param_ysize:
					param_distribution((*i).second, ly, distr_ly);
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_edge:
				case param_radius:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted values for shape type 'prism3x'" << std::endl;
					break;
				default:
					std::cerr << "warning: unknown parameters for shape type 'prism3x'. ignoring"
								<< std::endl;
			} // switch
		} // for

		// on cpu
		std::cout << "-- Computing prism3x FF on CPU ..." << std::endl;
		// initialize ff
		ff.clear();  ff.reserve(nqz_ * nqy_ * nqx_);
		for(unsigned int i = 0; i < nqz_ * nqy_ * nqx_; ++ i) ff.push_back(complex_t(0, 0));

		float_t d = 0.85;	// FIXME: hardcoded? variable?
		float_t gamma = 0.0;	// FIXME: hardcoded? variable?
		complex_t i(0.0, 1.0);

		#pragma omp parallel for collapse(3)
		for(unsigned int j_z = 0; j_z < nqz_; ++ j_z) {
			for(unsigned int j_y = 0; j_y < nqy_; ++ j_y) {
				for(unsigned int j_x = 0; j_x < nqx_; ++ j_x) {
					float_t temp_qx = QGrid::instance().qx(j_x);
					float_t temp_qy = QGrid::instance().qy(j_y);
					complex_t temp_qz = QGrid::instance().qz_extended(j_z);
					complex_t mqx, mqy, mqz;
					compute_meshpoints(temp_qx, temp_qy, temp_qz, rot_, mqx, mqy, mqz);
					float_t sg = sin(gamma);
					float_t cg = cos(gamma);
					float_t qx_rot = temp_qx * cg + temp_qy * sg;
					float_t qy_rot = temp_qy * cg - temp_qx * sg;
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {				// H
						for(unsigned int i_y = 0; i_y < ly.size(); ++ i_y) {		// L
							for(unsigned int i_x = 0; i_x < lx.size(); ++ i_x) {	// Lx
								float_t temp_lx = lx[i_x] * 2, temp_ly = ly[i_y] * 2;// multiply by 2 (why?)
								float_t a1 = h[i_h] / (d * temp_ly);
								float_t b1 = 0.0;
								float_t a2 = h[i_h] / ((d - 1) * temp_ly);
								float_t b2 = h[i_h] / (1 - d);
								float_t temp1 = qx_rot * temp_lx / 2;
								float_t fqx = temp_lx;
								if(boost::math::fpclassify(qx_rot) != FP_ZERO) {
									fqx *= sin(temp1) / temp1;
								} // if
								complex_t k1 = temp_qz * a1 + qy_rot;
								complex_t k2 = temp_qz * a2 + qy_rot;
								complex_t i1 = exp(i * temp_qz * b1) * integral_e(0, d * temp_ly, k1);
								complex_t i2 = exp(i * temp_qz * b2) * integral_e(d * temp_ly, temp_ly, k2);
								complex_t i3 = integral_e(0, temp_ly, qy_rot);
								complex_t iy;
								if(boost::math::fpclassify(temp_qz.real()) == FP_ZERO &&
										boost::math::fpclassify(temp_qz.imag()) == FP_ZERO) {
									if(boost::math::fpclassify(qy_rot) == FP_ZERO) {
										iy = h[i_h] * temp_ly / 2;
									} else {
										iy = integral_xe(0, d * temp_ly, a1, b1, qy_rot) +
												integral_xe(d * temp_ly, temp_ly, a2, b2, qy_rot);
									} // if-else
								} else {
									iy = (- i / temp_qz) * (i1 + i2 + i3);
								} // if-else
								temp_ff += fqx * iy;
							} // for i_x
						} // for i_y
					} // for i_h
					complex_t temp7 = (mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
					unsigned int curr_index = nqx_ * nqy_ * j_z + nqx_ * j_y + j_x;
					ff[curr_index] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_prism3x()

} // namespace hig
