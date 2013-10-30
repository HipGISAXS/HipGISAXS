/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_hcylinder.cpp
 *  Created: Jul 12, 2012
 *  Modified: Sat 28 Sep 2013 11:07:05 AM PDT
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
#include "../common/constants.hpp"
#include "../common/enums.hpp"
#include "../model/shape.hpp"
#include "../model/qgrid.hpp"
#include "../utils/utilities.hpp"
#include "../numerics/numeric_utils.hpp"

namespace hig {

	/**
	 * horizontal cylinder
	 */
	bool AnalyticFormFactor::compute_horizontal_cylinder(float_t tau, float_t eta, shape_param_list_t& params,
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

#ifdef TIME_DETAIL_2
		woo::BoostChronoTimer maintimer;
		maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
		// on gpu
		std::cout << "-- Computing hcylinder FF on GPU ..." << std::endl;

		std::vector<float_t> transvec_v;
		transvec_v.push_back(transvec[0]);
		transvec_v.push_back(transvec[1]);
		transvec_v.push_back(transvec[2]);

		gff_.compute_horizontal_cylinder(tau, eta, h, distr_h, r, distr_r, rot_, transvec_v, ff);
#else
		// on cpu
		std::cout << "-- Computing hcylinder FF on CPU ..." << std::endl;

		complex_t unitc(0, 1);

		ff.clear(); ff.reserve(nqx_ * nqy_ * nqy_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

		#pragma omp parallel for collapse(3)
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					complex_t mqx, mqy, mqz;
					compute_meshpoints(QGrid::instance().qx(x), QGrid::instance().qy(y),
										QGrid::instance().qz_extended(z), rot_, mqx, mqy, mqz);
					complex_t qpar = sqrt(mqz * mqz + mqy * mqy);
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							temp_ff += distr_r[i_r] * distr_h[i_h] * 2 * PI_ * r[i_r] * r[i_r] *
										(cbessj(qpar * r[i_r], 1) / (qpar * r[i_r])) *
										fq_inv(mqx, h[i_h]);
						} // for h
					} // for r
					complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					ff[index] = temp_ff * exp(complex_t(-temp1.imag(), temp1.real()));
				} // for x
			} // for y
		} // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
		maintimer.stop();
		std::cout << "**     HCylinder FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2
		return true;
	} // AnalyticFormFactor::compute_horizontal_cylinder()

} // namespace hig
