/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_prism6.cpp
  *  Created: Jul 12, 2012
  *  Modified: Thu 21 Feb 2013 04:54:06 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/math/special_functions/fpclassify.hpp>

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
		std::vector<float_t> l, distr_l;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
					param_distribution((*i).second, l, distr_l);
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_xsize:
				case param_ysize:
				case param_baseangle:
				case param_radius:
					std::cerr << "warning: ignoring unwanted parameters in prism6" << std::endl;
					break;
				default:
					std::cerr << "error: invalid parameter given for prism6" << std::endl;
					return false;
			} // switch
		} // for

#ifdef TIME_DETAIL_2
		woo::BoostChronoTimer maintimer;
		maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
		// on gpu
		std::cout << "-- Computing prism6 FF on GPU ..." << std::endl;

		std::vector<float_t> transvec_v;
		transvec_v.push_back(transvec[0]);
		transvec_v.push_back(transvec[1]); 
		transvec_v.push_back(transvec[2]);
		gff_.compute_prism6(tau, eta, l, distr_l, h, distr_h, rot_, transvec_v, ff);
#else
		// on cpu
		std::cout << "-- Computing prism6 FF on CPU ..." << std::endl;

		ff.clear(); ff.reserve(nqx_ * nqy_ * nqz_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

		float_t sqrt3 = sqrt(3.0);

		#pragma omp parallel for collapse(3)
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					complex_t mqx, mqy, mqz;
					compute_meshpoints(QGrid::instance().qx(x), QGrid::instance().qy(y),
										QGrid::instance().qz_extended(z), rot_, mqx, mqy, mqz);
					complex_t qm = tan(tau) * (mqx * sin(eta) + mqy * cos(eta));
					complex_t temp1 = ((float_t) 4.0 * sqrt3) / (3.0 * mqy * mqy - mqx * mqx);
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
						for(unsigned int i_l = 0; i_l < l.size(); ++ i_l) {
							complex_t rmqx = l[i_l] * mqx / sqrt3;
							complex_t rmqy = l[i_l] * mqy;
							complex_t temp2 = rmqy * rmqy *	sinc(rmqx) * sinc(rmqy);
							complex_t temp3 = cos(2.0 * rmqx);
							complex_t temp4 = cos(rmqy) * cos(rmqx);
							complex_t temp5 = temp1 * (temp2 + temp3 - temp4);
							complex_t temp6 = fq_inv(mqz + qm, h[i_h]);
							temp_ff += temp5 * temp6;
						} // for l
					} // for h
					complex_t temp7 = (mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					ff[index] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
				} // for x
			} // for y
		} // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
		maintimer.stop();
		std::cout << "**        Prism6 FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2
		return true;
	} // AnalyticFormFactor::compute_prism6()

} // namespace hig
