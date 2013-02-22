/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_cylinder.cpp
  *  Created: Jul 12, 2012
  *  Modified: Thu 21 Feb 2013 04:39:43 PM PST
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
	 * cylinder
	 */
	bool AnalyticFormFactor::compute_cylinder(shape_param_list_t& params, float_t tau, float_t eta,
			std::vector<complex_t>& ff, vector3_t transvec) {
		std::vector <float_t> h, distr_h;	// for h dimension: param_height
		std::vector <float_t> r, distr_r;	// for r dimension: param_radius
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			if(!(*i).second.isvalid()) {
				std::cerr << "warning: ignoring invalid shape parameter" << std::endl;
				continue;
			} // if
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted input parameters for 'cylinder'" << std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: unknown parameter type given" << std::endl;
					return false;
			} // switch
		} // for

		if(h.size() < 1 || r.size() < 1) {
			std::cerr << "error: missing parameters for cylinder" << std::endl;
			return false;
		} // if

#ifdef TIME_DETAIL_2
		woo::BoostChronoTimer maintimer;
		maintimer.start();
#endif // TIME_DETAIL_2
#ifdef FF_ANA_GPU
		// on gpu
		std::cout << "-- Computing cylinder FF on GPU ..." << std::endl;

		std::vector<float_t> transvec_v;
		transvec_v.push_back(transvec[0]);
		transvec_v.push_back(transvec[1]);
		transvec_v.push_back(transvec[2]);
		gff_.compute_cylinder(tau, eta, h, distr_h, r, distr_r, rot_, transvec_v, ff);
#else
		// on cpu
		std::cout << "-- Computing cylinder FF on CPU ..." << std::endl;

		ff.clear(); ff.reserve(nqx_ * nqy_ * nqz_);
		for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

		#pragma omp parallel for collapse(3)
		for(unsigned z = 0; z < nqz_; ++ z) {
			for(unsigned y = 0; y < nqy_; ++ y) {
				for(unsigned x = 0; x < nqx_; ++ x) {
					complex_t mqx, mqy, mqz;
					compute_meshpoints(QGrid::instance().qx(x), QGrid::instance().qy(y),
										QGrid::instance().qz_extended(z), rot_, mqx, mqy, mqz);
					complex_t qpar = sqrt(mqx * mqx + mqy * mqy);
					complex_t qm = mqz + (mqx * sin(eta) + mqy * cos(eta)) * tan(tau);
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							complex_t temp1 = fq_inv(qm, h[i_h]);
							complex_t temp2 = qpar * r[i_r];
							complex_t temp4 = (cbessj(temp2, 1) / temp2) * temp1;
							temp_ff += distr_r[i_r] * distr_h[i_h] * 2 * PI_ * r[i_r] * r[i_r] * temp4;
						} // for h
					} // for r
					complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
					unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
					ff[curr_index] = temp_ff * temp2;
				} // for x
			} // for y
		} // for z
#endif // FF_ANA_GPU
#ifdef TIME_DETAIL_2
		maintimer.stop();
		std::cout << "**      Cylinder FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
#endif // TIME_DETAIL_2

		return true;
	} // AnalyticFormFactor::compute_cylinder()

} // namespace hig
