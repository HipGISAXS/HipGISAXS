/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_rand_cylinder.cpp
  *  Created: Jul 12, 2012
  *  Modified: Thu 21 Feb 2013 01:02:52 PM PST
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
	 * random cylinders - for SAXS
	 */
	bool AnalyticFormFactor::compute_random_cylinders(shape_param_list_t& params,
			std::vector<complex_t>& ff,	float_t tau, float_t eta, vector3_t transvec) {
		std::vector<float_t> r, distr_r;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in random cylinders"
								<< std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: unknown or invalid parameter given for random cylinders"
								<< std::endl;
					return false;
			} // switch
		} // for

		if(r.size() < 1 || h.size() < 1) {
			std::cerr << "error: both radius and height parameters are required for random cylinders"
						<< std::endl;
			return false;
		} // if

#ifdef FF_ANA_GPU
		// on gpu
		std::cout << "-- Computing random cylinders FF on GPU ..." << std::endl;

		std::vector<float_t> transvec_v;
		transvec_v.push_back(transvec[0]);
		transvec_v.push_back(transvec[1]);
		transvec_v.push_back(transvec[2]);

		gff_.compute_random_cylinders(tau, eta, h, distr_h, r, distr_r, rot_, transvec_v, ff);
#else
		// on cpu
		std::cout << "-- Computing random cylinder FF on CPU ..." << std::endl;

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqy_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));
		complex_t unitc(0, 1);

		float_t dx = 0.001;		// FIXME: hard-coded ???
		unsigned int nx = (1.0 - dx) / dx + 1;

		for(unsigned int z = 0; z < nqz_; ++ z) {
			complex_t qz = QGrid::instance().qz_extended(z);
			complex_t temp_ff(0.0, 0.0);
			for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
				for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
					float_t x_val = 0.0;
					complex_t temp_ffx(0.0, 0.0);
					for(unsigned int i_x = 0; i_x < nx; ++ i_x, x_val += dx) {
						complex_t temp1 = sinc(qz * h[i_h] * x_val / (float_t) 2.0);
						complex_t temp2 = qz * r[i_r] * sqrt(1 - x_val * x_val);
						temp_ffx += temp1 * cbessj(temp2, 1) / temp2;
					} // for
					temp_ff += 4.0 * distr_r[i_r] * distr_h[i_h] * temp_ffx * temp_ffx;
				} // for h
			} // for r
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					ff[index] = temp_ff;
				} // for x
			} // for y
		} // for z
#endif // FF_ANA_GPU
		return true;

	} // AnalyticFormFactor::compute_random_cylinders()

} // namespace hig
