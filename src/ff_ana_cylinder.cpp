/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_cylinder.cpp
  *  Created: Jul 12, 2012
  *  Modified: Tue 19 Feb 2013 04:07:33 PM PST
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

#ifdef FF_ANA_GPU
		// on gpu
		std::cout << "-- Computing cylinder FF on GPU ..." << std::endl;

		std::vector<float_t> transvec_v;
		transvec_v.push_back(transvec[0]);
		transvec_v.push_back(transvec[1]);
		transvec_v.push_back(transvec[2]);
		// copy to gpu memory
		float_t *qx_h = new (std::nothrow) float_t[nqx_];
		float_t *qy_h = new (std::nothrow) float_t[nqy_];
		cucomplex_t *qz_h = new (std::nothrow) cucomplex_t[nqz_];
		if(qx_h == NULL || qy_h == NULL || qz_h == NULL) {
			std::cerr << "error: memory allocation for host mesh grid failed" << std::endl;
			return false;
		} // if
		for(unsigned int ix = 0; ix < nqx_; ++ ix) qx_h[ix] = QGrid::instance().qx(ix);
		for(unsigned int iy = 0; iy < nqy_; ++ iy) qy_h[iy] = QGrid::instance().qy(iy);
		for(unsigned int iz = 0; iz < nqz_; ++ iz) {
			qz_h[iz].x = QGrid::instance().qz_extended(iz).real();
			qz_h[iz].y = QGrid::instance().qz_extended(iz).imag();
		} // for qz

		ff_gpu_.compute_cylinder(tau, eta, h, distr_h, r, distr_r, qx_h, qy_h, qz_h, rot_, transvec_v, ff);
#else
		// on cpu
		std::cout << "-- Computing cylinder FF on CPU ..." << std::endl;

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqz_);
		for(int i_z = 0; i_z < nqz_; ++ i_z) {
			for(int i_y = 0; i_y < nqy_; ++ i_y) {
				for(int i_x = 0; i_x < nqx_; ++ i_x) {
					ff.push_back(complex_t(0.0, 0.0));
				} // for x
			} // for y
		} // for z

		for(unsigned z = 0; z < nqz_; ++ z) {
			for(unsigned y = 0; y < nqy_; ++ y) {
				for(unsigned x = 0; x < nqx_; ++ x) {
					complex_t mqx = QGrid::instance().qy(y) * rot_[0] +
									QGrid::instance().qx(x) * rot_[1] +
									QGrid::instance().qz_extended(z) * rot_[2];
					complex_t mqy = QGrid::instance().qy(y) * rot_[3] +
									QGrid::instance().qx(x) * rot_[4] +
									QGrid::instance().qz_extended(z) * rot_[5];
					complex_t mqz = QGrid::instance().qy(y) * rot_[6] +
									QGrid::instance().qx(x) * rot_[7] +
									QGrid::instance().qz_extended(z) * rot_[8];
					complex_t temp1 = mqx * mqx;
					complex_t temp2 = mqy * mqy;
					complex_t qpar = sqrt(temp1 + temp2);
					temp1 = mqx * sin(eta);
					temp2 = mqy * cos(eta);
					complex_t qm = (temp1 + temp2) * tan(tau);
					complex_t temp_fq = mqz + qm;
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							complex_t temp1 = fq_inv(temp_fq, h[i_h]);
							complex_t temp2 = qpar * r[i_r];
							complex_t temp3 = cbessj(temp2, 1);
							complex_t temp4 = (temp3 / temp2) * temp1;
							temp_ff += distr_r[i_r] * distr_h[i_h] * 2 * PI_ * r[i_r] * r[i_r] * temp4;
						} // for h
					} // for r
					temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
					unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
					ff[curr_index] = temp_ff * temp2;
				} // for x
			} // for y
		} // for z
#endif // FF_ANA_GPU

		return true;
	} // AnalyticFormFactor::compute_cylinder()

} // namespace hig
