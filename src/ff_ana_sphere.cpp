/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_sphere.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 12:53:04 PM PST
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
	 * sphere
	 */
	bool AnalyticFormFactor::compute_sphere(shape_param_list_t& params, std::vector<complex_t> &ff,
											vector3_t transvec) {
		std::vector<float_t> r, distr_r;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: unknown or invalid parameter given for sphere" << std::endl;
					return false;
			} // switch
		} // for

		if(r.size() < 1) {
			std::cerr << "error: radius parameter required for sphere" << std::endl;
			return false;
		} // if

		woo::BoostChronoTimer maintimer;
		maintimer.start();

#ifdef FF_ANA_GPU
		/* on gpu */
		std::cout << "-- Computing sphere FF on GPU ..." << std::endl;

		std::vector<float_t> transvec_v;
		transvec_v.push_back(transvec[0]);
		transvec_v.push_back(transvec[1]);
		transvec_v.push_back(transvec[2]);

		/*float_t *qx_h = new (std::nothrow) float_t[nqx_];
		float_t *qy_h = new (std::nothrow) float_t[nqy_];
		cucomplex_t *qz_h = new (std::nothrow) cucomplex_t[nqz_];
		if(qx_h == NULL || qy_h == NULL || qz_h == NULL) {
			std::cerr << "error: memory allocation for host mesh grid failed" << std::endl;
			return false;
		} // if
		for(unsigned int ix = 0; ix < nqx_; ++ ix) {
			qx_h[ix] = QGrid::instance().qx(ix);
		} // for qx
		for(unsigned int iy = 0; iy < nqy_; ++ iy) {
			qy_h[iy] = QGrid::instance().qy(iy);
		} // for qy
		for(unsigned int iz = 0; iz < nqz_; ++ iz) {
			qz_h[iz].x = QGrid::instance().qz_extended(iz).real();
			qz_h[iz].y = QGrid::instance().qz_extended(iz).imag();
		} // for qz*/

		gff_.compute_sphere(r, distr_r, /*qx_h, qy_h, qz_h,*/ rot_, transvec_v, ff);

		//std::cout << "** GPU analytic sphere computation time: " << gpu_time / 10e6 << " ms." << std::endl;
#else
		/* on cpu */
		std::cout << "-- Computing sphere FF on CPU ..." << std::endl;

		std::vector <complex_t> q;
		q.clear();
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					complex_t temp_qx = QGrid::instance().qy(y) * rot_[0] +
										QGrid::instance().qx(x) * rot_[1] +
										QGrid::instance().qz_extended(z) * rot_[2];
					complex_t temp_qy = QGrid::instance().qy(y) * rot_[3] +
										QGrid::instance().qx(x) * rot_[4] +
										QGrid::instance().qz_extended(z) * rot_[5];
					complex_t temp_qz = QGrid::instance().qy(y) * rot_[6] +
										QGrid::instance().qx(x) * rot_[7] +
										QGrid::instance().qz_extended(z) * rot_[8];
					temp_qx *= temp_qx;
					temp_qy *= temp_qy;
					temp_qz *= temp_qz;
					q.push_back(sqrt(temp_qx + temp_qy + temp_qz));
				} // for
			} // for
		} // for

		ff.clear();
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0, 0));

		std::vector<float_t>::iterator iter_r = r.begin();
		std::vector<float_t>::iterator iter_d = distr_r.begin();
		for(unsigned int i_r = 0; i_r < r.size(); ++ i_r, ++ iter_d) {
			for(unsigned int z = 0; z < nqz_; ++ z) {
				for(unsigned int y = 0; y < nqy_; ++ y) {
					for(unsigned int x = 0; x < nqx_; ++ x) {
						unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
						complex_t temp1 = q[index] * r[i_r];
						complex_t temp2 = sin(temp1) - temp1 * cos(temp1);
						complex_t temp3 = temp1 * temp1 * temp1;
						complex_t temp_mesh_qz = QGrid::instance().qy(y) * rot_[6] +
													QGrid::instance().qx(x) * rot_[7] +
													QGrid::instance().qz_extended(z) * rot_[8];
						ff[index] += distr_r[i_r] * 4 * PI_ * pow(r[i_r], 3) * (temp2 / temp3) *
										exp(complex_t(0, 1) * temp_mesh_qz * r[i_r]);
					} // for x
				} // for y
			} // for z
		} // for r

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
					ff[index] *= exp(complex_t(0, 1) * (mqx * transvec[0] +
									mqy * transvec[1] + mqz * transvec[2]));
				} // for x
			} // for y
		} // for z

		//boost::timer::cpu_times const elapsed_time(timer.elapsed());
		//boost::timer::nanosecond_type const elapsed(elapsed_time.system + elapsed_time.user);
		//double cpu_time = elapsed;
		//std::cout << "** CPU analytic sphere computation time: " << cpu_time / 10e6 << " ms." << std::endl;
		
#endif // FF_ANA_GPU

		maintimer.stop();
		std::cout << "**      Analytic FF compute time: " << maintimer.elapsed_msec() << " ms." << std::endl;
		
		return true;
	} // AnalyticFormFactor::compute_sphere()


} // namespace hig

