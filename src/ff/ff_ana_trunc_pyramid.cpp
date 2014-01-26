/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_trunc_pyramid.cpp
 *  Created: Jul 12, 2012
 *  Modified: Sun 26 Jan 2014 10:25:42 AM PST
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

#include <woo/timer/woo_boostchronotimers.hpp>

#include <ff/ff_ana.hpp>
#include <common/constants.hpp>
#include <common/enums.hpp>
#include <model/shape.hpp>
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#include <numerics/numeric_utils.hpp>

namespace hig {

	/**
	 * truncated pyramid
	 */
	bool AnalyticFormFactor::compute_truncated_pyramid(shape_param_list_t& params,
														std::vector<complex_t>& ff,
														vector3_t transvec) {
		std::vector<float_t> x, distr_x;
		std::vector<float_t> y, distr_y;
		std::vector<float_t> h, distr_h;
		std::vector<float_t> b, distr_b;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_xsize:
					param_distribution((*i).second, x, distr_x);
					break;
				case param_ysize:
					param_distribution((*i).second, y, distr_y);
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_baseangle:
					param_distribution((*i).second, b, distr_b);
					break;
				case param_edge:
				case param_radius:
					std::cerr << "warning: ignoring unwanted parameters in truncated pyramid" << std::endl;
					break;
				default:
					std::cerr << "error: unknown/invalid parameter given for truncated pyramid" << std::endl;
					return false;
			} // switch
		} // for

		float_t tr = 0.4;	// FIXME: hardcoded ... ???

		#ifdef TIME_DETAIL_2
			woo::BoostChronoTimer maintimer;
			maintimer.start();
		#endif
		//#ifdef FF_ANA_GPU
		//	std::cout << "-- Computing truncated pyramid FF on GPU ..." << std::endl;
		//	// ...
		//#else
			std::cout << "-- Computing truncated pyramid FF on CPU ..." << std::endl;
			ff.clear(); ff.reserve(nqx_ * nqy_ * nqz_);
			for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0, 0));

			#pragma omp parallel for collapse(3)
			for(int zi = 0; zi < nqz_; ++ zi) {
				for(int yi = 0; yi < nqy_; ++ yi) {
					for(int xi = 0; xi < nqx_; ++ xi) {
						complex_t mqx, mqy, mqz;
						compute_meshpoints(QGrid::instance().qx(xi), QGrid::instance().qy(yi),
											QGrid::instance().qz_extended(zi), rot_, mqx, mqy, mqz);

						complex_t temp_ff(0.0, 0.0);
						for(int i_x = 0; i_x < x.size(); ++ i_x) {
							for(int i_y = 0; i_y < y.size(); ++ i_y) {
								for(int i_h = 0; i_h < h.size(); ++ i_h) {
									for(int i_b = 0; i_b < b.size(); ++ i_b) {
										float_t xx = x[i_x] * (1 - tr);
										float_t hh = h[i_h] * (1 - tr);
										float_t h2 = h[i_h] * tr / 2;
										float_t yy = y[i_y] * (1 - tr);
										float_t bb = b[i_b] * PI_ / 180;
										float_t prob = distr_x[i_h] * distr_y[i_y] * distr_h[i_h];

										complex_t temp_val(0.0, 0.0);
										temp_val += 4 * yy * xx * sinc(mqx * xx) *
													sinc(mqy * yy) * sinc(mqz * hh / (float_t) 2.0);
										temp_val += truncated_pyramid_core(1, 0,
													mqx, mqy, mqz, xx, yy, h2, bb, vector3_t(0, 0, hh / 2));
										temp_val += truncated_pyramid_core(1, PI_ / 2,
													mqx, mqy, mqz, xx, yy, h2, bb, vector3_t(0, yy, 0));
										temp_val += truncated_pyramid_core(1, PI_,
													mqx, mqy, mqz, xx, yy, h2, bb, vector3_t(0, 0, - hh / 2));
										temp_val += truncated_pyramid_core(1, - PI_ / 2,
													mqx, mqy, mqz, xx, yy, h2, bb, vector3_t(0, - yy, 0));
										temp_val += truncated_pyramid_core(2, PI_ / 2,
													mqx, mqy, mqz, xx, yy, h2, bb, vector3_t(xx, 0, 0));
										temp_val += truncated_pyramid_core(2, - PI_ / 2,
													mqx, mqy, mqz, xx, yy, h2, bb, vector3_t(- xx, 0, 0));

										temp_ff += prob * temp_val;
									} // for b
								} // for h
							} // for y
						} // for x

						complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
						unsigned int index = nqx_ * nqy_ * zi + nqx_ * yi + xi;
						ff[index] = temp_ff * exp(complex_t(-temp1.imag(), temp1.real()));
					} // for x
				} // for y
			} // for z
		//#endif
		#ifdef TIME_DETAIL_2
			maintimer.stop();
			std::cout << "** Trunc Pyramid FF compute time: " << maintimer.elapsed_msec() << " ms."
						<< std::endl;
		#endif

		return true;
	} // AnalyticFormFactor::compute_truncated_pyramid()



	complex_t AnalyticFormFactor::truncated_pyramid_core(int axis, float_t rang,
									complex_t mqx, complex_t mqy, complex_t mqz,
									float_t x, float_t y, float_t h, float_t bang, vector3_t transvec) {
		vector3_t temp_rot1(0.0, 0.0, 0.0), temp_rot2(0.0, 0.0, 0.0), temp_rot3(0.0, 0.0, 0.0);
		compute_rotation_matrix(axis, - rang, temp_rot1, temp_rot2, temp_rot3);
		complex_t rqx = temp_rot1[0] * mqx + temp_rot2[0] * mqy + temp_rot3[0] * mqz;
		complex_t rqy = temp_rot1[1] * mqx + temp_rot2[1] * mqy + temp_rot3[1] * mqz;
		complex_t rqz = temp_rot1[2] * mqx + temp_rot2[2] * mqy + temp_rot3[2] * mqz;
		float_t tanbang = tan(bang);
		complex_t qxy_s = (rqx + rqy) / tanbang;
		complex_t qxy_d = (rqx - rqy) / tanbang;
		complex_t q1 = (rqz + qxy_d) / (float_t) 2.0;
		complex_t q2 = (-rqz + qxy_d) / (float_t) 2.0;
		complex_t q3 = (rqz + qxy_s) / (float_t) 2.0;
		complex_t q4 = (-rqz + qxy_s) / (float_t) 2.0;

		complex_t sinc_eiq1 = sinc(q1 * h) * complex_t(0, 1) * exp(q1 * h);
		complex_t sinc_eiq2 = sinc(q2 * h) * complex_t(0, 1) * exp(- q2 * h);
		complex_t sinc_eiq3 = sinc(q3 * h) * complex_t(0, 1) * exp(q3 * h);
		complex_t sinc_eiq4 = sinc(q4 * h) * complex_t(0, 1) * exp(- q4 * h);

		complex_t k1 = sinc_eiq1 + sinc_eiq2;
		complex_t k2 = - complex_t(0, 1) * (sinc_eiq1 - sinc_eiq2);
		complex_t k3 = sinc_eiq3 + sinc_eiq4;
		complex_t k4 = - complex_t(0, 1) * (sinc_eiq3 - sinc_eiq4);

		if((boost::math::fpclassify(rqx.real()) == FP_ZERO &&
				boost::math::fpclassify(rqx.imag()) == FP_ZERO) ||
				(boost::math::fpclassify(rqy.real()) == FP_ZERO &&
				 boost::math::fpclassify(rqy.imag()) == FP_ZERO)) {
			// FIXME: this is temporary fix !!!!! ...
			return complex_t(0.0, 0.0);
		} // if
		complex_t val = (h / (rqx * rqy)) * (cos(rqx * x - rqy * y) * k1 +
											sin(rqx * x - rqy * y) * k2 -
											cos(rqx * x + rqy * y) * k3 -
											sin(rqx * x + rqy * y) * k4) *
											exp(complex_t(0, 1) *
											(rqx * transvec[0] + rqy * transvec[1] + rqz * transvec[2]));

		return val;
	} // AnalyticFormFactor::truncated_pyramid_core()


	bool AnalyticFormFactor::compute_rotation_matrix(int axis, float_t rang,
													vector3_t& r1, vector3_t& r2, vector3_t& r3) {
		float_t s = sin(rang);
		float_t c = cos(rang);

		switch(axis) {
			case 1:
				r1[0] = 1; r1[1] = 0; r1[2] = 0;
				r2[0] = 0; r2[1] = c; r2[2] = - s;
				r3[0] = 0; r3[1] = s; r3[2] = c;
				break;

			case 2:
				r1[0] = c; r1[1] = 0; r1[2] = - s;
				r2[0] = 0; r2[1] = 1; r2[2] = 0;
				r3[0] = s; r3[1] = 0; r3[2] = c;
				break;

			case 3:
				r1[0] = c; r1[1] = - s; r1[2] = 0;
				r2[0] = s; r2[1] = c; r2[2] = 0;
				r3[0] = 0; r3[1] = 0; r3[2] = 1;
				break;

			default:
				std::cerr << "error: something went really bad in truncated pyramid FF" << std::endl;
				return false;
		} // switch

		return true;
	} // AnalyticFormFactor::compute_rotation_matrix()

} // namespace hig

