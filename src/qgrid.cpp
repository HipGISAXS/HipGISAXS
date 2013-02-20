/***
  *  $Id: qgrid.cpp 47 2012-08-23 21:05:16Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: qgrid.cpp
  *  Created: Jun 17, 2012
  *  Modified: Wed 20 Feb 2013 12:24:59 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "constants.hpp"
#include "qgrid.hpp"
#include "hig_input.hpp"
#include "utilities.hpp"

namespace hig {

	/**
	 * create Q-grid in reciprocal space
	 */
	bool QGrid::create(float_t freq, float_t alpha_i, float_t k0) {
										// x and y are reversed in slim's code ? ...
		vector2_t total_pixels = HiGInput::instance().detector_total_pixels();
		vector2_t min_pixel = HiGInput::instance().param_output_minpoint();
		vector2_t max_pixel = HiGInput::instance().param_output_maxpoint();
		OutputRegionType type = HiGInput::instance().param_output_type();

		// sanitize these values
		if(min_pixel[0] < 0 || min_pixel[0] >= max_pixel[0] || min_pixel[0] > total_pixels[0])
			min_pixel[0] = 0;
		if(min_pixel[1] < 0 || min_pixel[1] >= max_pixel[1] || min_pixel[1] > total_pixels[1])
			min_pixel[1] = 0;

		if(max_pixel[0] < 0 || max_pixel[0] <= min_pixel[0] || max_pixel[0] > total_pixels[0])
			max_pixel[0] = total_pixels[0];
		if(max_pixel[1] < 0 || max_pixel[1] <= min_pixel[1] || max_pixel[1] > total_pixels[1])
			max_pixel[1] = total_pixels[1];

		vector3_t qmax, qmin;
		vector3_t step;

		if(type == region_pixels) {
			vector2_t beam = HiGInput::instance().detector_direct_beam();
			float_t pixel_size = HiGInput::instance().detector_pixel_size();
			float_t sd_distance = HiGInput::instance().detector_sd_distance();
	
			vector3_t temp_qmin, temp_qmax;
			pixel_to_kspace(min_pixel, k0, alpha_i, pixel_size, sd_distance, beam, temp_qmin);
			pixel_to_kspace(max_pixel, k0, alpha_i, pixel_size, sd_distance, beam, temp_qmax);
	
			qmax[0] = max(temp_qmax[0], temp_qmin[0]);
			qmax[1] = max(fabs(temp_qmax[1]), fabs(temp_qmin[1]));
			qmax[2] = max(temp_qmax[2], temp_qmin[2]);
			qmin[0] = min(temp_qmax[0], temp_qmin[0]);
			qmin[1] = - qmax[1];
			qmin[2] = min(temp_qmax[2], temp_qmin[2]);
	
			//std::cout << "min q-point: " << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << std::endl;
			//std::cout << "max q-point: " << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << std::endl;
	
			/* one-pixel range in k-space
			 * in order to resolve a pixel, dq must be <= qpixel and lower to resolve subpixels */
			vector2_t pixel_res = HiGInput::instance().param_resolution();
			//int num_q_pixel = pixel_res[0] * pixel_res[1];	// number of q-points per pixel
	
			vector3_t q0, q1, q_pixel;
			pixel_to_kspace(vector2_t(0, 0), k0, alpha_i, pixel_size, sd_distance, beam, q0);
			pixel_to_kspace(vector2_t(1, 1), k0, alpha_i, pixel_size, sd_distance, beam, q1);
			//q_pixel = fabs(q0 - q1) / num_q_pixel;		// check this ... ??
			// or ... ??
			q_pixel[0] = (qmax[0] - qmin[0]) / 2;			// why ... ?
			q_pixel[1] = fabs(q0[1] - q1[1]) / pixel_res[0];
			q_pixel[2] = fabs(q0[2] - q1[2]) / pixel_res[1];
	
			//std::cout << "one pixel range in Q-grid: " << q_pixel[0] << ", "
			//			<< q_pixel[1] << ", " << q_pixel[2] << std::endl;
	
			// number of variables can be reduced here ... too many redundant ones ...
	
			vector3_t g1, g2, g3;
			if(HiGInput::instance().experiment() == "saxs") {
				g1[0] = g1[1] = g1[2] = 0;
				g2[0] = 0; g2[1] = 1; g2[2] = q_pixel[2];
				g3[0] = g3[1] = 0; g3[2] = qmax[2];
			} else if(HiGInput::instance().experiment() == "gisaxs") {
				g1 = qmin;
				g2 = q_pixel;
				g3 = qmax;
			} else {
				std::cerr << "error: unknown experiment type '" << HiGInput::instance().experiment()
						<< "'" << std::endl;
				return false;
			} // if-else
	
			// check what to do with qx ...
	
			step[0] = g2[0];	// either 0 or (max-min)/2 -> 0 or 2 points :-/ ... ???
			step[1] = g2[1];	// not dividing by pixel_res because already did before ... check ???
			step[2] = g2[2];	// "
	
			qmin[0] = g1[0]; qmin[1] = g1[1]; qmin[2] = g1[2];
			qmax[0] = g3[0]; qmax[1] = g3[1]; qmax[2] = g3[2];
	
			// temporary, for hard-coded q values
			qmin[0] = 0.0; qmax[0] = 0.0;
			/*qmin[1] = -PI_; qmax[1] = PI_;
			qmin[2] = 0.0122718; qmax[2] = PI_;
			step[0] = 0; step[1] = q_pixel[1]; step[2] = q_pixel[2]; */

		} else if(type == region_qspace) {

			vector2_t beam = HiGInput::instance().detector_direct_beam();
			float_t pixel_size = HiGInput::instance().detector_pixel_size();
			float_t sd_distance = HiGInput::instance().detector_sd_distance();
			vector2_t pixel_res = HiGInput::instance().param_resolution();

			qmin[0] = 0.0; qmax[0] = 0.0;
			qmax[1] = max(fabs(max_pixel[0]), fabs(min_pixel[0]));
			qmax[2] = max(max_pixel[1], min_pixel[1]);
			qmin[1] = - qmax[1];
			qmin[2] = min(max_pixel[1], min_pixel[1]);
	
			//std::cout << "min q-point: " << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << std::endl;
			//std::cout << "max q-point: " << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << std::endl;
	
			vector3_t q0, q1;
			pixel_to_kspace(vector2_t(0, 0), k0, alpha_i, pixel_size, sd_distance, beam, q0);
			pixel_to_kspace(vector2_t(1, 1), k0, alpha_i, pixel_size, sd_distance, beam, q1);
			step[0] = fabs(q1[0] - q0[0]) / 2;
			step[1] = fabs(q1[1] - q0[1]) / pixel_res[0];
			step[2] = fabs(q1[2] - q0[2]) / pixel_res[1];
		} else {
			std::cerr << "error: unknown output region type" << std::endl;
			return false;
		} // if-else

		//std::cout << "new min q-point: " << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << std::endl;
		//std::cout << "new max q-point: " << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << std::endl;
		//std::cout << "step: " << step[0] << ", " << step[1] << ", " << step[2] << std::endl;
		std::cout << "**                  Q-grid range: ("
					<< qmin[0] << ", " << qmin[1] << ", " << qmin[2] << ") x ("
					<< qmax[0] << ", " << qmax[1] << ", " << qmax[2] << ")" << std::endl;

		qx_.push_back(qmin[0]);
		qy_.push_back(qmin[1]);
		qz_.push_back(qmin[2]);
		for(float_t val = qmin[0] + step[0]; val < qmax[0]; val += step[0]) qx_.push_back(val);
		for(float_t val = qmin[1] + step[1]; val < qmax[1]; val += step[1]) qy_.push_back(val);
		for(float_t val = qmin[2] + step[2]; val < qmax[2]; val += step[2]) qz_.push_back(val);

		std::cout << "**               NQX x NQY x NQZ: " << qx_.size() << " x " << qy_.size()
					<< " x " << qz_.size() << std::endl;

		return true;
	} // QGrid::create()


	bool QGrid::create_qz_extended(float_t k0, float_t kzi_0, complex_t kzi_q, complex_t dnl_q) {
		cqvec_t qz_temp0, qz_temp1, qz_temp2, qz_temp3;
		qz_extended_.clear();
		for(qvec_iter_t q = qz_.begin(); q != qz_.end(); ++ q) {
			//std::complex<float_t> kz_q_temp = sqrt(pow((*q) + kzi_0, 2) -
			//									pow(k0, 2) * std::complex<float_t>(dnl_q.x, dnl_q.y)) -
			//									std::complex<float_t>(kzi_q.x, kzi_q.y);
			//complex_t kz_q = sqrt(pow((*q) + kzi_0, 2) - pow(k0, 2) * dnl_q) - kzi_q;
			float_t temp0 = (*q) + kzi_0;
			float_t temp1 = temp0 * temp0;
			complex_t temp2 = k0 * k0 * dnl_q;
			complex_t temp3 = sqrt(temp1 - temp2);

			complex_t temp4 = temp3 - kzi_q;
			qz_temp0.push_back(temp4);
			complex_t temp5 = - temp3 - kzi_q;
			qz_temp1.push_back(temp5);
			complex_t temp6 = temp3 + kzi_q;
			qz_temp2.push_back(temp6);
			complex_t temp7 = - temp3 + kzi_q;
			qz_temp3.push_back(temp7);
		} // for
		// can the 4 vectors be concatenated instead of copying ...
		for(cqvec_iter_t q = qz_temp0.begin(); q != qz_temp0.end(); ++ q)
			qz_extended_.push_back(*q);
		for(cqvec_iter_t q = qz_temp1.begin(); q != qz_temp1.end(); ++ q)
			qz_extended_.push_back(*q);
		for(cqvec_iter_t q = qz_temp2.begin(); q != qz_temp2.end(); ++ q)
			qz_extended_.push_back(*q);
		for(cqvec_iter_t q = qz_temp3.begin(); q != qz_temp3.end(); ++ q)
			qz_extended_.push_back(*q);

		return true;
	} // QGrid::create_qz_extended()


	bool QGrid::pixel_to_kspace(vector2_t pixel, float_t k0, float_t alpha_i,
					float_t rho, float_t distance, vector2_t beam,
					vector3_t& qpoint) {
		float_t alpha_f = atan(rho * (pixel[1] - 2 * beam[1] + 1043) / distance) - alpha_i;
		float_t theta_f = atan(rho * (pixel[0] - beam[0]) / distance);

		qpoint[0] = k0 * (cos(theta_f) * cos(alpha_f) - cos(alpha_i));	// x
		qpoint[1] = k0 * cos(alpha_f) * sin(theta_f);					// y
		qpoint[2] = k0 * (sin(alpha_f) + sin(alpha_i));					// z

		return true;
	} // QGrid::pixel_to_kspace()


	bool QGrid::kspace_to_pixel() {
		// needed ... ?
		// do it later ...
		return true;
	} // QGrid::kspace_to_pixel()


	// temporary for steepest descent fitting ...
	bool QGrid::create_z_cut(float_t freq, float_t alpha_i, float_t k0, float_t zcut) {
		vector2_t total_pixels = HiGInput::instance().detector_total_pixels();
		vector2_t min_pixel = HiGInput::instance().param_output_minpoint();
		vector2_t max_pixel = HiGInput::instance().param_output_maxpoint();
		OutputRegionType type = HiGInput::instance().param_output_type();

		// sanitize these values
		if(min_pixel[0] < 0 || min_pixel[0] >= max_pixel[0] || min_pixel[0] > total_pixels[0])
			min_pixel[0] = 0;
		if(min_pixel[1] < 0 || min_pixel[1] >= max_pixel[1] || min_pixel[1] > total_pixels[1])
			min_pixel[1] = 0;

		if(max_pixel[0] < 0 || max_pixel[0] <= min_pixel[0] || max_pixel[0] > total_pixels[0])
			max_pixel[0] = total_pixels[0];
		if(max_pixel[1] < 0 || max_pixel[1] <= min_pixel[1] || max_pixel[1] > total_pixels[1])
			max_pixel[1] = total_pixels[1];

		// create z just for the required zcut and not the whole thing
		min_pixel[1] = zcut;
		max_pixel[1] = zcut;

		vector3_t qmax, qmin;
		vector3_t step;

		if(type == region_pixels) {
			vector2_t beam = HiGInput::instance().detector_direct_beam();
			float_t pixel_size = HiGInput::instance().detector_pixel_size();
			float_t sd_distance = HiGInput::instance().detector_sd_distance();
	
			vector3_t temp_qmin, temp_qmax;
			pixel_to_kspace(min_pixel, k0, alpha_i, pixel_size, sd_distance, beam, temp_qmin);
			pixel_to_kspace(max_pixel, k0, alpha_i, pixel_size, sd_distance, beam, temp_qmax);
	
			qmax[0] = max(temp_qmax[0], temp_qmin[0]);
			qmax[1] = max(fabs(temp_qmax[1]), fabs(temp_qmin[1]));
			qmax[2] = max(temp_qmax[2], temp_qmin[2]);
			qmin[0] = min(temp_qmax[0], temp_qmin[0]);
			qmin[1] = - qmax[1];
			qmin[2] = min(temp_qmax[2], temp_qmin[2]);
	
			//std::cout << "min q-point: " << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << std::endl;
			//std::cout << "max q-point: " << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << std::endl;
	
			/* one-pixel range in k-space
			 * in order to resolve a pixel, dq must be <= qpixel and lower to resolve subpixels */
			vector2_t pixel_res = HiGInput::instance().param_resolution();
			//int num_q_pixel = pixel_res[0] * pixel_res[1];	// number of q-points per pixel
	
			vector3_t q0, q1, q_pixel;
			pixel_to_kspace(vector2_t(0, 0), k0, alpha_i, pixel_size, sd_distance, beam, q0);
			pixel_to_kspace(vector2_t(1, 1), k0, alpha_i, pixel_size, sd_distance, beam, q1);
			//q_pixel = fabs(q0 - q1) / num_q_pixel;		// check this ... ??
			// or ... ??
			q_pixel[0] = (qmax[0] - qmin[0]) / 2;			// why ... ?
			q_pixel[1] = fabs(q0[1] - q1[1]) / pixel_res[0];
			q_pixel[2] = fabs(q0[2] - q1[2]) / pixel_res[1];
	
			//std::cout << "one pixel range in Q-grid: " << q_pixel[0] << ", "
			//			<< q_pixel[1] << ", " << q_pixel[2] << std::endl;
	
			// number of variables can be reduced here ... too many redundant ones ...
	
			vector3_t g1, g2, g3;
			if(HiGInput::instance().experiment() == "saxs") {
				g1[0] = g1[1] = g1[2] = 0;
				g2[0] = 0; g2[1] = 1; g2[2] = q_pixel[2];
				g3[0] = g3[1] = 0; g3[2] = qmax[2];
			} else if(HiGInput::instance().experiment() == "gisaxs") {
				g1 = qmin;
				g2 = q_pixel;
				g3 = qmax;
			} else {
				std::cerr << "error: unknown experiment type '" << HiGInput::instance().experiment()
						<< "'" << std::endl;
				return false;
			} // if-else
	
			// check what to do with qx ...
	
			step[0] = g2[0];	// either 0 or (max-min)/2 -> 0 or 2 points :-/ ... ???
			step[1] = g2[1];	// not dividing by pixel_res because already did before ... check ???
			step[2] = g2[2];	// "
	
			qmin[0] = g1[0]; qmin[1] = g1[1]; qmin[2] = g1[2];
			qmax[0] = g3[0]; qmax[1] = g3[1]; qmax[2] = g3[2];
	
			// temporary, for hard-coded q values
			/*qmin[0] = 0.0; qmax[0] = 0.0;
			qmin[1] = -PI_; qmax[1] = PI_;
			qmin[2] = 0.0122718; qmax[2] = PI_;
			step[0] = 0; step[1] = q_pixel[1]; step[2] = q_pixel[2]; */

		} else if(type == region_qspace) {

			vector2_t beam = HiGInput::instance().detector_direct_beam();
			float_t pixel_size = HiGInput::instance().detector_pixel_size();
			float_t sd_distance = HiGInput::instance().detector_sd_distance();
			vector2_t pixel_res = HiGInput::instance().param_resolution();

			qmin[0] = 0.0; qmax[0] = 0.0;	// x is jsut 1D for now ...
			qmax[1] = max(fabs(max_pixel[0]), fabs(min_pixel[0]));
			qmax[2] = max(max_pixel[1], min_pixel[1]);
			qmin[1] = - qmax[1];
			qmin[2] = min(max_pixel[1], min_pixel[1]);
	
			//std::cout << "min q-point: " << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << std::endl;
			//std::cout << "max q-point: " << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << std::endl;
	
			vector3_t q0, q1;
			pixel_to_kspace(vector2_t(0, 0), k0, alpha_i, pixel_size, sd_distance, beam, q0);
			pixel_to_kspace(vector2_t(1, 1), k0, alpha_i, pixel_size, sd_distance, beam, q1);
			step[0] = fabs(q1[0] - q0[0]) / 2;
			step[1] = fabs(q1[1] - q0[1]) / pixel_res[0];
			step[2] = fabs(q1[2] - q0[2]) / pixel_res[1];
		} else {
			std::cerr << "error: unknown output region type" << std::endl;
			return false;
		} // if-else

		//std::cout << "new min q-point: " << qmin[0] << ", " << qmin[1] << ", " << qmin[2] << std::endl;
		//std::cout << "new max q-point: " << qmax[0] << ", " << qmax[1] << ", " << qmax[2] << std::endl;
		//std::cout << "step: " << step[0] << ", " << step[1] << ", " << step[2] << std::endl;
		std::cout << "**                  Q-grid range: ("
					<< qmin[0] << ", " << qmin[1] << ", " << qmin[2] << ") x ("
					<< qmax[0] << ", " << qmax[1] << ", " << qmax[2] << ")" << std::endl;

		qx_.push_back(qmin[0]);
		qy_.push_back(qmin[1]);
		qz_.push_back(qmin[2]);
		for(float_t val = qmin[0] + step[0]; val < qmax[0]; val += step[0]) qx_.push_back(val);
		for(float_t val = qmin[1] + step[1]; val < qmax[1]; val += step[1]) qy_.push_back(val);
		for(float_t val = qmin[2] + step[2]; val < qmax[2]; val += step[2]) qz_.push_back(val);

		std::cout << "**               NQX x NQY x NQZ: " << qx_.size() << " x " << qy_.size()
					<< " x " << qz_.size() << std::endl;

		return true;
	} // QGrid::create()



} // namespace hig
