/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: distributions.hpp
 *  Created: Jul 02, 2012
 *  Modified: Tue 16 Jul 2013 11:49:45 AM PDT
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

#ifndef _DISTRIBUTIONS_HPP_
#define _DISTRIBUTIONS_HPP_

#include <cmath>
#include <vector>

#include "globals.hpp"
#include "typedefs.hpp"

namespace hig {

	// this class is stateless
	class StatisticalDistributions {	// singleton

		private:
			StatisticalDistributions() { }
			StatisticalDistributions(const StatisticalDistributions&);
			StatisticalDistributions& operator=(const StatisticalDistributions&);

		public:
			static StatisticalDistributions& instance() {
				static StatisticalDistributions stat;
				return stat;
			} // instance()

			// use gnu gsl, or sgi scsl, or boost, or mkl ...
			// improve using mkl and other math libs ...
			// currently taking only integer mean
			float_t* gaussian_dist_2d(unsigned int x_size, unsigned int y_size,
										int x_mean, int y_mean,
										float_t x_sd, float_t y_sd,
										float_t* &out_matrix) {

				int x_min = x_mean - std::floor((float_t)x_size / 2.0);
				if((x_size & 0x01) == 0) x_min = x_min - 1;
				int x_max = x_mean + std::floor((float_t)x_size / 2.0);
				int y_min = y_mean - std::floor((float_t)y_size / 2.0);
				if((y_size & 0x01) == 0) y_min = y_min - 1;
				int y_max = y_mean + std::floor((float_t)y_size / 2.0);

				out_matrix = new (std::nothrow) float_t[x_size * y_size];

				int i = 0, j = 0;
				for(int y = y_min; y <= y_max; ++ y, ++ j) {
					i = 0;
					for(int x = x_min; x <= x_max; ++ x, ++ i) {
						float_t tempx = pow((x - x_mean), 2) / (2 * pow(x_sd, 2));
						float_t tempy = pow((y - y_mean), 2) / (2 * pow(y_sd, 2));
						out_matrix[x_size * j + i] = (float_t) 1.0 / (sqrt(pow(x_sd, 2) + pow(y_sd, 2)) *
															sqrt(2 * PI_)) *
															exp(-1.0 * (tempx + tempy));
					} // for x
				} // for y

				return out_matrix;
			} // gaussian_dist_2d()


			std::vector<float_t>& gaussian_dist_1d(float_t mean, float_t sd,
													std::vector<float_t> &in_arr,
													std::vector<float_t> &out_arr) {
				if(in_arr.size() == 0) return in_arr;
				if(!boost::math::isfinite(mean))
					mean = (in_arr[0] + in_arr[in_arr.size() - 1]) / 2;

				for(std::vector<float_t>::iterator i = in_arr.begin(); i != in_arr.end(); ++ i) {
					float_t temp = exp(pow((*i) - mean, 2) / (2 * sd * sd)) / sqrt(2 * PI_) * sd;
					out_arr.push_back(temp);
				} // for

				return out_arr;
			} // gaussian_dist_1d()


			std::vector<float_t>& uniform_dist() {
				// ...
			} // uniform_dist()


			std::vector<float_t>& uniform_dist_1d(unsigned int size, std::vector<float_t> &out_arr) {
				for(unsigned int i = 0; i < size; ++ i) out_arr.push_back(1.0);
				return out_arr;
			} // uniform_dist()


			// implement more distributions ... (for fun!)

	}; // class StatisticalDistributions

} // namespace

#endif /* _DISTRIBUTIONS_HPP_ */
