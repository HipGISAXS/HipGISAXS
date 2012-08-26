/***
  *  $Id: qgrid_test_create.cpp 38 2012-08-09 23:01:20Z asarje $
  *
  *  Project:
  *
  *  File: qgrid_test_create.cpp
  *  Created: Aug 01, 2012
  *  Modified: Thu 02 Aug 2012 05:22:25 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <complex>

#include "qgrid.hpp"

namespace hig {

	bool QGrid::create_test() {

		int nqx = 3;
		int nqy = 295;
		int nqz = 201;
		int nqz_ext = 4 * nqz;

		float_t xbegin = -1.9227496248816338;
		float_t xstep = -1.3420298154914343 + 1.9227496248816338;
		float_t ybegin = -295;
		float_t ystep = 2;
		//float_t zbegin = 141.0339451916621;
		float_t zbegin = 207.0339451916621;
		float_t zstep = 2;
		//std::complex<float_t> zebegin(35.381377193987376, -93.85377544760027);
		float_t qx = xbegin, qy = ybegin, qz = zbegin;
		for(int i = 0; i < nqx; ++ i, qx += xstep) qx_.push_back(qx);
		for(int i = 0; i < nqy; ++ i, qy += ystep) qy_.push_back(qy);
		for(int i = 0; i < nqz; ++ i, qz += zstep) qz_.push_back(qz);
		float_t k0 = 50642.25919495042;
		float_t alphai = 0.12 * 3.1419 / 180;
		float_t kzi_0 = -1.0 * k0 * sin(alphai);
		std::complex<float_t> dnl_q(2.0 * 1.95e-06, 2.0 * 2.15e-09);
		std::complex<float_t> kzi_q_temp(pow(sin(alphai), 2), 0);
		std::complex<float_t> kzi_q = -1.0 * k0 * sqrt(kzi_q_temp - dnl_q);
		create_qz_extended(k0, kzi_0, kzi_q, dnl_q);

		return true;
	} // QGrid::create_test()

} // namespace hig
