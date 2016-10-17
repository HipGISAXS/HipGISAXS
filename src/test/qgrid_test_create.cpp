/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: qgrid_test_create.cpp
 *  Created: Aug 01, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <complex>

#include <model/qgrid.hpp>

namespace hig {

  bool QGrid::create_test() {

    int nqx = 3;
    int nqy = 295;
    int nqz = 201;
    int nqz_ext = 4 * nqz;

    real_t xbegin = -1.9227496248816338;
    real_t xstep = -1.3420298154914343 + 1.9227496248816338;
    real_t ybegin = -295;
    real_t ystep = 2;
    //real_t zbegin = 141.0339451916621;
    real_t zbegin = 207.0339451916621;
    real_t zstep = 2;
    //std::complex<real_t> zebegin(35.381377193987376, -93.85377544760027);
    real_t qx = xbegin, qy = ybegin, qz = zbegin;
    for(int i = 0; i < nqx; ++ i, qx += xstep) qx_.push_back(qx);
    for(int i = 0; i < nqy; ++ i, qy += ystep) qy_.push_back(qy);
    for(int i = 0; i < nqz; ++ i, qz += zstep) qz_.push_back(qz);
    real_t k0 = 50642.25919495042;
    real_t alphai = 0.12 * 3.1419 / 180;
    real_t kzi_0 = -1.0 * k0 * sin(alphai);
    std::complex<real_t> dnl_q(2.0 * 1.95e-06, 2.0 * 2.15e-09);
    std::complex<real_t> kzi_q_temp(pow(sin(alphai), 2), 0);
    std::complex<real_t> kzi_q = (real_t)- k0 * sqrt(kzi_q_temp - dnl_q);
    create_qz_extended(k0, kzi_0, kzi_q, dnl_q);

    return true;
  } // QGrid::create_test()

} // namespace hig
/ namespace hig
