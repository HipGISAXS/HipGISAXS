/***
  *  $Id: sf.cpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: sf.cpp
  *  Created: Jun 18, 2012
  *  Modified: Fri 23 Nov 2012 12:16:35 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>

#include "woo/timer/woo_boostchronotimers.hpp"

#include "sf.hpp"
#include "qgrid.hpp"
#include "utilities.hpp"

namespace hig {

	StructureFactor::StructureFactor(): sf_(NULL), nx_(0), ny_(0), nz_(0) {
	} // StructureFactor::StructureFactor()


	StructureFactor::~StructureFactor() {
		if(sf_ != NULL) delete[] sf_;
		sf_ = NULL;
	} // StructureFactor::~StructureFactor()


	void StructureFactor::clear() {
		if(sf_ != NULL) delete[] sf_;
		sf_ = NULL;
		nx_ = ny_ = nz_ = 0;
	} // StructureFactor::clear()


	/**
	 * compute structure factor sequentially on cpu
	 */
	bool StructureFactor::compute_structure_factor(std::string expt, vector3_t center,
							Lattice* lattice, vector3_t repet,
							vector3_t rotation_1, vector3_t rotation_2, vector3_t rotation_3,
							MPI::Intracomm& world_comm) {
		int my_rank = world_comm.Get_rank();

		woo::BoostChronoTimer maintimer, computetimer;

		maintimer.start();

		nx_ = QGrid::instance().nqx();
		ny_ = QGrid::instance().nqy();
		if(expt == "saxs")
			nz_ = QGrid::instance().nqz();
		else if(expt == "gisaxs")
			nz_ = QGrid::instance().nqz_extended();
		else
			return false;

		if(repet[0] < 1) repet[0] = 1;
		if(repet[1] < 1) repet[1] = 1;
		if(repet[2] < 1) repet[2] = 1;

		vector3_t arot(0, 0, 0), brot(0, 0, 0), crot(0, 0, 0);
		vector3_t temp_la(lattice->a()), temp_lb(lattice->b()), temp_lc(lattice->c());
		temp_la[2] = 0; temp_lb[2] = 0;
		mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_la, arot);
		mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lb, brot);
		mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lc, crot);

		vector3_t l_t = lattice->t();

		sf_ = NULL;
		sf_ = new (std::nothrow) complex_t[nx_ * ny_ * nz_];
		if(sf_ == NULL) {
			if(my_rank == 0)
				std::cerr << "error: could not allocate memory for structure factor" << std::endl;
			return false;
		} // if

		if(my_rank == 0) std::cout << "-- Computing structure factor ... " << std::flush;

		std::complex<float_t> unit_c(1, 0);
		std::complex<float_t> unit_ci(0, 1);

		computetimer.start();

		// good for gpu ...
		for(unsigned int z = 0; z < nz_; ++ z) {
			for(unsigned int y = 0; y < ny_; ++ y) {
				for(unsigned int x = 0; x < nx_; ++ x) {
					complex_t temp1, temp_x2, temp_y3, temp_y4, temp_x5;
					float_t temp_f;
					complex_t sa, sb, sc;
					float_t qx = QGrid::instance().qx(x);
					float_t qy = QGrid::instance().qy(y);
					complex_t qz;
					if(expt == "saxs") qz = QGrid::instance().qz(z);
					else if(expt == "gisaxs") qz = QGrid::instance().qz_extended(z);

					temp1 = exp(unit_ci * (arot[0] * qx + arot[1] * qy + arot[2] * qz));
					temp_x2 = unit_c - pow(temp1, repet[0]);
					temp_y3 = unit_c / (unit_c - temp1);
					temp_f = (float_t)(!((boost::math::isfinite)(temp_y3.real()) &&
											(boost::math::isfinite)(temp_y3.imag())));
					temp_y4 = unit_c / (unit_c / temp_y3 + temp_f);
					temp_f = (float_t)(!((boost::math::isfinite)((unit_c / temp_x2).real()) &&
											(boost::math::isfinite)((unit_c / temp_x2).imag())));
					temp_x5 = temp_x2 + repet[0] * temp_f;
					sa = pow(temp1, ((float_t)1.0 - repet[0]) / (float_t)2.0) * temp_y4 * temp_x5;

					temp1 = exp(unit_ci * (brot[0] * qx + brot[1] * qy + brot[2] * qz));
					temp_x2 = unit_c - pow(temp1, repet[1]);
					temp_y3 = unit_c / (unit_c - temp1);
					temp_f = (float_t)(!((boost::math::isfinite)(temp_y3.real()) &&
											(boost::math::isfinite)(temp_y3.imag())));
					temp_y4 = unit_c / (unit_c / temp_y3 + temp_f);
					temp_f = (float_t)(!((boost::math::isfinite)((unit_c / temp_x2).real()) &&
											(boost::math::isfinite)((unit_c / temp_x2).imag())));
					temp_x5 = temp_x2 + repet[1] * temp_f;
					sb = pow(temp1, ((float_t)1.0 - repet[1]) / (float_t)2.0) * temp_y4 * temp_x5;

					temp1 = exp(unit_ci * (crot[0] * qx + crot[1] * qy + crot[2] * qz));
					temp_x2 = unit_c - pow(temp1, repet[2]);
					temp_y3 = unit_c / (unit_c - temp1);
					temp_f = (float_t)(!((boost::math::isfinite)(temp_y3.real()) &&
											(boost::math::isfinite)(temp_y3.imag())));
					temp_y4 = unit_c / (unit_c / temp_y3 + temp_f);
					temp_f = (float_t)(!((boost::math::isfinite)((unit_c / temp_x2).real()) &&
											(boost::math::isfinite)((unit_c / temp_x2).imag())));
					temp_x5 = temp_x2 + repet[2] * temp_f;
					sc = temp_y4 * temp_x5;

					/*if(!((boost::math::isfinite)(sa.real()) && (boost::math::isfinite)(sa.imag()))) {
						std::cout << "sa sa sa sa sa sa sa: " << x << ", " << y << ", " << z << std::endl; }
					if(!((boost::math::isfinite)(sb.real()) && (boost::math::isfinite)(sb.imag()))) {
						std::cout << "sb sb sb sb sb sb sb: " << x << ", " << y << ", " << z << std::endl; }
					if(!((boost::math::isfinite)(sc.real()) && (boost::math::isfinite)(sc.imag()))) {
						std::cout << "sc sc sc sc sc sc sc: " << x << ", " << y << ", " << z << std::endl; }*/

					sf_[nx_ * ny_ * z + nx_ * y + x] = exp(unit_ci *
									(center[0] * qx + center[1] * qy + center[2] * qz)) *
									sa * sb * sc *
									(unit_c + exp(unit_ci * (l_t[0] * qx + l_t[1] * qy + l_t[2] * qz)));

					/*if(!((boost::math::isfinite)(sf_[nx_ * ny_ * z + nx_ * y + x].real()) &&
								(boost::math::isfinite)(sf_[nx_ * ny_ * z + nx_ * y + x].imag()))) {
						std::cout << "sf sf sf sf sf sf sf: " << x << ", " << y << ", " << z << std::endl; }*/
				} // for x
			} // for y
		} // for z

		computetimer.stop();
		maintimer.stop();

		if(my_rank == 0) {
			std::cout << "done. " << std::endl;
			std::cout << "**               SF compute time: " << computetimer.elapsed_msec()  << " ms."
						<< std::endl
						<< "**                 Total SF time: " << maintimer.elapsed_msec() << " ms."
						<< std::endl;
			//int naninfs = count_naninfs(nx_, ny_, nz_, sf_);
			//std::cout << " ------- " << naninfs << " / " << nx_ * ny_ * nz_ << " nans or infs" << std::endl;
		} // if

		return true;
	} // StructureFactor::compute_structure_factor()


	void StructureFactor::save_sf(unsigned int nqx, unsigned int nqy, unsigned int nqz,
								const char* filename) {
		std::ofstream f(filename);
		for(unsigned int z = 0; z < nqz; ++ z) {
			for(unsigned int y = 0; y < nqy; ++ y) {
				for(unsigned int x = 0; x < nqx; ++ x) {
					unsigned int index = nqx * nqy * z + nqx * y + x;
					f << sf_[index].real() << "\t" << sf_[index].imag() << std::endl;
				} // for
				f << std::endl;
			} // for
			f << std::endl;
		} // for
		f.close();
	} // StructureFactor::save_sf()

} // namespace hig
