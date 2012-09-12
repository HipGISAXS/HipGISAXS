/***
  *  $Id: sf.cpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS - High Performance GISAXS
  *
  *  File: sf.cpp
  *  Created: Jun 18, 2012
  *  Modified: Wed 12 Sep 2012 02:48:49 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>

#include "sf.hpp"
#include "qgrid.hpp"
#include "utilities.hpp"

namespace hig {

	StructureFactor::StructureFactor(): sf_(NULL), nx_(0), ny_(0), nz_(0) {
	} // StructureFactor::StructureFactor()


	StructureFactor::~StructureFactor() {
		if(sf_ != NULL) delete[] sf_;
	} // StructureFactor::~StructureFactor()


	bool StructureFactor::compute_structure_factor(std::string expt, vector3_t center,
							Lattice* lattice, vector3_t repet,
							vector3_t rotation_1, vector3_t rotation_2, vector3_t rotation_3,
							MPI::Intracomm& world_comm) {
		int my_rank = world_comm.Get_rank();

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
		vector3_t temp_la(lattice->a()), temp_lb(lattice->b());
		temp_la[2] = 0; temp_lb[2] = 0;
		mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_la, arot);
		mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lb, brot);
		mat_mul_3x1(rotation_1, rotation_2, rotation_3, lattice->c(), crot);

		vector3_t l_t = lattice->t();

		sf_ = NULL;
		sf_ = new (std::nothrow) complex_t[nx_ * ny_ * nz_];
		if(sf_ == NULL) {
			if(my_rank == 0)
				std::cerr << "error: could not allocate memory for structure factor" << std::endl;
			return false;
		} // if

		if(my_rank == 0) {
			std::cout << "**        Using input shape file: data/opv_zyplane.hd5";
			std::cout << "-- Computing structure factor ... " << std::endl;
			std::cout << "**               nqx x nqy x nqz: " << nx_ << " x "
						<< ny_ << " x " << nz_ << std::endl;
		} // if

		if(my_rank == 0) {
			std::cout << "arot: " << arot[0] << ", " << arot[1] << ", " << arot[2] << std::endl;
			std::cout << "brot: " << brot[0] << ", " << brot[1] << ", " << brot[2] << std::endl;
			std::cout << "crot: " << crot[0] << ", " << crot[1] << ", " << crot[2] << std::endl;
		} // if

		std::complex<float_t> unit_c(1, 0);
		std::complex<float_t> unit_ci(0, 1);

		// good for gpu ...
		for(unsigned int z = 0; z < nz_; ++ z) {
			for(unsigned int y = 0; y < ny_; ++ y) {
				for(unsigned int x = 0; x < nx_; ++ x) {
					/*std::complex<float_t> temp1, temp_x2, temp_y3, temp_y4, temp_x5;
					std::complex<float_t> sa, sb, sc;

					temp1 = exp(std::complex<float_t>(0, arot[0] * x + arot[1] * y + arot[2] * z));
					temp_x2 = unit_c - pow(temp1, repet[0]);
					temp_y3 = unit_c / (unit_c - temp1);
					temp_y4 = unit_c / (unit_c / temp_y3 + (float_t)(!boost::math::isfinite(temp_y3)));
					temp_x5 = temp_x2 + repet[0] * (!boost::math::isfinite(unit_c / temp_x2));
					sa = pow(temp1, (1 - repet[0]) / 2) * temp_y4 * temp_x5;

					temp1 = exp(std::complex<float_t>(0, brot[0] * x + brot[1] * y + brot[2] * z));
					temp_x2 = unit_c - pow(temp1, repet[1]);
					temp_y3 = unit_c / (unit_c - temp1);
					temp_y4 = unit_c / (unit_c / temp_y3 + (float_t)(!boost::math::isfinite(temp_y3)));
					temp_x5 = temp_x2 + repet[1] * (!boost::math::isfinite(unit_c / temp_x2));
					sb = pow(temp1, (1 - repet[1]) / 2) * temp_y4 * temp_x5;

					temp1 = exp(std::complex<float_t>(0, crot[0] * x + crot[1] * y + crot[2] * z));
					temp_x2 = unit_c - pow(temp1, repet[2]);
					temp_y3 = unit_c / (unit_c - temp1);
					temp_y4 = unit_c / (unit_c / temp_y3 + (float_t)(!boost::math::isfinite(temp_y3)));
					temp_x5 = temp_x2 + repet[2] * (!boost::math::isfinite(unit_c / temp_x2));
					sc = temp_y4 * temp_x5;

					std::complex<float_t> temp6 = exp(std::complex<float_t>(0, 1) *
									(center[0] * x + center[1] * y + center[2] * z)) *
									sa * sb * sc *
									(unit_c + exp(std::complex<float_t>(0, 1) *
									(l_t[0] * x + l_t[1] * y + l_t[2] * z)));
					sf_[nx_ * ny_ * z + nx_ * y + x].x = temp6.real();
					sf_[nx_ * ny_ * z + nx_ * y + x].y = temp6.imag(); */

					complex_t temp1, temp_x2, temp_y3, temp_y4, temp_x5;
					complex_t sa, sb, sc;
					float_t qx = QGrid::instance().qx(x);
					float_t qy = QGrid::instance().qy(y);
					complex_t qz;
					if(expt == "saxs")
						qz = QGrid::instance().qz(z);
					else if(expt == "gisaxs")
						qz = QGrid::instance().qz_extended(z);

					temp1 = exp(unit_ci * (arot[0] * qx + arot[1] * qy + arot[2] * qz));
					temp_x2 = unit_c - pow(temp1, repet[0]);
					temp_y3 = unit_c / (unit_c - temp1);
					//std::cout << "temp_y3 = " << temp_y3 << ", isfinite = "
					//			<< boost::math::isfinite(temp_y3)
					//			<< ", unit_c / temp_x2 = " << unit_c / temp_x2 << ", isfinite = "
					//			<< boost::math::isfinite(unit_c / temp_x2)
					//			<< std::endl;
					temp_y4 = unit_c / (unit_c / temp_y3 + (float_t)(!boost::math::isfinite(temp_y3)));
					temp_x5 = temp_x2 + repet[0] * (!boost::math::isfinite(unit_c / temp_x2));
					sa = pow(temp1, (1 - repet[0]) / 2) * temp_y4 * temp_x5;

					temp1 = exp(unit_ci * (brot[0] * qx + brot[1] * qy + brot[2] * qz));
					temp_x2 = unit_c - pow(temp1, repet[1]);
					temp_y3 = unit_c / (unit_c - temp1);
					//std::cout << "temp_y3 = " << temp_y3 << ", isfinite = "
					//			<< boost::math::isfinite(temp_y3)
					//			<< ", unit_c / temp_x2 = " << unit_c / temp_x2 << ", isfinite = "
					//			<< boost::math::isfinite(unit_c / temp_x2)
					//			<< std::endl;
					temp_y4 = unit_c / (unit_c / temp_y3 + (float_t)(!boost::math::isfinite(temp_y3)));
					temp_x5 = temp_x2 + repet[1] * (!boost::math::isfinite(unit_c / temp_x2));
					sb = pow(temp1, (1 - repet[1]) / 2) * temp_y4 * temp_x5;

					temp1 = exp(unit_ci * (crot[0] * qx + crot[1] * qy + crot[2] * qz));
					temp_x2 = unit_c - pow(temp1, repet[2]);
					temp_y3 = unit_c / (unit_c - temp1);
					//std::cout << "temp_y3 = " << temp_y3 << ", isfinite = "
					//			<< boost::math::isfinite(temp_y3)
					//			<< ", unit_c / temp_x2 = " << unit_c / temp_x2 << ", isfinite = "
					//			<< boost::math::isfinite(unit_c / temp_x2)
					//			<< std::endl;
					temp_y4 = unit_c / (unit_c / temp_y3 + (float_t)(!boost::math::isfinite(temp_y3)));
					temp_x5 = temp_x2 + repet[2] * (!boost::math::isfinite(unit_c / temp_x2));
					sc = temp_y4 * temp_x5;

					sf_[nx_ * ny_ * z + nx_ * y + x] = exp(unit_ci *
									(center[0] * qx + center[1] * qy + center[2] * qz)) *
									sa * sb * sc *
									(unit_c + exp(unit_ci * (l_t[0] * qx + l_t[1] * qy + l_t[2] * qz)));

				} // for x
			} // for y
		} // for z

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
