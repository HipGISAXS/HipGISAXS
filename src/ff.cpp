/***
  *  $Id: ff.cpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff.cpp
  *  Created: Jul 17, 2012
  *  Modified: Fri 23 Nov 2012 01:23:54 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <fstream>

#include "ff.hpp"

namespace hig {

	void FormFactor::clear() {
		ff_.clear();
		analytic_ff_.clear();
		if(nff_ != NULL) delete[] nff_;
		nff_ = NULL;
		numeric_ff_.clear();
	} // FormFactor::clear()

	//template<typename float_t, typename complex_t>
	//bool FormFactor<float_t, complex_t>::compute_form_factor(ShapeName shape, std::string shape_filename,
	bool FormFactor::compute_form_factor(ShapeName shape, std::string shape_filename,
										shape_param_list_t& params, float_t single_thickness,
										vector3_t& transvec, float_t shp_tau, float_t shp_eta,
										vector3_t& rot1, vector3_t& rot2, vector3_t& rot3,
										MPI::Intracomm& world_comm) {
		if(shape == shape_custom) {
			/* compute numerically */
			is_analytic_ = false;
			numeric_ff_.init();
			numeric_ff_.compute(shape_filename.c_str(), nff_, world_comm);
		} else {
			/* compute analytically */
			is_analytic_ = true;
			analytic_ff_.init(rot1, rot2, rot3, ff_);
			analytic_ff_.compute(shape, shp_tau, shp_eta, transvec,
								ff_, params, single_thickness, rot1, rot2, rot3, world_comm);
		} // if-else

		return true;
	} // FormFactor::compute_form_factpr()


	/* temporaries */

	bool FormFactor::read_form_factor(const char* filename,
									unsigned int nqx, unsigned int nqy, unsigned int nqz) {
		std::ifstream f(filename);
		nff_ = new (std::nothrow) cucomplex_t[nqx * nqy * nqz];
		for(unsigned int z = 0; z < nqz; ++ z) {
			for(unsigned int y = 0; y < nqy; ++ y) {
				for(unsigned int x = 0; x < nqx; ++ x) {
					unsigned int index = nqx * nqy * z + nqx * y + x;
					f >> nff_[index].x;
					f >> nff_[index].y;
				} // for
			} // for
		} // for
		f.close();

		return true;
	} // FormFactor::read_form_factor()


	void FormFactor::print_ff(unsigned int nqx, unsigned int nqy, unsigned int nqz) {
		for(unsigned int z = 0; z < nqz; ++ z) {
			for(unsigned int y = 0; y < nqy; ++ y) {
				for(unsigned int x = 0; x < nqx; ++ x) {
					unsigned int index = nqx * nqy * z + nqx * y + x;
					std::cout << nff_[index].x << "," << nff_[index].y << "\t";
				} // for
				std::cout << std::endl;
			} // for
			std::cout << std::endl;
		} // for
	} // FormFactor::print_ff()

	void FormFactor::save_ff(unsigned int nqx, unsigned int nqy, unsigned int nqz, const char* filename) {
		std::ofstream f(filename);
		for(unsigned int z = 0; z < nqz; ++ z) {
			for(unsigned int y = 0; y < nqy; ++ y) {
				for(unsigned int x = 0; x < nqx; ++ x) {
					unsigned int index = nqx * nqy * z + nqx * y + x;
					if(is_analytic_) {
						f << ff_[index].real() << "\t" << ff_[index].imag() << std::endl;
					} else {
						f << nff_[index].x << "\t" << nff_[index].y << std::endl;
					}
				} // for
				f << std::endl;
			} // for
			f << std::endl;
		} // for
		f.close();
	} // FormFactor::save_ff()

} // namespace hig
