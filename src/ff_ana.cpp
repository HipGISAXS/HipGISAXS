/***
  *  $Id: ff_ana.cpp 37 2012-08-09 22:59:59Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana.cpp
  *  Created: Jul 12, 2012
  *  Modified: Tue 19 Feb 2013 11:25:34 AM PST
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

	bool AnalyticFormFactor::init(vector3_t &rot1, vector3_t &rot2, vector3_t &rot3,
				std::vector<complex_t> &ff) {
		nqx_ = QGrid::instance().nqx();
		nqy_ = QGrid::instance().nqy();
		nqz_ = QGrid::instance().nqz_extended();

		ff_gpu_.grid_size(nqx_, nqy_, nqz_);

		// first make sure there is no residue from any previous computations
		ff.clear();

		// construct mesh grid thingy
		/*mesh_qx_.clear();
		mesh_qy_.clear();
		mesh_qz_.clear();
		mesh_qx_.reserve(nqx_ * nqy_ * nqz_);
		mesh_qy_.reserve(nqx_ * nqy_ * nqz_);
		mesh_qz_.reserve(nqx_ * nqy_ * nqz_);
		// this is very inefficient! improve ...
		// x and y are swapped for qx, qy, qz 
		for(unsigned int i = 0; i < nqz_; ++ i) {
			for(unsigned int j = 0; j < nqy_; ++ j) {
				for(unsigned int k = 0; k < nqx_; ++ k) {
					mesh_qx_.push_back(rot1[0] * QGrid::instance().qy(j) +
										rot1[1] * QGrid::instance().qx(k) +
										rot1[2] * QGrid::instance().qz_extended(i));
					mesh_qy_.push_back(rot2[0] * QGrid::instance().qy(j) +
										rot2[1] * QGrid::instance().qx(k) +
										rot2[2] * QGrid::instance().qz_extended(i));
					mesh_qz_.push_back(rot3[0] * QGrid::instance().qy(j) +
										rot3[1] * QGrid::instance().qx(k) +
										rot3[2] * QGrid::instance().qz_extended(i));
				} // for k
			} // for j
		} // for i */

		rot_ = new (std::nothrow) float_t[9];
		rot_[0] = rot1[0]; rot_[1] = rot1[1]; rot_[2] = rot1[2];
		rot_[3] = rot2[0]; rot_[4] = rot2[1]; rot_[5] = rot2[2];
		rot_[6] = rot3[0]; rot_[7] = rot3[1]; rot_[8] = rot3[2];

		return true;
	} // AnalyticFormFactor::init()


	void AnalyticFormFactor::clear() {
		nqx_ = nqy_ = nqz_ = 0;
		/*mesh_qx_.clear();
		mesh_qy_.clear();
		mesh_qz_.clear();*/
	} // AnalyticFormFactor::clear()


	bool AnalyticFormFactor::compute(ShapeName shape, float_t tau, float_t eta, vector3_t transvec,
									std::vector<complex_t>& ff,
									shape_param_list_t& params, float_t single_layer_thickness,
									vector3_t rot1, vector3_t rot2, vector3_t rot3,
									MPI::Intracomm& world_comm) {

		std::cout << "-- Computing form factor analytically ... " << std::endl;
		switch(shape) {
			case shape_box:						// cube or box
				if(!compute_box(nqx_, nqy_, nqz_, ff, shape, params, tau, eta, transvec, rot1, rot2, rot3)) {
					std::cerr << "error: something went wrong while computing FF for a box"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_cylinder:				// standing cylinder
				if(!compute_cylinder(params, tau, eta, ff, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a cylinder"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_random_cylinders:		// randomly oriented cylinders, for saxs
				if(!compute_random_cylinders()) {
					std::cerr << "error: something went wrong while computing FF for random cylinders"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_horizontal_cylinder:		// horizontal cylinder
				if(!compute_horizontal_cylinder(params, transvec, ff)) {
					std::cerr << "error: something went wrong while computing FF for a horizontal cylinder"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sphere:					// simple sphere
				if(!compute_sphere(params, ff, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a sphere"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism3:					// triangular prism (prism with 3 sides)
				if(!compute_prism(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a prism3"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism6:					// hexagonal prism (prism with 6 sides)
				if(!compute_prism6(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a prism6"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sawtooth_down:			// downwards sawtooth
				if(!compute_sawtooth_down()) {
					std::cerr << "error: something went wrong while computing FF for a sawtooth down"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sawtooth_up:				// upwards sawtooth
				if(!compute_sawtooth_up()) {
					std::cerr << "error: something went wrong while computing FF for a sawtooth up"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism3x:					// triangular grating in x direction
				if(!compute_prism3x()) {
					std::cerr << "error: something went wrong while computing FF for a prism3x"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_pyramid:					// pyramid
				if(!compute_pyramid()) {
					std::cerr << "error: something went wrong while computing FF for a pyramid"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_trunccone:				// truncated cone
				if(!compute_truncated_cone()) {
					std::cerr << "error: something went wrong while computing FF for a truncated cone"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_truncpyr:				// truncated pyramid
				if(!compute_truncated_pyramid()) {
					std::cerr << "error: something went wrong while computing FF for a truncated pyramid"
								<< std::endl;
					return false;
				} // if
				break;
			default:
				std::cerr << "error: invalid shape. given shape is not supported" << std::endl;
				return false;
		} // switch

		/*int rank;
		MPI_Comm_rank(world_comm, &rank);
		if(rank == 0) {
			int naninfs = count_naninfs(nqx_, nqy_, nqz_, ff);
			std::cout << " ------ " << naninfs << " / " << nqx_ * nqy_ * nqz_ << " nans or infs" << std::endl;
		} // if
*/
		return true;
	} // AnalyticFormFactor::compute()


	/**
	 * box
	 */
/*	bool AnalyticFormFactor::compute_box(unsigned int nqx, unsigned int nqy, unsigned int nqz,
										std::vector<complex_t>& ff,
										ShapeName shape, shape_param_list_t& params,
										float_t tau, float_t eta, vector3_t &transvec,
										vector3_t &rot1, vector3_t &rot2, vector3_t &rot3) {
		std::vector <float_t> x, distr_x;	// for x dimension: param_xsize  param_edge
		std::vector <float_t> y, distr_y;	// for y dimension: param_ysize  param_edge
		std::vector <float_t> z, distr_z;	// for z dimension: param_height param_edge
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			if(!(*i).second.isvalid()) {
				std::cerr << "warning: invalid shape parameter found" << std::endl;
				continue;
			} // if
			switch((*i).second.type()) {
				case param_edge:
					param_distribution((*i).second, x, distr_x);	// x == RR, distr_x == RRD
					param_distribution((*i).second, y, distr_y);
					param_distribution((*i).second, z, distr_z);
					break;
				case param_xsize:
					param_distribution((*i).second, x, distr_x);
					break;
				case param_ysize:
					param_distribution((*i).second, y, distr_y);
					break;
				case param_height:
					param_distribution((*i).second, z, distr_z);
					break;
				case param_radius:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted values for shape type 'box'" << std::endl;
					break;
				default:
					std::cerr << "warning: unknown parameters for shape type 'box'. ignoring"
								<< std::endl;
			} // switch
		} // for

		// check if x y z etc are set or not
		if(x.size() < 1 || y.size() < 1 || z.size() < 1) {
			std::cerr << "error: invalid or not enough box parameters given" << std::endl;
			return false;
		} // if

		// initialize ff
		ff.clear();
		ff.reserve(nqz * nqy * nqx);
		for(unsigned int i = 0; i < nqz * nqy * nqx; ++ i) ff.push_back(complex_t(0, 0));

		// ff computation for a box
		for(unsigned int j_z = 0; j_z < nqz; ++ j_z) {
			for(unsigned int j_y = 0; j_y < nqy; ++ j_y) {
				for(unsigned int j_x = 0; j_x < nqx; ++ j_x) {
					complex_t mqx = QGrid::instance().qy(j_y) * rot_[0] +
									QGrid::instance().qx(j_x) * rot_[1] +
									QGrid::instance().qz_extended(j_z) * rot_[2];
					complex_t mqy = QGrid::instance().qy(j_y) * rot_[3] +
									QGrid::instance().qx(j_x) * rot_[4] +
									QGrid::instance().qz_extended(j_z) * rot_[5];
					complex_t mqz = QGrid::instance().qy(j_y) * rot_[6] +
									QGrid::instance().qx(j_x) * rot_[7] +
									QGrid::instance().qz_extended(j_z) * rot_[8];
					complex_t temp1 = sin(eta) * mqx;
					complex_t temp2 = cos(eta) * mqy;
					complex_t temp3 = temp1 + temp2;
					complex_t temp_qm = tan(tau) * temp3;
					unsigned int curr_index = nqx * nqy * j_z + nqx * j_y + j_x;
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_z = 0; i_z < z.size(); ++ i_z) {
						for(unsigned int i_y = 0; i_y < y.size(); ++ i_y) {
							for(unsigned int i_x = 0; i_x < x.size(); ++ i_x) {
								complex_t temp1 = mqz + temp_qm;
								complex_t temp2 = mqy * z[i_z];
								complex_t temp3 = mqx * x[i_x];
								complex_t temp4 = fq_inv(temp1, y[i_y]);
								complex_t temp5 = sinc(temp2);
								complex_t temp6 = sinc(temp3);
								complex_t temp7 = temp6 * temp5;
								complex_t temp8 = temp7 * temp4;
								complex_t temp9 = 4 * distr_x[i_x] * distr_y[i_y] * distr_z[i_z] *
													z[i_z] * x[i_x];
								complex_t temp10 = temp9 * temp8;
								temp_ff += temp10;
								if(!(boost::math::isfinite(temp10.real()) &&
											boost::math::isfinite(temp10.imag()))) {
									std::cerr << "+++++++++++++++ here it is +++++++ " << j_x << ", "
												<< j_y << ", " << j_z << std::endl;
									exit(1);
								} // if
							} // for i_x
						} // for i_y
					} // for i_z
					complex_t temp_e = exp(mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
					if(!(boost::math::isfinite(temp_e.real()) && boost::math::isfinite(temp_e.imag()))) {
						std::cerr << "---------------- here it is ------ " << j_x << ", "
									<< j_y << ", " << j_z << std::endl;
						exit(1);
					} // if
					ff[curr_index] = temp_ff * temp_e;
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_box()
*/

	/**
	 * cylinder
	 */
/*	bool AnalyticFormFactor::compute_cylinder(shape_param_list_t& params, float_t tau, float_t eta,
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
					std::cerr << "warning: ignoring unwanted input values for cylinder" << std::endl;
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

		return true;
	} // AnalyticFormFactor::compute_cylinder()
*/

	/**
	 * random cylinders
	 */
/*	bool AnalyticFormFactor::compute_random_cylinders() { // for saxs
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_random_cylinders"
					<< std::endl;
		return false;
*/		// ...
		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_random_cylinders()
*/

	/**
	 * horizontal cylinder
	 */
/*	bool AnalyticFormFactor::compute_horizontal_cylinder(shape_param_list_t& params,
														vector3_t transvec,
														std::vector<complex_t>& ff) {
		std::vector<float_t> r, distr_r;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in horizontal cylinder"
								<< std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: unknown or invalid parameter given for horizontal cylinder"
								<< std::endl;
					return false;
			} // switch
		} // for

		if(r.size() < 1 || h.size() < 1) {
			std::cerr << "error: both radius and height parameters are required for horizontal cylinder"
						<< std::endl;
			return false;
		} // if

		// in slims code, why not doing range of r and h ???

		complex_t unitc(0, 1);

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqy_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

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
					complex_t temp_qpar = sqrt(mqz * mqz + mqy * mqy);
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
								temp_ff += 2 * PI_ * h[i_h] * r[i_r] * r[i_r] *
											(cbessj(temp_qpar * r[i_r], 1) / (temp_qpar * r[i_r])) *
											exp(complex_t(-(mqz * r[i_r]).imag(), (mqz * r[i_r]).real())) *
											sinc(mqx * h[i_h] / (float_t)2.0);
						} // for h
					} // for r
					complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					ff[index] = temp_ff * exp(complex_t(-temp1.imag(), temp1.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_horizontal_cylinder()
*/

	/**
	 * prism - 3 face
	 */
/*	bool AnalyticFormFactor::compute_prism(shape_param_list_t& params, std::vector<complex_t>& ff,
											float_t tau, float_t eta, vector3_t transvec) {
		std::vector<float_t> r, distr_r;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: invalid parameter given for prism shape" << std::endl;
					return false;
			} // switch
		} // for

		if(h.size() < 1 || r.size() < 1) {
			std::cerr << "error: need radius and height for prism shape" << std::endl;
			return false;
		} // if

		// why not doing for a range of r, h for this? ... ???

		float_t sqrt3 = sqrt(3.0);
		complex_t unit(0, 1.0);

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqz_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

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
					complex_t temp1 = mqx * (mqx * mqx - 3.0 * mqy * mqy);
					complex_t temp2 = mqz + tan(tau) * (mqx * sin(eta) + mqy * cos(eta));
					complex_t temp_ff(0, 0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							complex_t temp3 = mqy * r[i_r] / sqrt3;
							complex_t temp4 = mqx * exp(unit * mqy * r[i_r] * sqrt3);
							complex_t temp5 = mqx * cos(mqx * r[i_r]);
							complex_t temp6 = unit * sqrt3 * mqy * sin(mqx * r[i_r]);
							complex_t temp7 = fq_inv(temp2, h[i_h]);
							complex_t temp8 = (temp4 - temp5 - temp6) * temp7;
							complex_t temp9 = 2.0 * sqrt3 * exp(complex_t(-temp3.imag(), temp3.real()));
							temp_ff += temp9 / temp1 * temp8;
						} // for h
					} // for r
					complex_t temp10 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					ff[index] = temp_ff * exp(complex_t(-temp10.imag(), temp10.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_prism()
*/

	/**
	 * six faceted prism
	 */
/*	bool AnalyticFormFactor::compute_prism6(shape_param_list_t& params, std::vector<complex_t>& ff,
											float_t tau, float_t eta, vector3_t transvec) {
		std::vector<float_t> r, distr_r;
		std::vector<float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: invalid parameter given for prism shape" << std::endl;
					return false;
			} // switch
		} // for

		float_t sqrt3 = sqrt(3.0);

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqz_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

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
					complex_t qm = tan(tau) * (mqx * sin(eta) + mqy * cos(eta));
					complex_t temp1 = ((float_t)4.0 * sqrt3) / (3.0 * mqy * mqy - mqx * mqx);
					complex_t temp_ff(0, 0);
					for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							complex_t rmqx = r[i_r] * mqx / sqrt3;
							complex_t rmqy = r[i_r] * mqy;
							complex_t temp2 = rmqy * rmqy *	sinc(rmqx) * sinc(rmqy);
							complex_t temp3 = cos(2.0 * rmqx);
							complex_t temp4 = cos(rmqy) * cos(rmqx);
							complex_t temp5 = temp1 * (temp2 + temp3 - temp4);
							complex_t temp6 = fq_inv(mqz + qm, h[i_h]);
							temp_ff += temp5 * temp6;
						} // for h
					} // for r
					complex_t temp7 = (mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2]);
					ff[index] = temp_ff * exp(complex_t(-temp7.imag(), temp7.real()));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_prism6()
*/

	/**
	 * triangular grating in the x-direction
	 */
/*	bool AnalyticFormFactor::compute_prism3x() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_prism3x"
					<< std::endl;
		return false;
*/		// ...
		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_prism3x()
*/

	/**
	 * upwards sawtooth
	 */
/*	bool AnalyticFormFactor::compute_sawtooth_up() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_up"
					<< std::endl;
		return false;
*/		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_sawtooth_up()
*/

	/**
	 * downwards sawtooth
	 */
/*	bool AnalyticFormFactor::compute_sawtooth_down() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_down"
					<< std::endl;
		return false;
*/		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_sawtooth_down()
*/

	/**
	 * pyramid
	 */
/*	bool AnalyticFormFactor::compute_pyramid() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_pyramid"
					<< std::endl;
		return false;
*/		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_pyramid()
*/

	/**
	 * truncated cone
	 */
/*	bool AnalyticFormFactor::compute_truncated_cone() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_cone"
					<< std::endl;
		return false;
*/		// ...
		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_pyramid()
*/

	/**
	 * truncated pyramid
	 */
/*	bool AnalyticFormFactor::compute_truncated_pyramid() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_pyramid"
					<< std::endl;
		return false;
		// ...
*/		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */
/*	} // AnalyticFormFactor::compute_truncated_pyramid()
*/

	/**
	 * matrix computation helpers
	 */

	bool AnalyticFormFactor::mat_fq_inv_in(unsigned int x_size,
											unsigned int y_size,
											unsigned int z_size,
											complex_vec_t& matrix, float_t y) {
		for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			*i = fq_inv(*i, y);
		} // for
		
		return true;
	} // AnalyticFormFactor::mat_fq_inv()

	bool AnalyticFormFactor::mat_fq_inv(unsigned int x_size, unsigned int y_size, unsigned int z_size,
										const complex_vec_t& matrix, float_t y, complex_vec_t& result) {
		result.clear();
		for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			result.push_back(fq_inv(*i, y));
		} // for
		return true;
	} // AnalyticFormFactor::mat_fq_inv()

	complex_t AnalyticFormFactor::fq_inv(complex_t value, float_t y) {
		complex_t unitc(0, 1.0);
		complex_t temp = 2.0 * exp(unitc * value * y / (float_t) 2.0) *
							sin(value * y / (float_t) 2.0) / value;
		if(fabs(temp) <= 1e-14) temp = y;
		return temp;
	} // AnalyticFormFactor::fq_inv()


	bool AnalyticFormFactor::mat_sinc(unsigned int x_size, unsigned int y_size, unsigned int z_size,
										const complex_vec_t& matrix, complex_vec_t& result) {
		result.clear();
		for(std::vector<complex_t>::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			result.push_back(sinc(*i));
		} // for
		return true;
	} // AnalyticFormFactor::mat_sinc()

	bool AnalyticFormFactor::mat_sinc_in(unsigned int x_size,
											unsigned int y_size,
											unsigned int z_size,
											std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			*i = sinc(*i);
		} // for
		return true;
	} // AnalyticFormFactor::mat_sinc()

	complex_t AnalyticFormFactor::sinc(complex_t value) {
		complex_t temp;
		if(fabs(value.real()) <= 1e-14 && fabs(value.imag()) <= 1e-14) temp = complex_t(1.0, 0.0);
		else temp = sin(value) / value;
		return temp;
	} // AnalyticFormFactor::sinc()


	/**
	 * generate parameter's distribution
	 */
	bool AnalyticFormFactor::param_distribution(ShapeParam& param, std::vector<float_t>& dim,
												std::vector<float_t>& dim_vals) {
		if(!param.isvalid()) {
			std::cerr << "error: invalid shape parameter encountered" << std::endl;
			return false;
		} // if

		if(param.nvalues() < 1) {
			std::cerr << "error: empty parameter found (nvalues = 0)" << std::endl;
			return false;
		} // if
		float_t pmax = param.max(), pmin = param.min();
		if(pmax < pmin) pmax = pmin;
		if(param.nvalues() > 1) {
			float_t step = fabs(pmax - pmin) / (param.nvalues() - 1);
			float_t curr = pmin;
			do {
				dim.push_back(curr);
				curr += step;
			} while(curr < param.max());	// assumes min < max ...
		} else {
			dim.push_back(pmin);
		} // if-else

		if(param.stat() == stat_none || param.stat() == stat_null) {	// just one value
			dim_vals.push_back(1.0);
		} else if(param.stat() == stat_uniform) {
			for(unsigned int i = 0; i < dim.size(); ++ i) {
				dim_vals.push_back(1.0);
			} // for
		} else if(param.stat() == stat_gaussian) {
			float_t mean = param.mean();
			if(!boost::math::isfinite(mean)) {
				mean = (dim[0] + dim[dim.size() - 1]) / 2;
			} // if
			for(unsigned int i = 0; i < dim.size(); ++ i) {
				dim_vals.push_back(exp(-1.0 * pow((dim[i] - mean), 2) / (2 * pow(param.deviation(), 2)))
									/ (sqrt(2 * PI_) * param.deviation()));
			} // for
		} else if(param.stat() == stat_random) {
			std::cerr << "uh-oh: random statistic has not been implemented yet" << std::endl;
			return false;
			// ...
		} else {
			std::cerr << "error: an invalid statistic value given for shape parameter" << std::endl;
			return false;
		} // if-else

		return true;
	} // AnalyticFormFactor::param_distribution()

} // namespace hig

