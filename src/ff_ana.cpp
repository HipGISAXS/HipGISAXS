/***
  *  $Id: ff_ana.cpp 37 2012-08-09 22:59:59Z asarje $
  *
  *  Project:
  *
  *  File: ff_ana.cpp
  *  Created: Jul 12, 2012
  *  Modified: Mon 27 Aug 2012 11:49:32 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/math/special_functions/fpclassify.hpp>

#include "ff_ana.hpp"
#include "shape.hpp"
#include "enums.hpp"
#include "qgrid.hpp"
#include "utilities.hpp"

namespace hig {

	bool AnalyticFormFactor::init(vector3_t &rot1, vector3_t &rot2, vector3_t &rot3,
				std::vector<complex_t> &ff) {
		nqx_ = QGrid::instance().nqx();
		nqy_ = QGrid::instance().nqy();
		nqz_ = QGrid::instance().nqz();

		//ff_ = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];	// change to vector? ...
		//ff = ff_;		// i hope this works
		ff.clear();

		// make these into vector ...
		//mesh_qx_ = new (std::nothrow) float_t[nqx_ * nqy_ * nqz_];
		//mesh_qy_ = new (std::nothrow) float_t[nqx_ * nqy_ * nqz_];
		//mesh_qz_ = new (std::nothrow) float_t[nqx_ * nqy_ * nqz_];
		mesh_qx_.clear();
		mesh_qy_.clear();
		mesh_qz_.clear();

		// this is very inefficient! improve ...
		for(unsigned int i = 0; i < nqz_; ++ i) {
			for(unsigned int j = 0; j < nqy_; ++ j) {
				for(unsigned int k = 0; k < nqx_; ++ k) {
					mesh_qx_.push_back(rot1[0] * QGrid::instance().qx(k) +
										rot1[1] * QGrid::instance().qy(j) +
										rot1[2] * QGrid::instance().qz(i));
					mesh_qy_.push_back(rot2[0] * QGrid::instance().qx(k) +
										rot2[1] * QGrid::instance().qy(j) +
										rot2[2] * QGrid::instance().qz(i));
					mesh_qz_.push_back(rot3[0] * QGrid::instance().qx(k) +
										rot3[1] * QGrid::instance().qy(j) +
										rot3[2] * QGrid::instance().qz(i));
					/*mesh_qx_[nqx_ * nqz_ * i + nqx_ * j + k] = rot1[0] * QGrid::instance().qx(k) +
																rot1[1] * QGrid::instance().qy(j) +
																rot1[2] * QGrid::instance().qz(i);
					mesh_qy_[nqx_ * nqz_ * i + nqx_ * j + k] = rot2[0] * QGrid::instance().qx(k) +
																rot2[1] * QGrid::instance().qy(j) +
																rot2[2] * QGrid::instance().qz(i);
					mesh_qz_[nqx_ * nqz_ * i + nqx_ * j + k] = rot3[0] * QGrid::instance().qx(k) +
																rot3[1] * QGrid::instance().qy(j) +
																rot3[2] * QGrid::instance().qz(i); */
				} // for k
			} // for j
		} // for i

		return true;
	} // AnalyticFormFactor::init()


	bool AnalyticFormFactor::compute(ShapeName shape, float_t tau, float_t eta, vector3_t transvec,
									std::vector<complex_t>& ff,
									shape_param_list_t& params, float_t single_layer_thickness,
									vector3_t rot1, vector3_t rot2, vector3_t rot3,
									MPI::Intracomm& world_comm) {
		switch(shape) {
			case shape_box:
				if(!compute_box(nqx_, nqy_, nqz_, ff, shape, params, tau, eta, transvec, rot1, rot2, rot3)) {
					std::cerr << "error: something went wrong while computing FF for a box"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_cylinder:
				if(!compute_cylinder(params, tau, eta, ff, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a cylinder"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_random_cylinders:		// for saxs
				if(!compute_random_cylinders()) {
					std::cerr << "error: something went wrong while computing FF for random cylinders"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_horizontal_cylinder:
				if(!compute_horizontal_cylinder(params, transvec, ff)) {
					std::cerr << "error: something went wrong while computing FF for a horizontal cylinder"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sphere:
				if(!compute_sphere(params, ff, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a sphere"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism3:
				if(!compute_prism(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a prism"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism6:
				if(!compute_prism6()) {
					std::cerr << "error: something went wrong while computing FF for a prism6"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sawtooth_down:
				if(!compute_sawtooth_down()) {
					std::cerr << "error: something went wrong while computing FF for a sawtooth down"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sawtooth_up:
				if(!compute_sawtooth_up()) {
					std::cerr << "error: something went wrong while computing FF for a sawtooth up"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism3x:
				if(!compute_prism3x()) {
					std::cerr << "error: something went wrong while computing FF for a prism3x"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_pyramid:
				if(!compute_pyramid()) {
					std::cerr << "error: something went wrong while computing FF for a pyramid"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_trunccone:
				if(!compute_truncated_cone()) {
					std::cerr << "error: something went wrong while computing FF for a truncated cone"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_truncpyr:
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

		return true;
	} // AnalyticFormFactor::compute()


	/**
	 * box
	 */
	bool AnalyticFormFactor::compute_box(unsigned int nqx, unsigned int nqy, unsigned int nqz,
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

		std::vector<complex_t> mesh_qm = mat_mul(tan(tau), //nqx, nqy, nqz,
												mat_add(nqx, nqy, nqz,
												mat_mul(/*nqx, nqy, nqz,*/ mesh_qx_, sin(eta)),
												nqx, nqy, nqz,
												mat_mul(/*nqx_, nqy_, nqz_,*/ mesh_qy_, cos(eta))));
		// ff computation for a box
		for(unsigned int i_z = 0; i_z < z.size(); ++ i_z) {
			for(unsigned int i_y = 0; i_y < y.size(); ++ i_y) {
				for(unsigned int i_x = 0; i_x < x.size(); ++ i_x) {
					// scalar matrix multiplication
					ff = mat_add(nqx_, nqy_, nqz_, ff, nqx_, nqy_, nqz_,
							mat_mul(distr_x[i_x] * distr_y[i_y] * distr_z[i_z] * 4 * z[i_z] * x[i_x],
							//nqx, nqy, nqz,
							mat_dot_prod(nqx, nqy, nqz,
							mat_dot_prod(nqx, nqy, nqz,
							mat_sinc(nqx, nqy, nqz, mat_mul(/*nqx, nqy, nqz,*/ mesh_qx_, x[i_x])),
							nqx, nqy, nqz, mat_sinc(nqx, nqy, nqz,
							mat_mul(/*nqx, nqy, nqz,*/ mesh_qz_, z[i_z]))),
							nqx, nqy, nqz,
							mat_fq_inv(nqx, nqy, nqz,
							mat_add(nqx, nqy, nqz, mesh_qz_, nqx, nqy, nqz, mesh_qm),
							y[i_y]))));
				} // for i_x
			} // for i_y
		} // for i_z

		std::vector<complex_t>::iterator i_qx = mesh_qx_.begin();
		std::vector<complex_t>::iterator i_qy = mesh_qy_.begin();
		std::vector<complex_t>::iterator i_qz = mesh_qz_.begin();
		for(int i_z = 0; i_qz != mesh_qz_.end(); ++ i_qz, ++ i_z) {
			for(int i_y = 0; i_qy != mesh_qy_.end(); ++ i_qy, ++ i_y) {
				for(int i_x = 0; i_qx != mesh_qx_.end(); ++ i_qx, ++ i_x) {
					complex_t temp = exp(mesh_qx_[i_x] * transvec[0] +
										mesh_qy_[i_y] * transvec[1] +
										mesh_qz_[i_z] * transvec[2]);
					ff[nqx * nqy * i_z + nqx * i_y + i_x] *= temp;
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_box()


	/**
	 * cylinder
	 */
	bool AnalyticFormFactor::compute_cylinder(shape_param_list_t& params, float_t tau, float_t eta,
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

		std::vector<complex_t> qpar = mat_sqrt(//nqx_, nqy_, nqz_,
										mat_add(nqx_, nqy_, nqz_,
										mat_sqr(/*nqx_, nqy_, nqz_,*/ mesh_qx_),
										nqx_, nqy_, nqz_,
										mat_sqr(/*nqx_, nqy_, nqz_,*/ mesh_qy_)));
		std::vector<complex_t> qm = mat_mul(/*nqx_, nqy_, nqz_,*/ tan(tau),
										mat_add(nqx_, nqy_, nqz_,
										mat_mul(/*nqx_, nqy_, nqz_,*/ mesh_qx_, sin(eta)),
										nqx_, nqy_, nqz_,
										mat_mul(/*nqx_, nqy_, nqz_,*/ mesh_qy_, cos(eta))));
		for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
			for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
				ff = mat_add(nqx_, nqy_, nqz_, ff, nqx_, nqy_, nqz_,
						mat_mul(distr_r[i_r] * distr_h[i_h] * 2.0 * PI_ * pow(r[i_r], 2), //nqx_, nqy_, nqz_,
						mat_dot_prod(nqx_, nqy_, nqz_,
						mat_dot_div(nqx_, nqy_, nqz_,
						mat_besselj(1, nqx_, nqy_, nqz_,
						mat_mul(/*nqx_, nqy_, nqz_,*/ qpar, r[i_r])), nqx_, nqy_, nqz_,
						mat_mul(/*nqx_, nqy_, nqz_,*/ qpar, r[i_r])), nqx_, nqy_, nqz_,
						mat_fq_inv(nqx_, nqy_, nqz_,
						mat_add(nqx_, nqy_, nqz_, mesh_qz_, nqx_, nqy_, nqz_, qm), h[i_h]))));
			} // for h
		} // for r

		//ff_ = mat_dot_prod(nqx_, nqy_, nqz_, ff_, nqx_, nqy_, nqz_,
		//		mat_exp(nqx_, nqy_, nqz_,
		//		mat_mul(complex_t(0, 1), nqx_, nqy_, nqz_,
		//		mat_add(nqx_, nqy_, nqz_,
		//		mat_mul(nqx_, nqy_, nqz_, mesh_qx_, transvec[0]), nqx_, nqy_, nqz_,
		//		mat_add(nqx_, nqy_, nqz_,
		//		mat_mul(nqx_, nqy_, nqz_, mesh_qy_, transvec[1]), nqx_, nqy_, nqz_,
		//		mat_mul(nqx_, nqy_, nqz_, mesh_qz_, transvec[2]))))));
		std::vector<complex_t>::iterator i_qx = mesh_qx_.begin();
		std::vector<complex_t>::iterator i_qy = mesh_qy_.begin();
		std::vector<complex_t>::iterator i_qz = mesh_qz_.begin();
		for(int i_z = 0; i_qz != mesh_qz_.end(); ++ i_qz, ++ i_z) {
			for(int i_y = 0; i_qy != mesh_qy_.end(); ++ i_qy, ++ i_y) {
				for(int i_x = 0; i_qx != mesh_qx_.end(); ++ i_qx, ++ i_x) {
					complex_t temp = exp(mesh_qx_[i_x] * transvec[0] +
										mesh_qy_[i_y] * transvec[1] +
										mesh_qz_[i_z] * transvec[2]);
					ff[nqx_ * nqy_ * i_z + nqx_ * i_y + i_x] *= temp;
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_cylinder()


	/**
	 * random cylinders
	 */
	bool AnalyticFormFactor::compute_random_cylinders() { // for saxs
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_random_cylinders"
					<< std::endl;
		return false;
		// ...
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
	} // AnalyticFormFactor::compute_random_cylinders()


	/**
	 * horizontal cylinder
	 */
	bool AnalyticFormFactor::compute_horizontal_cylinder(shape_param_list_t& params,
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
			std::cerr << "error: both radius and height parameters required for horizontal cylinder"
						<< std::endl;
			return false;
		} // if

		// why not doing range ??? ...

		complex_t unit(0, 1);

		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					complex_t temp_qpar = sqrt(mesh_qz_[index] * mesh_qz_[index] +
												mesh_qy_[index] * mesh_qy_[index]);
					ff.push_back(2 * PI_ * h[0] * r[0] * r[0] *
								(besselj(1, temp_qpar * r[0]) / (temp_qpar * r[0])) *
								exp(unit * mesh_qz_[index] * r[0]) *
								sinc(mesh_qx_[index] * h[0] / (float_t)2.0) *
								exp(unit * (mesh_qx_[index] * transvec[0] +
											mesh_qy_[index] * transvec[1] +
											mesh_qz_[index] * transvec[2])));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_horizontal_cylinder()


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

		//std::vector <float_t> q = mat_sqrt(nqx_, nqy_, nqz_,
		//							mat_add(nqx_, nqy_ nqz_,
		//							mat_add(nqx_, nqy_, nqz_,
		//							mat_sqr(mesh_qx_), nqz_, nqy_, nqz_,
		//							mat_sqr(nqx_, nqy_, nqz_, mesh_qy_)), nqx_, nqy_, nqz_,
		//							mat_sqr(nqx_, nqy_, nqz_, mesh_qz_)));
		std::vector <complex_t> q;
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					complex_t temp_qx = mesh_qx_[index] * mesh_qx_[index];
					complex_t temp_qy = mesh_qy_[index] * mesh_qy_[index];
					complex_t temp_qz = mesh_qz_[index] * mesh_qz_[index];
					q.push_back(sqrt(temp_qx + temp_qy + temp_qz));
				} // for
			} // for
		} // for

		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0, 0));

		//std::vector <float_t> iter_r = r.begin();
		//std::vector <float_t> iter_d = distr_r.begin();
		//for(; iter_r != r.end(); ++ iter_r, ++ iter_d) {
		//	ff_ = mat_add(nqx_, nqy_, nqz_, ff_, nqx_, nqy_, nqz_,
		//			mat_mul(
		//			mat_mul(((*iter_d) * 4 * PI_ * pow((*iter_r), 3)), nqx_, nqy_, nqz_, ((
		//			mat_sub(nqx_, nqy_, nqz_,
		//			mat_sin(nqx_, nqy_, nqz_,
		//			mat_mul(nqx_, nqy_, nqz_, q, (*iter_r))), nqx_, nqy_, nqz_,
		//			mat_dot_prod(nqx_, nqy_, nqz_,
		//			mat_mul(nqx_, nqy_, nqz_,
		//			mat_dot_prod(q, (*iter_r)), nqx_, nqy_, nqz_,
		//			mat_cos(nqx_, nqy_, nqz_,
		//			mat_mul(nqx_, nqy_, nqz_, q, (*iter_r))))))) /
		//			mat_pow(nqx_, nqy_, nqz_,
		//			mat_mul(nqx_, nqy_, nqz_, q, (*iter_r)), 3))), nqx_, nqy_, nqz_,
		//			mat_exp(nqx_, nqy_, nqz_,
		//			mat_complex(nqx_, nqy_, nqz_, 0,
		//			mat_mul(nqx_, nqy_, nqz_, mesh_qz_, (*iter_r))))));
		//} // for
		std::vector<float_t>::iterator iter_r = r.begin();
		std::vector<float_t>::iterator iter_d = distr_r.begin();
		for(; iter_r != r.end(); ++ iter_r, ++ iter_d) {
			for(unsigned int z = 0; z < nqz_; ++ z) {
				for(unsigned int y = 0; y < nqy_; ++ y) {
					for(unsigned int x = 0; x < nqx_; ++ x) {
						unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
						ff[index] += (*iter_d) * 4 * PI_ * pow((*iter_r), 3) *
										((sin(q[index] * (*iter_r)) -
										(*iter_r) * q[index] * cos(q[index] * (*iter_r))) /
										(complex_t)pow(q[index] * (*iter_r), 3)) *
										exp(complex_t(0, 1) * mesh_qz_[index] * (*iter_r));
					} // for x
				} // for y
			} // for z
			for(std::vector<complex_t>::iterator iter_f = ff.begin();
					iter_f != ff.end(); ++ iter_f) {
				if((*iter_f).real() <= 1e-14) (*iter_f) = complex_t(4 * PI_ * pow((*iter_r), 3), 0);
			} // for f
		} // for r

		std::vector<complex_t>::iterator iter_f = ff.begin();
		unsigned int i = 0;
		for(; iter_f != ff.end(); ++ iter_f, ++ i) {
			(*iter_f) *= exp(complex_t(0, 1) * (mesh_qx_[i] * transvec[0] +
							mesh_qy_[i] * transvec[1] + mesh_qz_[i] * transvec[2]));
		} // for

		return true;
	} // AnalyticFormFactor::compute_sphere()


	/**
	 * prism - 3 face
	 */
	bool AnalyticFormFactor::compute_prism(shape_param_list_t& params, std::vector<complex_t>& ff,
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
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
					//float_t qm = tan(tau) * (mesh_qx_[index] * sin(eta) +
					//				mesh_qy_[index] * cos(eta));
					ff.push_back(2.0 * sqrt3 * exp(unit * mesh_qy_[index] * r[0] / sqrt3) /
								(mesh_qx_[index] * (mesh_qx_[index] * mesh_qx_[index] -
												complex_t(3.0, 0.0) * mesh_qy_[index] * mesh_qy_[index])) *
								(mesh_qx_[index] * exp(unit * mesh_qy_[index] * r[0] * sqrt3) -
								 mesh_qx_[index] * cos(mesh_qx_[index] * r[0]) -
								 unit * sqrt3 * mesh_qy_[index] * sin(mesh_qx_[index] * r[0])) *
								 fq_inv(mesh_qz_[index] + tan(tau) *
								(mesh_qx_[index] * sin(eta) +
								 mesh_qy_[index] * cos(eta)), h[0]) *
								 exp(unit * (mesh_qx_[index] * transvec[0] +
											 mesh_qy_[index] * transvec[1] +
											 mesh_qz_[index] * transvec[2])));
				} // for x
			} // for y
		} // for z

		return true;
	} // AnalyticFormFactor::compute_prism()


	/**
	 * six faceted prism
	 */
	bool AnalyticFormFactor::compute_prism6() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_prism6"
					<< std::endl;
		return false;
		// ...
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
	} // AnalyticFormFactor::compute_prism6()


	/**
	 * triangular grating in the x-direction
	 */
	bool AnalyticFormFactor::compute_prism3x() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_prism3x"
					<< std::endl;
		return false;
		// ...
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
	} // AnalyticFormFactor::compute_prism3x()


	/**
	 * upwards sawtooth
	 */
	bool AnalyticFormFactor::compute_sawtooth_up() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_up"
					<< std::endl;
		return false;
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
	} // AnalyticFormFactor::compute_sawtooth_up()


	/**
	 * downwards sawtooth
	 */
	bool AnalyticFormFactor::compute_sawtooth_down() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_down"
					<< std::endl;
		return false;
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
	} // AnalyticFormFactor::compute_sawtooth_down()


	/**
	 * pyramid
	 */
	bool AnalyticFormFactor::compute_pyramid() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_pyramid"
					<< std::endl;
		return false;
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
	} // AnalyticFormFactor::compute_pyramid()


	/**
	 * truncated cone
	 */
	bool AnalyticFormFactor::compute_truncated_cone() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_cone"
					<< std::endl;
		return false;
		// ...
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
	} // AnalyticFormFactor::compute_pyramid()


	/**
	 * truncated pyramid
	 */
	bool AnalyticFormFactor::compute_truncated_pyramid() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_pyramid"
					<< std::endl;
		return false;
		// ...
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
	} // AnalyticFormFactor::compute_truncated_pyramid()


	/**
	 * matrix computation helpers
	 */

	std::vector<complex_t>& AnalyticFormFactor::mat_fq_inv(unsigned int x_size,
															unsigned int y_size,
															unsigned int z_size,
															std::vector<complex_t>& matrix, float_t y) {
		for(std::vector<complex_t>::iterator iter = matrix.begin(); iter != matrix.end(); ++ iter) {
			complex_t temp = 2.0 * exp((*iter) * y / (float_t)2.0) *
							sin((*iter) * y / (float_t)2.0) / (*iter);
			if(fabs(temp) <= 1e-14) temp = y;
		} // for
		
		return matrix;
	} // AnalyticFormFactor::mat_fq_inv()


	complex_t AnalyticFormFactor::fq_inv(complex_t value, float_t y) {
		complex_t temp = 2.0 * exp(value * y / (float_t)2.0) * sin(value * y / (float_t)2.0) / value;
		if(fabs(temp) <= 1e-14) temp = y;
		return temp;
	} // AnalyticFormFactor::fq_inv()


	std::vector<complex_t>& AnalyticFormFactor::mat_sinc(unsigned int x_size,
													unsigned int y_size,
													unsigned int z_size,
													std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator iter = matrix.begin(); iter != matrix.end(); ++ iter) {
			complex_t temp = sin((*iter).real()) / (*iter).real();
			if(fabs(temp) <= 1e-14) temp = 1.0;
			(*iter) = temp;
		} // for

		return matrix;
	} // AnalyticFormFactor::mat_sinc()


	float_t AnalyticFormFactor::sinc(complex_t value) {
		float_t temp = sin(value.real()) / value.real();
		if(fabs(temp) <= 1e-14) temp = 1.0;
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
		if(param.nvalues() > 1) {
			float_t step = fabs(param.max() - param.min()) / (param.nvalues() - 1);
			float_t curr = param.min();
			do {
				dim.push_back(curr);
				curr += step;
			} while(curr <= param.max());	// assumes min < max ...
		} else {
			dim.push_back(param.min());
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
				dim_vals.push_back(exp(-1.0 * pow((dim[i] - mean), 2) / 2 * pow(param.deviation(), 2))
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

