/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Tue 19 Feb 2013 11:45:01 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <complex>
#include <cuComplex.h>

#include "ff_ana_gpu.cuh"
#include "enums.hpp"

namespace hig {

	AnalyticFormFactorG::AnalyticFormFactorG(unsigned int nx, unsigned int ny, unsigned int nz) {
		nqx_ = nx; nqy_ = ny;  nqz_ = nz;
	} // AnalyticFormFactorG::AnalyticFormFactorG();


	AnalyticFormFactorG::AnalyticFormFactorG() { }


	AnalyticFormFactorG::~AnalyticFormFactorG() { }


	void AnalyticFormFactorG::grid_size(unsigned int nx, unsigned int ny, unsigned int nz) {
		nqx_ = nx; nqy_ = ny;  nqz_ = nz;
	} // AnalyticFormFactorG::grid_size()


//	/**
//	 * box on gpu
//	 */
//	bool AnalyticFormFactorG::compute_box_gpu(unsigned int nqx, unsigned int nqy, unsigned int nqz,
//										std::vector<complex_t>& ff,
//										ShapeName shape, shape_param_list_t& params,
//										float_t tau, float_t eta, std::vector<float_t> &transvec,
//										std::vector<float_t> &rot1, std::vector<float_t> &rot2, std::vector<float_t> &rot3) {
//		std::vector <float_t> x, distr_x;	// for x dimension: param_xsize  param_edge
//		std::vector <float_t> y, distr_y;	// for y dimension: param_ysize  param_edge
//		std::vector <float_t> z, distr_z;	// for z dimension: param_height param_edge
//		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			if(!(*i).second.isvalid()) {
//				std::cerr << "warning: invalid shape parameter found" << std::endl;
//				continue;
//			} // if
//			switch((*i).second.type()) {
//				case param_edge:
//					param_distribution((*i).second, x, distr_x);	// x == RR, distr_x == RRD
//					param_distribution((*i).second, y, distr_y);
//					param_distribution((*i).second, z, distr_z);
//					break;
//				case param_xsize:
//					param_distribution((*i).second, x, distr_x);
//					break;
//				case param_ysize:
//					param_distribution((*i).second, y, distr_y);
//					break;
//				case param_height:
//					param_distribution((*i).second, z, distr_z);
//					break;
//				case param_radius:
//				case param_baseangle:
//					std::cerr << "warning: ignoring unwanted values for shape type 'box'" << std::endl;
//					break;
//				default:
//					std::cerr << "warning: unknown parameters for shape type 'box'. ignoring"
//								<< std::endl;
//			} // switch
//		} // for
//
//		// check if x y z etc are set or not
//		if(x.size() < 1 || y.size() < 1 || z.size() < 1) {
//			std::cerr << "error: invalid or not enough box parameters given" << std::endl;
//			return false;
//		} // if
//
//		//std::vector<complex_t> mesh_qm = mat_mul(tan(tau),
//		//										mat_add(nqx, nqy, nqz,
//		//										mat_mul(mesh_qx_, sin(eta)),
//		//										nqx, nqy, nqz,
//		//										mat_mul(mesh_qy_, cos(eta))));
//		complex_vec_t mesh_qm;
//		// temp vars can be resuced here ...
//		complex_vec_t temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10;
//		mat_mul(mesh_qx_, sin(eta), temp1);
//		mat_mul(mesh_qy_, cos(eta), temp2);
//		mat_add(nqx, nqy, nqz, temp1, nqx, nqy, nqz, temp2, temp3);
//		mat_mul(tan(tau), temp3, mesh_qm);
//
//		// initialize ff
//		ff.clear();
//		ff.reserve(nqz * nqy * nqx);
//		for(unsigned int i = 0; i < nqz * nqy * nqx; ++ i) ff.push_back(complex_t(0, 0));
//		// ff computation for a box
//		for(unsigned int i_z = 0; i_z < z.size(); ++ i_z) {
//			for(unsigned int i_y = 0; i_y < y.size(); ++ i_y) {
//				for(unsigned int i_x = 0; i_x < x.size(); ++ i_x) {
//					/*mat_add(nqx, nqy, nqz, mesh_qz_, nqx, nqy, nqz, mesh_qm, temp1);
//					mat_mul(mesh_qy_, z[i_z], temp2);
//					mat_mul(mesh_qx_, x[i_x], temp3);
//					mat_fq_inv(nqx, nqy, nqz, temp1, y[i_y], temp4);
//					mat_sinc(nqx, nqy, nqz, temp2, temp5);
//					mat_sinc(nqx, nqy, nqz, temp3, temp6);
//					mat_dot_prod(nqx, nqy, nqz, temp6, nqx, nqy, nqz, temp5, temp7);
//					mat_dot_prod(nqx, nqy, nqz, temp7, nqx, nqy, nqz, temp4, temp8);
//					complex_t temp9 = distr_x[i_x] * distr_y[i_y] * distr_z[i_z] * 4 * z[i_z] * x[i_x];
//					mat_mul(temp9, temp8, temp10);
//					mat_add_in(nqx, nqy, nqz, ff, nqx, nqy, nqz, temp10); */
//					/*ff = mat_add(nqx, nqy, nqz, ff, nqx, nqy, nqz,
//							mat_mul(distr_x[i_x] * distr_y[i_y] * distr_z[i_z] * 4 * z[i_z] * x[i_x],
//							mat_dot_prod(nqx, nqy, nqz,
//							mat_dot_prod(nqx, nqy, nqz,
//							mat_sinc(nqx, nqy, nqz, mat_mul(mesh_qx_, x[i_x])),
//							nqx, nqy, nqz, mat_sinc(nqx, nqy, nqz,
//							mat_mul(mesh_qz_, z[i_z]))),
//							nqx, nqy, nqz,
//							mat_fq_inv(nqx, nqy, nqz,
//							mat_add(nqx, nqy, nqz, mesh_qz_, nqx, nqy, nqz, mesh_qm),
//							y[i_y])))); */
//					for(unsigned int j_z = 0; j_z < nqz; ++ j_z) {
//						for(unsigned int j_y = 0; j_y < nqy; ++ j_y) {
//							for(unsigned int j_x = 0; j_x < nqx; ++ j_x) {
//								unsigned int curr_index = nqx * nqy * j_z + nqx * j_y + j_x;
//								complex_t temp1 = mesh_qz_[curr_index] + mesh_qm[curr_index];
//								complex_t temp2 = mesh_qy_[curr_index] * z[i_z];
//								complex_t temp3 = mesh_qx_[curr_index] * x[i_x];
//								complex_t temp4 = fq_inv(temp1, y[i_y]);
//								complex_t temp5 = sinc(temp2);
//								complex_t temp6 = sinc(temp3);
//								complex_t temp7 = temp6 * temp5;
//								complex_t temp8 = temp7 * temp4;
//								complex_t temp9 = 4 * distr_x[i_x] * distr_y[i_y] * distr_z[i_z] *
//													z[i_z] * x[i_x];
//								complex_t temp10 = temp9 * temp8;
//								if(!(boost::math::isfinite(temp10.real()) &&
//											boost::math::isfinite(temp10.imag()))) {
//									std::cerr << "+++++++++++++++ here it is +++++++ " << j_x << ", "
//												<< j_y << ", " << j_z << std::endl;
//									exit(1);
//								} // if
//								ff[curr_index] += temp10;
//							} // for x
//						} // for y
//					} // for z
//				} // for i_x
//			} // for i_y
//		} // for i_z
//
//		for(unsigned int j_z = 0; j_z < nqz; ++ j_z) {
//			for(unsigned int j_y = 0; j_y < nqy; ++ j_y) {
//				for(unsigned int j_x = 0; j_x < nqx; ++ j_x) {
//					unsigned int curr_index = nqx * nqy * j_z + nqx * j_y + j_x;
//					complex_t temp = exp(mesh_qx_[curr_index] * transvec[0] +
//										mesh_qy_[curr_index] * transvec[1] +
//										mesh_qz_[curr_index] * transvec[2]);
//					if(!(boost::math::isfinite(temp.real()) &&
//								boost::math::isfinite(temp.imag()))) {
//						std::cerr << "---------------- here it is ------ " << j_x << ", "
//									<< j_y << ", " << j_z << std::endl;
//						exit(1);
//					} // if
//					ff[curr_index] *= temp;
//				} // for x
//			} // for y
//		} // for z
//
//		return true;
//	} // AnalyticFormFactorG::compute_box()
//
//
//	/**
//	 * cylinder
//	 */
//	bool AnalyticFormFactorG::compute_cylinder_gpu(shape_param_list_t& params, float_t tau, float_t eta,
//			std::vector<complex_t>& ff, std::vector<float_t> transvec) {
//		std::vector <float_t> h, distr_h;	// for h dimension: param_height
//		std::vector <float_t> r, distr_r;	// for r dimension: param_radius
//		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			if(!(*i).second.isvalid()) {
//				std::cerr << "warning: ignoring invalid shape parameter" << std::endl;
//				continue;
//			} // if
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_baseangle:
//					std::cerr << "warning: ignoring unwanted input values for cylinder" << std::endl;
//					break;
//				case param_height:
//					param_distribution((*i).second, h, distr_h);
//					break;
//				case param_radius:
//					param_distribution((*i).second, r, distr_r);
//					break;
//				default:
//					std::cerr << "error: unknown parameter type given" << std::endl;
//					return false;
//			} // switch
//		} // for
//
//		if(h.size() < 1 || r.size() < 1) {
//			std::cerr << "error: missing parameters for cylinder" << std::endl;
//			return false;
//		} // if
//
//		ff.clear();
//		ff.reserve(nqx_ * nqy_ * nqz_);
//		for(int i_z = 0; i_z < nqz_; ++ i_z) {
//			for(int i_y = 0; i_y < nqy_; ++ i_y) {
//				for(int i_x = 0; i_x < nqx_; ++ i_x) {
//					ff.push_back(complex_t(0.0, 0.0));
//				} // for x
//			} // for y
//		} // for z
//
//		std::vector<complex_t> qpar, qm;
//		std::vector<complex_t> temp1, temp2, temp3, temp4, temp5, temp6, temp7;
//
//		// unroll the following into combined loops
//		// ...
//
//		mat_sqr(mesh_qx_, temp1);
//		mat_sqr(mesh_qy_, temp2);
//		mat_add(nqx_, nqy_, nqz_, temp1, nqx_, nqy_, nqz_, temp2, qpar);
//		mat_sqrt_in(qpar);
//		mat_mul(mesh_qx_, sin(eta), temp3);
//		mat_mul(mesh_qy_, cos(eta), temp4);
//		mat_add(nqx_, nqy_, nqz_, temp3, nqx_, nqy_, nqz_, temp4, temp5);
//		mat_mul(tan(tau), temp5, qm);
//
//		temp1.clear(); temp2.clear(); temp3.clear(); temp4.clear(); temp5.clear();
//
//		for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
//			for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
//				mat_add(nqx_, nqy_, nqz_, mesh_qz_, nqx_, nqy_, nqz_, qm, temp1);
//				mat_fq_inv(nqx_, nqy_, nqz_, temp1, h[i_h], temp2);
//				mat_mul(qpar, r[i_r], temp3);
//				mat_besselj(1, nqx_, nqy_, nqz_, temp3, temp4);
//				mat_dot_div(nqx_, nqy_, nqz_, temp4, nqx_, nqy_, nqz_, temp3, temp5);
//				mat_dot_prod(nqx_, nqy_, nqz_, temp5, nqx_, nqy_, nqz_, temp2, temp6);
//				mat_mul(distr_r[i_r] * distr_h[i_h] * 2.0 * PI_ * pow(r[i_r], 2), temp6, temp7);
//				mat_add_in(nqx_, nqy_, nqz_, ff, nqx_, nqy_, nqz_, temp7);
//			} // for h
//		} // for r
//
//		temp1.clear(); temp2.clear(); temp3.clear(); temp4.clear();
//		temp5.clear(); temp6.clear(); temp7.clear();
//
//		mat_mul(mesh_qx_, transvec[0], temp1);
//		mat_mul(mesh_qy_, transvec[1], temp2);
//		mat_mul(mesh_qz_, transvec[2], temp3);
//		mat_add_in(nqx_, nqy_, nqz_, temp2, nqx_, nqy_, nqz_, temp3);
//		mat_add_in(nqx_, nqy_, nqz_, temp1, nqx_, nqy_, nqz_, temp2);
//		mat_mul_in(complex_t(0, 1), temp1);
//		mat_exp_in(temp1);
//		mat_dot_prod_in(nqx_, nqy_, nqz_, ff, nqx_, nqy_, nqz_, temp1);
//
//		return true;
//	} // AnalyticFormFactorG::compute_cylinder()
//
//
//	/**
//	 * random cylinders
//	 */
//	bool AnalyticFormFactorG::compute_random_cylinders_gpu() { // for saxs
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_random_cylinders"
//					<< std::endl;
//		return false;
//		// ...
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_random_cylinders()
//
//
//	/**
//	 * horizontal cylinder
//	 */
//	bool AnalyticFormFactorG::compute_horizontal_cylinder_gpu(shape_param_list_t& params,
//														std::vector<float_t> transvec,
//														std::vector<complex_t>& ff) {
//		std::vector<float_t> r, distr_r;
//		std::vector<float_t> h, distr_h;
//		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_baseangle:
//					std::cerr << "warning: ignoring unwanted parameters in horizontal cylinder"
//								<< std::endl;
//					break;
//				case param_height:
//					param_distribution((*i).second, h, distr_h);
//					break;
//				case param_radius:
//					param_distribution((*i).second, r, distr_r);
//					break;
//				default:
//					std::cerr << "error: unknown or invalid parameter given for horizontal cylinder"
//								<< std::endl;
//					return false;
//			} // switch
//		} // for
//
//		if(r.size() < 1 || h.size() < 1) {
//			std::cerr << "error: both radius and height parameters are required for horizontal cylinder"
//						<< std::endl;
//			return false;
//		} // if
//
//		// in slims code, why not doing range of r and h ???
//
//		complex_t unitc(0, 1);
//
//		ff.clear();
//		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));
//
//		for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
//			for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
//				for(unsigned int z = 0; z < nqz_; ++ z) {
//					for(unsigned int y = 0; y < nqy_; ++ y) {
//						for(unsigned int x = 0; x < nqx_; ++ x) {
//							unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
//							complex_t temp_qpar = sqrt(mesh_qz_[index] * mesh_qz_[index] +
//														mesh_qy_[index] * mesh_qy_[index]);
//							ff[index] += 2 * PI_ * h[i_h] * r[i_r] * r[i_r] *
//										(cbessj(temp_qpar * r[i_r], 1) / (temp_qpar * r[i_r])) *
//										exp(unitc * mesh_qz_[index] * r[i_r]) *
//										sinc(mesh_qx_[index] * h[i_h] / (float_t)2.0) *
//										exp(unitc * (mesh_qx_[index] * transvec[0] +
//													mesh_qy_[index] * transvec[1] +
//													mesh_qz_[index] * transvec[2]));
//						} // for x
//					} // for y
//				} // for z
//			} // for h
//		} // for r
//
//		return true;
//	} // AnalyticFormFactorG::compute_horizontal_cylinder()


//	/**
//	 * prism - 3 face
//	 */
//	bool AnalyticFormFactorG::compute_prism_gpu(shape_param_list_t& params, std::vector<complex_t>& ff,
//											float_t tau, float_t eta, std::vector<float_t> transvec) {
//		std::vector<float_t> r, distr_r;
//		std::vector<float_t> h, distr_h;
//		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_baseangle:
//					std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
//					break;
//				case param_height:
//					param_distribution((*i).second, h, distr_h);
//					break;
//				case param_radius:
//					param_distribution((*i).second, r, distr_r);
//					break;
//				default:
//					std::cerr << "error: invalid parameter given for prism shape" << std::endl;
//					return false;
//			} // switch
//		} // for
//
//		if(h.size() < 1 || r.size() < 1) {
//			std::cerr << "error: need radius and height for prism shape" << std::endl;
//			return false;
//		} // if
//
//		// why not doing for a range of r, h for this? ... ???
//
//		float_t sqrt3 = sqrt(3.0);
//		complex_t unit(0, 1.0);
//
//		ff.clear();
//		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));
//
//		for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
//			for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
//				for(unsigned int z = 0; z < nqz_; ++ z) {
//					for(unsigned int y = 0; y < nqy_; ++ y) {
//						for(unsigned int x = 0; x < nqx_; ++ x) {
//							unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
//							//float_t qm = tan(tau) * (mesh_qx_[index] * sin(eta) +
//							//				mesh_qy_[index] * cos(eta));
//							ff[index] += 2.0 * sqrt3 * exp(unit * mesh_qy_[index] * r[i_r] / sqrt3) /
//										(mesh_qx_[index] * (mesh_qx_[index] * mesh_qx_[index] -
//										 complex_t(3.0, 0.0) * mesh_qy_[index] * mesh_qy_[index])) *
//										(mesh_qx_[index] * exp(unit * mesh_qy_[index] * r[i_r] * sqrt3) -
//										 mesh_qx_[index] * cos(mesh_qx_[index] * r[i_r]) -
//										 unit * sqrt3 * mesh_qy_[index] * sin(mesh_qx_[index] * r[i_r])) *
//										 fq_inv(mesh_qz_[index] + tan(tau) *
//										(mesh_qx_[index] * sin(eta) +
//										 mesh_qy_[index] * cos(eta)), h[i_h]) *
//										 exp(unit * (mesh_qx_[index] * transvec[0] +
//													 mesh_qy_[index] * transvec[1] +
//													 mesh_qz_[index] * transvec[2]));
//						} // for x
//					} // for y
//				} // for z
//			} // for h
//		} // for r
//
//		return true;
//	} // AnalyticFormFactorG::compute_prism()
//
//
//	/**
//	 * six faceted prism
//	 */
//	bool AnalyticFormFactorG::compute_prism6_gpu(shape_param_list_t& params, std::vector<complex_t>& ff,
//											float_t tau, float_t eta, std::vector<float_t> transvec) {
//		std::vector<float_t> r, distr_r;
//		std::vector<float_t> h, distr_h;
//		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_baseangle:
//					std::cerr << "warning: ignoring unwanted parameters in sphere" << std::endl;
//					break;
//				case param_height:
//					param_distribution((*i).second, h, distr_h);
//					break;
//				case param_radius:
//					param_distribution((*i).second, r, distr_r);
//					break;
//				default:
//					std::cerr << "error: invalid parameter given for prism shape" << std::endl;
//					return false;
//			} // switch
//		} // for
//
//		// 4*sq3./( 3*qy.^2 - qx.^2 ) .* (R^2 * qy.^2 .*  SINC_Matrix(qx*R/sq3) .* SINC_Matrix(qy*R) + cos(2*qx*R/sq3) - cos(qy*R) .* cos(qx*R/sq3) )  .* Fq_Inv_Matrix(qz+qm, H) .* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3)))
//
//		float_t sqrt3 = sqrt(3.0);
//		complex_t unit(0, 1.0);
//
//		ff.clear();
//		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));
//
//		for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
//			for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
//				for(unsigned int z = 0; z < nqz_; ++ z) {
//					for(unsigned int y = 0; y < nqy_; ++ y) {
//						for(unsigned int x = 0; x < nqx_; ++ x) {
//							unsigned int index = nqx_ * nqy_ * z + nqx_ * y + x;
//							complex_t qm = tan(tau) * (mesh_qx_[index] * sin(eta) +
//														mesh_qy_[index] * cos(eta));
//							complex_t temp1 = ((float_t)4.0 * sqrt3) /
//												(3.0 * mesh_qy_[index] * mesh_qy_[index] -
//												 mesh_qx_[index] * mesh_qx_[index]);
//							complex_t temp2 = r[i_r] * r[i_r] * mesh_qy_[index] * mesh_qy_[index] *
//												sinc(mesh_qx_[index] * r[i_r] / sqrt3) *
//												sinc(mesh_qy_[index] * r[i_r]);
//							complex_t temp3 = cos(2.0 * mesh_qx_[index] * r[i_r] / sqrt3);
//							complex_t temp4 = cos(mesh_qy_[index] * r[i_r]) *
//												cos(mesh_qx_[index] * r[i_r] / sqrt3);
//							complex_t temp5 = temp1 * (temp2 + temp3 - temp4);
//							complex_t temp6 = fq_inv(mesh_qz_[index] + qm, h[i_h]);
//							complex_t temp7 = (mesh_qx_[index] * transvec[0] +
//												mesh_qy_[index] * transvec[1] +
//												mesh_qz_[index] * transvec[2]);
//							ff[index] += temp5 * temp6 * exp(unit * temp7);
//						} // for x
//					} // for y
//				} // for z
//			} // for h
//		} // for r
//
//		return true;
//	} // AnalyticFormFactorG::compute_prism6()
//
//
//	/**
//	 * triangular grating in the x-direction
//	 */
//	bool AnalyticFormFactorG::compute_prism3x_gpu() {
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_prism3x"
//					<< std::endl;
//		return false;
//		// ...
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_prism3x()
//
//
//	/**
//	 * upwards sawtooth
//	 */
//	bool AnalyticFormFactorG::compute_sawtooth_up_gpu() {
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_up"
//					<< std::endl;
//		return false;
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_sawtooth_up()
//
//
//	/**
//	 * downwards sawtooth
//	 */
//	bool AnalyticFormFactorG::compute_sawtooth_down_gpu() {
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_down"
//					<< std::endl;
//		return false;
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_sawtooth_down()
//
//
//	/**
//	 * pyramid
//	 */
//	bool AnalyticFormFactorG::compute_pyramid_gpu() {
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_pyramid"
//					<< std::endl;
//		return false;
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_pyramid()
//
//
//	/**
//	 * truncated cone
//	 */
//	bool AnalyticFormFactorG::compute_truncated_cone_gpu() {
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_cone"
//					<< std::endl;
//		return false;
//		// ...
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_pyramid()
//
//
//	/**
//	 * truncated pyramid
//	 */
//	bool AnalyticFormFactorG::compute_truncated_pyramid_gpu() {
//		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_pyramid"
//					<< std::endl;
//		return false;
//		// ...
//		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
//			switch((*i).second.type()) {
//				case param_edge:
//				case param_xsize:
//				case param_ysize:
//				case param_height:
//				case param_radius:
//				case param_baseangle:
//				default:
//			} // switch
//		} // for */
//	} // AnalyticFormFactorG::compute_truncated_pyramid()

} // namespace hig

