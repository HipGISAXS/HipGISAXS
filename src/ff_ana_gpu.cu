/***
  *  $Id: ff_ana_gpu.cu 37 2012-08-09 22:59:59Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Thu 06 Dec 2012 04:44:09 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <complex>
#include <cuComplex.h>

#include "ff_ana_gpu.cuh"
#include "enums.hpp"
#include "cu_complex_numeric.cuh"

namespace hig {

	// make gpu kernels members of the class ...

	// forward declarations of gpu kernels
	__global__ void form_factor_sphere_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
					//cucomplex_t* mesh_qx, cucomplex_t* mesh_qy, cucomplex_t* mesh_qz,
					float_t* qx, float_t* qy, cucomplex_t* mesh_qz, float_t* rot,
					unsigned int n_r, float_t* r, unsigned int n_distr_r, float_t* distr_r,
					unsigned int n_transvec, float_t* transvec, cucomplex_t* ff);
	__device__ void ff_sphere_kernel_compute_tff(float r, float distr_r,
										cuFloatComplex q, cuFloatComplex mqz, cuFloatComplex& f);
	__device__ void ff_sphere_kernel_compute_tff(double r, double distr_r,
										cuDoubleComplex q, cuDoubleComplex mqz, cuDoubleComplex& f);
	__device__ cuFloatComplex ff_sphere_kernel_compute_ff(cuFloatComplex temp_f,
								cuFloatComplex mqx, cuFloatComplex mqy, cuFloatComplex mqz,
								float tx, float ty, float tz);
	__device__ cuDoubleComplex ff_sphere_kernel_compute_ff(cuDoubleComplex temp_f,
								cuDoubleComplex mqx, cuDoubleComplex mqy, cuDoubleComplex mqz,
								double tx, double ty, double tz);


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


	/**
	 * sphere host function
	 */
	bool AnalyticFormFactorG::compute_sphere(const std::vector<float_t>& r,
											const std::vector<float_t>& distr_r,
											const float_t* qx_h,
											const float_t* qy_h,
											const cucomplex_t* qz_h,
											const float_t* rot_h,
											const std::vector<float_t>& transvec,
											std::vector<complex_t>& ff) {
		unsigned int n_r = r.size();
		unsigned int n_distr_r = distr_r.size();
		unsigned int n_transvec = transvec.size();		// this should be = 3
		const float_t *r_h = r.empty() ? NULL : &*r.begin();
		const float_t *distr_r_h = distr_r.empty() ? NULL : &*distr_r.begin();
		const float_t *transvec_h = transvec.empty() ? NULL : &*transvec.begin();

		unsigned int grid_size = nqx_ * nqy_ * nqz_;
		//std::cout << "nqx x nqy x nqz = " << nqx_ << " x " << nqy_ << " x " << nqz_ << std::endl;

		cudaEvent_t mem_begin_e, mem_end_e;
		cudaEventCreate(&mem_begin_e);
		cudaEventCreate(&mem_end_e);
		float mem_time = 0.0, temp_time = 0.0;

		cudaEventRecord(mem_begin_e, 0);

		// construct device buffers
		float_t *qx_d, *qy_d; cucomplex_t *qz_d, *ff_d;
		float_t *r_d, *distr_r_d, *transvec_d;
		float_t *rot_d;
		if(cudaMalloc((void **) &qx_d, nqx_ * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for mesh_qx_d" << std::endl;
			return false;
		} // if
		if(cudaMalloc((void **) &qy_d, nqy_ * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for mesh_qy_d" << std::endl;
			cudaFree(qx_d);
			return false;
		} // if
		if(cudaMalloc((void **) &qz_d, nqz_ * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for mesh_qz_d" << std::endl;
			cudaFree(qy_d);
			cudaFree(qx_d);
			return false;
		} // if
		if(cudaMalloc((void **) &ff_d, grid_size * sizeof(cucomplex_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for mesh_qz_d" << std::endl;
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			return false;
		} // if
		if(cudaMalloc((void **) &r_d, n_r * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for r_d" << std::endl;
			cudaFree(ff_d);
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			return false;
		} // if
		if(cudaMalloc((void **) &distr_r_d, n_distr_r * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for n_distr_r_d" << std::endl;
			cudaFree(r_d);
			cudaFree(ff_d);
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			return false;
		} // if
		if(cudaMalloc((void **) &transvec_d, n_transvec * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for transvec_d" << std::endl;
			cudaFree(distr_r_d);
			cudaFree(r_d);
			cudaFree(ff_d);
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			return false;
		} // if
		if(cudaMalloc((void **) &rot_d, 9 * sizeof(float_t)) != cudaSuccess) {
			std::cerr << "error: device memory allocation failed for rot_d" << std::endl;
			cudaFree(transvec_d);
			cudaFree(distr_r_d);
			cudaFree(r_d);
			cudaFree(ff_d);
			cudaFree(qz_d);
			cudaFree(qy_d);
			cudaFree(qx_d);
			return false;
		} // if

		// copy data to device buffers
		cudaMemcpy(qx_d, qx_h, nqx_ * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qy_d, qy_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qz_d, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);
		cudaMemcpy(r_d, r_h, n_r * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(distr_r_d, distr_r_h, n_distr_r * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(transvec_d, transvec_h, n_transvec * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(rot_d, rot_h, 9 * sizeof(float_t), cudaMemcpyHostToDevice);

		cudaEventRecord(mem_end_e, 0);
		cudaEventSynchronize(mem_end_e);
		cudaEventElapsedTime(&temp_time, mem_begin_e, mem_end_e);
		mem_time += temp_time;

		size_t device_mem_avail, device_mem_total, device_mem_used;
		cudaMemGetInfo(&device_mem_avail, &device_mem_total);
		device_mem_used = device_mem_total - device_mem_avail;
//		if(rank == 0) {
			std::cout << "++       Used device memory: " << (float) device_mem_used / 1024 / 1024
						<< " MB" << std::endl;
			std::cout << "++       Free device memory: " << (float) device_mem_avail / 1024 / 1024
						<< " MB" << std::endl;
//		}

		//for(int cby = 2; cby < 129; cby += 2) {
		//for(int cbz = 2; cbz < 129; cbz += 2) {
		cudaEvent_t begin_e, end_e;
		cudaEventCreate(&begin_e);
		cudaEventCreate(&end_e);
		cudaEventRecord(begin_e, 0);
		// decompose computations and construct and call the kernel
		// decomposing along y and z
		// note that (cuda x y z != dim x y z)
		unsigned int cuda_block_y = 16, cuda_block_z = 6;
		unsigned int cuda_num_blocks_y = (unsigned int) ceil((float_t) nqy_ / cuda_block_y);
		unsigned int cuda_num_blocks_z = (unsigned int) ceil((float_t) nqz_ / cuda_block_z);
		dim3 ff_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
		dim3 ff_block_size(cuda_block_y, cuda_block_z, 1);

		size_t shared_mem_size = (nqx_ + cuda_block_y) * sizeof(float_t) +
									cuda_block_z * sizeof(cucomplex_t);
		if(shared_mem_size > 49152) {
			std::cerr << "Too much shared memory requested!" << std::endl;
			return false;
		} // if

		// the kernel
		form_factor_sphere_kernel <<< ff_grid_size, ff_block_size, shared_mem_size >>> (
				nqx_, nqy_, nqz_, qx_d, qy_d, qz_d, rot_d,
				n_r, r_d, n_distr_r, distr_r_d, n_transvec, transvec_d,
				ff_d);

		cudaThreadSynchronize();
		cudaError_t err = cudaGetLastError();
		if(err != cudaSuccess) {
			std::cerr << "error: form factor kernel failed [" << __FILE__ << ":" << __LINE__ << "]: "
						<< cudaGetErrorString(err) << std::endl;
		} else {
			float kernel_time;
			cudaEventRecord(end_e, 0);
			cudaEventSynchronize(end_e);
			cudaEventElapsedTime(&kernel_time, begin_e, end_e);
			//std::cout << "block size: " << cby << " x " << cbz << ". ";
			std::cout << "Analytical Sphere Kernel completed in " << kernel_time << " ms." << std::endl;

			cudaEventRecord(mem_begin_e, 0);

			cucomplex_t* ff_h = new (std::nothrow) cucomplex_t[nqx_ * nqy_ * nqz_];
			// copy result to host
			cudaMemcpy(ff_h, ff_d, nqx_ * nqy_ * nqz_ * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
			for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
				ff.push_back(complex_t(ff_h[i].x, ff_h[i].y));
			} // for
			delete[] ff_h;

			cudaEventRecord(mem_end_e, 0);
			cudaEventSynchronize(mem_end_e);
			cudaEventElapsedTime(&temp_time, mem_begin_e, mem_end_e);
			mem_time += temp_time;
		} // if-else

		std::cout << "GPU memory time: " << mem_time << " ms." << std::endl;
		//} // for cbz
		//} // for cby

		cudaFree(rot_d);
		cudaFree(transvec_d);
		cudaFree(distr_r_d);
		cudaFree(r_d);
		cudaFree(ff_d);
		cudaFree(qz_d);
		cudaFree(qy_d);
		cudaFree(qx_d);

		return true;
	} // AnalyticFormFactorG::compute_sphere()


	// sphere gpu kernel
/*	__global__ void form_factor_sphere_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
					float_t* qx, float_t* qy, cucomplex_t* qz, float_t* rot,
					unsigned int n_r, float_t* r, unsigned int n_distr_r, float_t* distr_r,
					unsigned int n_transvec, float_t* transvec, cucomplex_t* ff) {
		// decomposition is along y and z
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;
		// compute
		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				unsigned int index = base_index + i_x;
				// computing mesh values on the fly instead of storing them
				cucomplex_t temp_mqx = make_cuC(qy[i_y] * rot[0] + qx[i_x] * rot[1] + qz[i_z].x * rot[2],
												qz[i_z].y * rot[2]);
				cucomplex_t temp_mqy = make_cuC(qy[i_y] * rot[3] + qx[i_x] * rot[4] + qz[i_z].x * rot[5],
												qz[i_z].y * rot[5]);
				cucomplex_t temp_mqz = make_cuC(qy[i_y] * rot[6] + qx[i_x] * rot[7] + qz[i_z].x * rot[8],
												qz[i_z].y * rot[8]);
				cucomplex_t q = cuCnorm3(temp_mqx, temp_mqy, temp_mqz);
				cucomplex_t temp_f = make_cuC((float_t)0.0, (float_t)0.0);
				for(unsigned int i_r = 0; i_r < n_r; ++ i_r) {
					float_t temp_r = r[i_r];
					ff_sphere_kernel_compute_tff(temp_r, distr_r[i_r], q, temp_mqz, temp_f);

				} // for i_r
				ff[index] = ff_sphere_kernel_compute_ff(temp_f,	temp_mqx, temp_mqy, temp_mqz,
													transvec[0], transvec[1], transvec[2]);
			} // for x
		} // if
	} // form_factor_sphere_kernel() 
*/
	extern __shared__ float_t dynamic_shared[];

	// sphere gpu kernel
	__global__ void form_factor_sphere_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
					float_t* qx, float_t* qy, cucomplex_t* qz, float_t* rot,
					unsigned int n_r, float_t* r, unsigned int n_distr_r, float_t* distr_r,
					unsigned int n_transvec, float_t* transvec, cucomplex_t* ff) {
		// decomposition is along y and z
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int base_index = nqx * nqy * i_z + nqx * i_y;

		// shared buffers
		//float_t* qx_s = (float_t*) dynamic_shared;
		//float_t* qy_s = (float_t*) &qx_s[nqx];
		//cucomplex_t* qz_s = (cucomplex_t*) &qy_s[blockDim.x];
		cucomplex_t* qz_s = (cucomplex_t*) dynamic_shared;
		float_t* qx_s = (float_t*) &qz_s[blockDim.y];
		float_t* qy_s = (float_t*) &qx_s[nqx];

		// load all qx
		unsigned int i_thread = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y;
		unsigned int num_loads = ceil((float_t) nqx / num_threads);
		for(int i = 0; i < num_loads; ++ i) {
			unsigned int index = i * num_threads + i_thread;
			if(index < nqx) qx_s[index] = qx[index];
			else ;	// nop
		} // for
		// load part of qy
		if(i_y < nqy && threadIdx.y == 0)	// first row of threads
			qy_s[threadIdx.x] = qy[i_y];
		// load part of qz
		if(i_z < nqz && threadIdx.x == 0)	// first column of threads
			qz_s[threadIdx.y] = qz[i_z];

		__syncthreads();

		// compute
		if(i_y < nqy && i_z < nqz) {
			for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
				// computing mesh values on the fly instead of storing them
				cucomplex_t temp_mqx = make_cuC(qy_s[threadIdx.x] * rot[0] + qx_s[i_x] * rot[1] +
												qz_s[threadIdx.y].x * rot[2],
												qz_s[threadIdx.y].y * rot[2]);
				cucomplex_t temp_mqy = make_cuC(qy_s[threadIdx.x] * rot[3] + qx_s[i_x] * rot[4] +
												qz_s[threadIdx.y].x * rot[5],
												qz_s[threadIdx.y].y * rot[5]);
				cucomplex_t temp_mqz = make_cuC(qy_s[threadIdx.x] * rot[6] + qx_s[i_x] * rot[7] +
												qz_s[threadIdx.y].x * rot[8],
												qz_s[threadIdx.y].y * rot[8]);
				cucomplex_t q = cuCnorm3(temp_mqx, temp_mqy, temp_mqz);
				cucomplex_t temp_f = make_cuC((float_t)0.0, (float_t)0.0);
				for(unsigned int i_r = 0; i_r < n_r; ++ i_r) {
					float_t temp_r = r[i_r];
					float_t temp_distr_r = distr_r[i_r];
					ff_sphere_kernel_compute_tff(temp_r, temp_distr_r, q, temp_mqz, temp_f);
				} // for i_r
				ff[base_index + i_x] = ff_sphere_kernel_compute_ff(temp_f,	temp_mqx, temp_mqy, temp_mqz,
													transvec[0], transvec[1], transvec[2]);
			} // for x
		} // if
	} // form_factor_sphere_kernel()

	__device__ void ff_sphere_kernel_compute_tff(float r, float distr_r,
										cuFloatComplex q, cuFloatComplex mqz, cuFloatComplex &f) {
		cuFloatComplex temp1 = make_cuFloatComplex(q.x * r, q.y * r);
		cuFloatComplex temp2 = cuCsubf(cuCsin(temp1), cuCmulf(temp1, cuCcos(temp1)));
		cuFloatComplex temp3 = cuCmulf(temp1, cuCmulf(temp1, temp1));
		float temp4 = distr_r * 4 * PI_ * r * r * r;
		cuFloatComplex temp5 = cuCmulf(cuCdivf(temp2, temp3),
									cuCexp(make_cuFloatComplex(-r * mqz.y, r * mqz.x)));
		cuFloatComplex tempff = make_cuFloatComplex(temp4 * temp5.x, temp4 * temp5.y);
		f = cuCaddf(f, tempff);
	} // ff_sphere_kernel_compute_tff()

	__device__ void ff_sphere_kernel_compute_tff(double r, double distr_r,
										cuDoubleComplex q, cuDoubleComplex mqz, cuDoubleComplex &f) {
		cuDoubleComplex temp1 = cuCmul(q, make_cuDoubleComplex(r, 0.0));
		cuDoubleComplex temp2 = cuCsub(cuCsin(temp1), cuCmul(temp1, cuCcos(temp1)));
		cuDoubleComplex temp3 = cuCmul(temp1, cuCmul(temp1, temp1));
		double temp4 = distr_r * 4 * PI_ * r * r * r;
		cuDoubleComplex temp5 = cuCmul(cuCdiv(temp2, temp3),
									cuCexp(make_cuDoubleComplex(-r * mqz.y, r * mqz.x)));
		cuDoubleComplex tempff = make_cuDoubleComplex(temp4 * temp5.x, temp4 * temp5.y);
		f = cuCadd(f, tempff);
	} // ff_sphere_kernel_compute_tff()

	__device__ cuFloatComplex ff_sphere_kernel_compute_ff(cuFloatComplex temp_f,
								cuFloatComplex mqx, cuFloatComplex mqy, cuFloatComplex mqz,
								float tx, float ty, float tz) {
			float rl = tx * mqx.x + ty * mqy.x + tz * mqz.x;
			float im = tx * mqx.y + ty * mqy.y + tz * mqz.y;
			cuFloatComplex temp1 = cuCexp(make_cuFloatComplex(-im, rl));
			return cuCmulf(temp_f, temp1);
	} // ff_sphere_kernel_compute_ff()

	__device__ cuDoubleComplex ff_sphere_kernel_compute_ff(cuDoubleComplex temp_f,
								cuDoubleComplex mqx, cuDoubleComplex mqy, cuDoubleComplex mqz,
								double tx, double ty, double tz) {
			double rl = tx * mqx.x + ty * mqy.x + tz * mqz.x;
			double im = tx * mqx.y + ty * mqy.y + tz * mqz.y;
			cuDoubleComplex temp1 = cuCexp(make_cuDoubleComplex(-im, rl));
			return cuCmul(temp_f, temp1);
	} // ff_sphere_kernel_compute_ff()


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

