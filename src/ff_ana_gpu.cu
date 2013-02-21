/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Thu 21 Feb 2013 10:37:17 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <complex>
#include <cuComplex.h>

#include "ff_ana_gpu.cuh"
#include "qgrid.hpp"
#include "enums.hpp"

namespace hig {

	AnalyticFormFactorG::AnalyticFormFactorG(unsigned int nx, unsigned int ny, unsigned int nz):
		nqx_(nx), nqy_(ny), nqz_(nz),
		qx_(NULL), qy_(NULL), qz_(NULL), ff_(NULL),
		transvec_(NULL), rot_(NULL)	{
	} // AnalyticFormFactorG::AnalyticFormFactorG();


	AnalyticFormFactorG::AnalyticFormFactorG():
		nqx_(0), nqy_(0), nqz_(0),
		qx_(NULL), qy_(NULL), qz_(NULL), ff_(NULL),
		transvec_(NULL), rot_(NULL)	{
	} // AnalyticFormFactorG::AnalyticFormFactorG()


	AnalyticFormFactorG::~AnalyticFormFactorG() {
		// release device memories
		destroy();
	} // AnalyticFormFactorG::~AnalyticFormFactorG()


	void AnalyticFormFactorG::grid_size(unsigned int nx, unsigned int ny, unsigned int nz) {
		nqx_ = nx; nqy_ = ny;  nqz_ = nz;
	} // AnalyticFormFactorG::grid_size()


	bool AnalyticFormFactorG::init(unsigned int nqx, unsigned int nqy, unsigned int nqz) {
		// this does the following:
		// 	+ allocate device buffers
		// 	+ copy qgrid to device memory
		nqx_ = nqx; nqy_ = nqy; nqz_ = nqz;
		cudaMalloc((void**) &qx_, nqx_ * sizeof(float_t));
		cudaMalloc((void**) &qy_, nqy_ * sizeof(float_t));
		cudaMalloc((void**) &qz_, nqz_ * sizeof(cucomplex_t));
		cudaMalloc((void**) &ff_, nqx_ * nqy_ * nqz_ * sizeof(cucomplex_t));
		cudaMalloc((void**) &transvec_, 3 * sizeof(float_t));	// check this ...
		cudaMalloc((void**) &rot_, 9 * sizeof(float_t));
		if(qx_ == NULL || qy_ == NULL || qz_ == NULL || ff_ == NULL || transvec_ == NULL || rot_ == NULL) {
			std::cerr << "error: device memory allocation failed" << std::endl;
			return false;
		} // if

		// first need to construct host buffers
		float_t* qx_h = new (std::nothrow) float_t[nqx_];
		float_t* qy_h = new (std::nothrow) float_t[nqy_];
		cucomplex_t* qz_h = new (std::nothrow) cucomplex_t[nqz_];
		if(qx_h == NULL || qy_h == NULL || qz_h == NULL) {
			std::cerr << "error: memory allocation for host mesh grid failed" << std::endl;
			return false;
		} // if
		for(unsigned int ix = 0; ix < nqx_; ++ ix) qx_h[ix] = QGrid::instance().qx(ix);
		for(unsigned int iy = 0; iy < nqy_; ++ iy) qy_h[iy] = QGrid::instance().qy(iy);
		for(unsigned int iz = 0; iz < nqz_; ++ iz) {
			qz_h[iz].x = QGrid::instance().qz_extended(iz).real();
			qz_h[iz].y = QGrid::instance().qz_extended(iz).imag();
		} // for qz

		cudaMemcpy(qx_, qx_h, nqx_ * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qy_, qy_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(qz_, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

		size_t device_mem_avail, device_mem_total, device_mem_used;
		cudaMemGetInfo(&device_mem_avail, &device_mem_total);
		device_mem_used = device_mem_total - device_mem_avail;
		std::cout << "++            Used device memory: " << (float) device_mem_used / 1024 / 1024
					<< " MB" << std::endl;
		std::cout << "++            Free device memory: " << (float) device_mem_avail / 1024 / 1024
					<< " MB" << std::endl;

		return true;
	} // AnalyticFormFactorG::init()


	bool AnalyticFormFactorG::run_init(const float_t* rot_h, const std::vector<float_t>& transvec) {
		// this does the following:
		// 	+ copy transvec and rotation matrices to device memory
		if(transvec_ == NULL || rot_ == NULL) {
			std::cerr << "error: AnalyticFormFactorG is not initialized" << std::endl;
			return false;
		} // if
		const float_t *transvec_h = transvec.empty() ? NULL : &*transvec.begin();

		cudaMemcpy(transvec_, transvec_h, 3 * sizeof(float_t), cudaMemcpyHostToDevice);
		cudaMemcpy(rot_, rot_h, 9 * sizeof(float_t), cudaMemcpyHostToDevice);

		return true;
	} // AnalyticFormFactorG::run_init()


	// release all resources and init to 0
	bool AnalyticFormFactorG::destroy() {
		if(rot_ != NULL) cudaFree(rot_);
		if(transvec_ != NULL) cudaFree(transvec_);
		if(ff_ != NULL) cudaFree(ff_);
		if(qz_ != NULL) cudaFree(qz_);
		if(qy_ != NULL) cudaFree(qy_);
		if(qx_ != NULL) cudaFree(qx_);

		rot_ = NULL; transvec_ = NULL; ff_ = NULL;
		qz_ = NULL; qy_ = NULL; qx_ = NULL;
		nqx_ = 0; nqy_ = 0; nqz_ = 0;

		return true;
	} // AnalyticFormFactorG::destroy()


	bool AnalyticFormFactorG::construct_output_ff(std::vector<complex_t>& ff) {
		unsigned int grid_size = nqx_ * nqy_ * nqz_;
		cucomplex_t* ff_h = new (std::nothrow) cucomplex_t[grid_size];
		cudaMemcpy(ff_h, ff_, grid_size * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
		ff.clear(); ff.reserve(grid_size);
		for(unsigned int i = 0; i < grid_size; ++ i) ff.push_back(complex_t(ff_h[i].x, ff_h[i].y));

		delete[] ff_h;
		return true;
	} // AnalyticFormFactorG::construct_output_ff()


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

