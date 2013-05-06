/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_gpu.cu
  *  Created: Oct 16, 2012
  *  Modified: Mon 06 May 2013 10:37:57 AM PDT
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
		cudaError_t err = cudaGetLastError();
		//std::cerr << "AFTER MALLOCS: " << cudaGetErrorString(err) << std::endl;

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


	bool AnalyticFormFactorG::clear() {
		return destroy();
	} // AnalyticFormFactorG::clear()


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


} // namespace hig

