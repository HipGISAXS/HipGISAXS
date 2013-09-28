/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_gpu.cu
 *  Created: Oct 16, 2012
 *  Modified: Thu 26 Sep 2013 10:22:14 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <complex>
#include <cuComplex.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include "ff_ana_gpu.cuh"
#include "../../model/qgrid.hpp"
#include "../../common/enums.hpp"

namespace hig {

	AnalyticFormFactorG::AnalyticFormFactorG(unsigned int nx, unsigned int ny, unsigned int nz):
		nqx_(nx), nqy_(ny), nqz_(nz),
		qx_(NULL), qy_(NULL), qz_(NULL), ff_(NULL),
		ff_buff_d_(NULL), ff_buff_h_(NULL),
		transvec_(NULL), rot_(NULL)	{
	} // AnalyticFormFactorG::AnalyticFormFactorG();


	AnalyticFormFactorG::AnalyticFormFactorG():
		nqx_(0), nqy_(0), nqz_(0),
		qx_(NULL), qy_(NULL), qz_(NULL), ff_(NULL),
		ff_buff_d_(NULL), ff_buff_h_(NULL),
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

		unsigned int est_device_mem_need = nqx_ * nqy_ * nqz_ * sizeof(cucomplex_t);
		compute_hyperblock_size(est_device_mem_need, device_mem_avail);

		nb_x_ = (unsigned int) ceil((float) nqx_ / b_nqx_);
		nb_y_ = (unsigned int) ceil((float) nqy_ / b_nqy_);
		nb_z_ = (unsigned int) ceil((float) nqz_ / b_nqz_);
		unsigned int num_hblocks = nb_x_ * nb_y_ * nb_z_;
		std::cout << "++         Number of hyperblocks: " << num_hblocks
					<< " [" << nb_x_ << " x " << nb_y_ << " x " << nb_z_ << "]" << std::endl;

		// allocate ff buffer memories
		cudaMalloc((void**) &ff_buff_d_, b_nqx_ * b_nqy_ * b_nqz_ * sizeof(cucomplex_t));
		cudaMemset(ff_buff_d_, 0, b_nqx_ * b_nqy_ * b_nqz_ * sizeof(cucomplex_t));
		ff_buff_h_ = new (std::nothrow) cucomplex_t[b_nqx_ * b_nqy_ * b_nqz_];

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


	bool AnalyticFormFactorG::run_init(const float_t* rot_h, const std::vector<float_t>& transvec,
										std::vector<complex_t>& ff) {
		complex_t zero(0, 0);
		ff.resize(nqx_ * nqy_ * nqz_);
		//std::fill(ff.begin(), ff.end(), zero);
		memset(&(ff[0]), 0, nqx_ * nqy_ * nqz_);

		return run_init(rot_h, transvec);
	} // AnalyticFormFactorG::run_init()


	bool AnalyticFormFactorG::clear() {
		return destroy();
	} // AnalyticFormFactorG::clear()


	// release all resources and init to 0
	bool AnalyticFormFactorG::destroy() {
		if(ff_buff_h_ != NULL) delete[] ff_buff_h_;
		if(ff_buff_d_ != NULL) cudaFree(ff_buff_d_);
		if(rot_ != NULL) cudaFree(rot_);
		if(transvec_ != NULL) cudaFree(transvec_);
		if(ff_ != NULL) cudaFree(ff_);
		if(qz_ != NULL) cudaFree(qz_);
		if(qy_ != NULL) cudaFree(qy_);
		if(qx_ != NULL) cudaFree(qx_);

		ff_buff_h_ = NULL; ff_buff_d_ = NULL;
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


	bool AnalyticFormFactorG::move_ff_buff_to_host_ff(std::vector<complex_t>& ff,
								unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
								unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
								int ib_x, int ib_y, int ib_z) {
		unsigned int size = curr_b_nqx * curr_b_nqy * curr_b_nqz;
		cudaMemcpy(ff_buff_h_, ff_buff_d_, size * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
		// make sure ff was already allocated and initialized to 0 during init ...
		unsigned long int base_i = nqx_ * nqy_ * ib_z * b_nqz + nqx_ * ib_y * b_nqy + ib_x * b_nqx;
		for(int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
			unsigned long int start_i = base_i + nqx_ * nqy_ * i_z;
			unsigned long int super_i = 0;
			for(int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
				super_i = start_i + nqx_ * i_y;
				for(int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
					unsigned long int final_i = super_i + i_x;
					unsigned long int block_i = curr_b_nqx * curr_b_nqy * i_z + curr_b_nqx * i_y + i_x;
					ff[final_i] += complex_t(ff_buff_h_[block_i].x, ff_buff_h_[block_i].y);
					/*if(!(boost::math::isnormal(ff_buff_h_[block_i].x) &&
								boost::math::isnormal(ff_buff_h_[block_i].y))) {
						std::cerr << "********************* oh here ..... "
									<< i_x << ", " << i_y << ", " << i_z
									<< ": " << ff_buff_h_[block_i].x << ", " << ff_buff_h_[block_i].y
									<< std::endl;
					}*/
				} // for x
			} // for y
		} // for z
		return true;
	} // AnalyticFormFactorG::move_ff_buff_to_host_ff()


	bool AnalyticFormFactorG::compute_hyperblock_size(unsigned int est_device_mem_need,
								unsigned int avail_device_mem) {
		unsigned int FF_ANA_HBLOCK_X_ = 100;
		unsigned int FF_ANA_HBLOCK_Y_ = 256;
		unsigned int FF_ANA_HBLOCK_Z_ = 256;
		b_nqx_ = std::min(nqx_, FF_ANA_HBLOCK_X_);
		b_nqy_ = std::min(nqy_, FF_ANA_HBLOCK_Y_);
		b_nqz_ = std::min(nqz_, FF_ANA_HBLOCK_Z_);
		b_nqx_ = nqx_;
		b_nqy_ = nqy_;
		b_nqz_ = nqz_;
		std::cout << "++               Hyperblock size: " << b_nqx_ << " x " << b_nqy_ << " x "
					<< b_nqz_ << std::endl;
		return true;
	} // AnalyticFormFactorG::compute_hyperblock_size()

} // namespace hig

