/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_gpu.cu
 *  Created: Oct 16, 2012
 *  Modified: Wed 08 Oct 2014 12:17:47 PM PDT
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

#include <ff/gpu/ff_ana_gpu.cuh>
#include <model/qgrid.hpp>
#include <common/enums.hpp>

namespace hig {

  __constant__ float_t tau_d;
  __constant__ float_t eta_d;
  __constant__ float_t transvec_d[3];
  __constant__ float_t rot_d[9];

  AnalyticFormFactorG::AnalyticFormFactorG(unsigned int ny, unsigned int nz):
    nqy_(ny), nqz_(nz),
    qx_(NULL), qy_(NULL), qz_(NULL), ff_(NULL) {
  } // AnalyticFormFactorG::AnalyticFormFactorG();


  AnalyticFormFactorG::AnalyticFormFactorG():
    nqy_(0), nqz_(0),
    qx_(NULL), qy_(NULL), qz_(NULL), ff_(NULL){
  } // AnalyticFormFactorG::AnalyticFormFactorG()


  AnalyticFormFactorG::~AnalyticFormFactorG() {
    // release device memories
    destroy();
  } // AnalyticFormFactorG::~AnalyticFormFactorG()


  void AnalyticFormFactorG::grid_size(unsigned int ny, unsigned int nz) {
    nqy_ = ny;  nqz_ = nz;
  } // AnalyticFormFactorG::grid_size()


  bool AnalyticFormFactorG::init(unsigned int nqy, unsigned int nqz) {
    // this does the following:
    //   + allocate device buffers
    //   + copy qgrid to device memory
    nqy_ = nqy; nqz_ = nqz;
    cudaMalloc((void**) &qx_, nqy_ * sizeof(float_t));
    cudaMalloc((void**) &qy_, nqy_ * sizeof(float_t));
    cudaMalloc((void**) &qz_, nqz_ * sizeof(cucomplex_t));
    cudaMalloc((void**) &ff_, nqz_ * sizeof(cucomplex_t));

    if(qx_ == NULL || qy_ == NULL || qz_ == NULL || ff_ == NULL) {
      std::cerr << "error: device memory allocation failed" << std::endl;
      return false;
    } // if

    cudaError_t err = cudaGetLastError();
    //std::cerr << "AFTER MALLOCS: " << cudaGetErrorString(err) << std::endl;

    // first need to construct host buffers
    float_t* qx_h = new (std::nothrow) float_t[nqy_];
    float_t* qy_h = new (std::nothrow) float_t[nqy_];
    cucomplex_t* qz_h = new (std::nothrow) cucomplex_t[nqz_];
    if(qx_h == NULL || qy_h == NULL || qz_h == NULL) {
      std::cerr << "error: memory allocation for host mesh grid failed" << std::endl;
      return false;
    } // if
    for(unsigned int ix = 0; ix < nqy_; ++ ix) qx_h[ix] = QGrid::instance().qx(ix);
    for(unsigned int iy = 0; iy < nqy_; ++ iy) qy_h[iy] = QGrid::instance().qy(iy);
    for(unsigned int iz = 0; iz < nqz_; ++ iz) {
      qz_h[iz].x = QGrid::instance().qz_extended(iz).real();
      qz_h[iz].y = QGrid::instance().qz_extended(iz).imag();
    } // for qz

    cudaMemcpy(qx_, qx_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_, qy_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

    delete [] qx_h;
    delete [] qy_h;
    delete [] qz_h;
    return true;
  } // AnalyticFormFactorG::init()

  bool AnalyticFormFactorG::clear() {
    return destroy();
  } // AnalyticFormFactorG::clear()


  // release all resources and init to 0
  bool AnalyticFormFactorG::destroy() {
    if(ff_ != NULL) cudaFree(ff_);
    if(qz_ != NULL) cudaFree(qz_);
    if(qy_ != NULL) cudaFree(qy_);
    if(qx_ != NULL) cudaFree(qx_);

    ff_ = NULL;
    qz_ = NULL; qy_ = NULL; qx_ = NULL;
    nqy_ = 0; nqz_ = 0;
    return true;
  } // AnalyticFormFactorG::destroy()


  bool AnalyticFormFactorG::construct_output_ff(std::vector<complex_t>& ff) {
    unsigned int grid_size = nqz_;
    cucomplex_t* ff_h = new (std::nothrow) cucomplex_t[grid_size];
    cudaMemcpy(ff_h, ff_, grid_size * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
    ff.clear(); ff.reserve(grid_size);
    for(unsigned int i = 0; i < grid_size; ++ i) ff.push_back(complex_t(ff_h[i].x, ff_h[i].y));
    delete[] ff_h;
    return true;
  } // AnalyticFormFactorG::construct_output_ff()

} // namespace hig

