/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf_gpu.cu
 *  Created: Oct 15, 2012
 *  Modified: Wed 08 Oct 2014 12:17:49 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Abhinav Sarje <asarje@lbl.gov>
 *              Dinesh Kumar <dkumar@lbl.gov>
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

#include <woo/timer/woo_boostchronotimers.hpp>
#include <sf/sf.hpp>
#include <model/qgrid.hpp>
#include <common/typedefs.hpp>
#include <common/constants.hpp>
#include <utils/utilities.hpp>
#include <utils/gpu/cu_utilities.cuh>

namespace hig {


#ifdef SF_GPU
  bool StructureFactor::compute_structure_factor_gpu(std::string expt, vector3_t center,
                                 Lattice* lattice, vector3_t repet, vector3_t scaling,
                                 RotMatrix_t & rot
                                 #ifdef USE_MPI
                                   , woo::MultiNode& world_comm, std::string comm_key
                                 #endif
                                ) {
    ny_ = QGrid::instance().nqy();
    if(expt == "saxs") nz_ = QGrid::instance().nqz();
    else if(expt == "gisaxs") nz_ = QGrid::instance().nqz_extended();
    sf_ = new (std::nothrow) complex_t[nz_];
    if(sf_ == NULL) return false;
    gsf_.init();
    bool ret = gsf_.compute(expt, center, lattice, repet, scaling, rot
                        #ifdef USE_MPI
                          , world_comm, comm_key
                        #endif
                       );
    gsf_.get_sf(sf_);
    gsf_.destroy();
    return ret;
  } // StructureFactor::compute_structure_factor_gpu()

  StructureFactorG::StructureFactorG():
                    inited_(false),
                    nqx_(0), nqy_(0), nqz_(0),
                    sf_(NULL), qx_(NULL), qy_(NULL), qz_(NULL),
                    repet_(NULL), rot_(NULL), center_(NULL), transvec_(NULL) {
  } // StructureFactorG::StructureFactorG()


  StructureFactorG::~StructureFactorG() {
    destroy();
  } // StructureFactorG::~StructureFactorG()


  bool StructureFactorG::init() {
    // allocate device buffers
    // copy qgrid to device memory

    // TODO: unify ff and sf stuff. put separate qgrid for gpu ...

    if(inited_) {
      std::cerr << "error: gpu sf already initialized" << std::endl;
      return false;
    } // if

    nqx_ = QGrid::instance().nqx();
    nqy_ = QGrid::instance().nqy();
    nqz_ = QGrid::instance().nqz_extended();
    cudaMalloc((void**) &qx_, nqx_ * sizeof(real_t));
    cudaMalloc((void**) &qy_, nqy_ * sizeof(real_t));
    cudaMalloc((void**) &qz_, nqz_ * sizeof(cucomplex_t));
    cudaMalloc((void**) &sf_, nqz_ * sizeof(cucomplex_t));
    cudaMalloc((void**) &repet_, 3 * sizeof(real_t));
    cudaMalloc((void**) &rot_, 9 * sizeof(real_t));
    cudaMalloc((void**) &center_, 3 * sizeof(real_t));
    cudaMalloc((void**) &transvec_, 3 * sizeof(real_t));
    if(qx_ == NULL || qy_ == NULL || qz_ == NULL || sf_ == NULL || repet_ == NULL || rot_ == NULL ||
        center_ == NULL || transvec_ == NULL) {
      std::cerr << "error: device memory allocation failed for structure factor" << std::endl;
      return false;
    } // if
    cudaError_t err = cudaGetLastError();

    // construct host buffers
    real_t* qx_h = new (std::nothrow) real_t[nqx_];
    real_t* qy_h = new (std::nothrow) real_t[nqy_];
    cucomplex_t* qz_h = new (std::nothrow) cucomplex_t[nqz_];
    if(qx_h == NULL || qy_h == NULL || qz_h == NULL) {
      std::cerr << "error: memory allocation for host grid in sf failed" << std::endl;
      return false;
    } // if
    for(unsigned int ix = 0; ix < nqx_; ++ ix) qx_h[ix] = QGrid::instance().qx(ix);
    for(unsigned int iy = 0; iy < nqy_; ++ iy) qy_h[iy] = QGrid::instance().qy(iy);
    for(unsigned int iz = 0; iz < nqz_; ++ iz) {
      qz_h[iz].x = QGrid::instance().qz_extended(iz).real();
      qz_h[iz].y = QGrid::instance().qz_extended(iz).imag();
    } // for

    // TODO: make them async if it helps ...
    cudaMemcpy(qx_, qx_h, nqx_ * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_, qy_h, nqy_ * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

    // TODO: do hyperblocking if needed just as in ff ...

    #ifdef MEM_DETAIL
      size_t device_mem_avail, device_mem_used, device_mem_total;
      cudaMemGetInfo(&device_mem_avail, &device_mem_total);
      device_mem_used = device_mem_total - device_mem_avail;
      std::cout << "++            Used device memory: " << (float) device_mem_used / 1024 / 1024
                << " MB" << std::endl;
      std::cout << "++            Free device memory: " << (float) device_mem_avail / 1024 / 1024
                << " MB" << std::endl;
    #endif

    inited_ = true;

    delete[] qz_h;
    delete[] qy_h;
    delete[] qx_h;

    return true;
  } // StructureFactorG::init()


  // TODO: avoid all this copying .........
  void StructureFactorG::get_sf(complex_t* &sf) {
    // copy sf to main memory from device mem
    cucomplex_t * temp_buf = new (std::nothrow) cucomplex_t[nqz_];
    if(temp_buf == NULL) {
      std::cerr << "error: memory allocation for temporary buffer failed" << std::endl;
      return;
    } // if

    // TODO: doing blocking copy for now. lots of room to improve ...
    cudaMemcpy(temp_buf, sf_, nqz_ * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);

    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "error: structure factor memory copy failed: " << cudaGetErrorString(err) << std::endl;
    } else {
      #pragma omp parallel for
      for(unsigned int i = 0; i < nqz_; ++ i)
        sf[i] = complex_t(temp_buf[i].x, temp_buf[i].y);
    } // if-else

    delete[] temp_buf;
  } // StructureFactorG::get_sf()


/*  bool StructureFactorG::run_init(const real_t* rot_h, const std::vector<real_t>& repet) {
    // copy repet and rotation to device memory
    if(repet_ == NULL || rot_ == NULL) {
      std::cerr << "error: StructureFactorG is not initialized" << std::endl;
      return false;
    } // if

    const real_t* repet_h = repet.empty() ? NULL : &*repet.begin();
    if(repet_h == NULL) return false;

    cudaMemcpy(repet_, repet_h, 3 * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(rot_, rot_h, 9 * sizeof(real_t), cudaMemcpyHostToDevice);

    return true;
  } // StructureFactorG::run_init()
*/

  bool StructureFactorG::clear() {
    return destroy();
  } // StructureFactorG::clear();


  bool StructureFactorG::destroy() {
    if(inited_) {
      if(transvec_ != NULL) cudaFree(transvec_); transvec_ = NULL;
      if(center_ != NULL) cudaFree(center_); center_ = NULL;
      if(rot_ != NULL) cudaFree(rot_); rot_ = NULL;
      if(repet_ != NULL) cudaFree(repet_); repet_ = NULL;
      if(sf_ != NULL) cudaFree(sf_); sf_ = NULL;
      if(qz_ != NULL) cudaFree(qz_); qz_ = NULL;
      if(qy_ != NULL) cudaFree(qy_); qy_ = NULL;
      if(qx_ != NULL) cudaFree(qx_); qx_ = NULL;
    } // if
    nqx_ = nqy_ = nqz_ = 0;
    inited_ = false;
    return true;
  } // StructureFactorG::destroy()


  bool StructureFactorG::compute(std::string expt, vector3_t center,
                                 Lattice* lattice, vector3_t repet, vector3_t scaling,
                                 RotMatrix_t rot
                                 #ifdef USE_MPI
                                   , woo::MultiNode& world_comm, std::string comm_key
                                 #endif
                                ) {
    #ifdef USE_MPI
      int my_rank = world_comm.rank(comm_key);
      bool master = world_comm.is_master(comm_key);
    #else
      int my_rank = 0;
      bool master = true;
    #endif

    #ifdef SF_VERBOSE
      if(master) std::cout << "-- Computing structure factor on GPU ... " << std::flush;
    #endif

    if(repet[0] < 1) repet[0] = 1;
    if(repet[1] < 1) repet[1] = 1;
    if(repet[2] < 1) repet[2] = 1;

    vector3_t temp_la(lattice->a() * scaling[0]),
              temp_lb(lattice->b() * scaling[1]),
              temp_lc(lattice->c() * scaling[2]);
    vector3_t l_t = lattice->t();
    vector3_t arot = rot * temp_la;
    vector3_t brot = rot * temp_lb;
    vector3_t crot = rot * temp_lc;

    real_t rot_h[9];
    rot_h[0] = arot[0]; rot_h[1] = arot[1]; rot_h[2] = arot[2];
    rot_h[3] = brot[0]; rot_h[4] = brot[1]; rot_h[5] = brot[2];
    rot_h[6] = crot[0]; rot_h[7] = crot[1]; rot_h[8] = crot[2];

    // copy current rot matrix and num repeats to device
    cudaMemcpy(rot_, rot_h, 9 * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(repet_, &repet[0], 3 * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(center_, &center[0], 3 * sizeof(real_t), cudaMemcpyHostToDevice);
    cudaMemcpy(transvec_, &l_t[0], 3 * sizeof(real_t), cudaMemcpyHostToDevice);

    //unsigned int cuda_block_y = 16, cuda_block_z = 8;
    //unsigned int cuda_num_blocks_y = ceil((real_t) nqy_ / cuda_block_y);
    //unsigned int cuda_num_blocks_z = ceil((real_t) nqz_ / cuda_block_z);
    //dim3 sf_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
    //dim3 sf_block_size(cuda_block_y, cuda_block_z, 1);
    unsigned int cuda_block = 256;
    unsigned int cuda_num_blocks = ceil((real_t) nqz_ / cuda_block);
    dim3 sf_grid_size(cuda_num_blocks, 1, 1);
    dim3 sf_block_size(cuda_block, 1, 1);

    // call the kernel
    structure_factor_kernel <<< sf_grid_size, sf_block_size >>>
                            (nqx_, nqy_, nqz_, qx_, qy_, qz_, rot_, repet_, center_, transvec_, sf_);
    
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "error: structure factor kernel failed: " << cudaGetErrorString(err) << std::endl;
      return false;
    } // if

    return true;
  } // StructureFactorG::compute()

/*    std::complex<real_t> unit_c(1, 0);
    std::complex<real_t> unit_ci(0, 1);

    // big data to transer to gpu:
    // qx, qy, qz
    // arot, brot, crot
    // repet

    // good for gpu ... TODO
    for(unsigned int z = 0; z < nz_; ++ z) {
      for(unsigned int y = 0; y < ny_; ++ y) {
        for(unsigned int x = 0; x < nx_; ++ x) {
          complex_t temp1, temp_x2, temp_y3, temp_y4, temp_x5;
          real_t temp_f;
          complex_t sa, sb, sc;
          real_t qx = QGrid::instance().qx(x);
          real_t qy = QGrid::instance().qy(y);
          complex_t qz;
          if(expt == "saxs")
            qz = QGrid::instance().qz(z);
          else if(expt == "gisaxs")
            qz = QGrid::instance().qz_extended(z);

          temp1 = exp(unit_ci * (arot[0] * qx + arot[1] * qy + arot[2] * qz));
          temp_x2 = unit_c - pow(temp1, repet[0]);
          temp_y3 = unit_c / (unit_c - temp1);
          temp_f = (real_t)(!((boost::math::isfinite)(temp_y3.real()) &&
                      (boost::math::isfinite)(temp_y3.imag())));
          temp_y4 = unit_c / (unit_c / temp_y3 + temp_f);
          temp_f = (real_t)(!((boost::math::isfinite)((unit_c / temp_x2).real()) &&
                      (boost::math::isfinite)((unit_c / temp_x2).imag())));
          temp_x5 = temp_x2 + repet[0] * temp_f;
          sa = pow(temp1, ((real_t)1.0 - repet[0]) / (real_t)2.0) * temp_y4 * temp_x5;

          temp1 = exp(unit_ci * (brot[0] * qx + brot[1] * qy + brot[2] * qz));
          temp_x2 = unit_c - pow(temp1, repet[1]);
          temp_y3 = unit_c / (unit_c - temp1);
          temp_f = (real_t)(!((boost::math::isfinite)(temp_y3.real()) &&
                      (boost::math::isfinite)(temp_y3.imag())));
          temp_y4 = unit_c / (unit_c / temp_y3 + temp_f);
          temp_f = (real_t)(!((boost::math::isfinite)((unit_c / temp_x2).real()) &&
                      (boost::math::isfinite)((unit_c / temp_x2).imag())));
          temp_x5 = temp_x2 + repet[1] * temp_f;
          sb = pow(temp1, ((real_t)1.0 - repet[1]) / (real_t)2.0) * temp_y4 * temp_x5;

          temp1 = exp(unit_ci * (crot[0] * qx + crot[1] * qy + crot[2] * qz));
          temp_x2 = unit_c - pow(temp1, repet[2]);
          temp_y3 = unit_c / (unit_c - temp1);
          temp_f = (real_t)(!((boost::math::isfinite)(temp_y3.real()) &&
                      (boost::math::isfinite)(temp_y3.imag())));
          temp_y4 = unit_c / (unit_c / temp_y3 + temp_f);
          temp_f = (real_t)(!((boost::math::isfinite)((unit_c / temp_x2).real()) &&
                      (boost::math::isfinite)((unit_c / temp_x2).imag())));
          temp_x5 = temp_x2 + repet[2] * temp_f;
          sc = temp_y4 * temp_x5;

*/          /*if(!((boost::math::isfinite)(sa.real()) && (boost::math::isfinite)(sa.imag()))) {
            std::cout << "sa sa sa sa sa sa sa: " << x << ", " << y << ", " << z << std::endl; }
          if(!((boost::math::isfinite)(sb.real()) && (boost::math::isfinite)(sb.imag()))) {
            std::cout << "sb sb sb sb sb sb sb: " << x << ", " << y << ", " << z << std::endl; }
          if(!((boost::math::isfinite)(sc.real()) && (boost::math::isfinite)(sc.imag()))) {
            std::cout << "sc sc sc sc sc sc sc: " << x << ", " << y << ", " << z << std::endl; }*/
/*
          sf_[nx_ * ny_ * z + nx_ * y + x] = exp(unit_ci *
                  (center[0] * qx + center[1] * qy + center[2] * qz)) *
                  sa * sb * sc *
                  (unit_c + exp(unit_ci * (l_t[0] * qx + l_t[1] * qy + l_t[2] * qz)));

  */        /*if(!((boost::math::isfinite)(sf_[nx_ * ny_ * z + nx_ * y + x].real()) &&
                (boost::math::isfinite)(sf_[nx_ * ny_ * z + nx_ * y + x].imag()))) {
            std::cout << "sf sf sf sf sf sf sf: " << x << ", " << y << ", " << z << std::endl; }*/
/*        } // for x
      } // for y
    } // for z

    if(my_rank == 0) {
      int naninfs = count_naninfs(nx_, ny_, nz_, sf_);
      std::cout << " ------- " << naninfs << " / " << nx_ * ny_ * nz_ << " nans or infs" << std::endl;
    } // if */

  /*__global__ void structure_factor_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                                          real_t* qx, real_t* qy, cucomplex_t* qz,
                                          real_t* rot, real_t* repet,
                                          real_t* center, real_t* transvec,
                                          cucomplex_t* sf) {
    unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
    unsigned int base_index = nqx * nqy * i_z + nqx * i_y;

    cucomplex_t unit_c = make_cuC((real_t) 1, (real_t) 0);

    if(i_y < nqy && i_z < nqz) {
      for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
        cucomplex_t sa, sb, sc;

        real_t q_x = qx[i_x], q_y = qy[i_y];
        cucomplex_t q_z = qz[i_z];

        cucomplex_t e_iqa = cuCexpi(rot[0] * q_x + rot[1] * q_y + rot[2] * q_z);
        cucomplex_t xa_0 = unit_c - cuCpow(e_iqa, repet[0]);
        cucomplex_t ya_0 = unit_c - e_iqa;
        if(fabs(ya_0.y) > REAL_ZERO_ || fabs(ya_0.x) > REAL_ZERO_) sa = xa_0 / ya_0;
        else sa = make_cuC(repet[0], 0);
        sa = cuCpow(e_iqa, ((real_t) 1.0 - repet[0]) / (real_t) 2.0) * sa;

        cucomplex_t e_iqb = cuCexpi(rot[3] * q_x + rot[4] * q_y + rot[5] * q_z);
        cucomplex_t xb_0 = unit_c - cuCpow(e_iqb, repet[1]);
        cucomplex_t yb_0 = unit_c - e_iqb;
        if(fabs(yb_0.y) > REAL_ZERO_ || fabs(yb_0.x) > REAL_ZERO_) sb = xb_0 / yb_0;
        else sb = make_cuC(repet[1], 0);
        sb = cuCpow(e_iqb, ((real_t) 1.0 - repet[1]) / (real_t) 2.0) * sb;

        cucomplex_t e_iqc = cuCexpi(rot[6] * q_x + rot[7] * q_y + rot[8] * q_z);
        cucomplex_t xc_0 = unit_c - cuCpow(e_iqc, repet[2]);
        cucomplex_t yc_0 = unit_c - e_iqc;
        if(fabs(yc_0.y) > REAL_ZERO_ || fabs(yc_0.x) > REAL_ZERO_) sc = xc_0 / yc_0;
        else sc = make_cuC(repet[2], 0);
        sc = cuCpow(e_iqc, ((real_t) 1.0 - repet[2]) / (real_t) 2.0) * sc;

*/
/*        cucomplex_t iqa = rot[0] * q_x + rot[1] * q_y + rot[2] * q_z;
        cucomplex_t iqna = repet[0] * iqa;
        cucomplex_t e_iqa = cuCexpi(iqa);
        cucomplex_t xa_0 = unit_c - cuCexpi(iqna);
        cucomplex_t ya_0 = unit_c - e_iqa;
        if(fabs(ya_0.y) > REAL_ZERO_ || fabs(ya_0.x) > REAL_ZERO_) sa = xa_0 / ya_0;
        else sa = make_cuC(repet[0], 0);
        sa = cuCpow(e_iqa, ((real_t) 1.0 - repet[0]) / (real_t) 2.0) * sa;

        cucomplex_t iqb = rot[3] * q_x + rot[4] * q_y + rot[5] * q_z;
        cucomplex_t iqnb = repet[1] * iqb;
        cucomplex_t e_iqb = cuCexpi(iqb);
        cucomplex_t xb_0 = unit_c - cuCexpi(iqnb);
        cucomplex_t yb_0 = unit_c - e_iqb;
        if(fabs(yb_0.y) > REAL_ZERO_ || fabs(yb_0.x) > REAL_ZERO_) sb = xb_0 / yb_0;
        else sb = make_cuC(repet[1], 0);
        sb = cuCpow(e_iqb, ((real_t) 1.0 - repet[1]) / (real_t) 2.0) * sb;

        cucomplex_t iqc = rot[6] * q_x + rot[7] * q_y + rot[8] * q_z;
        cucomplex_t iqnc = repet[2] * iqc;
        cucomplex_t e_iqc = cuCexpi(iqc);
        cucomplex_t xc_0 = unit_c - cuCexpi(iqnc);
        cucomplex_t yc_0 = unit_c - e_iqc;
        if(fabs(yc_0.y) > REAL_ZERO_ || fabs(yc_0.x) > REAL_ZERO_) sc = xc_0 / yc_0;
        else sc = make_cuC(repet[2], 0);
        sc = cuCpow(e_iqc, ((real_t) 1.0 - repet[2]) / (real_t) 2.0) * sc;
*/
/*        unsigned long int sf_i = base_index + i_x;
        cucomplex_t temp3 = cuCexpi(center[0] * q_x + center[1] * q_y + center[2] * q_z);
        cucomplex_t temp2 = unit_c + cuCexpi(transvec[0] * q_x + transvec[1] * q_y + transvec[2] * q_z);
        sf[sf_i] = temp3 * temp2 * sa * sb * sc;
      } // for x
    } // if
  } // StructureFactor::compute_structure_factor_gpu()
*/

  __global__ void structure_factor_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                                          real_t* qx, real_t* qy, cucomplex_t* qz,
                                          real_t* rot, real_t* repet,
                                          real_t* center, real_t* transvec,
                                          cucomplex_t* sf) {
    unsigned int i_z = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int i_y = i_z % nqy;
    unsigned int i_x = i_z % nqx;

    cucomplex_t unit_c = make_cuC(REAL_ONE_, REAL_ZERO_);

    if(i_z < nqz) {
      cucomplex_t sa, sb, sc;

      real_t q_x = qx[i_x], q_y = qy[i_y];
      cucomplex_t q_z = qz[i_z];

      cucomplex_t e_iqa = cuCexpi(rot[0] * q_x + rot[1] * q_y + rot[2] * q_z);
      cucomplex_t xa_0 = unit_c - cuCpow(e_iqa, repet[0]);
      cucomplex_t ya_0 = unit_c - e_iqa;
      if(fabs(ya_0.y) > REAL_ZERO_ || fabs(ya_0.x) > REAL_ZERO_) sa = xa_0 / ya_0;
      else sa = make_cuC(repet[0], 0);
      sa = cuCpow(e_iqa, ((real_t) 1.0 - repet[0]) / (real_t) 2.0) * sa;

      cucomplex_t e_iqb = cuCexpi(rot[3] * q_x + rot[4] * q_y + rot[5] * q_z);
      cucomplex_t xb_0 = unit_c - cuCpow(e_iqb, repet[1]);
      cucomplex_t yb_0 = unit_c - e_iqb;
      if(fabs(yb_0.y) > REAL_ZERO_ || fabs(yb_0.x) > REAL_ZERO_) sb = xb_0 / yb_0;
      else sb = make_cuC(repet[1], 0);
      sb = cuCpow(e_iqb, ((real_t) 1.0 - repet[1]) / (real_t) 2.0) * sb;

      cucomplex_t e_iqc = cuCexpi(rot[6] * q_x + rot[7] * q_y + rot[8] * q_z);
      cucomplex_t xc_0 = unit_c - cuCpow(e_iqc, repet[2]);
      cucomplex_t yc_0 = unit_c - e_iqc;
      if(fabs(yc_0.y) > REAL_ZERO_ || fabs(yc_0.x) > REAL_ZERO_) sc = xc_0 / yc_0;
      else sc = make_cuC(repet[2], 0);
      sc = cuCpow(e_iqc, ((real_t) 1.0 - repet[2]) / (real_t) 2.0) * sc;

      cucomplex_t temp3 = cuCexpi(center[0] * q_x + center[1] * q_y + center[2] * q_z);
      cucomplex_t temp2 = unit_c + cuCexpi(transvec[0] * q_x + transvec[1] * q_y + transvec[2] * q_z);
      sf[i_z] = temp3 * temp2 * sa * sb * sc;
    } // if
  } // StructureFactor::compute_structure_factor_gpu()

#endif // SF_GPU

} // namespace
