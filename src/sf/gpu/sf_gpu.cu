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
                                 vector3_t rotation_1, vector3_t rotation_2, vector3_t rotation_3
                                 #ifdef USE_MPI
                                   , woo::MultiNode& world_comm, std::string comm_key
                                 #endif
                                ) {
    nx_ = QGrid::instance().nqx();
    ny_ = QGrid::instance().nqy();
    if(expt == "saxs") nz_ = QGrid::instance().nqz();
    else if(expt == "gisaxs") nz_ = QGrid::instance().nqz_extended();
    else return false;
    sf_ = new (std::nothrow) complex_t[nx_ * ny_ * nz_];
    gsf_.init(nx_, ny_, nz_);
    bool ret = gsf_.compute(expt, center, lattice, repet, scaling, rotation_1, rotation_2, rotation_3
                        #ifdef USE_MPI
                          , world_comm, comm_key
                        #endif
                       );
    gsf_.get_sf(sf_);
    gsf_.destroy();
    return ret;
  } // StructureFactor::compute_structure_factor_gpu()

  StructureFactorG::StructureFactorG():
                    nqx_(0), nqy_(0), nqz_(0),
                    sf_(NULL), qx_(NULL), qy_(NULL), qz_(NULL),
                    repet_(NULL), rot_(NULL), center_(NULL), transvec_(NULL) {
  } // StructureFactorG::StructureFactorG()


  StructureFactorG::~StructureFactorG() {
    destroy();
  } // StructureFactorG::~StructureFactorG()


  bool StructureFactorG::init(unsigned int nqx, unsigned int nqy, unsigned int nqz) {
    // allocate device buffers
    // copy qgrid to device memory

    // TODO: unify ff and sf stuff. put separate qgrid for gpu ...

    nqx_ = nqx; nqy_ = nqy; nqz_ = nqz;
    cudaMalloc((void**) &qx_, nqx_ * sizeof(float_t));
    cudaMalloc((void**) &qy_, nqy_ * sizeof(float_t));
    cudaMalloc((void**) &qz_, nqz_ * sizeof(cucomplex_t));
    cudaMalloc((void**) &sf_, nqx_ * nqy_ * nqz_ * sizeof(cucomplex_t));
    cudaMalloc((void**) &repet_, 3 * sizeof(float_t));
    cudaMalloc((void**) &rot_, 9 * sizeof(float_t));
    cudaMalloc((void**) &center_, 3 * sizeof(float_t));
    cudaMalloc((void**) &transvec_, 3 * sizeof(float_t));
    if(qx_ == NULL || qy_ == NULL || qz_ == NULL || sf_ == NULL || repet_ == NULL || rot_ == NULL) {
      std::cerr << "error: device memory allocation failed for structure factor" << std::endl;
      return false;
    } // if
    cudaError_t err = cudaGetLastError();

    // construct host buffers
    float_t* qx_h = new (std::nothrow) float_t[nqx_];
    float_t* qy_h = new (std::nothrow) float_t[nqz_];
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
    cudaMemcpy(qx_, qx_h, nqx_ * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qy_, qy_h, nqy_ * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(qz_, qz_h, nqz_ * sizeof(cucomplex_t), cudaMemcpyHostToDevice);

    // TODO: do hyperblocking if needed just as in ff ...

    size_t device_mem_avail, device_mem_used, device_mem_total;
    cudaMemGetInfo(&device_mem_avail, &device_mem_total);
    device_mem_used = device_mem_total - device_mem_avail;
    #ifdef MEM_DETAIL
      std::cout << "++            Used device memory: " << (float) device_mem_used / 1024 / 1024
                << " MB" << std::endl;
      std::cout << "++            Free device memory: " << (float) device_mem_avail / 1024 / 1024
                << " MB" << std::endl;
    #endif

    return true;
  } // StructureFactorG::init()


  // TODO: avoid all this copying .........
  void StructureFactorG::get_sf(complex_t* sf) {
    // copy sf to main memory from device mem
    cucomplex_t * temp_buf = new (std::nothrow) cucomplex_t[nqx_ * nqy_ * nqz_];

    // TODO: doing blocking copy for now. lots of room to improve ...
    cudaMemcpy(temp_buf, sf_, nqx_ * nqy_ * nqz_ * sizeof(cucomplex_t), cudaMemcpyDeviceToHost);
    #pragma omp parallel for
    for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i)
      sf[i] = complex_t(temp_buf[i].x, temp_buf[i].y);

    delete[] temp_buf;
  } // StructureFactorG::get_sf()


/*  bool StructureFactorG::run_init(const float_t* rot_h, const std::vector<float_t>& repet) {
    // copy repet and rotation to device memory
    if(repet_ == NULL || rot_ == NULL) {
      std::cerr << "error: StructureFactorG is not initialized" << std::endl;
      return false;
    } // if

    const float_t* repet_h = repet.empty() ? NULL : &*repet.begin();
    if(repet_h == NULL) return false;

    cudaMemcpy(repet_, repet_h, 3 * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(rot_, rot_h, 9 * sizeof(float_t), cudaMemcpyHostToDevice);

    return true;
  } // StructureFactorG::run_init()
*/

  bool StructureFactorG::clear() {
    return destroy();
  } // StructureFactorG::clear();


  bool StructureFactorG::destroy() {
    if(transvec_ != NULL) cudaFree(transvec_); transvec_ = NULL;
    if(center_ != NULL) cudaFree(center_); center_ = NULL;
    if(rot_ != NULL) cudaFree(rot_); rot_ = NULL;
    if(repet_ != NULL) cudaFree(repet_); repet_ = NULL;
    if(sf_ != NULL) cudaFree(sf_); sf_ = NULL;
    if(qz_ != NULL) cudaFree(qz_); qz_ = NULL;
    if(qy_ != NULL) cudaFree(qy_); qy_ = NULL;
    if(qx_ != NULL) cudaFree(qx_); qx_ = NULL;
    nqx_ = nqy_ = nqz_ = 0;
    return true;
  } // StructureFactorG::destroy()


  bool StructureFactorG::compute(std::string expt, vector3_t center,
                                 Lattice* lattice, vector3_t repet, vector3_t scaling,
                                 vector3_t rotation_1, vector3_t rotation_2, vector3_t rotation_3
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

    vector3_t arot(0, 0, 0), brot(0, 0, 0), crot(0, 0, 0);
    vector3_t temp_la(lattice->a() * scaling[0]),
              temp_lb(lattice->b() * scaling[1]),
              temp_lc(lattice->c() * scaling[2]);
    mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_la, arot);
    mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lb, brot);
    mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lc, crot);
    vector3_t l_t = lattice->t();

    float_t rot_h[9];
    rot_h[0] = arot[0]; rot_h[1] = arot[1]; rot_h[2] = arot[2];
    rot_h[3] = brot[0]; rot_h[4] = brot[1]; rot_h[5] = brot[2];
    rot_h[6] = crot[0]; rot_h[7] = crot[1]; rot_h[8] = crot[2];

    // copy current rot matrix and num repeats to device
    cudaMemcpy(rot_, rot_h, 9 * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(repet_, &repet[0], 3 * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(center_, &center[0], 3 * sizeof(float_t), cudaMemcpyHostToDevice);
    cudaMemcpy(transvec_, &l_t[0], 3 * sizeof(float_t), cudaMemcpyHostToDevice);

    unsigned int cuda_block_y = 16, cuda_block_z = 8;
    unsigned int cuda_num_blocks_y = ceil((float_t) nqy_ / cuda_block_y);
    unsigned int cuda_num_blocks_z = ceil((float_t) nqz_ / cuda_block_z);
    dim3 sf_grid_size(cuda_num_blocks_y, cuda_num_blocks_z, 1);
    dim3 sf_block_size(cuda_block_y, cuda_block_z, 1);

    // call the kernel
    structure_factor_kernel <<< sf_grid_size, sf_block_size >>>
                            (nqx_, nqy_, nqz_, qx_, qy_, qz_, rot_, repet_, center_, transvec_, sf_);
    
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
      std::cerr << "error: structure factor kernel failed: " << cudaGetErrorString(err) << std::endl;
      return false;
    } else {
      // move sf to cpu
      // TODO: merge sf ff multiplication to a gpu routine, wont need to transfer this output ...
      //cudaMemcpy();
    } // if-else

    return true;
  } // StructureFactorG::compute()


  __global__ void structure_factor_kernel(unsigned int nqx, unsigned int nqy, unsigned int nqz,
                                          float_t* qx, float_t* qy, cucomplex_t* qz,
                                          float_t* rot, float_t* repet,
                                          float_t* center, float_t* transvec,
                                          cucomplex_t* sf) {
    unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
    unsigned int base_index = nqx * nqy * i_z + nqx * i_y;

    cucomplex_t unit_c = make_cuC((float_t) 1, (float_t) 0);

    if(i_y < nqy && i_z < nqz) {
      for(unsigned int i_x = 0; i_x < nqx; ++ i_x) {
        cucomplex_t sa, sb, sc;

        float_t q_x = qx[i_x], q_y = qy[i_y];
        cucomplex_t q_z; q_z.x = qz[i_z].x; q_z.y = qz[i_z].y;

        cucomplex_t e_iqa = cuCexpi(rot[0] * q_x + rot[1] * q_y + rot[2] * q_z);
        cucomplex_t xa_0 = unit_c - cuCpow(e_iqa, repet[0]);
        cucomplex_t ya_0 = unit_c - e_iqa;

        float_t tempya = sqrt(ya_0.x * ya_0.x + ya_0.y * ya_0.y);
        if(fabs(ya_0.x) > REAL_ZERO_ || fabs(ya_0.y) > REAL_ZERO_) sa = xa_0 / ya_0;
        else sa = make_cuC(repet[0], 0);
        sa = cuCpow(e_iqa, ((float_t) 1.0 - repet[0]) / (float_t) 2.0) * sa;

        cucomplex_t iqb = make_cuC(- rot[5] * q_z.y, rot[3] * q_x + rot[4] * q_y + rot[5] * q_z.x);
        cucomplex_t iqnb = repet[1] * iqb;

        cucomplex_t e_iqb = cuCexp(iqb);
        cucomplex_t xb_0 = unit_c - cuCexp(iqnb);
        cucomplex_t yb_0 = unit_c - cuCexp(iqb);

        float_t tempyb = sqrt(yb_0.x * yb_0.x + yb_0.y * yb_0.y);
        if(fabs(yb_0.y) > REAL_ZERO_ || fabs(yb_0.x) > REAL_ZERO_) sb = xb_0 / yb_0;
        else sb = make_cuC(repet[1], 0);
        sb = cuCpow(e_iqb, ((float_t) 1.0 - repet[1]) / (float_t) 2.0) * sb;

        cucomplex_t e_iqc = cuCexpi(rot[6] * q_x + rot[7] * q_y + rot[8] * q_z);
        cucomplex_t xc_0 = unit_c - cuCpow(e_iqc, repet[2]);
        cucomplex_t yc_0 = unit_c - e_iqc;

        float_t tempyc = sqrt(yc_0.x * yc_0.x + yc_0.y * yc_0.y);
        if(fabs(yc_0.y) > REAL_ZERO_ || fabs(yc_0.x) > REAL_ZERO_) sc = xc_0 / yc_0;
        else sc = make_cuC(repet[2], 0);
        sc = cuCpow(e_iqc, ((float_t) 1.0 - repet[2]) / (float_t) 2.0) * sc;

        unsigned long int sf_i = base_index + i_x;
        cucomplex_t temp3 = cuCexpi(center[0] * q_x + center[1] * q_y + center[2] * q_z);
        cucomplex_t temp2 = unit_c + cuCexpi(transvec[0] * q_x + transvec[1] * q_y + transvec[2] * q_z);
        sf[sf_i] = temp3 * temp2 * sa * sb * sc;
      } // for x
    } // if
  } // StructureFactor::compute_structure_factor_gpu()

#endif // SF_GPU

} // namespace
