/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_cylinder_opt.cpp
  *  Created: Aug 31, 2015
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  *
  *  Licensing: The HipGISAXS software is only available to be downloaded and
  *  used by employees of academic research institutions, not-for-profit
  *  research laboratories, or governmental research facilities. Please read the
  *  accompanying LICENSE file before downloading the software. By downloading
  *  the software, you are agreeing to be bound by the terms of this
  *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
  */

#include <cstring>

#ifdef PROFILE_PAPI
# include <papi.h>
#endif // PROFILE_PAPI

#ifdef _OPENMP
# include <omp.h>
#endif // _OPENMP

#include <woo/timer/woo_boostchronotimers.hpp>
#include <common/constants.hpp>
#include <numerics/numeric_utils.hpp>
#include <ff/ff_ana.hpp>
#include <model/qgrid.hpp>

#ifdef FF_CPU_OPT
# if defined FF_CPU_OPT_MKL         // using Intel MKL's VML and CBLAS for computations

#   ifdef DOUBLEP
#     define MKL_Complex16 hig::complex_t
#   else
#     define MKL_Complex8 hig::complex_t
#   endif // DOUBLEP
#   include <mkl.h>
#   include <mkl_cblas.h>

# elif defined FF_CPU_OPT_AVX       // using Intel AVX intrinsics for computations

#   include <numerics/cpu/avx_numerics.hpp>

# endif // FF_CPU_OPT_XXX
#endif // FF_CPU_OPT

namespace hig {

#ifdef FF_CPU_OPT

#if defined FF_CPU_OPT_AVX      // to use manual avx intrinsics vectorization

  // TODO: directly overload operations instead of all the avx_* functions ...

  // ////////////////////////////////////////////////////
  // intrinsic wrapper naming (only for floating-point)
  // avx_xxx_ab  => a = r|c, b = p|s
  // avx_xxx_abc => a = r|c, b = r|c
  // r = real,             c = complex,
  // p = packed (vector),  s = scalar,
  // ////////////////////////////////////////////////////

  inline avx_m256c_t AnalyticFormFactor::ff_cylinder_kernel_opt_vec(const avx_m256c_t& qpar,
                                                                    const avx_m256c_t& qz,
                                                                    real_t radius, real_t height) {
    real_t vol = 2. * PI_ * radius * radius * height;
    real_t half_height = 0.5 * height;
    avx_m256c_t temp1 = avx_mul_crp(qz, half_height);
    avx_m256c_t sinc_val = avx_sinc_cp(temp1);
    avx_m256c_t temp2 = avx_mul_ccp(CMPLX_ONE_, temp1);
    avx_m256c_t expt_val = avx_exp_cp(temp2);
    temp1 = avx_mul_crp(qpar, radius);
    avx_m256c_t bess_val = avx_cbessj_cp(temp1, 1);
    bess_val = avx_div_ccp(bess_val, temp1);
    temp1 = avx_mul_ccp(sinc_val, expt_val);
    temp2 = avx_mul_ccp(temp1, bess_val);
    return avx_mul_crp(temp2, vol);
  } // AnalyticFormFactor::ff_cylinder_kernel_opt_vec()


  inline bool AnalyticFormFactor::cylinder_kernel_opt(std::vector<real_t>& r, std::vector<real_t>& h,
                                                      vector3_t transvec, std::vector<complex_t>& ff_vec) {
    real_t* qx = QGrid::instance().qx();
    real_t* qy = QGrid::instance().qy();
    complex_t* qz_extended = QGrid::instance().qz_extended();
    complex_t* ff = &ff_vec[0];

    real_t *rp = &r[0], *hp = &h[0];
    int rsize = r.size(), hsize = h.size();

    unsigned int loop_limit = nqz_ / AVX_VEC_LEN_;
    int loop_rem = nqz_ % AVX_VEC_LEN_;  // TODO ...

    #pragma omp for schedule(runtime)
    for(unsigned int l = 0; l < loop_limit; ++ l) {
      unsigned int z = AVX_VEC_LEN_ * l;
      unsigned int y = z % nqy_;
      avx_m256c_t mqx, mqy, mqz;
      rot_.rotate_vec(qx + y, qy + y, qz_extended + z, mqx, mqy, mqz);
      avx_m256c_t qx2 = avx_mul_ccp(mqx, mqx);
      avx_m256c_t qy2 = avx_mul_ccp(mqy, mqy);
      avx_m256c_t qxy2 = avx_add_ccp(qx2, qy2);
      avx_m256c_t qpar = avx_sqrt_cp(qxy2);
      avx_m256c_t temp_ff = avx_setzero_cp();
      for(unsigned int i_r = 0; i_r < rsize; ++ i_r) {
        for(unsigned int i_h = 0; i_h < hsize; ++ i_h) {
          avx_m256c_t temp = ff_cylinder_kernel_opt_vec(qpar, mqz, rp[i_r], hp[i_h]);
          temp_ff = avx_add_ccp(temp_ff, temp);
        } // for h
      } // for r
      avx_m256c_t temp1 = avx_mul_crp(mqz, transvec[2]);
      avx_m256c_t temp2 = avx_fma_rccp(transvec[1], mqy, temp1);
      temp1 = avx_fma_rccp(transvec[0], mqx, temp2);
      temp2 = avx_mul_ccp(CMPLX_ONE_, temp1);
      temp1 = avx_exp_cp(temp2);
      temp2 = avx_mul_ccp(temp_ff, temp1);
      avx_store_cp(ff + z, temp2);
    } // for z

    return true;
  } // AnalyticFormFactor::cylinder_kernel_opt()

#elif defined FF_CPU_OPT_MKL    // to use Intel MKL vector functions (VML and CBLAS)

  inline void AnalyticFormFactor::ff_cylinder_kernel_opt_vec(complex_t *buf2,
                                                             const complex_t* qpar, const complex_t* qz,
                                                             real_t radius, real_t height, complex_t* ff) {
    real_t vol = 2. * PI_ * radius * radius * height;
    complex_t *temp0, *temp1, *temp2;
    complex_t *sinc_val, *exp_val, *bess_val;
    temp0 = buf2; temp1 = temp0 + VEC_LEN; temp2 = temp1 + VEC_LEN;
    sinc_val = temp2 + VEC_LEN; exp_val = sinc_val + VEC_LEN; bess_val = exp_val + VEC_LEN;
    cblas_zcopy(VEC_LEN, qz, 1, temp0, 1);
    cblas_zdscal(VEC_LEN, 0.5 * height, temp0, 1);
    sinc_vec(VEC_LEN, temp0, sinc_val);
    cblas_zscal(VEC_LEN, &CMPLX_ONE_, temp0, 1);
    vzExp(VEC_LEN, temp0, exp_val);
    cblas_zcopy(VEC_LEN, qpar, 1, temp1, 1);
    cblas_zdscal(VEC_LEN, radius, temp1, 1);
    cbessj_vec(VEC_LEN, temp1, 1, temp2);
    vzDiv(VEC_LEN, temp2, temp1, bess_val);
    vzMul(VEC_LEN, sinc_val, exp_val, temp0);
    vzMul(VEC_LEN, temp0, bess_val, temp1);
    cblas_zdscal(VEC_LEN, vol, bess_val, 1);
    vzAdd(VEC_LEN, ff, bess_val, temp0);
    cblas_zcopy(VEC_LEN, temp0, 1, ff, 1);
  } // ff_cylinder_kernel_opt_vec()


  inline bool AnalyticFormFactor::cylinder_kernel_opt(std::vector<real_t>& r, std::vector<real_t>& h,
                                                      vector3_t transvec, std::vector<complex_t>& ff_vec) {
    real_t* qx = QGrid::instance().qx();
    real_t* qy = QGrid::instance().qy();
    complex_t* qz_extended = QGrid::instance().qz_extended();
    complex_t* ff = &ff_vec[0];

    real_t *rp = &r[0], *hp = &h[0];
    int rsize = r.size(), hsize = h.size();

    complex_t *buf = new (std::nothrow) complex_t[VEC_LEN * 14];
    complex_t *temp0, *temp1, *temp2, *temp_ff, *qpar;
    temp0 = buf; temp1 = temp0 + VEC_LEN; temp2 = temp1 + VEC_LEN; temp_ff = temp2 + VEC_LEN;
    qpar = temp_ff + VEC_LEN;
    complex_t *mq, *mq0, *mq1, *mq2;
    mq = qpar + VEC_LEN; mq0 = mq; mq1 = mq + VEC_LEN; mq2 = mq + 2 * VEC_LEN;
    complex_t *buf2 = mq + 3 * VEC_LEN;

    unsigned int loop_limit = nqz_ / VEC_LEN;
    int loop_rem = nqz_ % VEC_LEN;  // TODO ...

    vmlSetMode(VML_EP);
    //vmlSetMode(VML_LA);

    #pragma omp for schedule(runtime)
    for(unsigned int l = 0; l < loop_limit; ++ l) {
      unsigned int z = VEC_LEN * l;
      unsigned int y = z % nqy_;
      rot_.rotate_vec(VEC_LEN, qx + y, qy + y, qz_extended + z, mq);
      vzMul(VEC_LEN, mq0, mq0, temp0);
      vzMul(VEC_LEN, mq1, mq1, temp1);
      vzAdd(VEC_LEN, temp0, temp1, temp2);
      vzSqrt(VEC_LEN, temp2, qpar);
      memset(temp_ff, 0, VEC_LEN * sizeof(complex_t));
      for(int i_r = 0; i_r < rsize; ++ i_r) {
        for(int i_h = 0; i_h < hsize; ++ i_h) {
          ff_cylinder_kernel_opt_vec(buf2, qpar, mq2, rp[i_r], hp[i_r], temp_ff);
        } // for h
      } // for r
      cblas_zdscal(VEC_LEN, transvec[0], mq0, 1);
      cblas_zdscal(VEC_LEN, transvec[1], mq1, 1);
      cblas_zdscal(VEC_LEN, transvec[2], mq2, 1);
      vzAdd(VEC_LEN, mq0, mq1, temp2);
      vzAdd(VEC_LEN, mq2, temp2, temp1);
      cblas_zscal(VEC_LEN, &CMPLX_ONE_, temp1, 1);
      vzExp(VEC_LEN, temp1, temp2);
      vzMul(VEC_LEN, temp_ff, temp2, ff + z);
    } // for z

    delete[] buf;

#if 0
    #pragma omp for schedule(runtime)
    for(unsigned int l = 0; l < loop_limit; ++ l) {
      for(int i = 0; i < VEC_LEN; ++ i) {
        unsigned int z = VEC_LEN * l + i;
        unsigned int y = z % nqy_;
        rot_.rotate_new(qx[y], qy[y], qz_extended[z], mq);
        complex_t qpar = sqrt(mq[0] * mq[0] + mq[1] * mq[1]);
        complex_t temp_ff(0.0, 0.0);
        for(int i_r = 0; i_r < rsize; ++ i_r) {
          for(int i_h = 0; i_h < hsize; ++ i_h) {
              temp_ff += ff_cylinder_kernel_opt(qpar, mq[2], rp[i_r], hp[i_h]);
          } // for h
        } // for r
        complex_t temp1 = mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2];
        complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
        ff[z] = temp_ff * temp2;
      } // for i
    } // for z
#endif // 0

    return true;
  } // AnalyticFormFactor::cylinder_kernel_opt()

#else   // default is the original case FF_CPU_OPT_ORIG

  inline complex_t AnalyticFormFactor::ff_cylinder_kernel_opt(
                   complex_t qpar, complex_t qz, real_t radius, real_t height) {
    real_t vol = 2. * PI_ * radius * radius * height;
    complex_t sinc_val = sinc(0.5 * qz * height);
    complex_t expt_val = std::exp(CMPLX_ONE_ * 0.5 * qz * height);
    complex_t t1 = qpar * radius;
    complex_t bess_val = cbesselj(t1, 1) * std::conj(t1) / std::norm(t1);
    return (vol * sinc_val * expt_val * bess_val);
  } // ff_cylinder_kernel_opt()


  inline bool AnalyticFormFactor::cylinder_kernel_opt(std::vector<real_t>& r, std::vector<real_t>& h,
                                                      vector3_t transvec, std::vector<complex_t>& ff) {
    #pragma omp for schedule(runtime)
    for(unsigned z = 0; z < nqz_; ++ z) {
      unsigned int y = z % nqy_;
      std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), QGrid::instance().qy(y),
                                              QGrid::instance().qz_extended(z));
      complex_t qpar = sqrt(mq[0] * mq[0] + mq[1] * mq[1]);
      complex_t temp_ff(0.0, 0.0);
      for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
        for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
          temp_ff += ff_cylinder_kernel_opt(qpar, mq[2], r[i_r], h[i_h]);
        } // for h
      } // for r
      complex_t temp1 = mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2];
      complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
      ff[z] = temp_ff * temp2;
    } // for z

#if 0
    #pragma omp shared(transvec, r, h, ff)
    {
    for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
      for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {

        #pragma omp for schedule(runtime)
        for(unsigned z = 0; z < nqz_; ++ z) {
          unsigned int y = z % nqy_;
          std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), QGrid::instance().qy(y),
                  QGrid::instance().qz_extended(z));
          complex_t qpar = sqrt(mq[0] * mq[0] + mq[1] * mq[1]);
          complex_t temp_ff = ff_cylinder_kernel(qpar, mq[2], r[i_r], h[i_h]);
          complex_t temp1 = mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2];
          complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
          ff[z] += temp_ff * temp2;
        } // for z

      } // for h
    } // for r
    }
#endif

    return true;
  } // AnalyticFormFactor::cylinder_kernel_opt()

#endif // FF_CPU_OPT_XXX


  /**
   * The main wrapper function for cylinder form factor
   */
  bool AnalyticFormFactor::cylinder_opt(std::vector<real_t>& r, std::vector<real_t>& h, vector3_t transvec,
                                   std::vector<complex_t>& ff) {

    woo::BoostChronoTimer timer;

    int num_threads = 1, thread_id = 0;
    #ifdef _OPENMP
      #pragma omp parallel
      {
      if(omp_get_thread_num() == 0) num_threads = omp_get_num_threads();
      }
    #endif

    #ifdef PROFILE_PAPI
    double overall_ipc = 0.0;
    int papi_events[3] = { PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_DP_OPS };
    long long papi_counter_values[num_threads][3];
    if(PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
      std::cerr << "error: failed in PAPI_library_init()" << std::endl;
      return false;
    } // if
    if(PAPI_thread_init((unsigned long (*)(void))(omp_get_thread_num)) != PAPI_OK) {
      std::cerr << "error: failed in PAPI_thread_init()" << std::endl;
      return false;
    } // if
    #endif  // PROFILE_PAPI

    timer.start();

    #pragma omp parallel shared(transvec, r, h, ff) private(thread_id)
    {

    #ifdef _OPENMP
      thread_id = omp_get_thread_num();
    #endif

    #ifdef PROFILE_PAPI
    PAPI_start_counters(papi_events, 3);
    #endif // PROFILE_PAPI

    cylinder_kernel_opt(r, h, transvec, ff);

    #ifdef PROFILE_PAPI
    PAPI_stop_counters(papi_counter_values[thread_id], 3);
    #endif

    } // omp parallel

    timer.stop();
    std::cout << "++                  Compute time: " << timer.elapsed_msec() << " ms." << std::endl;

    #ifdef PROFILE_PAPI
    long long int papi_cycles = 0, papi_inst = 0, papi_flop = 0;
    long long int papi_total_cycles = 0, papi_total_inst = 0, papi_total_flop = 0;
    for(int i = 0; i < num_threads; ++ i) {
      papi_cycles = papi_counter_values[i][0];
      papi_inst = papi_counter_values[i][1];
      papi_flop = papi_counter_values[i][2];
      papi_total_cycles += papi_cycles;
      papi_total_inst += papi_inst;
      papi_total_flop += papi_flop;
      std::cout << "++                  PAPI_TOT_CYC: " << papi_cycles << std::endl;
      std::cout << "++                  PAPI_TOT_INS: " << papi_inst << std::endl;
      std::cout << "++                   PAPI_DP_OPS: " << papi_flop << std::endl;
      std::cout << "++                           IPC: "
                << (double) papi_inst / papi_cycles << std::endl;
      std::cout << "++                        MFLOPS: "
                << papi_flop / timer.elapsed_msec() / 1000. << std::endl;
    } // for
    std::cout << "++                  Total GFLOPS: "
              << papi_total_flop / timer.elapsed_msec() / 1000. / 1000. << std::endl;
    #endif // PROFILE_PAPI

    return true;
  } // AnalyticFormFactor::cylinder_opt()

#endif // FF_CPU_OPT

} // namespace hig
