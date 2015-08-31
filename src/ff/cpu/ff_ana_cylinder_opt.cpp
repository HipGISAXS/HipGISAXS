/***
  *  Project:
  *
  *  File:
  *  Created: Aug 31, 2015
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */


#ifdef PROFILE_PAPI
#include <papi.h>
#endif // PROFILE_PAPI

#include <ff/ff_ana.hpp>
#include <model/qgrid.hpp>

namespace hig {

#ifdef FF_CPU_OPT

  bool AnalyticFormFactor::cylinder_opt(std::vector<real_t>& r, std::vector<real_t>& h, vector3_t transvec,
                                   std::vector<complex_t>& ff) {

    #ifdef PROFILE_PAPI
    long long int papi_total_cycles = 0, papi_total_inst = 0, papi_total_flop = 0;
    double overall_ipc = 0.0;
    int papi_events[3] = { PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_OPS };
    long long papi_counter_values[3];
    #endif  // PROFILE_PAPI

    #ifdef PROFILE_PAPI
    PAPI_start_counters(papi_events, 3);
    #endif // PROFILE_PAPI

    #pragma omp parallel for schedule(runtime) shared(transvec, r, h, ff)
    for(unsigned z = 0; z < nqz_; ++ z) {
      unsigned int y = z % nqy_;
      std::vector<complex_t> mq = rot_.rotate(QGrid::instance().qx(y), QGrid::instance().qy(y),
              QGrid::instance().qz_extended(z));
      complex_t qpar = sqrt(mq[0] * mq[0] + mq[1] * mq[1]);
      complex_t temp_ff(0.0, 0.0);
      for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
        for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
          temp_ff += ff_cylinder_kernel(qpar, mq[2], r[i_r], h[i_h]);
        } // for h
      } // for r
      complex_t temp1 = mq[0] * transvec[0] + mq[1] * transvec[1] + mq[2] * transvec[2];
      complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
      ff[z] = temp_ff * temp2;
    } // for z

    #ifdef PROFILE_PAPI
    PAPI_stop_counters(papi_counter_values, 3);
    papi_total_cycles += papi_counter_values[0];
    papi_total_inst += papi_counter_values[1];
    papi_total_flop += papi_counter_values[2];
    std::cout << "++                  PAPI_TOT_CYC: " << papi_total_cycles << std::endl;
    std::cout << "++                  PAPI_TOT_INS: " << papi_total_inst << std::endl;
    std::cout << "++                   PAPI_FP_OPS: " << papi_total_flop << std::endl;
    std::cout << "++                           IPC: "
              << (double) papi_total_inst / papi_total_cycles << std::endl;
    #endif // PROFILE_PAPI

    return true;
  } // AnalyticFormFactor::cylinder_opt()

#endif // FF_CPU_OPT

} // namespace hig
