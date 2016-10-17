/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff.hpp
 *  Created: Jul 18, 2012
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

#ifndef _FF_HPP_
#define _FF_HPP_

#include <string>

#ifdef USE_MPI
#include <woo/comm/multi_node_comm.hpp>
#endif

#include <common/typedefs.hpp>
#include <numerics/matrix.hpp>
#include <ff/ff_ana.hpp>
#include <ff/ff_num.hpp>


namespace hig {

  /***
   * The main Form Factor class which stores the computed form factor data
   */
  class FormFactor {
    private:
      // TODO: fix the numeric and non-numeric version ... unify them ...
      bool is_analytic_;
      AnalyticFormFactor analytic_ff_;
      NumericFormFactor numeric_ff_;             // same as MFormFactor
      std::vector <complex_t> ff_;            /* the form factor data */

    public:
      // TODO: clean/unify and improve constructors for gpu/cpu ana/num etc ...
      #ifdef FF_NUM_GPU
        #ifdef FF_NUM_GPU_FUSED
          FormFactor(): numeric_ff_(64, 8), is_analytic_(false) { }
          FormFactor(int a, int b): numeric_ff_(a, b), is_analytic_(false) { }
        #elif defined KERNEL2
          FormFactor(): numeric_ff_(2, 4, 4), is_analytic_(false) { } // default cuda block size
          FormFactor(int a, int b, int c): numeric_ff_(a, b, c), is_analytic_(false) { }
        #else
          FormFactor(): numeric_ff_(64), is_analytic_(false) { } // default cuda block size
          FormFactor(int s): numeric_ff_(s), is_analytic_(false) { }
        #endif // KERNEL2
      #else  // use CPU or MIC
        FormFactor(): numeric_ff_(), is_analytic_(false) { }
      #endif

      ~FormFactor() { }

      void clear(void);

      bool compute_form_factor(ShapeName shape, std::string shape_filename,
                  shape_param_list_t& params,
                  real_t single_thickness,
                  vector3_t& transvec, real_t shp_tau, real_t shp_eta,
                  RotMatrix_t &
                  #ifdef USE_MPI
                    , woo::MultiNode&, std::string
                  #endif
                  );

      complex_t operator[](unsigned int i) const { return ff_[i]; }

      // for testing only ... remove ...
      bool read_form_factor(const char* filename,
                  unsigned int nqx, unsigned int nqy, unsigned int nqz);
      void print_ff(unsigned int nqx, unsigned int nqy, unsigned int nqz);
      void save_ff(unsigned int nqz, const char* filename);
      void save (unsigned, unsigned, const char *);
      void printff(unsigned int nqx, unsigned int nqy, unsigned int nqz) {
        std::cout << "ff:" << std::endl;
        for(unsigned int i = 0; i < nqx * nqy * nqz; ++ i)
          std::cout << ff_[i].real() << "," << ff_[i].imag() << "\t";
        std::cout << std::endl;
      } // printff()

      complex_t* ff(void) { return &ff_[0]; }
  }; // class FormFactor

} // namespace


#endif // _FF_HPP_
