/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff.cpp
 *  Created: Jul 17, 2012
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

#include <iostream>
#include <fstream>

#include <ff/ff.hpp>

namespace hig {

  void FormFactor::clear() {
    ff_.clear();
    analytic_ff_.clear();
    numeric_ff_.clear();
    is_analytic_ = false;
  } // FormFactor::clear()


  bool FormFactor::compute_form_factor(ShapeName shape, std::string shape_filename,
                    shape_param_list_t& params, real_t single_thickness,
                    vector3_t& transvec, real_t shp_tau, real_t shp_eta,
                    RotMatrix_t & rot
                    #ifdef USE_MPI
                      , woo::MultiNode& multi_node, std::string comm_key
                    #endif
                    ) {
    if(shape == shape_custom) {
      /* compute numerically */
      is_analytic_ = false;
      numeric_ff_.init(rot, ff_);
      numeric_ff_.compute(shape_filename.c_str(), ff_, rot
                #ifdef USE_MPI
                  , multi_node, comm_key
                #endif
                );
    } else {
      /* compute analytically */
      is_analytic_ = true;
      analytic_ff_.init(rot, ff_);
      analytic_ff_.compute(shape, shp_tau, shp_eta, transvec,
                  ff_, params, single_thickness, rot
                  #ifdef USE_MPI
                    , multi_node, comm_key
                  #endif
                  );
    } // if-else
    return true;
  } // FormFactor::compute_form_factor()


  /* temporaries */

  bool FormFactor::read_form_factor(const char* filename,
                    unsigned int nqx, unsigned int nqy, unsigned int nqz) {
    std::ifstream f(filename);
    ff_.clear();
    for(unsigned int z = 0; z < nqz; ++ z) {
      for(unsigned int y = 0; y < nqy; ++ y) {
        for(unsigned int x = 0; x < nqx; ++ x) {
          real_t tempr, tempi;
          f >> tempr; f >> tempi;
          ff_.push_back(complex_t(tempr, tempi));
        } // for
      } // for
    } // for
    f.close();

    return true;
  } // FormFactor::read_form_factor()


  void FormFactor::print_ff(unsigned int nqx, unsigned int nqy, unsigned int nqz) {
    for(unsigned int z = 0; z < nqz; ++ z) {
      for(unsigned int y = 0; y < nqy; ++ y) {
        for(unsigned int x = 0; x < nqx; ++ x) {
          unsigned int index = nqx * nqy * z + nqx * y + x;
          std::cout << ff_[index].real() << "," << ff_[index].imag() << "\t";
        } // for
        std::cout << std::endl;
      } // for
      std::cout << std::endl;
    } // for
  } // FormFactor::print_ff()


  void FormFactor::save (unsigned nrow, unsigned ncol, const char * filename) {
    int size = 4 * ncol * nrow;
    std::ofstream out (filename);
    for (unsigned i = 0; i < size; i++){
      out << ff_[i].real() << ", " << ff_[i].imag() << std::endl;
      //if ((i+1)%ncol == 0)
      //  out << std::endl;
    }
    out.close();
  }

  void FormFactor::save_ff(unsigned int nqz, const char* filename) {
    std::ofstream f(filename);
    for(unsigned int z = 0; z < nqz; ++ z) {
      f << std::abs(ff_[z]) << std::endl;
    } // for
    f.close();
  } // FormFactor::save_ff()

} // namespace hig
