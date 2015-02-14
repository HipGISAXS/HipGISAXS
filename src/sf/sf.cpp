/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf.cpp
 *  Created: Jun 18, 2012
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
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>
#include <sf/sf.hpp>
#include <model/qgrid.hpp> 
#include <utils/utilities.hpp>
#include <common/constants.hpp>

namespace hig {

  StructureFactor::StructureFactor(): sf_(NULL), nrow_(0), ncol_(0) {
  } // StructureFactor::StructureFactor()


/*  StructureFactor::StructureFactor(int nqx, int nqy, int nqz) {
    nx_ = nqx;
    ny_ = nqy;
    nz_ = nqz;
    sf_ = new (std::nothrow) complex_t[nx_ * ny_ * nz_];
    #ifdef SF_GPU
      gsf_.init(nqx, nqy, nqz);
    #endif
  } // StructureFactor::StructureFactor()
*/

  StructureFactor::StructureFactor (int nqx, int nqy) {
    nrow_ = nqx;
    ncol_ = nqy;
    sf_ = new (std::nothrow) complex_t[nrow_ * ncol_];
    #ifdef SF_GPU
      gsf_.init(1, nrow_, ncol_);
    #endif
  } // StructureFactor::StructureFactor()

  StructureFactor & StructureFactor::operator=(const StructureFactor & rhs) {
    /*nx_ = rhs.nx_;
    ny_ = rhs.ny_;
    nz_ = rhs.nz_;
    size_t size = nx_ * ny_ * nz_;*/
    nrow_ = rhs.nrow_;
    ncol_ = rhs.ncol_;
    size_t size = nrow_ * ncol_;
    if(sf_ != NULL) delete[] sf_;
    sf_ = new (std::nothrow) complex_t[size];
    if(sf_ == NULL) {
      std::cerr << "error: could not allocate memeory" << std::endl;
      // TODO should call MPI_Abort
      std::exit(1);
    }
    memcpy(sf_, rhs.sf_, size * sizeof(complex_t));
    #ifdef SF_GPU
      gsf_.init(1, nrow_, ncol_);
    #endif
    return *this;
  } // StructureFactor::operator=()

  StructureFactor & StructureFactor::operator+= (const StructureFactor & rhs) {
    if(nrow_ != rhs.nrow_ || ncol_ != rhs.ncol_ ) {
      std::cerr << "error: cannot add non-similar shaped arrays" << std::endl;
      std::exit(1);
    }
    //int num_of_el = nx_ * ny_ * nz_;
    int num_of_el = nrow_ * ncol_;
    #pragma omp parallel for
    for(int i = 0; i < num_of_el; ++ i) sf_[i] += rhs.sf_[i];
    return *this;
  } // StructureFactor::operator+=()

  StructureFactor::~StructureFactor() {
    if(sf_ != NULL) delete[] sf_;
    sf_ = NULL;
  } // StructureFactor::~StructureFactor()


  void StructureFactor::clear() {
    if(sf_ != NULL) delete[] sf_;
    sf_ = NULL;
    //nx_ = ny_ = nz_ = 0;
    nrow_ = ncol_ = 0;
    #ifdef SF_GPU
      gsf_.destroy();
    #endif
  } // StructureFactor::clear()

  /**
   * compute structure factor on cpu
   */
  bool StructureFactor::compute_structure_factor(std::string expt, vector3_t center,
               Lattice* lattice, vector3_t repet, vector3_t scaling,
               vector3_t rotation_1, vector3_t rotation_2, vector3_t rotation_3
               #ifdef USE_MPI
                 , woo::MultiNode& world_comm, std::string comm_key
               #endif
               ) {
    #ifdef USE_MPI
      bool master = world_comm.is_master(comm_key);
    #else
      bool master = true;
    #endif

    woo::BoostChronoTimer maintimer, computetimer;

    maintimer.start();

    nrow_ = QGrid::instance().nrows();
    ncol_ = QGrid::instance().ncols();
    ny_ = QGrid::instance().nqy();
    if(expt == "saxs") nz_ = QGrid::instance().nqz();
    else if(expt == "gisaxs") nz_ = QGrid::instance().nqz_extended();
    else return false;

    if(repet[0] < 1) repet[0] = 1;
    if(repet[1] < 1) repet[1] = 1;
    if(repet[2] < 1) repet[2] = 1;

    vector3_t arot(0, 0, 0), brot(0, 0, 0), crot(0, 0, 0);
    vector3_t temp_la(lattice->a() * scaling[0]),
              temp_lb(lattice->b() * scaling[1]),
              temp_lc(lattice->c() * scaling[2]);
    //temp_la[2] = 0; temp_lb[2] = 0;
    mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_la, arot);
    mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lb, brot);
    mat_mul_3x1(rotation_1, rotation_2, rotation_3, temp_lc, crot);

    //vector3_t l_t = lattice->t() * scaling;
    vector3_t l_t = lattice->t();

    sf_ = NULL;
    sf_ = new (std::nothrow) complex_t[nz_];
    if(sf_ == NULL) {
      if(master)
        std::cerr << "error: could not allocate memory for structure factor" << std::endl;
      return false;
    } // if

    #ifdef SF_VERBOSE
      if(master) std::cout << "-- Computing structure factor on CPU ... " << std::flush;
    #endif

    std::complex<real_t> unit_c(1, 0);
    std::complex<real_t> unit_ci(0, 1);

    /*#pragma omp parallel for collapse(3)
    for(unsigned int z = 0; z < nz_; ++ z) {
      for(unsigned int y = 0; y < ny_; ++ y) {
        for(unsigned int x = 0; x < nx_; ++ x) {
          complex_t temp1, temp_x2, temp_y3, temp_y4, temp_x5;
          float_t temp_f;
          complex_t sa, sb, sc;
          float_t qx = QGrid::instance().qx(x);
          float_t qy = QGrid::instance().qy(y);
          complex_t qz;
          if(expt == "saxs") qz = QGrid::instance().qz(z);
          else if(expt == "gisaxs") qz = QGrid::instance().qz_extended(z);

          complex_t e_iqa = exp(unit_ci * (arot[0] * qx + arot[1] * qy + arot[2] * qz));
          complex_t Xa_0 = unit_c - pow(e_iqa, repet[0]);
          complex_t Ya_0 = unit_c - e_iqa;
          if(fabs(Ya_0.imag()) > REAL_ZERO_ || fabs(Ya_0.real()) > REAL_ZERO_) sa = Xa_0 / Ya_0;
          else sa = repet[0];
          sa = pow(e_iqa, ((float_t) 1.0 - repet[0]) / (float_t) 2.0) * sa;

          complex_t e_iqb = exp(unit_ci * (brot[0] * qx + brot[1] * qy + brot[2] * qz));
          complex_t Xb_0 = unit_c - pow(e_iqb, repet[1]);
          complex_t Yb_0 = unit_c - e_iqb;
          if(fabs(Yb_0.imag()) > REAL_ZERO_ || fabs(Yb_0.real()) > REAL_ZERO_) sb = Xb_0 / Yb_0;
          else sb = repet[1];
          sb = pow(e_iqb, ((float_t) 1.0 - repet[1]) / (float_t) 2.0) * sb;

          complex_t e_iqc = exp(unit_ci * (crot[0] * qx + crot[1] * qy + crot[2] * qz));
          complex_t Xc_0 = unit_c - pow(e_iqc, repet[2]);
          complex_t Yc_0 = unit_c - e_iqc;
          if(fabs(Yc_0.imag()) > REAL_ZERO_ || fabs(Yc_0.real()) > REAL_ZERO_) sc = Xc_0 / Yc_0;
          else sc = repet[2];
          sc = pow(e_iqc, ((float_t) 1.0 - repet[2]) / (float_t) 2.0) * sc;

          if(!((boost::math::isfinite)(sa.real()) && (boost::math::isfinite)(sa.imag()))) {
            std::cout << "sa sa sa sa sa sa sa: " << x << ", " << y << ", " << z << std::endl; }
          if(!((boost::math::isfinite)(sb.real()) && (boost::math::isfinite)(sb.imag()))) {
            std::cout << "sb sb sb sb sb sb sb: " << x << ", " << y << ", " << z << std::endl; }
          if(!((boost::math::isfinite)(sc.real()) && (boost::math::isfinite)(sc.imag()))) {
            std::cout << "sc sc sc sc sc sc sc: " << x << ", " << y << ", " << z << std::endl; }

          unsigned long int sf_i = nx_ * ny_ * z + nx_ * y + x;
          complex_t temp3 = center[0] * qx + center[1] * qy + center[2] * qz;
          temp3 = exp(complex_t(-temp3.imag(), temp3.real()));
          complex_t temp2 = l_t[0] * qx + l_t[1] * qy + l_t[2] * qz;
          temp2 = unit_c + exp(complex_t(-temp2.imag(), temp2.real()));
          sf_[sf_i] = temp3 * temp2 * sa * sb * sc;

          if(!((boost::math::isfinite)(temp3.real()) && (boost::math::isfinite)(temp3.imag()))) {
            std::cerr << "error: here it is not finite (666) " << x << ", " << y << ", " << z
                  << ": here it is finite (444) " << center[0] << ", "
                  << center[1] << ", " << center[2] << std::endl;
          } // if

          if(!((boost::math::isfinite)(temp2.real()) && (boost::math::isfinite)(temp2.imag()))) {
            std::cerr << "error: here it is not finite (888) " << x << ", " << y << ", " << z
                  << std::endl;
          } // if
        } // for x
      } // for y
      */

    real_t mach_eps = std::numeric_limits<real_t>::epsilon() ; //move to constants.hpp if not there

    computetimer.start();

    /* TODO:
     *  Rotating q-vector is done for every single q-vector.
     *  This can be avoided by rotation of lattice vectors which 
     *  will not change for the grain.
     */

    #pragma omp parallel for 
    for(unsigned int i = 0; i < nz_; ++ i) {
      complex_t temp1, temp_x2, temp_y3, temp_y4, temp_x5;
      unsigned j = i % ny_;
      real_t temp_f;
      complex_t sa, sb, sc;
      real_t qx = QGrid::instance().qx(j);
      real_t qy = QGrid::instance().qy(j);
      complex_t qz;
      if(expt == "saxs") qz = QGrid::instance().qz(i);
      else if(expt == "gisaxs") qz = QGrid::instance().qz_extended(i);

      complex_t e_iqa = exp(unit_ci * (arot[0] * qx + arot[1] * qy + arot[2] * qz));
      complex_t Xa_0 = unit_c - pow(e_iqa, repet[0]);
      complex_t Ya_0 = unit_c - e_iqa;

      real_t tempya = sqrt(Ya_0.real() * Ya_0.real() + Ya_0.imag() * Ya_0.imag());
      if(fabs(Ya_0.imag()) > mach_eps || fabs(Ya_0.real()) > mach_eps) sa = Xa_0 / Ya_0;
      else sa = repet[0];
      sa = pow(e_iqa, ((real_t) 1.0 - repet[0]) / (real_t) 2.0) * sa;

      complex_t iqb = unit_ci *  (brot[0] * qx + brot[1] * qy + brot[2] * qz) ;
      complex_t iqNb =  repet[1] * iqb;

      complex_t e_iqb = exp(iqb);
      complex_t Xb_0 = unit_c - exp(iqNb);
      complex_t Yb_0 = unit_c - exp(iqb);

      real_t tempyb = sqrt(Yb_0.real() * Yb_0.real() + Yb_0.imag() * Yb_0.imag());
      if(fabs(Yb_0.imag()) > mach_eps || fabs(Yb_0.real()) > mach_eps) sb = Xb_0 / Yb_0;
      else sb = repet[1];
      sb = pow(e_iqb, ((real_t) 1.0 - repet[1]) / (real_t) 2.0) * sb;


      complex_t e_iqc = exp(unit_ci * (crot[0] * qx + crot[1] * qy + crot[2] * qz));
      complex_t Xc_0 = unit_c - pow(e_iqc, repet[2]);
      complex_t Yc_0 = unit_c - e_iqc;

      real_t tempyc = sqrt(Yc_0.real() * Yc_0.real() + Yc_0.imag() * Yc_0.imag());
      if(fabs(Yc_0.imag()) > mach_eps || fabs(Yc_0.real()) > mach_eps) sc = Xc_0 / Yc_0;
      else sc = repet[2];
      sc = pow(e_iqc, ((real_t) 1.0 - repet[2]) / (real_t) 2.0) * sc;

      if(!((boost::math::isfinite)(sa.real()) && (boost::math::isfinite)(sa.imag()))) {
        std::cout << "sa sa sa sa sa sa sa: " << i << ", " << j << std::endl; }
      if(!((boost::math::isfinite)(sb.real()) && (boost::math::isfinite)(sb.imag()))) {
        std::cout << "sb sb sb sb sb sb sb: " << i << ", " << j << std::endl; }
      if(!((boost::math::isfinite)(sc.real()) && (boost::math::isfinite)(sc.imag()))) {
        std::cout << "sc sc sc sc sc sc sc: " << i << ", " << j << std::endl; }

      complex_t temp3 = center[0] * qx + center[1] * qy + center[2] * qz;
      temp3 = exp(complex_t(-temp3.imag(), temp3.real()));
      complex_t temp2 = l_t[0] * qx + l_t[1] * qy + l_t[2] * qz;
      temp2 = unit_c + exp(complex_t(-temp2.imag(), temp2.real()));
      sf_[i] = temp3 * temp2 * sa * sb * sc;

      if(!((boost::math::isfinite)(temp3.real()) && (boost::math::isfinite)(temp3.imag()))) {
        std::cerr << "error: here it is not finite (666) " << i << ", " << j
              << ": here it is finite (444) " << center[0] << ", "
              << center[1] << ", " << center[2] << std::endl;
      } // if

      if(!((boost::math::isfinite)(temp2.real()) && (boost::math::isfinite)(temp2.imag()))) {
        std::cerr << "error: here it is not finite (888) " << i << ", " << j << std::endl;
      } // if
    } // for i

    computetimer.stop();
    maintimer.stop();

    #ifdef SF_VERBOSE
      if(master) std::cout << "done. " << std::endl;
    #endif

    if(master) {
      #ifdef TIME_DETAIL_1
        std::cout << "**               SF compute time: "
                << computetimer.elapsed_msec()  << " ms." << std::endl;
        //std::cout << "**                 Total SF time: "
        //          << maintimer.elapsed_msec() << " ms." << std::endl;


        //  save_sf(  QGrid::instance().nqx(),   QGrid::instance().nqy(),   QGrid::instance().nqz(), "/home/stchourou/sf.dat");

        //int naninfs = count_naninfs(nx_, ny_, nz_, sf_);
        //std::cout << " ------- " << naninfs << " / " << nx_ * ny_ * nz_ << " nans or infs" << std::endl;
      #endif // TIME_DETAIL_1
    } // if

    return true;
  } // StructureFactor::compute_structure_factor()


  void StructureFactor::save (const char* filename) {
    std::ofstream f(filename);
    for(unsigned int z = 0; z < nrow_; ++ z) {
      for(unsigned int y = 0; y < ncol_; ++ y) {
          unsigned int index = nrow_ * y + ncol_;
          f << std::norm (sf_[index]) << " ";
      } // for
      f << std::endl;
    } // for
    f.close();
  } // StructureFactor::save_sf()

} // namespace hig
