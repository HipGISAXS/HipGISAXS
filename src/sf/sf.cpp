/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf.cpp
 *  Created: Jun 18, 2012
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
#include <cmath>
#include <boost/math/special_functions/fpclassify.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>
#include <sf/sf.hpp>
#include <model/qgrid.hpp> 
#include <utils/utilities.hpp>
#include <common/constants.hpp>
#include <config/hig_input.hpp>

namespace hig {
  StructureFactor::StructureFactor(){
    ny_ = 0;
    nz_ = 0;
    type_ = default_type;
    sf_ = nullptr;
    // #ifdef SF_GPU
    //  gsf_.init(HiGInput::instance().experiment());
    //#endif
  } // StructureFactor::StructureFactor()

  StructureFactor::~StructureFactor() {
    if(sf_ != NULL) delete[] sf_;
    sf_ = NULL;
  } // StructureFactor::~StructureFactor()

  StructureFactor & StructureFactor::operator=(const StructureFactor & rhs) {
    ny_ = rhs.ny_;
    nz_ = rhs.nz_;
    type_ = rhs.type_;
    if(sf_ != NULL) delete[] sf_;
    sf_ = new (std::nothrow) complex_t[nz_];
    if(sf_ == NULL) {
      std::cerr << "error: could not allocate memeory" << std::endl;
      // TODO should call MPI_Abort or raise exception
      std::exit(1);
    }
    memcpy(sf_, rhs.sf_, nz_ * sizeof(complex_t));
    #ifdef SF_GPU
      gsf_.init(ny_, ny_, nz_);
    #endif
    return *this;
  } // StructureFactor::operator=()

  StructureFactor & StructureFactor::operator+= (const StructureFactor & rhs) {
    if(nz_ != rhs.nz_){
      std::cerr << "error: cannot add non-similar shaped arrays" << std::endl;
      std::exit(1);
    }
    #pragma omp parallel for
    for(int i = 0; i < nz_; ++ i) sf_[i] += rhs.sf_[i];
    return *this;
  } // StructureFactor::operator+=()


  void StructureFactor::clear() {
    if(sf_ != NULL) delete[] sf_;
    sf_ = NULL;
    #ifdef SF_GPU
      gsf_.destroy();
    #endif
  } // StructureFactor::clear()

  /**
   * compute structure factor on cpu
   */
  bool StructureFactor::compute_structure_factor(std::string expt, vector3_t center,
               Lattice* lattice, vector3_t repet, vector3_t scaling,
               RotMatrix_t & rot,
               std::shared_ptr<Paracrystal> pc, std::shared_ptr<PercusYevick> py
               #ifdef USE_MPI
                 , woo::MultiNode& world_comm, std::string comm_key
               #endif
               ) {
    #ifdef USE_MPI
      bool master = world_comm.is_master(comm_key);
    #else
      bool master = true;
    #endif

    ny_ = QGrid::instance().nqy();
    if(expt == "saxs") nz_ = QGrid::instance().nqz();
    else if(expt == "gisaxs") nz_ = QGrid::instance().nqz_extended();
    else return false;

    woo::BoostChronoTimer maintimer, computetimer;
    maintimer.start();

    if(sf_ == NULL) {
      sf_ = new (std::nothrow) complex_t[nz_];
      if(sf_ == NULL) {
        if(master)
          std::cerr << "error: could not allocate memory for structure factor" << std::endl;
        return false;
      }
    } // if

    if (type_ == paracrystal_type){
      if (pc == nullptr){
        std::cerr <<"Error: no paracrystal data found" << std::endl;
        return false;
      }
      if (pc->getDims() == 1){
        if(!paracrystal1d(pc->getDistYMean(), pc->getDistYStddev(), pc->getDomSize())){
          std::cerr << "Error: failed to calculate 1-D Paracrystal" << std::endl;
          return false;
        }
        return true;
      } else if (pc->getDims() == 2){
        switch(lattice->type()){
          case lattice_cubic:
            if(!paracrystal2d_cubic(pc->getDistXMean(), pc->getDistYMean(), 
                  pc->getDistXStddev(), pc->getDistYStddev(), pc->getDomSize())){
              std::cerr<<"Error: falied to calculate 2-D Paracrystal" << std::endl;
              return false;
            }
            break;
          case lattice_hex:
            if(!paracrystal2d_hex(pc->getDistXMean(), pc->getDistYMean(),
                  pc->getDistXStddev(), pc->getDistYStddev(), pc->getDomSize())){
              std::cerr<< "Error: failed to calculate 2-D Paracrystal" << std::endl; 
              return false; 
            }
            break;
          default:
            std::cerr << "Error: Paracrystal is defined for \"cubic\" and \"hex\" lattice types only." << std::endl;
            return false;
        }
        return true;
      } else {
        std::cerr <<"Error: something wrong with paracrystal data" << std::endl;
        return false;
      }
    }
    if (type_ == percusyevick_type){
      if(py == nullptr){
        std::cerr << "Error: no Percus-Yevick data found" << std::endl;
        return false;
      }
      if(!percus_yevik(py->getDiameter(), py->getVolfract(), py->getDims())){
        std::cerr << "Error: failed to calculate Percus-Yevick" << std::endl;
        return false;
      }
      return true;
    }

    if(repet[0] < 1) repet[0] = 1;
    if(repet[1] < 1) repet[1] = 1;
    if(repet[2] < 1) repet[2] = 1;

    vector3_t la(lattice->a() * scaling[0]),
              lb(lattice->b() * scaling[1]),
              lc(lattice->c() * scaling[2]);

    #ifdef SF_VERBOSE
      if(master) std::cerr << "-- Computing structure factor on CPU ... " << std::flush;
    #endif

    std::complex<real_t> unit_c(1, 0);
    std::complex<real_t> unit_ci(0, 1);

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
      std::vector<complex_t> mq = rot.rotate(qx, qy, qz);

      complex_t e_iqa = exp(unit_ci * (la[0] * mq[0] + la[1] * mq[1] + la[2] * mq[2]));
      complex_t Xa_0 = unit_c - pow(e_iqa, repet[0]);
      complex_t Ya_0 = unit_c - e_iqa;

      real_t tempya = sqrt(Ya_0.real() * Ya_0.real() + Ya_0.imag() * Ya_0.imag());
      if(fabs(Ya_0.imag()) > mach_eps || fabs(Ya_0.real()) > mach_eps) sa = Xa_0 / Ya_0;
      else sa = repet[0];
      sa = pow(e_iqa, ((real_t) 1.0 - repet[0]) / (real_t) 2.0) * sa;

      complex_t iqb = unit_ci *  (lb[0] * mq[0] + lb[1] * mq[1] + lb[2] * lb[2]) ;
      complex_t iqNb =  repet[1] * iqb;

      complex_t e_iqb = exp(iqb);
      complex_t Xb_0 = unit_c - exp(iqNb);
      complex_t Yb_0 = unit_c - exp(iqb);

      real_t tempyb = sqrt(Yb_0.real() * Yb_0.real() + Yb_0.imag() * Yb_0.imag());
      if(fabs(Yb_0.imag()) > mach_eps || fabs(Yb_0.real()) > mach_eps) sb = Xb_0 / Yb_0;
      else sb = repet[1];
      sb = pow(e_iqb, ((real_t) 1.0 - repet[1]) / (real_t) 2.0) * sb;


      complex_t e_iqc = exp(unit_ci * (lc[0] * mq[0] + lc[1] * mq[1] + lc[2] * mq[2]));
      complex_t Xc_0 = unit_c - pow(e_iqc, repet[2]);
      complex_t Yc_0 = unit_c - e_iqc;

      real_t tempyc = sqrt(Yc_0.real() * Yc_0.real() + Yc_0.imag() * Yc_0.imag());
      if(fabs(Yc_0.imag()) > mach_eps || fabs(Yc_0.real()) > mach_eps) sc = Xc_0 / Yc_0;
      else sc = repet[2];
      sc = pow(e_iqc, ((real_t) 1.0 - repet[2]) / (real_t) 2.0) * sc;

      if(!((boost::math::isfinite)(sa.real()) && (boost::math::isfinite)(sa.imag()))) {
        std::cerr << "sa sa sa sa sa sa sa: " << i << ", " << j << std::endl; }
      if(!((boost::math::isfinite)(sb.real()) && (boost::math::isfinite)(sb.imag()))) {
        std::cerr << "sb sb sb sb sb sb sb: " << i << ", " << j << std::endl; }
      if(!((boost::math::isfinite)(sc.real()) && (boost::math::isfinite)(sc.imag()))) {
        std::cerr << "sc sc sc sc sc sc sc: " << i << ", " << j << std::endl; }

      complex_t temp3 = center[0] * mq[0] + center[1] * mq[1] + center[2] * mq[2];
      temp3 = exp(complex_t(-temp3.imag(), temp3.real()));
      //complex_t temp2 = l_t[0] * qx + l_t[1] * qy + l_t[2] * qz;
      //temp2 = unit_c + exp(complex_t(-temp2.imag(), temp2.real()));
      sf_[i] = temp3 * sa * sb * sc;

      /************
      if(!((boost::math::isfinite)(temp3.real()) && (boost::math::isfinite)(temp3.imag()))) {
        std::cerr << "error: here it is not finite (666) " << i << ", " << j
              << ": here it is finite (444) " << center[0] << ", "
              << center[1] << ", " << center[2] << std::endl;
      } // if

      if(!((boost::math::isfinite)(temp2.real()) && (boost::math::isfinite)(temp2.imag()))) {
        std::cerr << "error: here it is not finite (888) " << i << ", " << j << std::endl;
      } // if
      ************/

    } // for i

    computetimer.stop();
    maintimer.stop();

    #ifdef SF_VERBOSE
      if(master) std::cerr << "done. " << std::endl;
    #endif

    if(master) {
      #ifdef TIME_DETAIL_1
        std::cerr << "**               SF compute time: "
                << computetimer.elapsed_msec()  << " ms." << std::endl;
        //std::cerr << "**                 Total SF time: "
        //          << maintimer.elapsed_msec() << " ms." << std::endl;


        //  save_sf(  QGrid::instance().nqx(),   QGrid::instance().nqy(),   QGrid::instance().nqz(), "/home/stchourou/sf.dat");

        //int naninfs = count_naninfs(nx_, ny_, nz_, sf_);
        //std::cerr << " ------- " << naninfs << " / " << nx_ * ny_ * nz_ << " nans or infs" << std::endl;
      #endif // TIME_DETAIL_1
    } // if

    return true;
  } // StructureFactor::compute_structure_factor()


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
    gsf_.init(ny_, ny_, nz_);
    bool ret = gsf_.compute(expt, center, lattice, repet, scaling, rot
                        //#ifdef USE_MPI
                        //  , world_comm, comm_key
                        //#endif
                       );
    gsf_.get_sf(sf_);
    gsf_.destroy();
    return ret;
  } // StructureFactor::compute_structure_factor_gpu()
#endif // SF_GPU


  void StructureFactor::save (const char* filename) {
    std::ofstream f(filename);
    for(int i = 0; i < nz_; i++){
      f << sf_[i].real() << ", " << sf_[i].imag() << std::endl;
    } // for
    f.close();
  } // StructureFactor::save_sf()

  void StructureFactor::save_sf(const std::string& filename) {
    std::ofstream f(filename.c_str());
    for(unsigned int z = 0; z < nz_; ++ z) {
      f << std::abs(sf_[z]) << std::endl;
    } // for
    f.close();
  } // StructureFactor::save_sf()
} // namespace hig
