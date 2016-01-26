/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_main.cpp
 *  Created: Jun 14, 2012
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

#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iomanip>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <woo/timer/woo_boostchronotimers.hpp>
#include <woo/random/woo_mtrandom.hpp>

#include <sim/hipgisaxs_main.hpp>
#include <common/typedefs.hpp>
#include <utils/utilities.hpp>
#include <numerics/matrix.hpp>
#include <numerics/numeric_utils.hpp>

#if defined USE_GPU || defined FF_ANA_GPU || defined FF_NUM_GPU
  #include <init/gpu/init_gpu.cuh>
#elif defined USE_MIC
  #include <init/mic/init_mic.hpp>
#endif


namespace hig {

  HipGISAXS::HipGISAXS(int narg, char** args): freq_(0.0), k0_(0.0),
        num_layers_(0), num_structures_(0),
        nqx_(0), nqy_(0), nqz_(0), nqz_extended_(0)
        #ifdef USE_MPI
          , multi_node_(narg, args)
        #endif
        //#ifdef FF_NUM_GPU   // use GPU
        //  #ifdef FF_NUM_GPU_FUSED
        //    ff_(64, 8)
        //  #elif defined KERNEL2
        //    ff_(2, 4, 4)
        //  #else
        //    ff_(64)
        //  #endif
        //#else   // use CPU
        //  ff_()
        //#endif
          {
    single_layer_refindex_.delta(0.0);
    single_layer_refindex_.beta(0.0);
    single_layer_thickness_ = 0.0;
    // construct hig input instance
    HiGInput::instance();
    // construct the qgrid -- to be changed to be not singleton ...
    QGrid::instance();
  } // HipGISAXS::HipGISAXS()


  HipGISAXS::~HipGISAXS() {
    // nothing to do here yet ...
  } // HipGISAXS::~HipGISAXS()


  bool HipGISAXS::init() {
            // is called at the beginning of the runs (after input is read)
            // it does the following:
            //   + set detector/system stuff
            //   + initialize output dir
            //   + create q-grid
            //   + construct layer profile
            //   + get layer profile data
            //   + get structure info
            //   + construct lattice vectors
            //   + compute cell size
    // TODO first check if the input has been constructed ...

    #ifdef USE_MPI
      root_comm_ = multi_node_.universe_key();
      int mpi_rank = multi_node_.rank(root_comm_);
      bool master = multi_node_.is_master(root_comm_);
    #else
      root_comm_ = "world";    // doesnt matter in this case
      int mpi_rank = 0;
      bool master = true;
    #endif
    // initially, the simulator communicator is same as the universe
    sim_comm_ = root_comm_;

    if(master) {
      std::cout << std::endl
          << "*******************************************************************" << std::endl
          << "******************** HipGISAXS v1.0-beta Testing ******************" << std::endl
          << "*******************************************************************" << std::endl
          << std::endl;
    } // if

    //photon conversion
    real_t photon = 0.0;
    real_t lambda = 0.0;
    std::string unit;
    freq_ = 0; k0_ = 0;
    HiGInput::instance().photon_energy(photon, unit);
    if(unit == "ev") {
    //  freq_ = photon * 1.60217646e9 / 6.626068;
      photon = photon / 1000;    // in keV
      lambda = 1.23984 / photon;
      freq_ = 1e-9 * photon * 1.60217646e-19 * 1000 / 6.626068e-34;
    } else if (unit == "kev") {
      lambda = 1.23984 / photon;
      freq_ = 1e-9 * photon * 1.60217646e-19 * 1000 / 6.626068e-34;
    } else { /* do something else ? ... */
      if(master) std::cerr << "error: photon energy is not given in 'ev' or 'kev'" << std::endl;
      return false;
    } // if-else

    //k0_ = 2 * PI_ * freq_ / LIGHT_SPEED_;
    std::cout << "**                    Wavelength: " << lambda << std::endl;
    k0_ = 2 * PI_ / lambda;

    // create output directory
    if(master) {    // this is not quite good for mpi ... improve ...
      const std::string p = HiGInput::instance().path() + "/" + HiGInput::instance().runname();
      if(!boost::filesystem::create_directory(p)) {
        std::cerr << "error: could not create output directory " << p << std::endl;
        return false;
      } // if
    } // if

    #ifdef USE_MPI
      multi_node_.barrier(root_comm_);
    #endif

    // create Q-grid
    real_t min_alphai = HiGInput::instance().scattering_min_alpha_i() * PI_ / 180;
    if(!QGrid::instance().create(freq_, min_alphai, k0_, mpi_rank)) {
      if(master) std::cerr << "error: could not create Q-grid" << std::endl;
      return false;
    } // if

    nrow_ = QGrid::instance().nrows();
    ncol_ = QGrid::instance().ncols();
    nqx_ = QGrid::instance().nqx();
    nqy_ = QGrid::instance().nqy();
    nqz_ = QGrid::instance().nqz();
    nqz_extended_ = QGrid::instance().nqz_extended();

    /* construct layer profile */
    if(!HiGInput::instance().construct_layer_profile()) {  // also can be done at input reading ...
      if(master) std::cerr << "error: could not construct layer profile" << std::endl;
      return false;
    } // if

    /* get initialization data from layers */

    if(HiGInput::instance().has_substrate_layer()) {  //
      substrate_refindex_ = HiGInput::instance().substrate_layer().refindex();
    } else {
      substrate_refindex_.delta(0); substrate_refindex_.beta(0);
    } // if-else


    num_layers_ = HiGInput::instance().num_layers();  // this excludes substrate layer
    if(num_layers_ == 1) {    // is this really needed? ...
      single_layer_refindex_ = HiGInput::instance().single_layer().refindex();
      single_layer_thickness_ = HiGInput::instance().single_layer().thickness();
    } else {
      // initialize
      single_layer_refindex_.delta(0); single_layer_refindex_.beta(0);
      single_layer_thickness_ = 0.0;
    } // if-else

    multilayer_.init();
    complex_t temp(substrate_refindex_.delta(), substrate_refindex_.beta());
    dns2_ = ((real_t) 2.0) * temp - (complex_t) pow(temp, 2);

    /* get initialization data from structures */
    
    num_structures_ = HiGInput::instance().num_structures();
    /* construct lattice vectors in each structure, if needed */
    // TODO: this can also be done at input reading ...
    if(!HiGInput::instance().construct_lattice_vectors()) {
      if(master) std::cerr << "error: could not construct lattice vectors" << std::endl;
      return false;
    } // if

    /* domain size */
    vector3_t min_vec(0.0, 0.0, 0.0), max_vec(0.0, 0.0, 0.0);
    real_t z_min_0 = 0.0, z_max_0 = 0.0;
    if(!HiGInput::instance().compute_domain_size(min_vec, max_vec, z_min_0, z_max_0)) {
      if(master) std::cerr << "error: could not compute domain size" << std::endl;
      return false;
    } // if

    // cell size ... check ... FIXME
    cell_[0] = fabs(max_vec[0] - min_vec[0]);
    cell_[1] = fabs(max_vec[1] - min_vec[1]);
    cell_[2] = fabs(z_max_0 - z_min_0);


    #ifdef _OPENMP
      if(master)
        std::cout << "++      Number of OpenMP threads: "
              << omp_get_max_threads() << std::endl;
    #endif

    #if defined USE_GPU || defined FF_ANA_GPU || defined FF_NUM_GPU
      if(master) std::cout << "-- Waking up GPU(s) ..." << std::flush;
      init_gpu();
      if(master) std::cout << " it woke up!" << std::endl;
    #elif defined USE_MIC
      if(master) std::cout << "-- Waking up MIC(s) ..." << std::flush;
      init_mic();
      if(master) std::cout << " done." << std::endl;
    #else
      if(master) std::cout << "-- Not set up to use any accelerator!" << std::endl;
    #endif

    return true;
  } // HipGISAXS::init()


  bool HipGISAXS::run_init(real_t alphai, real_t phi, real_t tilt, SampleRotation& rot_matrix) {
          // this is called for each config-run during the main run
          // it does the following:
          //   + construct the illuminated volume
          //  + construct rotation matrices

    #ifdef USE_MPI
      bool master = multi_node_.is_master(sim_comm_);
    #else
      bool master = true;
    #endif

    /*
    if(!illuminated_volume(alphai, HiGInput::instance().scattering_spot_area(),
        HiGInput::instance().min_layer_order(), HiGInput::instance().substrate_refindex())) {
      if(master) std::cerr << "error: something went wrong in illuminated_volume()" << std::endl;
      return false;
    } // if
    */

    /* rotation matrices */
    // compute r_phi = sample rotation by phi
    vector3_t rp1(0.0, 0.0, 0.0), rp2(0.0, 0.0, 0.0), rp3(0.0, 0.0, 0.0);
    vector3_t rt1(0.0, 0.0, 0.0), rt2(0.0, 0.0, 0.0), rt3(0.0, 0.0, 0.0);
    compute_rotation_matrix_z(phi, rp1, rp2, rp3);
    compute_rotation_matrix_x(tilt, rt1, rt2, rt3);
    mat_mul_3x3(rp1, rp2, rp3, rt1, rt2, rt3, rot_matrix.r1_, rot_matrix.r2_, rot_matrix.r3_);

    return true;
  } // HipGISAXS::run_init()


  // for fitting, update the qgrid when they are different
  bool HipGISAXS::override_qregion(unsigned int ny, unsigned int nz, unsigned int i) {

    OutputRegionType type = HiGInput::instance().reference_region_type(i);
    real_t miny = HiGInput::instance().reference_region_min_x(i);
    real_t minz = HiGInput::instance().reference_region_min_y(i);
    real_t maxy = HiGInput::instance().reference_region_max_x(i);
    real_t maxz = HiGInput::instance().reference_region_max_y(i);

    #ifdef USE_MPI
      // this is done at the universal level
      int mpi_rank = multi_node_.rank(root_comm_);
      bool master = multi_node_.is_master(root_comm_);
    #else
      int mpi_rank = 0;
      bool master = true;
    #endif

    if(type == region_qspace) {
      // update Q-grid
      real_t min_alphai = HiGInput::instance().scattering_min_alpha_i() * PI_ / 180;
      if(!QGrid::instance().update(ny, nz, miny, minz, maxy, maxz,
                                   freq_, min_alphai, k0_, mpi_rank)) {
        if(master) std::cerr << "error: could not update Q-grid" << std::endl;
        return false;
      } // if

      nrow_ = QGrid::instance().nrows();
      ncol_ = QGrid::instance().ncols();
      nqx_ = QGrid::instance().nqx();
      nqy_ = QGrid::instance().nqy();
      nqz_ = QGrid::instance().nqz();
      nqz_extended_ = QGrid::instance().nqz_extended();

    } else if(type == region_pixels) {
      std::cerr << "uh-oh: override option for pixels has not yet been implemented" << std::endl;
      return false;
    } else if(type == region_angles) {
      std::cerr << "uh-oh: override option for angles has not yet been implemented" << std::endl;
      return false;
    } // if-else
    return true;
  } // HipGISAXS::override_qregion()

  /**
   * This is the main function called from outside
   * It loops over all configurations and calls the simulation routine
   */
  bool HipGISAXS::run_all_gisaxs(int x_min, int x_max, int x_step) {
    if(!init()) return false;

    #ifdef USE_MPI
      // this is for the comm world for this simulation
      bool master = multi_node_.is_master(sim_comm_);
    #else
      bool master = true;
    #endif

    woo::BoostChronoTimer sim_timer;
    sim_timer.start();

    int num_alphai = 0, num_phi = 0, num_tilt = 0;;

    real_t alphai_min, alphai_max, alphai_step;
    HiGInput::instance().scattering_alphai(alphai_min, alphai_max, alphai_step);
    if(alphai_max < alphai_min) alphai_max = alphai_min;
    if(alphai_min == alphai_max || alphai_step == 0) num_alphai = 1;
    else num_alphai = (alphai_max - alphai_min) / alphai_step + 1;

    real_t phi_min, phi_max, phi_step;
    HiGInput::instance().scattering_inplanerot(phi_min, phi_max, phi_step);
    if(phi_step == 0) num_phi = 1;
    else num_phi = (phi_max - phi_min) / phi_step + 1;

    real_t tilt_min, tilt_max, tilt_step;
    HiGInput::instance().scattering_tilt(tilt_min, tilt_max, tilt_step);
    if(tilt_step == 0) num_tilt = 1;
    else num_tilt = (tilt_max - tilt_min) / tilt_step + 1;

    if(master) {
      std::cout << "**                  Num alphai: " << num_alphai << std::endl
            << "**                     Num phi: " << num_phi << std::endl
            << "**                    Num tilt: " << num_tilt << std::endl;
    } // if

    // loop over all alphai, phi, and tilt

    #ifdef USE_MPI
      // divide among processors
      int num_procs = multi_node_.size(sim_comm_);
      int rank = multi_node_.rank(sim_comm_);
      int alphai_color = 0;
      if(num_procs > num_alphai) {
        alphai_color = rank % num_alphai;
        alphai_min = alphai_min + alphai_step * alphai_color;
        num_alphai = 1;
      } else {
        alphai_color = rank;
        alphai_min = alphai_min + alphai_step *
              ((num_alphai / num_procs) * rank + min(rank, num_alphai % num_procs));
        num_alphai = (num_alphai / num_procs) + (rank < num_alphai % num_procs);
      } // if-else
      woo::comm_t alphai_comm = "alphai";
      multi_node_.split(alphai_comm, sim_comm_, alphai_color);

      bool amaster = multi_node_.is_master(alphai_comm);
      int temp_amaster = amaster;
      int *amasters = new (std::nothrow) int[multi_node_.size(sim_comm_)];
      // all alphai masters tell the world master about who they are
      multi_node_.allgather(sim_comm_, &temp_amaster, 1, amasters, 1);
    #else
      bool amaster = true;
    #endif // USE_MPI

    // for each incidence angle
    real_t alpha_i = alphai_min;
    for(int i = 0; i < num_alphai; i ++, alpha_i += alphai_step) {
      real_t alphai = alpha_i * PI_ / 180;

      real_t* averaged_data = NULL;    // to hold summation of all (if needed)

      #ifdef USE_MPI
        // divide among processors
        int num_procs = multi_node_.size(alphai_comm);
        int rank = multi_node_.rank(alphai_comm);
        int phi_color = 0;
        if(num_procs > num_phi) {
          phi_color = rank % num_phi;
          phi_min = phi_min + phi_step * phi_color;
          num_phi = 1;
        } else {
          phi_color = rank;
          phi_min = phi_min + phi_step *
                ((num_phi / num_procs) * rank + min(rank, num_phi % num_procs));
          num_phi = (num_phi / num_procs) + (rank < num_phi % num_procs);
        } // if-else
        woo::comm_t phi_comm = "phi";
        multi_node_.split(phi_comm, alphai_comm, phi_color);

        bool pmaster = multi_node_.is_master(phi_comm);
        int temp_pmaster = pmaster;
        int *pmasters = new (std::nothrow) int[multi_node_.size(alphai_comm)];
        // all phi masters tell the alphai master about who they are
        multi_node_.allgather(alphai_comm, &temp_pmaster, 1, pmasters, 1);
      #else
        bool pmaster = true;
      #endif // USE_MPI

      // for each inplane rotation angle
      real_t phi = phi_min;
      for(int j = 0; j < num_phi; j ++, phi += phi_step) {
        real_t phi_rad = phi * PI_ / 180;

        #ifdef USE_MPI
          // divide among processors
          int num_procs = multi_node_.size(phi_comm);
          int rank = multi_node_.rank(phi_comm);
          int tilt_color = 0;
          if(num_procs > num_tilt) {
            tilt_color = rank % num_tilt;
            tilt_min = tilt_min + tilt_step * tilt_color;
            num_tilt = 1;
          } else {
            tilt_color = rank;
            tilt_min = tilt_min + tilt_step *
                ((num_tilt / num_procs) * rank + min(rank, num_tilt % num_procs));
            num_tilt = (num_tilt / num_procs) + (rank < num_tilt % num_procs);
          } // if-else
          woo::comm_t tilt_comm = "tilt";
          multi_node_.split(tilt_comm, phi_comm, tilt_color);

          bool tmaster = multi_node_.is_master(tilt_comm);
          int temp_tmaster = tmaster;
          int *tmasters = new (std::nothrow) int[multi_node_.size(phi_comm)];
          // all tilt masters tell the phi master about who they are
          multi_node_.allgather(phi_comm, &temp_tmaster, 1, tmasters, 1);
        #else
          bool tmaster = true;
        #endif // USE_MPI

        // for each tilt angle
        real_t tilt = tilt_min;
        for(int k = 0; k < num_tilt; k ++, tilt += tilt_step) {
          real_t tilt_rad = tilt * PI_ / 180;

          if(tmaster) {
            std::cout << "-- Computing GISAXS "
                  << i * num_phi * num_tilt + j * num_tilt + k + 1 << " / "
                  << num_alphai * num_phi * num_tilt
                  << " [alphai = " << alpha_i << ", phi = " << phi
                  << ", tilt = " << tilt << "] ..." << std::endl << std::flush;
          } // if

          /* run a gisaxs simulation */

          real_t* final_data = NULL;
          if(!run_gisaxs(alpha_i, alphai, phi_rad, tilt_rad, final_data,
                #ifdef USE_MPI
                  tilt_comm,
                #endif
                0)) {
            if(tmaster)
              std::cerr << "error: could not finish successfully" << std::endl;
            return false;
          } // if

          if(tmaster) {
            std::cout << "-- Constructing GISAXS image ... " << std::flush;
            Image img(ncol_, nrow_, HiGInput::instance().palette());
            img.construct_image(final_data, 0); // merge this into the contructor ...
            std::cout << "done." << std::endl;

            if(x_max < x_min) x_max = x_min;
            // define output filename
            std::stringstream alphai_b, phi_b, tilt_b;
            std::string alphai_s, phi_s, tilt_s;
            alphai_b << alpha_i; alphai_s = alphai_b.str();
            phi_b << phi; phi_s = phi_b.str();
            tilt_b << tilt; tilt_s = tilt_b.str();
            std::string output(HiGInput::instance().param_pathprefix() +
                      "/" + HiGInput::instance().runname() +
                      "/img_ai=" + alphai_s + "_rot=" + phi_s +
                      "_tilt=" + tilt_s + ".tif");

            //img.save(output, x_min, x_max);
            std::cout << "**                    Image size: " << ncol_  << " x " << nrow_
                  << std::endl;
            std::cout << "-- Saving image in " << output << " ... " << std::flush;
            img.save(output);
            std::cout << "done." << std::endl;

            // save the actual data into a file also
            std::string data_file(HiGInput::instance().param_pathprefix() +
                    "/" + HiGInput::instance().runname() +
                    "/gisaxs_ai=" + alphai_s + "_rot=" + phi_s +
                    "_tilt=" + tilt_s + ".out");
            std::cout << "-- Saving raw data in " << data_file << " ... "
                << std::flush;
            save_gisaxs(final_data, data_file);
            std::cout << "done." << std::endl;
          } // if

          // also compute averaged values over phi and tilt
          if(num_phi > 1 || num_tilt > 1) {
            if(tmaster) {
              if(averaged_data == NULL) {
                averaged_data = new (std::nothrow) real_t[nrow_ * ncol_];
                memset(averaged_data, 0, nrow_ * ncol_ * sizeof(real_t));
              } // if
              int i = 0;
              add_data_elements(averaged_data, final_data, averaged_data, nrow_ * ncol_);
            } // if
            delete[] final_data;
          } else {
            averaged_data = final_data;
          } // if-else

          #ifdef USE_MPI
            multi_node_.barrier(tilt_comm);
          #endif

        } // for tilt
        #ifdef USE_MPI
          multi_node_.free(tilt_comm);

//          if(num_phi > 1 || num_tilt > 1) {
            // get data from all other phi_comm processors which were tmasters
            int psize = multi_node_.size(phi_comm);
            real_t* temp_data = NULL;
            if(psize > 1) {
              int *proc_sizes = new (std::nothrow) int[psize];
              int *proc_displacements = new (std::nothrow) int[psize];
              if(pmaster) {
                temp_data = new (std::nothrow) real_t[psize * nrow_ * ncol_];
              } // if
              for(int i = 0; i < psize; ++ i) {
                proc_sizes[i] = tmasters[i] * nrow_ * ncol_;
              } // for
              proc_displacements[0] = 0;
              for(int i = 1; i < psize; ++ i) {
                proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
              } // for
              int prank = multi_node_.rank(phi_comm);
              multi_node_.gatherv(phi_comm, averaged_data, proc_sizes[prank],
                        temp_data, proc_sizes, proc_displacements);
              if(pmaster) {
                for(int i = 1; i < psize; ++ i) {
                  add_data_elements(averaged_data, temp_data + proc_displacements[i],
                            averaged_data, proc_sizes[i]);
                } // for
                delete[] temp_data;
              } // if
              delete[] proc_displacements;
              delete[] proc_sizes;
            } // if
//          } // if

          multi_node_.barrier(phi_comm);

        #endif
      } // for phi
      #ifdef USE_MPI
        multi_node_.free(phi_comm);

//        if(num_phi > 1 || num_tilt > 1) {
          // get data from all other phi_comm processors which were tmasters
          int asize = multi_node_.size(alphai_comm);
          real_t* temp_data = NULL;
          if(asize > 1) {
            int *proc_sizes = new (std::nothrow) int[asize];
            int *proc_displacements = new (std::nothrow) int[asize];
            if(amaster) {
              temp_data = new (std::nothrow) real_t[asize * nrow_ * ncol_];
            } // if
            for(int i = 0; i < asize; ++ i) {
              proc_sizes[i] = pmasters[i] * nrow_ * ncol_;
            } // for
            proc_displacements[0] = 0;
            for(int i = 1; i < asize; ++ i) {
              proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
            } // for
            int arank = multi_node_.rank(alphai_comm);
            multi_node_.gatherv(alphai_comm, averaged_data, proc_sizes[arank],
                      temp_data, proc_sizes, proc_displacements);
            if(amaster) {
              for(int i = 1; i < asize; ++ i) {
                add_data_elements(averaged_data, temp_data + proc_displacements[i],
                          averaged_data, nrow_ * ncol_);
              } // for
              delete[] temp_data;
            } // if
            delete[] proc_displacements;
            delete[] proc_sizes;
          } // if
//        } // if

        multi_node_.barrier(alphai_comm);

      #endif

      if(amaster && (num_phi > 1 || num_tilt > 1)) {
        if(averaged_data != NULL) {
          Image img(ncol_, nrow_);
          img.construct_image(averaged_data, 0); // slice x = 0

          // define output filename
          std::stringstream alphai_b;
          std::string alphai_s;
          alphai_b << alpha_i; alphai_s = alphai_b.str();
          std::string output(HiGInput::instance().param_pathprefix() +
                    "/" + HiGInput::instance().runname() +
                    "/img_ai=" + alphai_s + "_averaged.tif");
          std::cout << "-- Saving averaged image in " << output << " ... " << std::flush;
          img.save(output);
          std::cout << "done." << std::endl;

          // save the actual data into a file also
          std::string data_file(HiGInput::instance().param_pathprefix() +
                  "/" + HiGInput::instance().runname() +
                  "/gisaxs_ai=" + alphai_s + "_averaged.out");
          std::cout << "-- Saving averaged raw data in " << data_file << " ... " << std::flush;
          save_gisaxs(averaged_data, data_file);
          std::cout << "done." << std::endl;

          delete[] averaged_data;
        } // if
      } // if

    } // for alphai
    #ifdef USE_MPI
      multi_node_.free(alphai_comm);
    #endif

    sim_timer.stop();
    if(master) {
      std::cout << "**         Total simulation time: " << sim_timer.elapsed_msec() << " ms."
            << std::endl;
    } // if

    return true;
  } // HipGISAXS::run_all_gisaxs()


  /**
   * used in fitting
   */

  bool HipGISAXS::fit_init() { return init(); }


  bool HipGISAXS::compute_gisaxs(real_t* &final_data, woo::comm_t comm_key) {
    if(!comm_key.empty()) sim_comm_ = comm_key;        // communicator for this simulation
    #ifdef USE_MPI
      bool master = multi_node_.is_master(sim_comm_);
    #else
      bool master = true;
    #endif

    int num_alphai = 0, num_phi = 0, num_tilt = 0;
    real_t alphai_min, alphai_max, alphai_step;
    HiGInput::instance().scattering_alphai(alphai_min, alphai_max, alphai_step);
    if(alphai_max < alphai_min) alphai_max = alphai_min;
    if(alphai_min == alphai_max || alphai_step == 0) num_alphai = 1;
    else num_alphai = (alphai_max - alphai_min) / alphai_step + 1;
    real_t phi_min, phi_max, phi_step;
    HiGInput::instance().scattering_inplanerot(phi_min, phi_max, phi_step);
    if(phi_step == 0) num_phi = 1;
    else num_phi = (phi_max - phi_min) / phi_step + 1;
    real_t tilt_min, tilt_max, tilt_step;
    HiGInput::instance().scattering_tilt(tilt_min, tilt_max, tilt_step);
    if(tilt_step == 0) num_tilt = 1;
    else num_tilt = (tilt_max - tilt_min) / tilt_step + 1;
    if(num_alphai > 1 || num_phi > 1 || num_tilt > 1) {
      if(master)
        std::cerr << "error: currently you can simulate only for single "
              << "alpha_i, phi and tilt angles"
              << std::endl;
      // TODO ...
      return -1.0;
    } // if

    woo::BoostChronoTimer sim_timer;
    sim_timer.start();
    real_t alpha_i = alphai_min;
    real_t alphai = alpha_i * PI_ / 180;
    real_t phi_rad = phi_min * PI_ / 180;
    real_t tilt_rad = tilt_min * PI_ / 180;
    if(master) std::cout << "-- Computing GISAXS ... " << std::endl << std::flush;
    /* run a gisaxs simulation */
    if(!run_gisaxs(alpha_i, alphai, phi_rad, tilt_rad, final_data,
          #ifdef USE_MPI
            sim_comm_,
          #endif
          0)) {
      if(master) std::cerr << "error: could not finish successfully" << std::endl;
      return -1.0;
    } // if
    sim_timer.stop();
    if(master)
      std::cout << "**        Total Simulation time: " << sim_timer.elapsed_msec()
            << " ms." << std::endl;

    // temporary for testing ...
    /*if(master) {
      Image img(ncol_, nrow_);
      img.construct_image(final_data, 0); // slice x = 0
      // define output filename
      std::stringstream alphai_b;
      std::string alphai_s;
      alphai_b << alpha_i; alphai_s = alphai_b.str();
      std::string output(HiGInput::instance().param_pathprefix() +
                "/" + HiGInput::instance().runname() +
                "/img_ai=" + alphai_s + ".tif");
      std::cout << "-- Saving image in " << output << " ... " << std::flush;
      img.save(output);
      std::cout << "done." << std::endl;
    } // if*/

    return true;
  } // HipGISAXS::compute_gisaxs()

  
  /**
   * run an experiment configuration
   * called for each configuration
   */

  /* all the real juice is here */
  bool HipGISAXS::run_gisaxs(real_t alpha_i, real_t alphai, real_t phi, real_t tilt,
                real_t* &img3d,
                #ifdef USE_MPI
                  woo::comm_t comm_key,
                #endif
                int corr_doms) {

    //SampleRotation rotation_matrix;
    // if(!run_init(alphai, phi, tilt, rotation_matrix)) return false;
    rot_ = RotMatrix_t(2, phi);

    //QGrid::instance().save ("qgrid.out");
    #ifdef USE_MPI
      bool master = multi_node_.is_master(comm_key);
      int ss = multi_node_.size(comm_key);
      //std::cout << "****************** MPI size for this simulation: " << ss << std::endl;
    #else
      bool master = true;
    #endif

    /**************************
    // compute propagation coefficients/fresnel coefficients
    // this can also go into run_init() ...
    complex_t *amm = NULL, *apm = NULL, *amp = NULL, *app = NULL;
    complex_t *rk1 = NULL, *rk2 = NULL, *rk1rk2 = NULL, *tk1tk2 = NULL, *h0 = NULL;
    // TODO : where is fc used ?????????? ....................
    complex_t* fc = NULL;
    if(!compute_propagation_coefficients(alphai, amm, apm, amp, app,
                        rk1, rk2, rk1rk2, tk1tk2, h0, fc)) {
      if(master) std::cerr << "error: failed to compute propogation coefficients" << std::endl;
      return false;
    } // if
    ****************************/

    /* loop over all structures and grains/grains */
    structure_iterator_t s = HiGInput::instance().structure_begin();
    int num_structs = num_structures_;
    #ifdef USE_MPI
      // divide among processors
      int num_procs = multi_node_.size(comm_key);
      int rank = multi_node_.rank(comm_key);
      int struct_color = 0;
      int soffset = 0;
      if(num_procs > num_structs) {
        struct_color = rank % num_structs;
        soffset = struct_color;
        num_structs = 1;
      } else {
        struct_color = rank;
        soffset = ((num_structs / num_procs) * rank + min(rank, num_structs % num_procs));
        num_structs = (num_structs / num_procs) + (rank < num_structs % num_procs);
      } // if-else
      woo::comm_t struct_comm = "structure";
      multi_node_.split(struct_comm, comm_key, struct_color);
      // goto structures i am responsible for
      for(int i = 0; i < soffset; ++ i) ++ s;

      bool smaster = multi_node_.is_master(struct_comm);
      int temp_smaster = smaster;
      int *smasters = new (std::nothrow) int[multi_node_.size(comm_key)];
      // all grain masters tell the structure master about who they are
      multi_node_.gather(comm_key, &temp_smaster, 1, smasters, 1);
    #else
      bool smaster = true;
    #endif // USE_MPI

    // initialize memory for struct_intensity
    unsigned int size = nrow_ * ncol_;
    real_t* struct_intensity = NULL;
    complex_t* c_struct_intensity = NULL;
    if(smaster) {
      struct_intensity = new (std::nothrow) real_t[num_structs * size];
      c_struct_intensity = new (std::nothrow) complex_t[num_structs * size];
    } // master

    // loop over structures
    for(int s_num = 0; s_num < num_structs && s != HiGInput::instance().structure_end();
        ++ s, ++ s_num) {
      // get all repetitions of the structure in the volume
      // with position dr and orientation (tau, eta)
      // i.e. compute matrix that defines all num_grains grains inside the volume
      // distr = [ drx_i dry_i drz_i  tau_i eta_i ] (NDISTR x 5 )
      // get params of the shape repeated in this structure: [dn2, id, dims, invar, t, seta, stau]
      // allocate id(nqy, nqz, num_grains)

      // get the shape
      // compute t, lattice, ndoms, dd, nn, id etc.

      if(smaster) {
        std::cout << "-- Processing structure " << s_num + 1 << " ..." << std::endl;
      } // if

      Structure *curr_struct = &((*s).second);
      Lattice *curr_lattice = (Lattice*) HiGInput::instance().lattice(*curr_struct);
      Unitcell *curr_unitcell = (Unitcell*) HiGInput::instance().unitcell(*curr_struct);
      bool struct_in_layer = (*s).second.grain_in_layer();

      /* calulate propagation coefficients for current layer*/
      int order = HiGInput::instance().structure_layer_order(*curr_struct);
      layer_qgrid_qz(alphai, multilayer_[order].one_minus_n2());
      complex_vec_t fc; 
      if (!multilayer_.propagation_coeffs(fc, k0_, alphai, order)){
        //TODO call mpi abort
        std::exit(1);
      }

      real_t *dd = NULL, *nn = NULL, *wght = NULL;    // come back to this ...
                        // these structures can be improved ...
      real_t tz = 0;
      int num_dimen = 3;
      int ndx = 0, ndy = 0;
      // compute dd and nn
      spatial_distribution(s, tz, num_dimen, ndx, ndy, dd);
      orientation_distribution(s, dd, ndx, ndy, nn, wght);
      std::string struct_dist = (*s).second.grain_orientation();
      int num_grains = ndx;

      /*for(int i = 0; i < ndx; ++ i) std::cout << nn[i] << " ";
      std::cout << std::endl;
      for(int i = 0; i < ndx; ++ i) std::cout << nn[ndx + i] << " ";
      std::cout << std::endl;
      for(int i = 0; i < ndx; ++ i) std::cout << nn[2 * ndx + i] << " ";
      std::cout << std::endl;
      std::exit(1); */

      int r1axis, r2axis, r3axis;
      if ( struct_dist == "bragg" ){
        r1axis = 2;
        r2axis = 0;
        r3axis = 1;
        // TODO this is a hack
        if (dd) delete [] dd;
        dd = new real_t[ndx * 3];
        for (int i = 0; i < ndx * 3; i++) dd[i] = REAL_ZERO_;
      } else {
        r1axis = (int) (*s).second.rotation_rot1()[0];
        r2axis = (int) (*s).second.rotation_rot2()[0];
        r3axis = (int) (*s).second.rotation_rot3()[0];
      }

      if(smaster) {
        std::cout << "-- Grains: " << num_grains << std::endl;
      } // if

      //std::cout << "DD: " << std::endl;
      //for(unsigned int d = 0; d < num_grains; ++ d) {
      //  std::cerr << dd[d] << "\t" << dd[num_grains + d] << "\t" << dd[num_grains * 2 + d]
      //        << std::endl;
      //} // for*/

      /* grain scalings */
      std::vector<vector3_t> scaling_samples;
      std::vector<real_t> scaling_weights;
      bool is_scaling_distribution = false;
      if (s->second.grain_scaling_is_dist()) {
        is_scaling_distribution = true;
        std::vector<StatisticType> dist = s->second.grain_scaling_dist();
        vector3_t mean = s->second.grain_scaling();
        vector3_t stddev = s->second.grain_scaling_stddev();
        std::vector<int> scaling_nvals = s->second.grain_scaling_nvals();

        // if sigma is zeros set sampling count to 1
        for (int i = 0; i < 3; i++)
          if (stddev[i] == 0)
            scaling_nvals[i] = 1;
        construct_scaling_distribution(dist, mean, stddev, scaling_nvals, scaling_samples,
                     scaling_weights);
      } else {
        scaling_samples.push_back(s->second.grain_scaling());
        scaling_weights.push_back(1.);
      } // if-else

      /*
         // for testing
         if(master) {
         std::cerr << "************ grain scaling distribution **************" << std::endl;
         std::cerr << " LENGTH = " << scaling_samples.size() << std::endl;
         for(std::vector<vector3_t>::iterator i = scaling_samples.begin();
         i != scaling_samples.end(); ++i) {
         std::cerr << "scaling = [ " << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << " ]"
         << std::endl;
         } // for
         } // if
       */

      /* grain repetitions */
      bool is_grain_repetition_dist = false;
      std::vector<vector3_t> all_grains_repeats;
      if((*s).second.grain_is_repetition_dist()) {
        is_grain_repetition_dist = true;
        // get nvalues from scaling distribution
        int num_repeats;
        if (scaling_samples.size() > 1) 
                    num_repeats = scaling_samples.size();    // CHECK ...
        else
          num_repeats = num_grains;          // CHECK ...
        construct_repetition_distribution((*s).second.grain_repetitiondist(), 
                        num_repeats, all_grains_repeats);
      } else {
        vector3_t grain_repeats = (*s).second.grain_repetition();
        all_grains_repeats.push_back(grain_repeats);
      } // if-else

      // for testing
      //if(master) {
      //  std::cout << "************ grain repetition distribution **************" << std::endl;
      //  std::cout << " LENGTH = " << all_grains_repeats.size() << std::endl;
      //  for(std::vector<vector3_t>::iterator i = all_grains_repeats.begin();
      //      i != all_grains_repeats.end(); ++ i) {
      //    std::cout << "[ " << (*i)[0] << " " << (*i)[1] << " " << (*i)[2] << " ]"
      //        << std::endl;
      //  } // for
      //} // if

      int num_repeat_scaling;
      if (scaling_samples.size() == all_grains_repeats.size())
        num_repeat_scaling = scaling_samples.size();
      else if ((scaling_samples.size() == 1) && (all_grains_repeats.size() > 1))
        num_repeat_scaling = all_grains_repeats.size();
      else if ((scaling_samples.size() > 1) && (all_grains_repeats.size() == 1))
        num_repeat_scaling = scaling_samples.size();
      else {
        std::cerr << "error: scaling and repetition distributions are not aligned." << std::endl;
        return false;
      }

      vector3_t curr_transvec = (*s).second.grain_transvec();
      if(HiGInput::instance().param_nslices() <= 1) {
        curr_transvec[2] = curr_transvec[2] - single_layer_thickness_;
      } else {
        curr_transvec[2] = curr_transvec[2]; // TODO/FIXME... for more than 1 layers ... 
                           // ... incomplete
      } // if-else
      //ShapeName shape_name = HiGInput::instance().shape_name((*s).second);
      //std::string shape_file = HiGInput::instance().shape_filename((*s).second);
      //shape_param_list_t shape_params = HiGInput::instance().shape_params((*s).second);

      /* computing dwba ff for each grain in structure (*s) with
       * ensemble containing num_grains grains */

      int grain_min = 0;
      int num_gr = num_grains;
      int grain_max = num_gr;

      #ifdef USE_MPI
        // divide among processors
        int num_procs = multi_node_.size(struct_comm);
        int rank = multi_node_.rank(struct_comm);
        int grain_color = 0;
        if(num_procs > num_gr) {
          grain_color = rank % num_gr;
          grain_min = grain_color;
          num_gr = 1;
        } else {
          grain_color = rank;
          grain_min = ((num_gr / num_procs) * rank + min(rank, num_gr % num_procs));
          num_gr = (num_gr / num_procs) + (rank < num_gr % num_procs);
        } // if-else
        grain_max = grain_min + num_gr;
        woo::comm_t grain_comm = "grain";
        multi_node_.split(grain_comm, struct_comm, grain_color);

        bool gmaster = multi_node_.is_master(grain_comm);
        int temp_gmaster = gmaster;
        int *gmasters = new (std::nothrow) int[multi_node_.size(struct_comm)];
        // all grain masters tell the structure master about who they are
        multi_node_.gather(struct_comm, &temp_gmaster, 1, gmasters, 1);
      #else
        bool gmaster = true;
      #endif // USE_MPI

      complex_t *grain_ids = NULL;
      if(gmaster) {    // only the structure master needs this
        grain_ids = new (std::nothrow) complex_t[num_gr * nrow_ * ncol_];
        if(grain_ids == NULL) {
          std::cerr << "error: could not allocate memory for 'id'" << std::endl;
          return false;
        } // if
        // initialize to 0
        memset(grain_ids, 0 , num_gr * nrow_ * ncol_ * sizeof(complex_t));
      } // if

      // loop over grains - each process processes num_gr grains
      for(int grain_i = grain_min; grain_i < grain_max; grain_i ++) {  // or distributions

        if(gmaster) {
          std::cout << "-- Processing grain " << grain_i + 1 << " / " << num_grains << " ..."
                << std::endl;
        } // if

        // define r_norm (grain orientation by tau and eta)
        // define full grain rotation matrix r_total = r_phi * r_norm
        // TODO: ... i think these tau eta zeta can be computed on the fly to save memory ...
        real_t rot1 = nn[0 * num_grains + grain_i];
        real_t rot2 = nn[1 * num_grains + grain_i];
        real_t rot3 = nn[2 * num_grains + grain_i];
        real_t gauss_weight = 1.0;
        if(struct_dist == "gaussian") {
          gauss_weight = wght[grain_i] * wght[num_grains + grain_i] *
                  wght[2 * num_grains + grain_i];
        } // if
        RotMatrix_t r1(r1axis, rot1);
        RotMatrix_t r2(r2axis, rot2);
        RotMatrix_t r3(r3axis, rot3); 

        // order of multiplication is important
        RotMatrix_t rot = rot_ * r3 * r2 * r1;

        /* center of unit cell replica */
        vector3_t curr_dd_vec(dd[grain_i + 0],
                    dd[grain_i + num_grains],
                    dd[grain_i + 2 * num_grains]);
        vector3_t center = rot * curr_dd_vec + curr_transvec;

        /* compute structure factor and form factor */

        // TODO: parallel tasks ...

        woo::BoostChronoTimer sftimer, fftimer;
        fftimer.start(); fftimer.pause();

        // temporarily do this ...
        std::vector<complex_t> ff;
        unsigned int sz = nqz_;
        if(HiGInput::instance().experiment() == "gisaxs") sz = nqz_extended_;
        ff.clear();
        ff.resize(sz, CMPLX_ZERO_);
        /*std::cout << "** ARRAY CHECK ** (sz=" << sz << ")" << std::endl;
        if(!check_finite(QGrid::instance().qx(), QGrid::instance().nqx())) {
          std::cerr << "** ARRAY CHECK ** qx failed check (sz=" << QGrid::instance().nqx() << ")" << std::endl;
        } // if
        if(!check_finite(QGrid::instance().qy(), QGrid::instance().nqy())) {
          std::cerr << "** ARRAY CHECK ** qy failed check (sz=" << QGrid::instance().nqy() << ")" << std::endl;
        } // if
        if(!check_finite(QGrid::instance().qz_extended(), QGrid::instance().nqz_extended())) {
          std::cerr << "** ARRAY CHECK ** qz failed check (sz=" << QGrid::instance().nqz_extended() << ")" << std::endl;
        } // if*/

        // loop over all elements in the unit cell
        for(Unitcell::element_iterator_t e = curr_unitcell->element_begin();
            e != curr_unitcell->element_end(); ++ e) {

          ShapeName shape_name = HiGInput::instance().shape_name((*e).first);
          real_t zrot = HiGInput::instance().shape_zrot((*e).first);
          real_t yrot = HiGInput::instance().shape_yrot((*e).first);
          real_t xrot = HiGInput::instance().shape_xrot((*e).first);
          std::string shape_file = HiGInput::instance().shape_filename((*e).first);
          shape_param_list_t shape_params = HiGInput::instance().shape_params((*e).first);

          // update rotation matrix
          rot = rot *  RotMatrix_t(0, xrot) * RotMatrix_t(1, yrot) * RotMatrix_t(2, zrot);

          for(Unitcell::location_iterator_t l = (*e).second.begin(); l != (*e).second.end(); ++ l) {
            vector3_t transvec = (*l);
            #ifdef FF_NUM_GPU   // use GPU
              #ifdef FF_NUM_GPU_FUSED
                FormFactor eff(64, 8);
              #elif defined KERNEL2
                FormFactor eff(2, 4, 4);
              #else
                FormFactor eff(64);
              #endif
            #else   // use CPU or MIC
              FormFactor eff;
            #endif

            // TODO remove these later
            real_t shape_tau = 0., shape_eta = 0.;
            fftimer.resume();
            //read_form_factor("curr_ff.out");
            form_factor(eff, shape_name, shape_file, shape_params, transvec,
                  shape_tau, shape_eta, rot
                  #ifdef USE_MPI
                    , grain_comm
                  #endif
                  );
            fftimer.pause();
            /*if(!check_finite(eff.ff(), sz)) {
              std::cerr << "** ARRAY CHECK ** eff failed check" << std::endl;
            } // if*/

            // for each location, add the FFs
            for(unsigned int i = 0; i < sz; ++ i) ff[i] += eff[i]; 
          } // for l
        } // for e

        fftimer.stop();
#ifndef FF_VERBOSE
          std::cout << "**               FF compute time: "
                    << fftimer.elapsed_msec() << " ms." << std::endl;
#endif
        /*if(!check_finite(&ff[0], ff.size())) {
          std::cerr << "** ARRAY CHECK ** ff failed check" << std::endl;
        } // if*/


        sftimer.start(); sftimer.pause();
        for (int i_scale = 0; i_scale < num_repeat_scaling; i_scale++) {

          /* set the scalig for this grain */
          vector3_t grain_scaling;
          real_t scaling_wght = 1.0;
          if (is_scaling_distribution) {
            grain_scaling = scaling_samples[i_scale];
            scaling_wght = scaling_weights[i_scale];
          } else {
            grain_scaling = scaling_samples[0];
            scaling_wght = scaling_weights[0];
          } // if-else

          /* set the repetition for this grain */
          vector3_t grain_repeats;
          if (is_grain_repetition_dist)
            grain_repeats = all_grains_repeats[grain_i % all_grains_repeats.size()];
          else
            grain_repeats = all_grains_repeats[0]; // same repeat for all grains

#ifdef SF_VERBOSE
          std::cout << "-- Distribution sample " << i_scale + 1 << " / "
              << scaling_samples.size() << ".\n";
 #endif

          /* calulate structure factor for the grain */
          real_t weight = gauss_weight * scaling_wght;
          StructureFactor sf;
          sf.putStructureType(curr_struct->getStructureType());
          std::shared_ptr<Paracrystal> pc = curr_struct->paracrystal();
          std::shared_ptr<PercusYevick> py = curr_struct->percusyevick();
          sftimer.resume();
          if(!structure_factor(sf, HiGInput::instance().experiment(), center, curr_lattice,
                  grain_repeats, grain_scaling, rot, pc, py
                  #ifdef USE_MPI
                    , grain_comm
                  #endif
                  )){
              std::cerr << "Error: aborting run due to previous errors" << std::endl;
              std::exit(1);
          }
          sftimer.pause();

          /* compute intensities using sf and ff */
          // TODO: parallelize ...
          if(gmaster) {  // grain master
            complex_t* base_id = grain_ids + (grain_i - grain_min) * nrow_ * ncol_;
            unsigned int nslices = HiGInput::instance().param_nslices();
            unsigned int imsize = nrow_ * ncol_;
            if(nslices <= 1) {
              /* without slicing */
              //if(single_layer_refindex_.delta() < 0 || single_layer_refindex_.beta() < 0) {
              //  // this should never happen
              //  std::cerr << "error: single layer information not correctly set"
              //            << std::endl;
              //  return false;
              //} // if
              complex_t dn2 = multilayer_[order].one_minus_n2() - curr_struct->one_minus_n2();
              if(HiGInput::instance().experiment() == "gisaxs") {  // GISAXS
                for(unsigned int i = 0; i < imsize; ++ i) {
                  unsigned int curr_index   = i;
                  unsigned int curr_index_0 = i; 
                  unsigned int curr_index_1 = 1 * imsize + i;
                  unsigned int curr_index_2 = 2 * imsize + i;
                  unsigned int curr_index_3 = 3 * imsize + i;
                  base_id[curr_index] = dn2 * weight * 
                                       (fc[curr_index_0] * sf[curr_index_0] * ff[curr_index_0] +
                                        fc[curr_index_1] * sf[curr_index_1] * ff[curr_index_1] +
                                        fc[curr_index_2] * sf[curr_index_2] * ff[curr_index_2] +
                                        fc[curr_index_3] * sf[curr_index_3] * ff[curr_index_3]);
                } // for i
              } else { // SAXS
                for(unsigned int i = 0; i < imsize; ++ i) {
                  unsigned int curr_index   = i;
                  unsigned int curr_index_0 = i; 
                  base_id[curr_index] = dn2 * weight * 
                                       (fc[curr_index_0] * sf[curr_index_0] * ff[curr_index_0]);
                } // for i
              } // if-else
              if(HiGInput::instance().saveff()) {
                std::string ffoutput(HiGInput::instance().param_pathprefix() +
                                     "/" + HiGInput::instance().runname() +
                                     "/ff.out");
                //ff.save_ff(4 * imsize, ffoutput);
                std::ofstream fout(ffoutput, std::ios::out);
                for (int i = 0; i < nqz_extended_; i++)
                  fout << std::abs(ff[i]) << std::endl;
                fout.close();
              } // if
              if(HiGInput::instance().savesf()) {
                std::string sfoutput(HiGInput::instance().param_pathprefix() +
                                     "/" + HiGInput::instance().runname() +
                                     "/sf.out");
                sf.save_sf(sfoutput);
              } // if
            } else {
              /* perform slicing */
              // not yet implemented ...
              std::cout << "uh-oh: ever thought about implementing the slicing scheme?"
                        << std::endl;
              //return false;
            } // if-else
          } // if gmaster
          //sf.clear(); // clean structure factor
        } // for i_scale
        sftimer.stop();
        #ifndef SF_VERBOSE
          std::cout << "**               SF compute time: "
                    << sftimer.elapsed_msec() << " ms." << std::endl;
        #else
          std::cout << sftimer.elapsed_msec() << " ms." << std::endl;
        #endif

        // clean form factor before going to next
        ff.clear();

      } // for num_gr

      complex_t* id = NULL;
      #ifdef USE_MPI
        if(multi_node_.size(struct_comm) > 1) {
          // collect grain_ids from all procs in struct_comm
          if(smaster) {
            id = new (std::nothrow) complex_t[num_grains * nrow_ * ncol_];
          } // if
          int *proc_sizes = new (std::nothrow) int[multi_node_.size(struct_comm)];
          int *proc_displacements = new (std::nothrow) int[multi_node_.size(struct_comm)];
          multi_node_.gather(struct_comm, &num_gr, 1, proc_sizes, 1);
          if(smaster) {
            // make sure you get data only from the gmasters
            for(int i = 0; i < multi_node_.size(struct_comm); ++ i)
              proc_sizes[i] *= (gmasters[i] * nrow_ * ncol_);
            proc_displacements[0] = 0;
            for(int i = 1; i < multi_node_.size(struct_comm); ++ i)
              proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
          } // if
          multi_node_.gatherv(struct_comm, grain_ids, gmaster * num_gr * nrow_ * ncol_,
                    id, proc_sizes, proc_displacements);
          delete[] proc_displacements;
          delete[] proc_sizes;
        } else {
          id = grain_ids;
        } // if-else

        delete[] gmasters;
        multi_node_.free(grain_comm);
      #else
        id = grain_ids;
      #endif

      delete[] nn;
      delete[] dd;
      delete[] wght;

      if(smaster) {
        // FIXME: double check the following correlation stuff ...
        // new stuff for grain/ensemble correlation
        unsigned int soffset = 0;
        switch(HiGInput::instance().param_structcorrelation()) {
          case structcorr_null:  // default
          case structcorr_nGnE:  // no correlation
            // struct_intensity = sum_grain(abs(grain_intensity)^2)
            // intensity = sum_struct(struct_intensity)
            soffset = s_num * nrow_ * ncol_;
            for(unsigned int i = 0; i < nrow_ * ncol_ ; ++ i) {
              unsigned int curr_index = i;
              real_t sum = 0.0;
              for(int d = 0; d < num_grains; ++ d) {
                unsigned int id_index = d * nrow_ * ncol_ + curr_index;
                sum += id[id_index].real() * id[id_index].real() + 
                  id[id_index].imag() * id[id_index].imag();
              } // for d
              struct_intensity[soffset + curr_index] = sum;
            } // for i
            break;

          case structcorr_nGE:  // non corr grains, corr ensemble
            // grain_intensity = abs(sum_struct(struct_intensity))^2
            // intensity = sum_grain(grain_itensity)
            // TODO ...
            std::cerr << "uh-oh: this nGE correlation is not yet implemented" << std::endl;
            return false;
            break;

          case structcorr_GnE:  // corr grains, non corr ensemble
            // struct_intensity = abs(sum_grain(grain_intensity))^2
            // intensty = sum_struct(struct_intensity)
            soffset = s_num * nrow_ * ncol_;
            for(unsigned int i = 0; i < nrow_ * ncol_; ++ i) {
              unsigned int curr_index = i;
              complex_t sum(0.0, 0.0);
              for(int d = 0; d < num_grains; ++ d) {
                unsigned int id_index = d * nrow_ * ncol_ + curr_index;
                sum += id[id_index];
              } // for d
              struct_intensity[soffset + curr_index] = sum.real() * sum.real() +
                                    sum.imag() * sum.imag();
            } // for i
            break;

          case structcorr_GE:    // both correlated
            // struct_intensity = sum_grain(grain_intensity)
            // intensity = abs(sum_struct(struct_intensity))^2
            soffset = s_num * nrow_ * ncol_;
            for(unsigned int i = 0; i < nqz_; ++ i) {
              unsigned int curr_index = i;
              complex_t sum(0.0, 0.0);
              for(int d = 0; d < num_grains; ++ d) {
                unsigned int id_index = d * nrow_ * ncol_ + curr_index;
                sum += id[id_index];
              } // for d
              c_struct_intensity[soffset + curr_index] = sum;
            } // for i
            break;

          default:        // error
            std::cerr << "error: unknown correlation type." << std::endl;
            return false;
        } // switch

        delete[] id;
      } // if smaster

    } // for num_structs

    std::vector<real_t> iratios;
    real_t iratios_sum = 0.0;
    for(structure_iterator_t s = HiGInput::instance().structure_begin(); s != HiGInput::instance().structure_end(); ++ s) {
      iratios.push_back((*s).second.iratio());
      iratios_sum += (*s).second.iratio();
    } // for

    /*if(iratios_sum != 1.0) {
      std::cerr << "error: iratios of all structures must add to 1.0 ("
            << iratios_sum << ")" << std::endl;
      return false;
    } // if*/

    real_t* all_struct_intensity = NULL;
    complex_t* all_c_struct_intensity = NULL;
    #ifdef USE_MPI
      if(multi_node_.size(comm_key) > 1) {
        // collect struct_intensity from all procs in comm_key
        if(master) {
          if(HiGInput::instance().param_structcorrelation() == structcorr_GE) {
            all_c_struct_intensity =
                new (std::nothrow) complex_t[num_structures_ * nrow_ * ncol_];
          } else {
            all_struct_intensity =
                new (std::nothrow) real_t[num_structures_ * nrow_ * ncol_];
          } // if-else
        } // if
        int *proc_sizes = new (std::nothrow) int[multi_node_.size(comm_key)];
        int *proc_displacements = new (std::nothrow) int[multi_node_.size(comm_key)];
        multi_node_.gather(comm_key, &num_structs, 1, proc_sizes, 1);
        if(master) {
          // make sure you get data only from the smasters
          for(int i = 0; i < multi_node_.size(comm_key); ++ i)
            proc_sizes[i] *= (smasters[i] * nrow_ * ncol_);
          proc_displacements[0] = 0;
          for(int i = 1; i < multi_node_.size(comm_key); ++ i)
            proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
        } // if
        if(HiGInput::instance().param_structcorrelation() == structcorr_GE) {
          multi_node_.gatherv(comm_key, c_struct_intensity,
                    smaster * num_structs * nrow_ * ncol_,
                    all_c_struct_intensity, proc_sizes, proc_displacements);
        } else {
          multi_node_.gatherv(comm_key, struct_intensity,
                    smaster * num_structs * nrow_ * ncol_,
                    all_struct_intensity, proc_sizes, proc_displacements);
        } // if-else
        delete[] proc_displacements;
        delete[] proc_sizes;
      } else {
        all_struct_intensity = struct_intensity;
        all_c_struct_intensity = c_struct_intensity;
      } // if-else
    #else
      all_struct_intensity = struct_intensity;
      all_c_struct_intensity = c_struct_intensity;
    #endif

    #ifdef USE_MPI
      multi_node_.free(struct_comm);
      multi_node_.barrier(comm_key);
    #endif

    if(master) {
      // normalize
      if(all_struct_intensity != NULL)
        normalize(all_struct_intensity, nrow_ * ncol_);

      img3d = new (std::nothrow) real_t[nrow_ * ncol_];
      if (img3d == nullptr) {
        std::cerr << "error: unable to allocate memeory." << std::endl;
        std::exit(1);
      }
      // sum of all struct_intensity into intensity

      // new stuff for correlation
      switch(HiGInput::instance().param_structcorrelation()) {
        case structcorr_null:  // default
        case structcorr_nGnE:  // no correlation
          // struct_intensity = sum_grain(abs(grain_intensity)^2)
          // intensity = sum_struct(struct_intensity)
          for(unsigned int z = 0; z < nrow_ * ncol_; ++ z) {
            real_t sum = 0.0;
            for(int s = 0; s < num_structures_; ++ s) {
              unsigned int index = s * nrow_ * ncol_ + z;
              sum += all_struct_intensity[index] * iratios[s];
            } // for d
            img3d[z] = sum;
          } // for z
          break;

        case structcorr_nGE:  // non corr grains, corr ensemble
          // grain_intensity = abs(sum_struct(struct_intensity))^2
          // intensity = sum_grain(grain_itensity)
          // TODO ...
          std::cerr << "uh-oh: this nGE correlation is not yet implemented" << std::endl;
          return false;
          break;

        case structcorr_GnE:  // corr grains, non corr ensemble
          // struct_intensity = abs(sum_grain(grain_intensity))^2
          // intensty = sum_struct(struct_intensity)
          for(unsigned int z = 0; z < nrow_ * ncol_; ++ z) {
            real_t sum = 0.0;
            for(int s = 0; s < num_structures_; ++ s) {
              unsigned int index = s * nrow_ * ncol_ + z;
              sum += all_struct_intensity[index] * iratios[s];
            } // for d
            img3d[z] = sum;
          } // for z
          break;

        case structcorr_GE:    // both correlated
          // struct_intensity = sum_grain(grain_intensity)
          // intensity = abs(sum_struct(struct_intensity))^2
          for(unsigned int z = 0; z < nrow_ * ncol_; ++ z) {
            complex_t sum = 0.0;
            for(int s = 0; s < num_structures_; ++ s) {
              unsigned int index = s * nrow_ * ncol_ + z;
              sum += all_c_struct_intensity[index] * iratios[s];
            } // for d
            img3d[z] = sum.real() * sum.real() + sum.imag() * sum.imag();
          } // for z
          break;

        default:        // error
          std::cerr << "error: unknown correlation type." << std::endl;
          return false;
      } // switch

      delete[] all_c_struct_intensity;
      delete[] all_struct_intensity;
    } // if master

    if(master) {
      // convolute/smear the computed intensities
      real_t sigma = HiGInput::instance().scattering_smearing();
      if(sigma > 0.0) {
        woo::BoostChronoTimer smear_timer;
        std::cout << "-- Smearing the result with sigma = " << sigma << " ... " << std::flush;
        if(img3d == NULL) {
          std::cerr << "error: there is no img3d. you are so in the dumps!" << std::endl;
        } // if
        smear_timer.start();
        gaussian_smearing(img3d, sigma);
        smear_timer.stop();
        std::cout << "done." << std::endl;
        std::cout << "**                 Smearing time: "
              << smear_timer.elapsed_msec() << " ms." << std::endl;
      } // if
    } // if master

    return true;
  } // HipGISAXS::run_gisaxs()


  bool HipGISAXS::normalize(real_t*& data, unsigned int size) {
    real_t min_val = data[0], max_val = data[0];
    for(unsigned int i = 1; i < size; ++ i) {
      min_val = std::min(min_val, data[i]);
      max_val = std::max(max_val, data[i]);
    } // for
    real_t d = max_val - min_val;
    if(d < 1e-30) return false;
    #pragma omp parallel for
    for(unsigned int i = 1; i < size; ++ i) data[i] = (data[i] - min_val) / d;
    return true;
  } // HipGISAXS::normalize()


  bool HipGISAXS::structure_factor(StructureFactor& sf,
                  std::string expt, vector3_t& center, Lattice* &curr_lattice,
                  vector3_t& grain_repeats, vector3_t& grain_scaling,
                  RotMatrix_t & rot, 
                  std::shared_ptr<Paracrystal> pc, std::shared_ptr<PercusYevick> py
                  #ifdef USE_MPI
                    , woo::comm_t comm_key
                  #endif
                  ) {
    #ifndef SF_GPU
      return sf.compute_structure_factor(expt, center, curr_lattice, grain_repeats,
                      grain_scaling, rot, pc, py
                      #ifdef USE_MPI
                        , multi_node_, comm_key 
                      #endif
                      );
    #else
      if (pc == nullptr && py == nullptr)
        return sf.compute_structure_factor_gpu(expt, center, curr_lattice, grain_repeats,
                      grain_scaling, rot
                      #ifdef USE_MPI
                        , multi_node_, comm_key
                      #endif
                      );
      else
        return sf.compute_structure_factor(expt, center, curr_lattice, grain_repeats,
                      grain_scaling, rot, pc, py
                      #ifdef USE_MPI
                        , multi_node_, comm_key 
                      #endif
                      );
    #endif
  } // HipGISAXS::structure_factor()


  bool HipGISAXS::form_factor(FormFactor& ff,
                ShapeName shape_name, std::string shape_file,
                shape_param_list_t& shape_params, vector3_t &curr_transvec,
                real_t shp_tau, real_t shp_eta,
                RotMatrix_t & rot
                #ifdef USE_MPI
                  , woo::comm_t comm_key
                #endif
                ) {
    return ff.compute_form_factor(shape_name, shape_file, shape_params,
                      single_layer_thickness_,
                      curr_transvec, shp_tau, shp_eta, rot
                      #ifdef USE_MPI
                        , multi_node_, comm_key
                      #endif
                      );
  } // HipGISAXS::form_factor()


  bool HipGISAXS::compute_rotation_matrix_z(real_t angle,
                        vector3_t& r1, vector3_t& r2, vector3_t& r3) {
    real_t s = sin(angle);
    real_t c = cos(angle);

    r1[0] = c; r1[1] = -s; r1[2] = 0.0;
    r2[0] = s; r2[1] = c; r2[2]  = 0.0;
    r3[0] = 0.0; r3[1] = 0.0; r3[2]  = 1.0;

    return true;
  } // HipGISAXS::comptue_rotation_matrix_z()


  bool HipGISAXS::compute_rotation_matrix_y(real_t angle,
                        vector3_t& r1, vector3_t& r2, vector3_t& r3) {
    real_t s = sin(angle);
    real_t c = cos(angle);

    //r1[0] = c; r1[1] = 0.0; r1[2] = -s;
    r1[0] = c; r1[1] = 0.0; r1[2] = s;
    r2[0] = 0.0; r2[1] = 1.0; r2[2]  = 0.0;
    //r3[0] = s; r3[1] = 0.0; r3[2]  = c;
    r3[0] = -s; r3[1] = 0.0; r3[2]  = c;

    return true;
  } // HipGISAXS::compute_rotation_matrix_y()


  bool HipGISAXS::compute_rotation_matrix_x(real_t angle,
                        vector3_t& r1, vector3_t& r2, vector3_t& r3) {
    real_t s = sin(angle);
    real_t c = cos(angle);

    r1[0] = 1.0; r1[1] = 0.0; r1[2] = 0.0;
    r2[0] = 0.0; r2[1] = c; r2[2]  = -s;
    r3[0] = 0.0; r3[1] = s; r3[2]  = c;

    return true;
  } // HipGISAXS::compute_rotation_matrix_x()


  bool HipGISAXS::illuminated_volume(real_t alpha_i, real_t spot_area, int min_layer_order,
                    RefractiveIndex substrate_refindex) {
    real_t spot_diameter = 2.0 * sqrt(spot_area / PI_) * 1e6;  // in nm
    real_t substr_delta = substrate_refindex.delta();
    real_t substr_beta = substrate_refindex.beta();

    if(HiGInput::instance().experiment() == "saxs") {  // SAXS
      vol_[0] = vol_[1] = vol_[2] = spot_diameter;
    } else if(HiGInput::instance().experiment() == "gisaxs") {    // GISAXS
      complex_t c_max_depth = complex_t(MAX_DEPTH_, 0);
      complex_t penetration_depth_layer;
      if(HiGInput::instance().is_single_layer()) {
        RefractiveIndex r1 = HiGInput::instance().single_layer().refindex();
        real_t alpha_c = sqrt(2.0 * r1.delta());
        penetration_depth_layer = -1.0 /
                      (2.0 * k0_ * (sqrt(alpha_i * alpha_i - alpha_c * alpha_c -
                      complex_t(0, 2.0 * r1.beta()))).imag());
      } else if(HiGInput::instance().num_layers() == 0) {
        real_t alpha_c = sqrt(2.0 * substr_delta);
        penetration_depth_layer = -1.0 /
                      (2.0 * k0_ * (sqrt(alpha_i * alpha_i - alpha_c * alpha_c -
                      complex_t(0, 2.0 * substr_beta))).imag());
      } else {
        // the sample is described by 2 or more layers, slicing scheme will be applied
        // NOTE: the case where a structure is implicitly in 2 different layers is
        // TODO: not currently handled ...
        std::cerr << "uh-oh: this case (num_layers > 1) has not yet been implemented yet"
              << std::endl;
        std::cerr << "go get yourself a nice cup of yummmy hot chocolate instead!" << std::endl;
        return false;
      } // if-else
      vol_[0] = vol_[1] = spot_diameter;
      vol_[2] = (penetration_depth_layer.real() < c_max_depth.real()) ?
                penetration_depth_layer.real() : c_max_depth.real();
    } else {
      std::cerr << "error: experiment type '" << HiGInput::instance().experiment()
            << "' is either unknown or has not been implemented." << std::endl;
      return false;
    } // if-else

    return true;
  } // HipGISAXS::illuminated_volume()


  bool HipGISAXS::compute_propagation_coefficients(real_t alpha_i,
          complex_t* &amm, complex_t* &apm, complex_t* &amp, complex_t* &app,
          complex_t* &rk1, complex_t* &rk2, complex_t* &rk1rk2, complex_t* &tk1tk2,
          complex_t* &h0, complex_t* &fc) {
    amm = apm = amp = app = NULL;
    rk1 = rk2 = rk1rk2 = tk1tk2 = NULL;
    h0 = NULL;

    // doesnt this also depend on the type of polarization of the light? ... ?

    if(HiGInput::instance().param_nslices() <= 1) {  /* computing without sample slicing */
      complex_t dnl2 = 2.0 * complex_t(single_layer_refindex_.delta(),
                      single_layer_refindex_.beta());
      size_t imsize = nrow_ * ncol_;

      if(!layer_qgrid_qz(alpha_i, dnl2)) {
        std::cerr << "error: could not compute extended qz" << std::endl;
        return false;
      } // if

      if(HiGInput::instance().struct_in_layer() && HiGInput::instance().num_layers() == 1) {
        /* compute fresnel coefficients for a structure
         * embedded inside a layer on a substrate */
        if(!compute_fresnel_coefficients_embedded(alpha_i, fc)) {
          std::cerr << "error: could not compute fresnel coefficients" << std::endl;
          return false;
        } // if
        // set aliases

        amm = fc;
        apm = fc + imsize;
        amp = fc + 2 * imsize;
        app = fc + 3 * imsize;
        h0 = fc + 4 * imsize;
      } else if(HiGInput::instance().num_layers() == 1 || HiGInput::instance().num_layers() == 0) {
        /* compute fresnel coefficients for a structure
         * on top of or buried inside a substrate */
        if(!compute_fresnel_coefficients_top_buried(alpha_i, fc)) {
          std::cerr << "error: could not compute fresnel coefficients" << std::endl;
          return false;
        } // if
        // set aliases
        h0 = fc;              // re-check these ...
        rk2 = fc + imsize;
        rk1 = fc + 2 * imsize;
        rk1rk2 = fc + 3 * imsize;
      } else {
        std::cerr << "error: invalid number of layers in non-slicing scheme" << std::endl;
        return false;
      } // if-else
    } else {
      /* computing using the sample slicing scheme */
      // TODO: not yet implemented ...
      std::cerr << "uh-oh: you hit a slicing part that has not been implemented yet" << std::endl;
      return false;
    } // if

    return true;
  } // HipGISAXS::compute_propagation_coefficients()


  bool HipGISAXS::layer_qgrid_qz(real_t alpha_i, complex_t dnl_j) {

    if(!QGrid::instance().create_qz_extended(k0_, alpha_i, dnl_j)){
      std::cerr << "error: something went wrong while creating qz_extended" << std::endl;
      return false;
    } // if
    nqz_extended_ = QGrid::instance().nqz_extended();
    return true;
  } // HipGISAXS::layer_qgrid_qz()


  // TODO optimize this later ...
  bool HipGISAXS::compute_fresnel_coefficients_embedded(real_t alpha_i, complex_t* &fc) {
    RefractiveIndex nl = single_layer_refindex_;
    RefractiveIndex ns = substrate_refindex_;
    real_t lt = single_layer_thickness_;
    size_t imsize = nrow_ * ncol_;

    complex_t dnl2(2.0 * nl.delta(), 2.0 * nl.beta());
    complex_t dns2(2.0 * ns.delta(), 2.0 * ns.beta());

    real_t sinai = sin(alpha_i);
    real_t kiz0 = -1.0 * k0_ * sinai;
    complex_t kiz1 = -1.0 * k0_ * sqrt(sinai * sinai - dnl2);
    complex_t kiz2 = -1.0 * k0_ * sqrt(sinai * sinai - dns2);

    complex_t r01_kiz1 = (kiz0 - kiz1) / (kiz0 + kiz1);
    complex_t r12_kiz1 = (kiz1 - kiz2) / (kiz1 + kiz2);
    complex_t t01_kiz1 = 2.0 * (kiz0 / (kiz0 + kiz1));

    complex_t a1m_kiz1 = t01_kiz1 /
              ((real_t) 1.0 + r01_kiz1 * r12_kiz1 * exp(complex_t(0, 2) * kiz1 * lt));
    complex_t a1p_kiz1 = a1m_kiz1 * r12_kiz1 * exp(complex_t(0, 2) * kiz1 * lt);

    complex_t *a1mi = NULL, *a1pi = NULL;
    a1mi = new (std::nothrow) complex_t[imsize];
    a1pi = new (std::nothrow) complex_t[imsize];
    if(a1mi == NULL || a1pi == NULL) {
      std::cerr << "error: failed to allocate memory for a1mi, a1pi" << std::endl;
      return false;
    } // if
    for(unsigned int i = 0; i < imsize; ++ i) {
      a1mi[i] = a1m_kiz1; a1pi[i] = a1p_kiz1;
    } // for

    complex_t *a1mf = NULL, *a1pf = NULL;
    a1mf = new (std::nothrow) complex_t[imsize];
    a1pf = new (std::nothrow) complex_t[imsize];
    if(a1mf == NULL || a1pf == NULL) {
      std::cerr << "error: failed to allocate memory for a1mf, a1pf" << std::endl;
      return false;
    } // if
    for(unsigned int i = 0; i < imsize; ++ i) {
      a1mf[i] = a1pf[i] = complex_t(0.0, 0.0);
    } // for

    // allocate fc memory
    fc = new (std::nothrow) complex_t[5 * imsize];  // 5 sets
    if(fc == NULL) {
      std::cerr << "error: failed to allocate memory for fc" << std::endl;
      return false;
    } // if


    for(unsigned int z = 0; z < nqz_; ++ z) {
      complex_t a1m_nkfz1, a1p_nkfz1;
      real_t kfz0 = QGrid::instance().qz(z) + kiz0;
      unsigned int idx = 4 * imsize + z;

      if(kfz0 < 0) {
        a1m_nkfz1 = complex_t(0.0, 0.0);
        a1p_nkfz1 = complex_t(0.0, 0.0);
        fc[idx] = complex_t(0.0, 0.0);
      } else {  // kfz0 >= 0
        real_t nkfz0 = -1.0 * kfz0;
        complex_t nkfz1 = -1.0 * sqrt(kfz0 * kfz0 - k0_ * k0_ * dnl2);
        complex_t nkfz2 = -1.0 * sqrt(kfz0 * kfz0 - k0_ * k0_ * dns2);

        complex_t r01_nkfz1 = (nkfz0 - nkfz1) / (nkfz0 + nkfz1);
        complex_t r12_nkfz1 = (nkfz1 - nkfz2) / (nkfz1 + nkfz2);
        complex_t t01_nkfz1 = 2.0 * (nkfz0 / (nkfz0 + nkfz1));

        complex_t uniti = complex_t(0, 1);
        complex_t temp4 = 2.0 * nkfz1 * lt;
        real_t temp3 = exp(-1.0 * temp4.imag());
        complex_t temp5 = exp(uniti * temp4.real());
                 
        complex_t temp6 = temp3 * temp5;
        a1m_nkfz1 = t01_nkfz1 /
              ((real_t) 1.0 + r01_nkfz1 * r12_nkfz1 * temp6);
        a1p_nkfz1 = a1m_nkfz1 * r12_nkfz1 * temp6;

        fc[idx] = complex_t(1.0, 0.0);
      } // if-else

      a1mf[z] = a1m_nkfz1;
      a1pf[z] = a1p_nkfz1;
    } // for z

    // the element-element products
    for(unsigned int z = 0; z < nqz_; z ++) {
      fc[0 * imsize + z ] = a1mi[z] * a1mf[z];
      fc[1 * imsize + z ] = a1pi[z] * a1mf[z];
      fc[2 * imsize + z ] = a1mi[z] * a1pf[z];
      fc[3 * imsize + z ] = a1pi[z] * a1pf[z];
    } // for z

    delete[] a1pf;
    delete[] a1mf;
    delete[] a1pi;
    delete[] a1mi;

    return true;
  } // HipGISAXS::compute_fresnel_coefficients_embedded()


  // optimize ...
  bool HipGISAXS::compute_fresnel_coefficients_top_buried(real_t alpha_i, complex_t* &fc) {

    RefractiveIndex ns = substrate_refindex_;
    complex_t ns2 = 2 * complex_t(ns.delta(), ns.beta());

    // incoming wave
    real_t sin_ai = std::sin(alpha_i);
    real_t kzi  = -1 * k0_ * sin_ai;
    complex_t tmp = std::sqrt(sin_ai * sin_ai - ns2);
    complex_t rki = (sin_ai - tmp)/(sin_ai + tmp);

    size_t imsize = nrow_ * ncol_;
    fc = new (std::nothrow) complex_t[5 * imsize];
    if(fc == NULL) {
      std::cerr << "error: failed to allocate memory for fc" << std::endl;
      return false;
    } // if

    for(unsigned int z = 0; z < nqz_; ++ z) {
      real_t kzf = QGrid::instance().qz(z) + kzi;
      if(kzf < 0) {
        fc[z] = fc[imsize + z] = fc[2*imsize + z] = fc[3*imsize + z] = CMPLX_ZERO_;
      } else {
        real_t sin_af = kzf / k0_;
        tmp = std::sqrt(sin_af * sin_af - ns2);
        complex_t rkf = (sin_af - tmp)/(sin_af + tmp);
        fc[z] = complex_t(1.0, 0);
        fc[imsize + z] = rkf;
        fc[2 * imsize + z] = rki;
        fc[3 * imsize + z] = rki * rkf;
      }
    } // for z

    return true;
  } // HipGISAXS::compute_fresnel_coefficients_top_buried()


  bool HipGISAXS::spatial_distribution(structure_iterator_t s, real_t tz, int dim,
                    int& rand_dim_x, int& rand_dim_y, real_t* &d) {
    vector3_t spacing = (*s).second.ensemble_spacing();
    vector3_t maxgrains = (*s).second.ensemble_maxgrains();
    std::string distribution = (*s).second.ensemble_distribution();
    // vol_, cell_

    vector3_t spaced_cell = cell_ + spacing;

    if(distribution == "random") {
      srand(time(NULL));
      if(dim == 3) {
        // find max density - number of grains in vol
        //vector3_t max_density = min(floor(vol_ / spaced_cell) + 1, maxgrains);
        //rand_dim_x = max(1,
        //      (int)std::floor(max_density[0] * max_density[1] * max_density[2] / 4.0));
        // TEMPORARY HACK ... FIXME ...
        rand_dim_x = maxgrains[0] * maxgrains[1] * maxgrains[2];
        rand_dim_y = 3;

        // construct random matrix
        real_t *d_rand = new (std::nothrow) real_t[rand_dim_x * rand_dim_y];
        srand(time(NULL));
        for(int i = 0; i < rand_dim_x * rand_dim_y; ++ i)
          d_rand[i] = ((real_t)rand() / RAND_MAX);

        d = new (std::nothrow) real_t[rand_dim_x * rand_dim_y * 4];

        //int base_index = 0;
        real_t mul_val1 = 1.; //vol_[0] / 2;
        real_t mul_val2 = 1.; //vol_[1] / 2;
        real_t mul_val3 = 1.; //vol_[2];

        // D0
        if(maxgrains[0] == 1) {
          for(int x = 0; x < rand_dim_x; ++ x) d[x] = 0;
        } else {
          for(int x = 0; x < rand_dim_x; ++ x)
            d[x] = d_rand[x] * mul_val1;
        } // if-else
        if(maxgrains[1] == 1) {
          for(int x = 0; x < rand_dim_x; ++ x) d[4 * rand_dim_x + x] = 0;
        } else {
          for(int x = 0; x < rand_dim_x; ++ x)
            d[4 * rand_dim_x + x] = d_rand[rand_dim_x + x] * mul_val2;
        } // if-else
        if(maxgrains[2] == 1) {
          for(int x = 0; x < rand_dim_x; ++ x) d[8 * rand_dim_x + x] = tz;
        } else {
          for(int x = 0; x < rand_dim_x; ++ x)
            d[8 * rand_dim_x + x] = d_rand[2 * rand_dim_x + x] * mul_val3 + tz;
        } // if-else
        // D1
        for(int x = 0; x < rand_dim_x; ++ x)
          d[rand_dim_x + x] = - d[x];
        for(int x = 0; x < rand_dim_x; ++ x)
          d[4 * rand_dim_x + rand_dim_x + x] = d[4 * rand_dim_x + x];
        for(int x = 0; x < rand_dim_x; ++ x)
          d[8 * rand_dim_x + rand_dim_x + x] = d[8 * rand_dim_x + x];
        // D2
        for(int x = 0; x < rand_dim_x; ++ x)
          d[2 * rand_dim_x + x] = d[x];
        for(int x = 0; x < rand_dim_x; ++ x)
          d[4 * rand_dim_x + 2 * rand_dim_x + x] = - d[4 * rand_dim_x + x];
        for(int x = 0; x < rand_dim_x; ++ x)
          d[8 * rand_dim_x + 2 * rand_dim_x + x] = d[8 * rand_dim_x + x];
        // D3
        for(int x = 0; x < rand_dim_x; ++ x)
          d[3 * rand_dim_x + x] = - d[x];
        for(int x = 0; x < rand_dim_x; ++ x)
          d[4 * rand_dim_x + 3 * rand_dim_x + x] = - d[4 * rand_dim_x + x];
        for(int x = 0; x < rand_dim_x; ++ x)
          d[8 * rand_dim_x + 3 * rand_dim_x + x] = d[8 * rand_dim_x + x];

        // fix the rand_dim_x
        // rand_dim_x *= 4; ///////////////////////////////////////////////

      } else if(dim == 2) {
        std::cerr << "error: dim == 2 case not implemented" << std::endl;
        // ...
        return false;
      } else if(dim == 1) {
        std::cerr << "error: dim == 1 case not implemented" << std::endl;
        // ...
        return false;
      } else {
        std::cerr << "error: invalid dim size " << dim << std::endl;
        return false;
      } // if-else

    } else if(distribution == "regular") {
      if(dim == 3) {
        //vector3_t nd = min(floor(vol_ / spaced_cell) + 1, maxgrains);
        // TEMPORARY HACK .......... FIXME ...
        vector3_t nd = maxgrains;

        int size = nd[0] * nd[1] * nd[2];
        d = new (std::nothrow) real_t[3 * size];
        real_t* d1 = d;
        real_t* d2 = d + size;
        real_t* d3 = d + 2 * size;

        rand_dim_x = size;    // total number of rows in d
        rand_dim_y = 3;      // 3 dims, number of cols in d

        // this is more correct I think:
        real_t* dx = new (std::nothrow) real_t[(int)nd[0]];
        real_t* dy = new (std::nothrow) real_t[(int)nd[1]];
        real_t* dz = new (std::nothrow) real_t[(int)nd[2]];
        real_t val_dx = 0, step_dx = spaced_cell[0], max_dx = spaced_cell[0] * (nd[0] - 1);
        for(int i = 0; i < (int)nd[0]; ++ i, val_dx += step_dx) {
          dx[i] = val_dx - max_dx / 2;
        } // for
        real_t val_dy = 0, step_dy = spaced_cell[1], max_dy = spaced_cell[1] * (nd[1] - 1);
        for(int i = 0; i < (int)nd[1]; ++ i, val_dy += step_dy) {
          dy[i] = val_dy - max_dy / 2;
        } // for
        real_t val_dz = 0, step_dz = spaced_cell[2];
        for(int i = 0; i < (int)nd[2]; ++ i, val_dz += step_dz) {
          dz[i] = val_dz + tz;
        } // for

        for(int x = 0; x < (int)nd[0]; ++ x) {
          for(int yz = 0; yz < (int)nd[1] * (int)nd[2]; ++ yz) {
            d1[(int)nd[1] * (int)nd[2] * x + yz] = dx[x];
          } // for
        } // for
        for(int x = 0; x < (int)nd[0]; ++ x) {
          for(int y = 0; y < (int)nd[1]; ++ y) {
            for(int z = 0; z < (int)nd[2]; ++ z) {
              d2[(int)nd[1] * (int)nd[2] * x + (int)nd[2] * y + z] = dy[y];
            } // for
          } // for
        } // for
        for(int xy = 0; xy < (int)nd[0] * (int)nd[1]; ++ xy) {
          for(int z = 0; z < nd[2]; ++ z) {
            d3[(int)nd[2] * xy + z] = dz[z];
          } // for z
        } // for xy
        delete[] dz;
        delete[] dy;
        delete[] dx;
      } else if(dim == 2) {
        std::cerr << "error: dim == 2 case not implemented" << std::endl;
        // ...
        return false;
      } else if(dim == 1) {
        std::cerr << "error: dim == 1 case not implemented" << std::endl;
        // ...
        return false;
      } else {
        std::cerr << "error: invalid dim size " << dim << std::endl;
        return false;
      } // if-else
    } else if(distribution == "gaussian" || distribution == "normal") {
      // gaussian distribution
/*      real_t mean = (*s).second.ensemble_distribution_p1();
      real_t sd = (*s).second.ensemble_distribution_p2();
      woo::MTNormalRandomNumberGenerator rgen(mean, sd);
      if(dim == 3) {
        vector3_t nd = maxgrains;
        unsigned int size = nd[0] * nd[1] * nd[2];
        d = new (std::nothrow) real_t[dim * size];
        real_t* d1 = d;
        real_t* d2 = d + size;
        real_t* d3 = d2 + size;

        //real_t max_x = vol_[0];
        //real_t max_y = vol_[1];
        //real_t max_z = vol_[2];
        real_t max_x = spaced_cell[0] * (nd[0] - 1);
        real_t max_y = spaced_cell[1] * (nd[1] - 1);
        real_t max_z = spaced_cell[2] * (nd[2] - 1);

        // compute x for each grain
        for(unsigned int i = 0; i < size; ++ i);

      } else if(dim == 2) {
        std::cerr << "error: dim == 2 case not implemented" << std::endl;
        // ...
        return false;
      } else if(dim == 1) {
        std::cerr << "error: dim == 1 case not implemented" << std::endl;
        // ...
        return false;
      } else {
        std::cerr << "error: invalid dim size " << dim << std::endl;
        return false;
      } // if-else
*/    } else {
      // read .spa file ...
      std::cerr << "uh-oh: seems like you wanted to read distribution from a file" << std::endl;
      std::cerr << "sorry dear, this has not been implemented yet" << std::endl;
      return false;
    } // if-else

    return true;
  } // HipGISAXS::spatial_distribution()


  bool HipGISAXS::orientation_distribution(structure_iterator_t s, real_t* dd, int & ndx, int ndy, 
                        real_t* &nn, real_t* &wght) {
    std::string distribution = (*s).second.grain_orientation();
    //vector2_t tau = (*s).second.rotation_tau();
    //vector2_t eta = (*s).second.rotation_eta();
    //vector2_t zeta = (*s).second.rotation_zeta();
    vector3_t rot1 = (*s).second.rotation_rot1();
    vector3_t rot2 = (*s).second.rotation_rot2();
    vector3_t rot3 = (*s).second.rotation_rot3();

    if (distribution == "bragg"){
      real_vec_t angles;
      Lattice * lattice = (Lattice *) HiGInput::instance().lattice(s->second);
      vector3_t gr_scaling = s->second.grain_scaling();
      vector3_t gr_repetitions = s->second.grain_repetition();
      lattice->bragg_angles(gr_repetitions, gr_scaling, k0_, angles);
      if (angles.size() > 0 ) {
        ndx = angles.size();
        nn = new (std::nothrow) real_t[ndx * 3];
        for (int i = 0; i < ndx * 3; i++) nn[i] = REAL_ZERO_;
        for (int i = 0; i < ndx; i++) nn[i] = angles[i];
      } else {
        std::cerr << "Failed to calculate orientations for Bragg condition" << std::endl;
        std::cerr << "This is a bug, please report it." << std::endl;
        return false;
      }
      return true;
    }
    nn = new (std::nothrow) real_t[ndx * ndy];
    if (nn == NULL) {
      std::cerr << "error: could not allocate memory" << std::endl;
      return false;
    }
    wght = new (std::nothrow) real_t[ndx * ndy];
    if (wght == NULL) {
      std::cerr << "error: could not allocate memory" << std::endl;
      return false;
    }
    // TODO i believe constructing nn may not be needed ...
    if(distribution == "single") {    // single
      for(int x = 0; x < ndx; ++ x) {
        //nn[x] = tau[0] * PI_ / 180;
        nn[x] = rot1[1] * PI_ / 180;
      } // for x
      for(int x = 0; x < ndx; ++ x) {
        //nn[ndx + x] = eta[0] * PI_ / 180;
        nn[ndx + x] = rot2[1] * PI_ / 180;
      } // for x
      for(int x = 0; x < ndx; ++ x) {
        //nn[2 * ndx + x] = zeta[0] * PI_ / 180;
        nn[2 * ndx + x] = rot3[1] * PI_ / 180;
      } // for x
    } else if(distribution == "random") {  // random
      for(int x = 0; x < 3 * ndx; ++ x) {
        nn[x] = (real_t(rand()) / RAND_MAX) * 2 * PI_;
      } // for x
    } else if(distribution == "range") {  // range
      //real_t dtau = fabs(tau[1] - tau[0]);
      //real_t deta = fabs(eta[1] - eta[0]);
      //real_t dzeta = fabs(zeta[1] - zeta[0]);
      real_t drot1 = fabs(rot1[2] - rot1[1]);
      real_t drot2 = fabs(rot2[2] - rot2[1]);
      real_t drot3 = fabs(rot3[2] - rot3[1]);
      for(int x = 0; x < ndx; ++ x) {
        //nn[x] = (tau[0] + (real_t(rand()) / RAND_MAX) * dtau) * PI_ / 180;
        nn[x] = (rot1[1] + (real_t(rand()) / RAND_MAX) * drot1) * PI_ / 180;
      } // for x
      for(int x = 0; x < ndx; ++ x) {
        //nn[ndx + x] = (eta[0] + (real_t(rand()) / RAND_MAX) * deta) * PI_ / 180;
        nn[ndx + x] = (rot2[1] + (real_t(rand()) / RAND_MAX) * drot2) * PI_ / 180;
      } // for x
      for(int x = 0; x < ndx; ++ x) {
        //nn[2 * ndx + x] = (zeta[0] + (real_t(rand()) / RAND_MAX) * dzeta) * PI_ / 180;
        nn[2 * ndx + x] = (rot3[1] + (real_t(rand()) / RAND_MAX) * drot3) * PI_ / 180;
      } // for x
      /*real_t da = PI_ / (2 * (ndx - 1));
      for(int x = 0; x < ndx; ++ x) {
        nn[x] = x * da;
      } // for x
      for(int x = 0; x < 2 * ndx; ++ x) nn[ndx + x] = 0;*/
    } else if(distribution == "gaussian" || distribution == "normal") {  // gaussian
      real_t mean1 = (*s).second.rotation_rot1_anglemean();
      real_t sd1 = (*s).second.rotation_rot1_anglesd();
      real_t mean2 = (*s).second.rotation_rot2_anglemean();
      real_t sd2 = (*s).second.rotation_rot2_anglesd();
      real_t mean3 = (*s).second.rotation_rot3_anglemean();
      real_t sd3 = (*s).second.rotation_rot3_anglesd();
      /*woo::MTNormalRandomNumberGenerator rgen1(mean1, sd1);
      woo::MTNormalRandomNumberGenerator rgen2(mean1, sd2);
      woo::MTNormalRandomNumberGenerator rgen3(mean1, sd3);
      real_t drot1 = fabs(rot1[2] - rot1[1]);
      real_t drot2 = fabs(rot2[2] - rot2[1]);
      real_t drot3 = fabs(rot3[2] - rot3[1]);
      if(drot1 > 1e-20) {
        for(int x = 0; x < ndx; ++ x) {
          real_t temp;
          while(1) {
            temp = rgen1.rand();
            if(temp > rot1[1] && temp < rot1[2]) break;
          } // while
          nn[x] = temp * PI_ / 180;
        } // for x
      } else {
        for(int x = 0; x < ndx; ++ x) nn[x] = rot1[1] * PI_ / 180;
      } // if
      if(drot2 > 1e-20) {
        for(int x = 0; x < ndx; ++ x) {
          real_t temp;
          while(1) {
            temp = rgen2.rand();
            if(temp > rot2[1] && temp < rot2[2]) break;
          } // while
          nn[ndx + x] = temp * PI_ / 180;
        } // for x
      } else {
        for(int x = 0; x < ndx; ++ x) nn[ndx + x] = rot2[1] * PI_ / 180;
      } // if
      if(drot3 > 1e-20) {
        for(int x = 0; x < ndx; ++ x) {
          real_t temp;
          while(1) {
            temp = rgen3.rand();
            if(temp > rot3[1] && temp < rot3[2]) break;
          } // while
          nn[2 * ndx + x] = temp * PI_ / 180;
        } // for x
      } else {
        for(int x = 0; x < ndx; ++ x) nn[2 * ndx + x] = rot3[1] * PI_ / 180;
      } // if*/

      // for the gaussian weighting with regular instead of gaussian sampling
      real_t drot1 = rot1[2] - rot1[1];
      real_t drot2 = rot2[2] - rot2[1];
      real_t drot3 = rot3[2] - rot3[1];
      if(fabs(drot1) > 1e-20) {
        drot1 /= ndx;
        for(int x = 0; x < ndx; ++ x) {
          real_t temp = rot1[1] + x * drot1;
          //real_t temp = rot1[1] + (real_t(rand()) / RAND_MAX) * drot1;
          nn[x] = temp * PI_ / 180;
          wght[x] = gaussian(temp, mean1, sd1);
        } // for x
      } else {
        for(int x = 0; x < ndx; ++ x) { 
          nn[x] = rot1[1] * PI_ / 180;
          wght[x] = 1;
        }
      } // if
      if(fabs(drot2) > 1e-20) {
        drot2 /= ndx;
        for(int x = 0; x < ndx; ++ x) {
          real_t temp = rot2[1] + x * drot2;
          nn[ndx + x] =  temp * PI_ / 180;
          wght[ndx + x] = gaussian(temp, mean2, sd2);
        } // for x
      } else {
        for(int x = 0; x < ndx; ++ x) {
          nn[ndx + x] = rot2[1] * PI_ / 180;
          wght[ndx + x] = 1;
        }
      } // if
      if(fabs(drot3) > 1e-20) {
        drot3 /= ndx;
        for(int x = 0; x < ndx; ++ x) {
          real_t temp = rot3[1] + x * drot3;
          nn[2 * ndx + x] = temp * PI_ / 180;
          wght[2 * ndx + x] = gaussian(temp, mean3, sd3);
        } // for x
      } else {
        for(int x = 0; x < ndx; ++ x) { 
          nn[2 * ndx + x] = rot3[1] * PI_ / 180;
          wght[2 * ndx + x] = 1;
        }
      } // if
    } else {
      // TODO read .ori file ...
      std::cerr << "uh-oh: I guess you wanted to read orientations from a file" << std::endl;
      std::cerr << "too bad, its not implemented yet" << std::endl;
      return false;
    } // if-else

    return true;
  } // HipGISAXS::orientation_distribution()


  bool HipGISAXS::generate_repetition_range(unsigned int min_val, unsigned int max_val,
                        int num_vals, std::vector<unsigned int>& vals) {
    // set up random number generator
    woo::MTRandomNumberGenerator rng;
    for(int i = 0; i < num_vals; ++ i) {
      unsigned int val = min_val + rng.rand() * (max_val - min_val);  // floor
      vals.push_back(val);
    } // for
    return true;
  } // HipGISAXS::generate_repetition_range()

  bool HipGISAXS::construct_repetition_distribution(const GrainRepetitions& repetition,
                            int num_grains,
                            std::vector<vector3_t>& all_grains_repeats) {
    unsigned int min_x = 0, max_x = 0;
    unsigned int min_y = 0, max_y = 0;
    unsigned int min_z = 0, max_z = 0;
    real_t mean_x = 0, sd_x = 0;  // for gaussian
    real_t mean_y = 0, sd_y = 0;  // for gaussian
    real_t mean_z = 0, sd_z = 0;  // for gaussian
    // vectors for values
    std::vector<unsigned int> vals_x, vals_y, vals_z;

    // for x
    switch(repetition.xstat()) {
      case stat_none:
        vals_x.push_back(repetition.xmin());
        break;
      case stat_range:
        min_x = repetition.xmin(); max_x = repetition.xmax();
        generate_repetition_range(min_x, max_x, num_grains, vals_x);
        break;
      case stat_gaussian:
        std::cout << "uh-oh: gaussian distribution for grain repetitions is not yet supported!"
              << std::endl;
        return false;
        break;
      default:
        std::cerr << "error: invalid grain repetition distribution token encountered"
              << std::endl;
        return false;
    } // switch
    // for y
    switch(repetition.ystat()) {
      case stat_none:
        vals_y.push_back(repetition.ymin());
        break;
      case stat_range:
        min_y = repetition.ymin(); max_y = repetition.ymax();
        generate_repetition_range(min_y, max_y, num_grains, vals_y);
        break;
      case stat_gaussian:
        std::cout << "uh-oh: gaussian distribution for grain repetitions is not yet supported!"
              << std::endl;
        return false;
        break;
      default:
        std::cerr << "error: invalid grain repetition distribution token encountered"
              << std::endl;
        return false;
    } // switch
    // for z
    switch(repetition.zstat()) {
      case stat_none:
        vals_z.push_back(repetition.zmin());
        break;
      case stat_range:
        min_z = repetition.zmin(); max_z = repetition.zmax();
        generate_repetition_range(min_z, max_z, num_grains, vals_z);
        break;
      case stat_gaussian:
        std::cout << "uh-oh: gaussian distribution for grain repetitions is not yet supported!"
              << std::endl;
        return false;
        break;
      default:
        std::cerr << "error: invalid grain repetition distribution token encountered"
              << std::endl;
        return false;
    } // switch

    // construct the vectors by zipping the x, y and z vals
    int sizes = std::min(std::min(vals_x.size(), vals_y.size()), vals_z.size());
    for(int i = 0; i < sizes; ++ i)
      all_grains_repeats.push_back(vector3_t(vals_x[i], vals_y[i], vals_z[i]));

    return true;
  } // HipGISAXS::construct_repetition_distribution()

  bool HipGISAXS::construct_scaling_distribution(std::vector<StatisticType> dist,
                          vector3_t mean, vector3_t stddev,
                          std::vector<int> nvals,
                          std::vector<vector3_t> &samples,
                           std::vector<real_t> &weights) {
        woo::MTRandomNumberGenerator rng;
    real_t width = 4;
    for (int i = 0; i < 3; ++ i) {
      if (dist[i] != stat_gaussian) {
        std::cerr << "error: only Gaussian distribution for scaling is implemented."
            << std::endl;
        return false;
      } // if
    } // for
    vector3_t temp;
    vector3_t min, dx;
    for (int i = 0; i < 3; ++ i) {
      min[i] = mean[i] - width * stddev[i];
      dx[i] = (2 * width) * stddev[i] / (nvals[i] + 1);
    } // for
    for (int i = 0; i < nvals[0]; ++ i) {
      for (int j = 0; j < nvals[1]; ++ j) {
        for (int k = 0; k < nvals[2]; ++ k) {
          temp[0] = min[0] + i * dx[0] + rng.rand() * dx[0];
          temp[1] = min[1] + j * dx[1] + rng.rand() * dx[1];
          temp[2] = min[2] + k * dx[2] + rng.rand() * dx[2];
          samples.push_back(temp);
          weights.push_back(gaussian3d(temp, mean, stddev));
        } // for
      } // for
    } // for
    return true;
  } // HipGISAXS::construct_scaling_distribution()

  /**
   * fitting related functions
   */

  bool HipGISAXS::update_params(const map_t& params) {
    std::cout << "Updating HipGISAXS parameters: ";
    for(map_t::const_iterator i = params.begin(); i != params.end(); ++ i)
      std::cout << (*i).first << " = " << (*i).second << "; ";
    std::cout << std::endl;
    //HiGInput::instance().print_all();
    return HiGInput::instance().update_params(params);
  } // HipGISAXS::update_params()


  /**
   * miscellaneous functions
   */

  void HipGISAXS::save_gisaxs(real_t *final_data, std::string output) {
    std::ofstream f(output.c_str());
    for(unsigned int z = 0; z < nrow_; ++ z) {
      for(unsigned int y = 0; y < ncol_; ++ y) {
        unsigned int index = ncol_ * z + y;
        f << final_data[index] << "\t";
      } // for
      f << std::endl;
    } // for
    f.close();
  } // HipGISAXS::save_gisaxs()
  
  
  void HipGISAXS::printfr(const char* name, real_t* arr, unsigned int size) {
    std::cout << name << ":" << std::endl;
    if(arr == NULL) { std::cout << "NULL" << std::endl; return; }
    for(unsigned int i = 0; i < size; ++ i) {
      std::cout << arr[i] << "\t";
    } // for
    std::cout << std::endl;
  } // HipGISAXS::printfr()
 

  void HipGISAXS::printfc(const char* name, complex_t* arr, unsigned int size) {
    std::cout << name << ":" << std::endl;
    if(arr == NULL) { std::cout << "NULL" << std::endl; return; }
    for(unsigned int i = 0; i < size; ++ i) {
      std::cout << arr[i].real() << "," << arr[i].imag() << "\t";
    } // for
    std::cout << std::endl;
  } // HipGISAXS::printfc()
 

  // for testing
  bool HipGISAXS::write_qgrid(char* filename) {
    std::ofstream qout(filename);
 
    qout << nqx_ << " " << nqy_ << " " << nqz_extended_ << std::endl;
    for(unsigned int i = 0; i < nqx_; ++ i) {
      qout << QGrid::instance().qx(i) << " ";
    } // for
    qout << std::endl;
    for(unsigned int i = 0; i < nqy_; ++ i) {
      qout << QGrid::instance().qy(i) << " ";
    } // for
    qout << std::endl;
    for(unsigned int i = 0; i < nqz_extended_; ++ i) {
      qout << QGrid::instance().qz_extended(i).real() << " "
          << QGrid::instance().qz_extended(i).imag();
    }
    qout.close();
    return true;
  } // HipGISAXS::write_qgrid()


  bool HipGISAXS::read_form_factor(FormFactor& ff, const char* filename) {
    return ff.read_form_factor(filename, nqx_, nqy_, nqz_extended_);
  } // HipGISAXS::read_form_factor()
} // namespace hig
