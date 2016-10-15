/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_fit.cpp
 *  Created: Jun 14, 2012
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

#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
//#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>

#include "woo/timer/woo_boostchronotimers.hpp"

#include "hipgisaxs_main.hpp"
#include "typedefs.hpp"
#include "utilities.hpp"


/* The main for HipGISAXS
 */
int main(int narg, char** args) {

  if(narg != 10) {
    std::cout << "usage: hipgisaxs-fit <input_config> <z-cut> <param1 min> <max> <num> "
          << "<param2 min> <max> <num> <num_dimensions>" << std::endl;
    return 1;
  } // if

  /* initialize MPI */
  MPI::Init(narg, args);
  int mpi_rank = MPI::COMM_WORLD.Get_rank();
  int mpi_num_procs = MPI::COMM_WORLD.Get_size();

  if(mpi_rank == 0) {
    std::cout << std::endl
          << "********************************************************************************"
          << std::endl
          << "***************************** HipGISAXS v0.01-alpha ****************************"
          << std::endl
          << "********************************************************************************"
          << std::endl << std::endl;
  } // if

  woo::BoostChronoTimer maintimer, readtimer, computetimer;
  maintimer.start();
  readtimer.start();
  /* read input file and construct input structures */
  hig::HipGISAXS my_gisaxs;
  if(mpi_rank == 0) std::cout << "**                HiG input file: " << args[1] << std::endl;
  if(!my_gisaxs.construct_input(args[1])) {
    if(mpi_rank == 0) std::cerr << "error: failed to construct input containers" << std::endl;
    MPI::Finalize();
    return 1;
  } // if
  //hig::HiGInput::instance().print_all();  // for testing
  readtimer.stop();
  std::cout << "**       Input construction time: " << readtimer.elapsed_msec() << " ms." << std::endl;

  real_t zcut = std::atof(args[2]);
  real_t radius_min = std::atof(args[3]);
  real_t radius_max = std::atof(args[4]);
  real_t radius_num = std::atof(args[5]);
  real_t sd_min = std::atof(args[6]);
  real_t sd_max = std::atof(args[7]);
  real_t sd_num = std::atof(args[8]);
  unsigned int dim = std::atoi(args[9]);

  computetimer.start();
  /* run fitting stuff */
  if(!my_gisaxs.fit_steepest_descent(zcut, radius_min, radius_max, radius_num,
                    sd_min, sd_max, sd_num, dim, MPI::COMM_WORLD)) {
    std::cerr << "error: could not complete the fitting run; some error occured" << std::endl;
    computetimer.stop();
    maintimer.stop();
    MPI::Finalize();
    return 1;
  } // if
  computetimer.stop();
  maintimer.stop();
  std::cout << "**            Total fitting time: " << computetimer.elapsed_msec() << " ms." << std::endl;
  std::cout << "**                    Total time: " << maintimer.elapsed_msec() << " ms." << std::endl;
  std::cout << std::endl;

  MPI::Finalize();
  return 0;
} // main()
