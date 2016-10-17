/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_sim.cpp
 *  Created: Dec 06, 2012
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
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cmath>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <woo/timer/woo_boostchronotimers.hpp>

#include <sim/hipgisaxs_main.hpp>
#include <common/typedefs.hpp>
#include <utils/utilities.hpp>


/* The main for HipGISAXS Simulation
 */
int main(int narg, char** args) {

  if(narg != 2) {
    std::cout << "usage: hipgisaxs <input_config>" << std::endl;
    return 1;
  } // if

  woo::BoostChronoTimer maintimer, readtimer;
  
  maintimer.start();
  readtimer.start();

  /* read input file and construct input structures */
  hig::HipGISAXS my_gisaxs(narg, args);
  
  if(!my_gisaxs.construct_input(args[1])) {
    std::cerr << "error: failed to construct input containers" << std::endl;
    return 1;
  } // if

  readtimer.stop();

  /* run the simulation */
  if(!my_gisaxs.run_all_gisaxs()) {
    std::cerr << "error: could not run the simulation - some error occured" << std::endl;
    maintimer.stop();
    return 1;
  } // if

  maintimer.stop();

  return 0;
} // main()
