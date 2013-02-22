/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: hipgisaxs_sim.cpp
  *  Created: Dec 06, 2012
  *  Modified: Fri 22 Feb 2013 01:29:58 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>

#include "woo/timer/woo_boostchronotimers.hpp"

#include "hipgisaxs_main.hpp"
#include "typedefs.hpp"
#include "utilities.hpp"


/* The main for HipGISAXS Simulation
 */
int main(int narg, char** args) {

	if(narg != 2) {
		std::cout << "usage: hipgisaxs <input_config>" << std::endl;
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

	woo::BoostChronoTimer maintimer, readtimer; //, computetimer;
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
	//hig::HiGInput::instance().print_all();	// for testing
	readtimer.stop();
	std::cout << "**       Input construction time: " << readtimer.elapsed_msec() << " ms." << std::endl;

	//computetimer.start();
	/* run the simulation */
	if(!my_gisaxs.run_all_gisaxs(MPI::COMM_WORLD)) {
		std::cerr << "error: could not run the simulation - some error occured" << std::endl;
		//computetimer.stop();
		maintimer.stop();
		MPI::Finalize();
		return 1;
	} // if
	//computetimer.stop();
	maintimer.stop();
	//std::cout << "**         Total simulation time: " << computetimer.elapsed_msec() << " ms." << std::endl;
	std::cout << "**                    Total time: " << maintimer.elapsed_msec() << " ms." << std::endl;
	std::cout << std::endl;

	MPI::Finalize();
	return 0;
} // main()
