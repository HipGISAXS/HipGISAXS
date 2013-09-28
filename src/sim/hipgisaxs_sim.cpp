/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_sim.cpp
 *  Created: Dec 06, 2012
 *  Modified: Sat 28 Sep 2013 10:01:19 AM PDT
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
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "../woo/timer/woo_boostchronotimers.hpp"

#include "hipgisaxs_main.hpp"
#include "../common/typedefs.hpp"
#include "../utils/utilities.hpp"


/* The main for HipGISAXS Simulation
 */
int main(int narg, char** args) {

	if(narg != 2) {
		std::cout << "usage: hipgisaxs <input_config>" << std::endl;
		return 1;
	} // if

	/* initialize MPI */
	/*MPI::Init(narg, args);
	int mpi_rank = MPI::COMM_WORLD.Get_rank();
	int mpi_num_procs = MPI::COMM_WORLD.Get_size();

	if(mpi_rank == 0) {
		std::cout << std::endl
					<< "********************************************************************************"
					<< std::endl
					<< "***************************** HipGISAXS v0.8 ***********************************"
					<< std::endl
					<< "********************************************************************************"
					<< std::endl << std::endl;
	} // if*/

	woo::BoostChronoTimer maintimer, readtimer; //, computetimer;
	
	maintimer.start();
	readtimer.start();

	/* read input file and construct input structures */
	hig::HipGISAXS my_gisaxs(narg, args);
	
	//if(mpi_rank == 0) std::cout << "**                HiG input file: " << args[1] << std::endl;
	if(!my_gisaxs.construct_input(args[1])) {
		//if(mpi_rank == 0)
		std::cerr << "error: failed to construct input containers" << std::endl;
		//MPI::Finalize();
		return 1;
	} // if
	//hig::HiGInput::instance().print_all();	// for testing
	readtimer.stop();
	//if(mpi_rank == 0)
	//	std::cout << "**       Input construction time: " << readtimer.elapsed_msec() << " ms."
	//				<< std::endl;

	//computetimer.start();
	/* run the simulation */
	if(!my_gisaxs.run_all_gisaxs()) {
		std::cerr << "error: could not run the simulation - some error occured" << std::endl;
		//computetimer.stop();
		maintimer.stop();
		//MPI::Finalize();
		return 1;
	} // if
	//computetimer.stop();
	maintimer.stop();
	//if(mpi_rank == 0) {
		//std::cout << "**         Total simulation time: " << computetimer.elapsed_msec() << " ms."
		//				<< std::endl;
		//std::cout << "**                    Total time: " << maintimer.elapsed_msec() << " ms."
		//			<< std::endl;
		//std::cout << std::endl;
	//} // if

	//MPI::Finalize();

	return 0;
} // main()
