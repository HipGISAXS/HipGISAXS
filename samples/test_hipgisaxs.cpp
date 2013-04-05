/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: hipgisaxs_sim.cpp
  *  Created: Dec 06, 2012
  *  Modified: Fri 05 Apr 2013 03:31:43 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  *
  *  Description: An example file demonstrating the use of hipgisaxs library.
  */

#include <mpi.h>

#include <hipgisaxs_main.hpp>


/* The main for HipGISAXS Simulation
 */
int main(int narg, char** args) {

	if(narg != 2) {
		std::cout << "usage: test_hipgisaxs <input_config>" << std::endl;
		return 1;
	} // if

	// initialize MPI (required for now)
	MPI::Init(narg, args);

	// create hipgisaxs library
	hig::HipGISAXS my_gisaxs;

	// fill the input object in hipgisaxs object with the input file
	my_gisaxs.construct_input(args[1]);

	// run the simulation
	my_gisaxs.run_all_gisaxs(MPI::COMM_WORLD);

	// finalize MPI
	MPI::Finalize();

	return 0;
} // main()
