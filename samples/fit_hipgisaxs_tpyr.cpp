/***
  *  Project:
  *
  *  File: fit_hipgisaxs.cpp
  *  Created: May 17, 2013
  *  Modified: Sun 23 Jun 2013 10:26:58 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <mpi.h>
#include <hipgisaxs_main.hpp>

#include "error_function.hpp"

int main(int narg, char** args) {
	if(narg != 3) {
		std::cout << "usage: fit_hipgisaxs <ref_data> <input_config>" << std::endl;
		return 1;
	} // if

	MPI::Init(narg, args);

	ErrorFunction my_err_func;		// the error functor

	float new_xsize = 40.0f;
	float new_ysize = 50.0f;
	float new_height = 100.0f;
	float new_baseangle = 60.0f;

	/* construct hipgisaxs object */
	hig::HipGISAXS my_gisaxs;
	/* construct the input config from the config file */
	my_gisaxs.construct_input(args[2]);			// input file contains the initial values of fit parameters

	/* read in reference gisaxs data from file */
	my_gisaxs.read_reference_gisaxs_data(args[1]);	// make sure the dimensions are same!
	std::cout << "-- Reference data reading done." << std::endl;

	/* compute the gisaxs pattern and error wrt reference */
	std::cout << "-- Computing GISAXS error ... " << std::endl << std::flush;
	double error = my_gisaxs.compute_gisaxs_error(my_err_func, MPI::COMM_WORLD);
	std::cout << "-- Computed error = " << error << std::endl;

	/* do some looping, updating parameter values etc */
	for(int i = 0; error > 1000 && i < 10; ++ i) {
		/* get new parameter values */
		new_xsize += 1.0f;
		new_ysize += 1.0f;
		new_height += 1.0f;
		new_baseangle = (float) (((int) new_baseangle + 1) % 90);
		/* update these paramaters in the hipgisaxs object */
		my_gisaxs.update_fit_param("xsize", new_xsize, "ysize", new_ysize, "height", new_height,
									"baseangle", new_baseangle);
		std::cout << "-- Updated xsize to " << new_xsize << ", ysize to " << new_ysize
					<< ", height to " << new_height << ", and baseangle to " << new_baseangle << std::endl;

		/* compute the gisaxs pattern and new error value */
		std::cout << "-- Computing GISAXS error ... " << std::endl << std::flush;
		error = my_gisaxs.compute_gisaxs_error(my_err_func, MPI::COMM_WORLD);
		std::cout << "-- Computed error = " << error << std::endl;
	} // if
	// etc ....

	std::cout << "Finished with xsize = " << new_xsize << ", ysize = " << new_ysize
				<< ", height = " << new_height << ", baseangle = " << new_baseangle << std::endl;
	
	MPI::Finalize();

	return 0;
} // main()
