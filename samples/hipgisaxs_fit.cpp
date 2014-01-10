/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_sim.cpp
 *  Created: Dec 06, 2012
 *  Modified: Fri 10 Jan 2014 10:24:09 AM PST
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
#include <mpi.h>
#include <hipgisaxs_main.hpp>

#include "error_function.hpp"


bool read_reference_gisaxs_data(char* filename, int nqx, int nqy, int nqz, float* &data) {
	data = new (std::nothrow) float[nqx * nqy * nqz];
	std::ifstream fin(filename);
	if(!fin.is_open()) {
		std::cerr << "error: failed to open reference data file" << std::endl;
		return false;
	} // if
	unsigned int i = 0;
	while(true) {
		double temp = -1.0;
		fin >> temp;
		if(fin.eof() || !fin.good()) break;
		data[i ++] = temp;
	} // while
	if(i != nqx * nqy * nqz) {
		std::cerr << "error: mismatch in reference data size" << std::endl;
		return false;
	} // if
	fin.close();
	return true;
} // read_reference_gisaxs_data()


/* The main for HipGISAXS testing
 */
int main(int narg, char** args) {

	if(narg != 6) {
		std::cout << "usage: hipgisaxs <input_config> <ref_data> <nqx> <nqy> <nqz>" << std::endl;
		return 1;
	} // if

	/* initialize hipgisaxs simulation object */
	hig::HipGISAXS my_gisaxs(narg, args);
	
	/* read input file and construct input structures */
	if(!my_gisaxs.construct_input(args[1])) {
		std::cerr << "error: failed to construct input containers" << std::endl;
		return 1;
	} // if

	my_gisaxs.fit_init();

	int nqx = atoi(args[3]), nqy = atoi(args[4]), nqz = atoi(args[5]);

	/* read the reference raw data */
	float* ref_data = NULL;
	if(!read_reference_gisaxs_data(args[2], nqx, nqy, nqz, ref_data)) {
		std::cerr << "error: failed to read reference data" << std::endl;
		return 1;
	} // if

	std::map <std::string, float> my_params;
	my_params["p1"] = 1.0;
	my_params["p2"] = 1.0;
	my_gisaxs.update_params(my_params);

	float* gisaxs_data = NULL;
	my_gisaxs.compute_gisaxs(gisaxs_data);

	ErrorFunction my_err_function;
	double error = my_err_function(ref_data, gisaxs_data, nqx, nqy, nqz);
	delete[] gisaxs_data;
	std::cout << "initial error = " << error << std::endl;

	/* fitting iterations */
	for(int i = 1; i <= 15; ++ i) {
		for(int j = 1; j <= 10; ++ j) {
			my_params["p1"] = i;
			my_params["p2"] = j;
			my_gisaxs.update_params(my_params);
			my_gisaxs.compute_gisaxs(gisaxs_data);
			error = my_err_function(ref_data, gisaxs_data, nqx, nqy, nqz);
			std::cout << i << " " << j << " error: " << error << std::endl;
			delete[] gisaxs_data;
		} // for
	} // for

	delete[] ref_data;

	return 0;
} // main()
