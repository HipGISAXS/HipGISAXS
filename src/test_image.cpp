/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: test_image.cpp
 *  Created: Aug 06, 2012
 *  Modified: Tue 16 Jul 2013 11:52:22 AM PDT
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
#include <ctime>
#include <cstdlib>
#include <cmath>

#include "image.hpp"


bool make_data(float* &data, int nx, int ny) {
	if(nx < 1 || ny < 1) return false;

	srand(time(NULL));
	data = new (std::nothrow) float[nx * ny];
	for(int i = 0; i < ny; ++ i) {
		for(int j = 0 ; j < nx; ++ j) {
			data[nx * i + j] = ((float)rand() / RAND_MAX) / (pow(10, rand() % 5));
		} // for
	} // for

	return true;
} // make_data()

int main(int narg, char** args) {
	int nx = 10, ny = 10;
	float *data = NULL;
	if(!make_data(data, nx, ny)) return -1;

	for(int i = 0; i < ny; ++ i) {
		for(int j = 0; j < nx; ++ j) {
			std::cout << data[nx * i + j] << "\t";
		} // for
		std::cout << std::endl;
	} // for

	hig::Image image(nx, ny);
	image.construct_image(data);
	image.save(std::string("image_file.tif"));

	return 0;
} // main()
