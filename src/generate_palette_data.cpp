/***
  *  $Id$
  *
  *  Project:
  *
  *  File: generate_palettes.cpp
  *  Created: Oct 12, 2012
  *  Modified: Fri 12 Oct 2012 09:49:22 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <fstream>
#include <cstdlib>


int main(int narg, char** args) {
	if(narg != 3) {
		std::cout << "usage: generate_palettes <width> <height>" << std::endl;
		return 0;
	} // if

	int width = std::atoi(args[1]);
	int height = std::atoi(args[2]);

	double start = 0.0, end = 1.0;
	double step = (end - start) / height;

	double value = start;
	std::ofstream fdata("palette_data.dat");
	for(int y = 0; y < height; ++ y) {
		for(int x = 0; x < width; ++ x) {
			fdata << value << "\t";
		} // for
		fdata << std::endl;
		value += step;
	} // for
	fdata.close();

	return 0;
} // main()
