/***
  *  $Id$
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: palette_generator.cpp
  *  Created: Aug 22, 2012
  *  Modified: Fri 12 Oct 2012 10:12:43 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "generate_palette.hpp"

int main(int narg, char** args) {

	if(narg != 9) {
		std::cout << "usage: generate_palette nx ny nz infile outfile r g b" << std::endl;
		exit(1);
	} // if

	unsigned int nx = std::atoi(args[1]);
	unsigned int ny = std::atoi(args[2]);
	unsigned int nz = std::atoi(args[3]);
	char* infile = args[4];
	char* outfile = args[5];
	unsigned int r = std::atoi(args[6]);
	unsigned int g = std::atoi(args[7]);
	unsigned int b = std::atoi(args[8]);

	hig::PaletteGenerator new_plot(nx, ny, nz, r, g, b);
	new_plot.generate(infile, outfile);

	return 0;
} // main()
