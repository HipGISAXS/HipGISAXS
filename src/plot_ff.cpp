/***
  *  $Id$
  *
  *  Project:
  *
  *  File: plot_ff.cpp
  *  Created: Aug 22, 2012
  *  Modified: Wed 22 Aug 2012 06:02:37 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "plot_ff.hpp"

int main(int narg, char** args) {

	if(narg != 6) {
		std::cout << "usage: plot_ff nx ny nz infile outfile" << std::endl;
		exit(1);
	} // if

	unsigned int nx = std::atoi(args[1]);
	unsigned int ny = std::atoi(args[2]);
	unsigned int nz = std::atoi(args[3]);
	char* infile = args[4];
	char* outfile = args[5];

	hig::Plot new_plot(nx, ny, nz);
	new_plot.plot(infile, outfile);

	return 0;
} // main()
