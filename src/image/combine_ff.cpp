/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: combine_ff.cpp
 *  Created: Aug 22, 2012
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

#include <image/combine_ff.hpp>

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

  hig::PlotCombined new_plot(nx, ny, nz);
  new_plot.plot(infile, outfile);

  return 0;
} // main()
