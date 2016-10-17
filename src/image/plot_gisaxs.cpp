/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: plot_gisaxs.cpp
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

#include <image/plot_gisaxs.hpp>

int main(int narg, char** args) {

  if(narg != 9) {
    //std::cout << "usage: plot_ff nx ny nz infile outfile" << std::endl;
    std::cout << "usage: plot_ff nx ny nz infile outfile r g b" << std::endl;
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

  //hig::Plot new_plot(nx, ny, nz);
  hig::PlotGISAXS new_plot(nx, ny, nz, r, g, b);
  new_plot.plot(infile, outfile);

  return 0;
} // main()
