/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: generate_palettes.cpp
 *  Created: Oct 12, 2012
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
