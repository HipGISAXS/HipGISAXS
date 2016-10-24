/**
 *  Project:
 *
 *  File: edf2ascii.cpp
 *  Created: Oct 21, 2016
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <file/edf_reader.hpp>

int main(int narg, char** args) {
  std::string infile(args[1]);
  std::string outfile(args[2]);

  hig::EDFReader edf(infile.c_str());
  hig::real_t *data = NULL;
  unsigned int nx = 0, ny = 0;
  edf.get_data(data, nx, ny);
  std::cout << "Data size: " << nx << " x " << ny << std::endl;

  std::ofstream out(outfile);
  for(auto i = 0; i < ny; ++ i) {
    for(auto j = 0; j < nx; ++ j) {
      out << data[i * nx + j] << "\t";
    } // for
    out << std::endl;
  } // for
  out.close();

  return 0;
} // main()
