/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: shape2hdf5.cpp
 *  Created: Aug 25, 2012
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

#include <iomanip>
#include <file/shape2hdf5.hpp>


shape2hdf5_converter::shape2hdf5_converter(char* filename, char* outfilename, MPI_Comm comm) {
  filename_ = NULL; outfilename_ = NULL; shape_def_ = NULL;

  filename_ = new std::string(filename);
  outfilename_ = new std::string(outfilename);
  comm_ = comm;

  std::vector<double> shape_def;

  int num = load_shape(filename, shape_def);
  std::cout << "NUM: " << num << std::endl;
  convert(outfilename, shape_def);
} // o2s_converter()


int shape2hdf5_converter::load_shape(char* filename, std::vector<double> &shape_def) {
    std::ifstream f(filename);
    if(!f.is_open()) {
        std::cout << "Cannot open file " << filename << std::endl;
        return -1;
    } // if
    double s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;

    while(true) {
        f >> s;
        if(f.eof() || !f.good()) break;
        f >> nx; f >> ny; f >> nz;
        f >> cx; f >> cy; f >> cz;
        shape_def.push_back(s);
        shape_def.push_back(nx);
        shape_def.push_back(ny);
        shape_def.push_back(nz);
        shape_def.push_back(cx);
        shape_def.push_back(cy);
        shape_def.push_back(cz);
    } // while

    f.close();
    return shape_def.size() / 7;
} // load_shape()


void shape2hdf5_converter::convert(char* outfilename, std::vector<double> shape_def) {
  shape_def_ = new double[shape_def.size()];

  int count = 0;
  for(std::vector<double>::iterator i = shape_def.begin(); i != shape_def.end(); ++ i)
    shape_def_[count ++] = *i;

  // call the C HDF5 function
  s2h_converter(&shape_def_, shape_def.size() / 7, outfilename, comm_);
} // convert()
