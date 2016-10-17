/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: rawshape_reader.cpp
 *  Created: Aug 25, 2013
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

#include <file/rawshape_reader.hpp>

namespace hig {

  /**
   * constructor: read raw shape file
   */
  RawShapeReader::RawShapeReader(const char* filename, double* &shape_def,
                  unsigned int& num_triangles) {
    std::vector<real_t> temp_shape_def;
    load_raw(filename, temp_shape_def);
    num_triangles = temp_shape_def.size() / 7;
    shape_def = new (std::nothrow) double[temp_shape_def.size()];
    for(int i = 0; i < temp_shape_def.size(); ++ i) shape_def[i] = temp_shape_def[i];
  } // o2s_converter()


  /**
   * load input raw file into local data structures in memory
   */
  bool RawShapeReader::load_raw(const char* filename, std::vector<real_t> &shape_def) {
    std::ifstream input(filename);
    if(!input.is_open()) {
      std::cout << "Unable to open raw shape definition file " << filename << std::endl;
      return false;
    } // if

      double s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;
    while(true) {
      input >> s;
      if(input.eof() || !input.good()) break;
      input >> nx; input >> ny; input >> nz;
      input >> cx; input >> cy; input >> cz;
      shape_def.push_back(s);
      shape_def.push_back(nx);
      shape_def.push_back(ny);
      shape_def.push_back(nz);
      shape_def.push_back(cx);
      shape_def.push_back(cy);
      shape_def.push_back(cz);
    } // while
    input.close();
    return true;
  } // RawShapeReader::load_raw()


} // namespace hig
