/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: orient_sample.cpp
 *  Created: Sep 06, 2012
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
#include <cstdlib>

#include <utils/orient_sample.hpp>

namespace hig {

  OrientSample::OrientSample(std::string infile, std::string outfile, int rotate_around_axis):
      infile_(infile), outfile_(outfile) {
    load_shape_def();
    orient_shape_def(rotate_around_axis);
  } // OrientSample::OrientSample()


  bool OrientSample::load_shape_def() {
    std::ifstream f(infile_.c_str());
    if(!f.is_open()) {
      std::cout << "Cannot open file " << infile_ << std::endl;
      return false;
    } // if
    double s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;
    shape_def_.clear();
    while(true) {
      f >> s;
      if(f.eof() || !f.good()) break;
      f >> nx; f >> ny; f >> nz;
      f >> cx; f >> cy; f >> cz;
      shape_def_.push_back(s);
      shape_def_.push_back(nx);
      shape_def_.push_back(ny);
      shape_def_.push_back(nz);
      shape_def_.push_back(cx);
      shape_def_.push_back(cy);
      shape_def_.push_back(cz);
    } // while
    f.close();
    num_triangles_ = (unsigned int)((double)shape_def_.size() / 7.0);
    return true;
  } // OrientSample::load_shape_def()


  bool OrientSample::orient_shape_def(int rotate_around_axis) {
    // rotate_around_axis  = 0, if x
    //             = 1, if y
    //             = 2, if z
    std::vector<int> rot_matrix;
    get_rotation_matrix(rotate_around_axis, rot_matrix);
    std::ofstream f(outfile_.c_str());
    if(!f.is_open()) {
      std::cout << "Cannot open file " << outfile_ << std::endl;
      return false;
    } // if
    for(unsigned int i = 0; i < num_triangles_; ++ i) {
      std::vector<double> orig_vector, orig_normal;
      std::vector<double> new_vector, new_normal;
      double area = shape_def_[7 * i];
      orig_normal.push_back(shape_def_[7 * i + 1]);
      orig_normal.push_back(shape_def_[7 * i + 2]);
      orig_normal.push_back(shape_def_[7 * i + 3]);
      orig_vector.push_back(shape_def_[7 * i + 4]);
      orig_vector.push_back(shape_def_[7 * i + 5]);
      orig_vector.push_back(shape_def_[7 * i + 6]);

      mat_vec_mul3x1(rot_matrix, orig_vector, new_vector);
      mat_vec_mul3x1(rot_matrix, orig_normal, new_normal);

      f << area << "\t" << new_normal[0] << "\t" << new_normal[1] << "\t" << new_normal[2] << "\t"
        << new_vector[0] << "\t"<< new_vector[1] << "\t" << new_vector[2] << std::endl;
    } // for
    f.close();
  } // OrientSample::orient_shape_def()


  bool OrientSample::get_rotation_matrix(int axis, std::vector<int> &rot_matrix) {
    // precomputed rotation matrices for 90 degree counter-clockwise rotations
    if(axis == 0) {
      // around x:
      // 1 0 0
      // 0 0 -1
      // 0 1 0
      rot_matrix.push_back(1); rot_matrix.push_back(0); rot_matrix.push_back(0);
      rot_matrix.push_back(0); rot_matrix.push_back(0); rot_matrix.push_back(-1);
      rot_matrix.push_back(0); rot_matrix.push_back(1); rot_matrix.push_back(0);
    } else if(axis == 1) {
      // around y:
      // 0 0 1
      // 0 1 0
      // -1 0 0
      rot_matrix.push_back(0); rot_matrix.push_back(0); rot_matrix.push_back(1);
      rot_matrix.push_back(0); rot_matrix.push_back(1); rot_matrix.push_back(0);
      rot_matrix.push_back(-1); rot_matrix.push_back(0); rot_matrix.push_back(0);
    } else if(axis == 2) {
      // around z:
      // 0 -1 0
      // 1 0 0
      // 0 0 1
      rot_matrix.push_back(0); rot_matrix.push_back(-1); rot_matrix.push_back(0);
      rot_matrix.push_back(1); rot_matrix.push_back(0); rot_matrix.push_back(0);
      rot_matrix.push_back(0); rot_matrix.push_back(0); rot_matrix.push_back(1);
    } else { // error
      // set all to 0 and notify error
      rot_matrix.push_back(0); rot_matrix.push_back(0); rot_matrix.push_back(0);
      rot_matrix.push_back(0); rot_matrix.push_back(0); rot_matrix.push_back(0);
      rot_matrix.push_back(0); rot_matrix.push_back(0); rot_matrix.push_back(0);
      return false;
    } // if-else
    return true;
  } // OrientSample::get_rotation_matrix()


  bool OrientSample::mat_vec_mul3x1(std::vector<int> rot, std::vector<double> vec,
                    std::vector<double> &out_vec) {
    out_vec.push_back(rot[0] * vec[0] + rot[1] * vec[1] + rot[2] * vec[2]);
    out_vec.push_back(rot[3] * vec[0] + rot[4] * vec[1] + rot[5] * vec[2]);
    out_vec.push_back(rot[6] * vec[0] + rot[7] * vec[1] + rot[8] * vec[2]);

    return true;
  } // OrientSample::mat_vec_mul3x1()

} // namespace hig


int main(int narg, char** args) {
  if(narg != 4) {
    std::cout << "usage: orient_sample <infile> <outfile> <axis>" << std::endl;
    return 0;
  } // if

  std::string infile(args[1]);
  std::string outfile(args[2]);
  int axis = atoi(args[3]);

  hig::OrientSample orienter(infile, outfile, axis);
//    std::cerr << "error: could not orient, junk the output if any" << std::endl;
//    return -1;
//  } // if
  return 0;
} // main()
