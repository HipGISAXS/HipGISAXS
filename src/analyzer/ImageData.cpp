/**
 *  Project: HipGISAXS
 *
 *  File: ImageData.cpp
 *  Created: Dec 26, 2013
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *          Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <analyzer/ImageData.hpp>

namespace hig {

  real_t ImageData::operator()(int i, int j) const {
    if(i < 0 || i >= n_par_ || j < 0 || j >= n_ver_) return 0;
    return data_[j * n_par_ + i];
  } // ImageData::operator()()


  real_t ImageData::operator()(unsigned int i, unsigned int j) const {
    if(i >= n_par_ || j >= n_ver_) return 0;
    return data_[j * n_par_ + i];
  } // ImageData::operator()()


  real_t ImageData::operator()(real_t qi, real_t qj) const {
    // TODO ...
    return 0;
  } // ImageData::operator()()


  void  ImageData::print() const {
    std::cout << filename_ << std::endl;
    for(int iv=0; iv<n_ver_; iv++)
      {
  for(int ip=0; ip<n_par_; ip++)
    {
      std::cout << data_[iv * n_par_ + ip] << "  ";
    }
  std::cout <<std::endl;
      }
  }

  void ImageData::save(string_t filename) const {
    std::ofstream file;
    file.open(filename.c_str());

    for(int iv=0; iv<n_ver_; iv++)
      {
        for(int ip=0; ip<n_par_; ip++)
          {
      file << data_[iv * n_par_ + ip] << "  ";
          }
  file << "\n";
      }
    file.close();
  }

  real_vec_t ImageData::read_string_values(string_t line){
    real_vec_t array;
    std::stringstream ssin(line);
    std::copy(std::istream_iterator<real_t>(ssin),
          std::istream_iterator<real_t>(),
          std::back_inserter(array));
    return array;
  } // ImageData::read_string_values()

  bool ImageData::read(string_t filename) {
    int nv = -1;
    int np = 0;
    data_.clear();

    std::ifstream file(filename);
    if(!file.is_open()) {
      std::cerr << "error: unable to open file " << filename << std::endl;
      return false;
    } // if

    std::string line;
    while(!file.eof()) {
      getline(file,line);
      ++ nv;
      real_vec_t img_z = read_string_values(line);
      if(nv == 0) np = img_z.size();
      data_.insert(data_.end(), img_z.begin(), img_z.end());
    } // while

    file.close();
    n_par_ = np;
    n_ver_ = nv;

    return true;
  } // ImageData::read()

} // namespace hig
