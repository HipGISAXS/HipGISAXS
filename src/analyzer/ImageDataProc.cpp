/**
 *  Project: HipGISAXS
 *
 *  File: ImageDataProc.cpp
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

#include <iostream>
#include <analyzer/ImageDataProc.hpp>

namespace hig{

  ImageDataProc& ImageDataProc::operator=(const ImageDataProc& idp){
    img_data_set_ = idp.get_data();
  }

  void ImageDataProc::push(string_t filename){
    ImageData img_data(filename);
    img_data_set_.push_back(img_data);
  }

  void ImageDataProc::push(ImageData img_data){
    img_data_set_.push_back(img_data);
  }

  void ImageDataProc::print(){
    std::cout<< "Data set size: " << size()<< std::endl;
    int cnt=1;
    img_data_vec_cit it;

    for(it=img_data_set_.begin(); it!=img_data_set_.end(); it++){
      std::cout << cnt << ": ";
      it->print();
      cnt++;
    }

  }
}
