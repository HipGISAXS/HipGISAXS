/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: ImageDataProc.hpp
 *  Created: Dec 26, 2013
 *  Description: Stores and processes datasets of gisaxs images
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _IMAGEDATAPROC_HPP_
#define _IMAGEDATAPROC_HPP_


#include <analyzer/ImageData.hpp>
//#include "Dist.hpp"

namespace hig{

  typedef  std::vector<ImageData> img_data_vec_t;
  typedef  std::vector<ImageData>::const_iterator img_data_vec_cit;

  class ImageDataProc {
  private :
    img_data_vec_t img_data_set_;
    //    Dist* pdist_;
  public:
    ImageDataProc(){}
    ImageDataProc(string_t path){}
    img_data_vec_t get_data() const {return img_data_set_;}
    //  ImageDataProc(img_data_vec_t img_data_set, Dist* pdist){}
    ~ImageDataProc(){};

    /* copy */
    ImageDataProc& operator=(const ImageDataProc& idp);

    void clear(){img_data_set_.clear();}
    void push(string_t filename);
    void push(ImageData img_data);
    int size(){return img_data_set_.size();}

    virtual bool load(string_t path){return false;}
    virtual bool init(){return false;}

    void print() ;

  }; /* class ImageDataProc */

} /* namespace hig */

#endif /* IMAGEDATAPROC_HPP_ */
