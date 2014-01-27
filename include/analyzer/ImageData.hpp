/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: ImageData.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:06 AM PST
 *  Description: gisaxs image data and basic operations for img processing/analysis
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _IMAGEDATA_HPP_
#define _IMAGEDATA_HPP_

#include <analyzer/enums.hpp>
#include <analyzer/Data.hpp>

namespace hig{

  class ImageData : public Data{

  private :
    float_mat_t img_;
    float_vec_t axis_par_;
    float_vec_t axis_ver_;
    frame_t frm_;
    ImgMode mode_;
    int n_par_;
    int n_ver_;
    bool is_valid_;

  public:
    ImageData(){}
    ImageData(string_t filename){ filename_ = filename; read(filename_); /* img_ =new float*[1]; axis_par_=new float[1]; axis_ver_=new float[1];*/ }
    ImageData(float_mat_t img, float_vec_t axis_par, float_vec_t axis_ver,frame_t frm):img_(img), axis_par_(axis_par), axis_ver_(axis_ver), frm_(frm){
      n_par_ = axis_par_.size();
      n_ver_ = axis_ver_.size();
    }
    ~ImageData(){ /* delete [] img_; delete [] axis_par_; delete [] axis_ver_; */}
    virtual bool load(string_t filename){return false;}
    virtual bool init(){return false;}
    /* setters */
    void set_img(float_mat_t img) { img_= img; }

    /* getters */
    float img_p (int iv, int ip) const;
    float img_q  (float qv, float qp) const;
    float_mat_t img() const {return img_;}
    int get_n_par() const { return n_par_;}
    int get_n_ver() const { return n_ver_;}

    void print() const;
    void save(string_t filename) const;
    void read(string_t filename);
    float_vec_t read_string_values(string_t line);

  }; /* class ImageData */

} /* namespace hig */

#endif /* IMAGEDATA_HPP_ */
