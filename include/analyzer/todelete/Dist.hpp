/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Dist.hpp
 *  Created: Dec 26, 2013
 *  Modified: Sun 02 Feb 2014 09:13:37 AM PST
 *  Description: Base class for computing residuals & distances between data, in various senses
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

#ifndef _DIST_HPP_
#define _DIST_HPP_

#include <analyzer/typedefs.hpp>

namespace hig{

  class ImageData;

  class Dist{
  public:
    Dist(){}
    ~Dist(){}
    //    virtual float_t dist(const float_t** mat1, const float_t** mat2, frame_t frm);
    //virtual float_mat_t dist(const ImageData& img1, const ImageData& img2);
    virtual float_vec_t dist(float* img1, float* img2, unsigned int);
    //float_t** residual(const float_t** mat1, const float_t** mat2);
  }; /* class Dist */

} /* namespace hig */

#endif /* DIST_HPP_ */
