/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Data.cpp
 *  Created: Dec 26, 2013
 *  Modified: Wed 08 Oct 2014 12:17:47 PM PDT
 *  Description: Abstract class for gisaxs data
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

#include <iostream>
#include <analyzer/Data.hpp>

namespace hig{

  float_vec_t Data::create_vec(float v0, float vN, float& stp, int& nv){
    float_vec_t v;
    float min = (v0<=vN ? v0 : vN);
    float max = (v0<vN  ? vN : v0);
    if(stp>0){
      nv =1+ (max-min)/stp;
      for(int i = 0; i<nv;i++){
  v.push_back(min+stp*i);
      }
    }
    else if(nv>1){
      stp = (max-min)/(nv-1);
      for(int i = 0; i<nv;i++){
  //std::cout << min+stp*i << "- ";
        v.push_back(min+stp*i);
      }
    }
    else{
      v.push_back(min);
      v.push_back(max);
    }
    return v;
  }

}
