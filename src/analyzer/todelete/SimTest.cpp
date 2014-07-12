/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: SimTest.cpp
 *  Created: Dec 26, 2013
 *  Modified: Sun 02 Feb 2014 09:07:40 AM PST
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
#include <math.h>
#include <analyzer/SimTest.hpp>

namespace hig{

  float SimTest::sincard(float x){
    if( fabs(x) < 1e-10)
      return 1;
    else return sin(x)/x;
  }

  SimTest::SimTest(string_t inp){


  }

  /******************************************************/

  int SimTest::update_vars(float_vec_t X){
    fit_vars_.clear();
    for(float_vec_t::iterator it=X.begin(); it!=X.end(); it++){
      fit_vars_.push_back( (*it) );
    }
  }

  float SimTest::model(float qy,float qz, float_vec_t x){
    float f=10;
    int n=0;
    for(float_vec_t::iterator it=x.begin(); it!=x.end();it++){
      float xi = (*it);
      f+=  xi  *  ( pow(qy,n) + pow(qz,n) );
      //      f+=  10* exp(-( (qy-xi)*(qy-xi) + (qz-xi)*(qz-xi) )/ 20 );
      n++;
    }

    /*
    if (x.size() == 2){
      float x0= x[0];
      float x1= x[1];
      //std::cout << qy << ", " << qz << " = " << fabs( sincard(x0*qy) ) << std::endl;
      //return x0 * qy +  x1 * qz + 10 ;
      return fabs( sincard(x0*qy) *  sincard(x1*qz)  );
      //      return 10* exp(-( (qy-x0)*(qy-x0) + (qz-x1)*(qz-x1) )/ 20 );
    }
    */

    return f;
  }

  int SimTest::compute(){
    float_mat_t img;
    int nqy = qy_.size();
    int nqz = qz_.size();

    for(int iz=0;iz<nqz; iz++)
      {
	//float_vec_t img_z;
	for(int iy=0;iy<nqy; iy++){
	  img.push_back( model(qy_[iy],qz_[iz], fit_vars_) );
	}
	//img.push_back(img_z);
      }

    data_sim_ =     ImageData(img, qy_, qz_, region_qspace);
    //data_sim_.set_img(img);
    return 0;
  }

}
