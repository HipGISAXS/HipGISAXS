/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: Dist.cpp
 *  Created: Dec 26, 2013
 *  Modified: Sun 02 Feb 2014 08:52:43 AM PST
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
#include <analyzer/Dist.hpp>
#include <analyzer/ImageData.hpp>

namespace hig{
  /*
  float_t Dist::dist(const float_t** mat1, const float_t** mat2, frame_t frm){
    return 0;
  }
  */
  float_mat_t Dist::dist(float_t* img1, float_t* img2, unsigned int np1){
    float_mat_t dm ;
	for(int ip=0; ip<np1; ip++)
	  {
	    dm.push_back( img1[ip] -  img2[ip] );
	  }
    return dm;
  }

/*  float_mat_t Dist::dist(const ImageData& img1, const ImageData& img2){
    float_mat_t dm ;
  */  /* TODO: check img1 and img2 sizes are consistent - otherwise interpolate  */
    /*int np1 = img1.get_n_par();
    int np2 = img2.get_n_par();
    if (np1!=np2) return dm;

    int nv1 = img1.get_n_ver();
    int nv2 = img2.get_n_ver();
    if (nv1!=nv2) return dm;

    for(int iv=0; iv<nv1; iv++)
      {
	float_vec_t dist_z ;
	for(int ip=0; ip<np1; ip++)
	  {
	    dist_z.push_back( img1.img_p(iv,ip) -  img2.img_p(iv,ip) );
	  }
	dm.push_back(dist_z);
      }
    return dm;
  }*/

  /*
  float_t** Dist::residual(const float_t** mat1, const float_t** mat2){

    return NULL;
  }
  */
}
