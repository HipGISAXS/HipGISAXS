/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: SimTest.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:18 AM PST
 *  Description: Class that computes a forward simulation test model.
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

#ifndef _SIMTEST_HPP_
#define _SIMTEST_HPP_

#include <iostream>
#include <analyzer/SimModel.hpp>

namespace hig{

  class SimTest : public SimModel{

  private :
    /*
      float_vec_t fit_vars_;*/
    int dim_;

    float_vec_t qy_;
    float_vec_t qz_;
    frame_t frm_;

  public:
    SimTest() {}
    SimTest(int dim, float_vec_t qy, float_vec_t qz): dim_(dim), qy_(qy), qz_(qz){}
    SimTest(string_t inp);

    ~SimTest(){}
    virtual  bool init(){is_valid_=true; std::cout << "Test Simulation Initialized!\n" ; return is_valid_;}
    /*
    virtual SimTest& operator=(const SimTest &rhs){
      std::cout << "copying...\n";
      this->pdata_sim_= rhs.pdata_sim_;
      this->fit_vars_ = rhs.fit_vars_;
      this->dim_= rhs.dim_;
      this->qy_= rhs.qy_;
      this->qz_=rhs.qz_;
      this->frm_= rhs.frm_;
    }
    */

    /* computer  */
    float model(float qy,float qz, float_vec_t x);
    virtual int compute();

    /* accessors  */
    virtual int update_vars(float_vec_t X);

    float sincard(float x);

  }; /* class SimTest */

} /* namespace hig */

#endif /* SIMTEST_HPP_ */
