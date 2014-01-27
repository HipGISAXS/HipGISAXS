/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: SimModel.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:15 AM PST
 *  Description: Base class that computes a forward simulation model. Maybe HipGISAXS or other models.
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

#ifndef _SIMMODEL_HPP_
#define _SIMMODEL_HPP_

#include <iostream>
#include <analyzer/ImageData.hpp>
#include <analyzer/typedefs.hpp>

namespace hig{

  class SimModel{

  protected :
    ImageData data_sim_; /* stores the computed sim. data  */
    bool is_valid_;

    float_vec_t fit_vars_;
    //    int dim_;

    /* TODO: get these vars from ImageData
    float_vec_t qy_;
    float_vec_t qz_;
    frame_t frm_;
    */

  public:
    SimModel(){}
    ~SimModel(){}
    virtual bool init(){is_valid_=false; return is_valid_;}
    /*
    virtual SimModel& operator=(const SimModel& rhs){
      this->pdata_sim_= rhs.pdata_sim_;
      std::cout<< "copied base!\n";
    }
    */

    /*  getters */
    ImageData get_data()  { return data_sim_; }

    virtual int update_vars(float_vec_t X);
    void    generate_rand_vars(float min, float max, int dim);

    /* computer  */
    virtual int compute();

    void print();
    void save(string_t filename);

  }; /* class SimModel */

} /* namespace hig */

#endif /* SIMMODEL_HPP_ */
