/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: AnaAlgorithm.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:21:44 AM PST
 *  Description: Defines the analysis algorithm to be used in workflow with required parameters
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
 *  NON-COMMERCIAL EN USER LICENSE AGREEMENT.
 */

#ifndef _ANAALGORITHM_HPP_
#define _ANAALGORITHM_HPP_

#include <analyzer/typedefs.hpp>
#include <analyzer/enums.hpp>
#include <analyzer/Analysis.hpp>
#include <analyzer/ObjFct.hpp>
#include <analyzer/AnaInput.hpp>

/*
f(X) - f(X*) (estimated)            <= fatol
|f(X) - f(X*)| (estimated) / |f(X)| <= frtol
||g(X)||                            <= gatol
||g(X)|| / |f(X)|                   <= grtol
||g(X)|| / ||g(X0)||                <= gttol
*/

namespace hig{

  class AnaAlgorithm{
  protected :
    AlgoType type_;
    bool is_valid;
    ObjFct* pobj_fct_;
    float_t* params_;
    float tol_; // error tolerance, maybe included in params_
    int max_iter_;
    int max_hist_; /* max iteration history saved  */
    float_vec_t X0_;/* initial vector */
    float_vec_t XN_;/* final vector */

  public:
    AnaAlgorithm(){ max_iter_=200; max_hist_=100; tol_=1e-4; }
    AnaAlgorithm(AlgoType type){type_= type; max_iter_=200; max_hist_=100; tol_=1e-4; }
    AnaAlgorithm(string_t algo_str){ max_iter_=200; max_hist_=100; tol_=1e-4; }
    AnaAlgorithm(AlgoType type, float_t* params, int  max_iter, int max_hist, int tol){ type_= type; params_= params;  max_iter_=max_iter; max_hist_=max_hist; tol_=tol;}

    static AlgoType get_type(string_t algo_str);
    string_t get_type_string();

    ~AnaAlgorithm(){/* delete[] params_; */}

    /*  setters- NOTE: must check for consistency of dim_ between X0_ and ObjFct */
    void set_initial_vec(float_vec_t X0){  X0_ = X0; if(pobj_fct_) pobj_fct_->set_dim(X0_.size()); }
    void set_initial_vec(AnaInput X0);
    void set_obj_fct(ObjFct* pobj_fct) { pobj_fct_= pobj_fct;  }

    /* getters  */
    float_vec_t get_final_vec() { return XN_; }

    void generate_rand_vars(float min, float max, int dim) ;
    //    Vec tao_get_vec( float_vec_t x , int dim );

    virtual bool init(){return false;}
    virtual bool parse_params(){return false;}
    virtual Analysis run(int argc,char **argv);
    virtual void print();

  }; /* class AnaAlgorithm */

} /* namespace hig */

#endif /* ANAALGORITHM_HPP_ */
