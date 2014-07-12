/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: StructVar.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:21 AM PST
 *  Description: Holds the structural variable parameters
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

#ifndef _STRUCTVAR_HPP_
#define _STRUCTVAR_HPP_

#include <analyzer/typedefs.hpp>

namespace hig{

  enum VarType{
    FLOAT=0,
    INT,
    BOOL,
    STRING
  };

  enum PhysType{
    length_nm=0,
    angle_deg,
    refindex_x
  };

  typedef struct{
    float magn_;
    VarType var_type_;
    PhysType phys_type_;
  } var_attr_t;

  typedef struct{
    bool to_fit_;
    float min_;
    float max_;
    int n_vals_;
    float tol_;
  } var_fit_params_t;


  class StructVar{

  private :
    float  val_;
    float init_val_;
    float fin_val_;
    string_t uid_;
    var_attr_t attr_;
    var_fit_params_t fit_params_;

  public:
    StructVar(){}
    StructVar(float init_val, string_t uid) : init_val_(init_val), uid_(uid) {}
    StructVar(float init_val, string_t uid, var_fit_params_t params);
    StructVar(float init_val, string_t uid, var_attr_t attr , var_fit_params_t params);
    ~StructVar(){}

    /* accessors  */
    void set_val(float val) { val_= val;}
    void set_final_val(float val) { fin_val_= val;}
    void set_attr(float magn, VarType vt, PhysType pt){ attr_.magn_=magn; attr_.var_type_=vt; attr_.phys_type_ = pt ;}
    void set_attr(var_attr_t attr){ attr_.magn_=attr.magn_; attr_.var_type_=attr.var_type_; attr_.phys_type_ = attr.phys_type_ ;}
    void set_fit_params(bool bf, float min, float max, int n_vals, float tol){
      fit_params_.to_fit_=bf; fit_params_.min_=min; fit_params_.max_=max ; fit_params_.n_vals_=n_vals; fit_params_.tol_=tol;}
    void set_fit_params(  var_fit_params_t params ){
      fit_params_.to_fit_=params.to_fit_; fit_params_.min_=params.min_; fit_params_.max_=params.max_ ; fit_params_.n_vals_=params.n_vals_; fit_params_.tol_=params.tol_;}

    float get_init_val(){return init_val_;}
    float get_val(){return val_;}
    float get_fin_val() {return fin_val_;}
    string_t get_uid(){return uid_;}

    static StructVar dummy() { return  StructVar(0, "-1"); }
    void print();

  }; /* class StructVar */


} /* namespace hig */

#endif /* STRUCTVAR_HPP_ */
