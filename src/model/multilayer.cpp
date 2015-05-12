/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: multilayer.cpp
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */


#include <model/layer.hpp>
#include <model/qgrid.hpp>
#include <config/hig_input.hpp>


namespace hig {

  MultiLayer::MultiLayer(){
    layers_.resize(0);
  }
  void MultiLayer::init(){
    // get number of layer including substrate
    int num_layers = HiGInput::instance().num_of_layers();

    // first layer is always the air/vacuum
    Layer substr;
    layer_citerator_t curr_;
    for (curr_ = HiGInput::instance().layers_begin(); 
            curr_ != HiGInput::instance().layers_end(); curr_++){
      int order = curr_->second.order();
      if (order == -1)
        substr = curr_->second;
      else
        layers_.push_back(curr_->second);
    }
    layers_.push_back(substr);

    // z value is measured from the top
    real_t z = 0;
    for (int i = 0; i < layers_.size(); i++){
      z += layers_[i].thickness();
      layers_[i].z_val(z);
    }
  }

  void MultiLayer::clear(){
    layers_.clear();
  }

  complex_vec_t MultiLayer::parratt_recursion(real_t alpha, real_t k0, 
          int order){

    real_t sina = std::sin(alpha);
    real_t kz0 = k0 * sina;
    complex_vec_t T; T.resize(layers_.size());
    complex_vec_t R; R.resize(layers_.size());

    // calc k-value
    complex_vec_t kz; kz.resize(layers_.size());
    for (int i = 0; i < layers_.size(); i++)
      kz[i] = k0 * std::sqrt(sina * sina - layers_[i].one_minus_n2());


    // Fresnel coefficents
    complex_vec_t r; r.resize(layers_.size()-1);
    complex_vec_t t; t.resize(layers_.size()-1); 
    for (int i = 0; i < layers_.size()-1; i++){
      r[i] = (kz[i] - kz[i+1]) / (kz[i] + kz[i+1]);
      t[i] = 1. + r[i];
    }

    // calc tilde{R}s
    complex_vec_t tmpR; tmpR.resize(layers_.size());
    //tmpR[layers_.size()-1] = 0;
    tmpR.back() = 0;
    for (int i = layers_.size()-2; i >= 0; i--){
      complex_t expv = std::exp(CMPLX_ONE_ * 2. * kz[i+1] * layers_[i].thickness());
      tmpR[i] = (r[i] + tmpR[i+1] * expv) / (1. + r[i] * tmpR[i+1] * expv);
    }

    complex_vec_t tmpT; tmpT.resize(layers_.size());
    tmpT[0] = 1;
    // calc tilde{T}s
    for (int i = 1; i < layers_.size()-1; i++){
      complex_t tmp1 = t[i-1] * tmpT[i-1] * std::exp(CMPLX_ONE_ * kz[i] * layers_[i-1].thickness());
      complex_t tmp2 = 1. + r[i-1] * tmpR[i] * std::exp(CMPLX_ONE_ * 2. * kz[i] * layers_[i-1].thickness());
      tmpT[i] = tmp1 / tmp2;
    }
    tmpT[layers_.size()-1] = t.back() * tmpT[layers_.size()-2];

    for (int i = 0; i < layers_.size()-1; i++){
      T[i] = tmpT[i] * std::exp(CMPLX_MINUS_ONE_ * kz[i] * layers_[i].z_val());
      R[i] = tmpR[i] * std::exp(CMPLX_ONE_ * kz[i] * layers_[i].z_val());
    }
    complex_vec_t coef;
    coef.push_back(T[order]); coef.push_back(R[order]);
    return coef;
  }

  bool MultiLayer::propagation_coeffs(complex_vec_t & coeff, real_t k0, real_t alpha_i, int order){
    if (order == -1){
      std::cerr << "Error: shapes can't be buried inside the substrate." << std::endl;
      std::cerr  << "***** if your know what your are doing," << std::endl;
      std::cerr  << "*****  add another layer with refrective index same as substrate." << std::endl;
      return false;
    }
    coeff.clear();
    int nqz = QGrid::instance().nqz_extended();
    int nqy = QGrid::instance().nqy();
    int ncol= QGrid::instance().ncols();
    coeff.resize(nqz);
 
    size_t nalpha = QGrid::instance().nalpha();
    complex_t Ti, Ri;
    complex_vec_t Tf, Rf;

    // vacuum and substrate only
    if ((order == 0) && (layers_[1].order() == -1)){
      Ti = complex_t(1., 0);
      Tf.resize(nalpha);

      // R(alpha_i)
      real_t sin_ai = std::sin(alpha_i);
      real_t kzi = -1. * k0 * sin_ai;
      complex_t tmp = std::sqrt(sin_ai * sin_ai - layers_[1].one_minus_n2());
      Ri = (sin_ai - tmp) / (sin_ai + tmp);

      // R(alpha_f)
      Rf.resize(nalpha);
      for (int i = 0; i < nalpha; i++){
        real_t sin_af = std::sin(QGrid::instance().alpha(i));
        real_t kzf = k0 * sin_af;
        tmp = std::sqrt(sin_af * sin_af - layers_[1].one_minus_n2());
        if ( kzf > 0 ){
          Rf[i] = (sin_af - tmp) / (sin_af + tmp);
          Tf[i] = complex_t(1., 0);
        } else {
          Rf[i] = Tf[i] = CMPLX_ZERO_;
        }
      }
    } 
    // one layer " This should work with multilayer "
    //else if ((layers_[order-1] == 0) && (layers_[order+1].order() == -1)){ } 
    // compute multilayers 
    else { 
      // Rs and Ts for incoming
      complex_vec_t coef_in = parratt_recursion(alpha_i, k0, order);
      Ti = coef_in[0];
      Ri = coef_in[1];

      // Rs and Ts for outgoing
      Tf.resize(nalpha);
      Rf.resize(nalpha);
#pragma omp parallel for
      for (int i = 0; i < nalpha; i++){
        real_t alpha = QGrid::instance().alpha(i);
        if (alpha > 0){
          complex_vec_t coef_out = parratt_recursion(alpha, k0, order);
          Tf[i] = coef_out[0];
          Rf[i] = coef_out[1];
        } else {
          Tf[i] = Rf[i] = CMPLX_ZERO_;
        }
      }
    }

    // fill in the Coefficients
#pragma omp parallel for
    for (int i = 0; i < nqy; i++){
      int j = i / ncol;
      coeff[i          ] = Ti * Tf[j];
      coeff[i +     nqy] = Ti * Rf[j];
      coeff[i + 2 * nqy] = Ri * Tf[j];
      coeff[i + 3 * nqy] = Ri * Rf[j];
    }
    return true;
  }
} // namespace hig
