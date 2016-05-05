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

  // TODO check errors and return false on error
  bool MultiLayer::init(layer_list_t & layers){
    // get number of layer including substrate
    int num_layers = layers.size();

    // first layer is always the air/vacuum
    Layer substr;
    substr.thickness(1.0E+10);
    layer_citerator_t curr;
    for (curr = layers.begin(); curr != layers.end(); curr++){
      int order = curr->second.order();
      if (order == -1)
        substr = curr->second;
      else
        layers_.push_back(curr->second);
    }
    layers_.push_back(substr);

    // z value is measured from the top
    real_t z = 0;
    for (int i = 1; i < layers_.size(); i++){
      z -= layers_[i].thickness();
      layers_[i].z_val(z);
    }
    return true;
  }

  void MultiLayer::clear(){
    layers_.clear();
  }

  complex_vec_t MultiLayer::parratt_recursion(real_t alpha, real_t k0, 
          int order){

    const int NL = layers_.size();
    real_t sina = std::sin(alpha);
    complex_vec_t T; T.resize(NL, CMPLX_ZERO_);
    complex_vec_t R; R.resize(NL, CMPLX_ZERO_);
    T.back() = complex_t(1., 0);

    // calc k-value
    complex_vec_t kz; kz.resize(NL, CMPLX_ZERO_);
    for (int i = 0; i < NL; i++)
      kz[i] = - k0 * std::sqrt(sina * sina - layers_[i].one_minus_n2());

    for (int i = NL-2; i > -1; i--){
      complex_t pij = (kz[i] + kz[i+1])/(2. * kz[i]);
      complex_t mij = (kz[i] - kz[i+1])/(2. * kz[i]);
      real_t z = layers_[i].z_val();
      complex_t exp_p = std::exp(CMPLX_MINUS_ONE_ * (kz[i+1] + kz[i]) * z);
      complex_t exp_m = std::exp(CMPLX_MINUS_ONE_ * (kz[i+1] - kz[i]) * z);
      complex_t a00 = pij * exp_m;
      complex_t a01 = mij * std::conj(exp_p);
      complex_t a10 = mij * exp_p;
      complex_t a11 = pij * std::conj(exp_m);
      T[i] = a00 * T[i+1] + a01 * R[i+1];
      R[i] = a10 * T[i+1] + a11 * R[i+1];

    }
      
    complex_t t0 = T[0];
    for (int i = 0; i < NL; i++){
      T[i] /= t0;
      R[i] /= t0;
    }
    complex_vec_t coef;
    coef.push_back(T[order]); 
    coef.push_back(R[order]);
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
    coeff.resize(nqz, CMPLX_ZERO_);
 
    size_t nalpha = QGrid::instance().nalpha();
    complex_t Ti, Ri;
    complex_vec_t Tf, Rf;
    Tf.resize(nalpha, CMPLX_ZERO_);
    Rf.resize(nalpha, CMPLX_ZERO_);

    // Rs and Ts for incoming
    complex_vec_t coef_in = parratt_recursion(alpha_i, k0, order);
    Ti = coef_in[0];
    Ri = coef_in[1];

    // Rs and Ts for outgoing
#pragma omp parallel for
    for (int i = 0; i < nalpha; i++){
      real_t alpha = QGrid::instance().alpha(i);
      if (alpha > 0){
        complex_vec_t coef_out = parratt_recursion(alpha, k0, order);
        Tf[i] = coef_out[0];
        Rf[i] = coef_out[1];
      }
    }

    // fill in the Coefficients
#pragma omp parallel for
    for (int i = 0; i < nqy; i++){
      int j = i / ncol;
      coeff[i          ] = Ti * std::conj(Tf[j]);
      coeff[i +     nqy] = Ti * std::conj(Rf[j]);
      coeff[i + 2 * nqy] = Ri * std::conj(Tf[j]);
      coeff[i + 3 * nqy] = Ri * std::conj(Rf[j]);
    }
    return true;
  }

  void MultiLayer::debug_multilayer(){
    const int NQZ = 200;
    real_t k0 = 50.679;
    real_t alpha_i = 0.00261799;
    const int order = 1; 
    complex_vec_t ci = parratt_recursion(alpha_i, k0, order);
    real_t sin_ai = std::sin(alpha_i);

    real_t dq = 2.0/(NQZ-1);
    for (int i = 0; i < 200; i++){
      real_t qz = i * dq;
      real_t alpha = std::asin(qz/k0 - sin_ai);
      complex_vec_t cf = parratt_recursion(alpha, k0, order);
      if (alpha > 0 )
        std::cerr << ci[0] * cf[0] << std::endl;
      else
        std::cerr << CMPLX_ZERO_ << std::endl;
    }
  }

} // namespace hig
