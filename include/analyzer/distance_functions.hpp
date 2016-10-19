/***
  *  Project:
  *
  *  File: distance_functions.hpp
  *  Created: May 17, 2013
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __DISTANCE_FUNCTIONS_HPP__
#define __DISTANCE_FUNCTIONS_HPP__

#include <vector>
#include <cmath>

#include <common/typedefs.hpp>
#include <common/constants.hpp>


/**
 * The base class used to define distance functors
 */

class DistanceMeasure {
  public:
    virtual bool operator()(hig::real_t*& ref, hig::real_t*& data,
                            unsigned int*& mask, unsigned int size,
                            std::vector<hig::real_t>& dist,
                            hig::real_t* mean = NULL) const = 0;
  protected:

    // L1 norm
    hig::real_t norm_l1(hig::real_t* arr, unsigned int*& mask, unsigned int size) const {
      hig::real_t sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) sum += mask[i] * fabs(arr[i]);
      return sum;
    } // norm_l1()

    // L2 norm
    hig::real_t norm_l2(hig::real_t* arr, unsigned int*& mask, unsigned int size) const {
      hig::real_t sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t temp = arr[i] * arr[i];
        sum += mask[i] * (std::isnan(temp) ? 0 : temp);
      } // for
      return sqrt(sum);
    } // norm_l2()

    // element-wise product
    hig::real_t vec_dot(hig::real_t* a, hig::real_t* b, unsigned int*& mask, unsigned int size) const {
      hig::real_t sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t temp = fabs(a[i] * b[i]);
        sum += mask[i] * (std::isnan(temp) ? 0 : temp);
      } // for
      return sum;
    } // vec_dot()

}; // class DistanceMeasure


/**
 * Various Distance Functors to be used with fitting algorithms
 */


// sum of absolute differences
class AbsoluteDifferenceError : public DistanceMeasure {
  public:
    AbsoluteDifferenceError() { }
    ~AbsoluteDifferenceError() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || data == NULL) return false;
      double dist_sum = 0.0;
      for(int i = 0; i < size; ++ i) {
        dist_sum += mask[i] * fabs(ref[i] - data[i]);
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class AbsoluteDifferenceError


// residual vector of differences
class ResidualVector : public DistanceMeasure {
  public:
    ResidualVector() { }
    ~ResidualVector() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      dist.clear();
      for(int i = 0; i < size; ++ i) {
        dist.push_back(mask[i] * (ref[i] - data[i]));
      } // for
      return true;
    } // operator()
}; // class ResidualVector


class RelativeResidualVector : public DistanceMeasure {
  public:
    RelativeResidualVector() { }
    ~RelativeResidualVector() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || data == NULL) return false;
      dist.clear();
      for(int i = 0; i < size; ++ i) {
        dist.push_back(mask[i] * (data[i] - ref[i]) / fabs(ref[i]));
      } // for
      return true;
    } // operator()
}; // class RelativeResidualVector


// sum of squares of absolute differences
class AbsoluteDifferenceSquare : public DistanceMeasure {
  public:
    AbsoluteDifferenceSquare() { }
    ~AbsoluteDifferenceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || data == NULL || mask == NULL) return false;
      double dist_sum = 0.0;
      for(int i = 0; i < size; ++ i) {
        double temp = mask[i] * fabs(ref[i] - data[i]);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class AbsoluteDifferenceNorm


// sum of squares of relative absolute differences
class RelativeAbsoluteDifferenceSquare : public DistanceMeasure {
  public:
    RelativeAbsoluteDifferenceSquare() { }
    ~RelativeAbsoluteDifferenceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || data == NULL) return false;
      double dist_sum = 0.0;
      for(int i = 0; i < size; ++ i) {
        double temp = mask[i] * fabs((ref[i] - data[i]) / ref[i]);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class RelativeAbsoluteDifferenceSquare


// sum of squares of scaled relative absolute differences
class ScaledRelativeAbsoluteDifferenceSquare : public DistanceMeasure {
  private:
    bool find_minmax(hig::real_t*& arr, unsigned int*& mask, unsigned int size,
                     hig::real_t& arr_min, hig::real_t& arr_max) const {
      arr_min = arr[0]; arr_max = arr[0];
      for(auto i = 0; i < size; ++ i) {
        if(!mask[i]) continue;
        arr_min = std::min(arr_min, arr[i]);
        arr_max = std::max(arr_max, arr[i]);
      } // for
      return true;
    } // find_minmax()

  public:
    ScaledRelativeAbsoluteDifferenceSquare() { }
    ~ScaledRelativeAbsoluteDifferenceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || dat == NULL) return false;
      // find min and max of masked ref and data in order to normalize each
      hig::real_t ref_min, ref_max, dat_min, dat_max;
      find_minmax(ref, mask, size, ref_min, ref_max);
      find_minmax(dat, mask, size, dat_min, dat_max);
      hig::real_t ref_range = ref_max - ref_min;
      hig::real_t dat_range = dat_max - dat_min;
      if(ref_range < hig::TINY_) ref_range = 1.0;
      if(dat_range < hig::TINY_) dat_range = 1.0;
      double dist_sum = 0.0;
      for(int i = 0; i < size; ++ i) {
        hig::real_t scaled_ref = (ref[i] - ref_min) / ref_range;
        hig::real_t scaled_dat = (dat[i] - dat_min) / dat_range;
        double temp = mask[i] * fabs(scaled_ref - scaled_dat);
        dist_sum += temp * temp;
      } // for
      //std::cout << "DISTANCE: " << dist_sum << std::endl;
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class ScaledRelativeAbsoluteDifferenceSquare


// normalized sum of squares of absolute differences
class AbsoluteDifferenceSquareNorm : public DistanceMeasure {
  public:
    AbsoluteDifferenceSquareNorm() { }
    ~AbsoluteDifferenceSquareNorm() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || data == NULL) return false;
      double dist_sum = 0.0;
      double ref_sum = 0.0;
      for(int i = 0; i < size; ++ i) {
        double temp = mask[i] * fabs(ref[i] - data[i]);
        dist_sum += temp * temp;
        ref_sum += mask[i] * ref[i] * ref[i];
      } // for
      dist_sum /= ref_sum;
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class AbsoluteDifferenceNorm

//! L-2 norm of difference in logrithms \f$ d = \| \ln R - \ln S \|_2 \f$
class LogDifferenceNorm2 : public DistanceMeasure {
  public:
    bool operator()(hig::real_t *& ref, hig::real_t *& data, unsigned *& mask, unsigned size,
      std::vector<hig::real_t> & dist) const {
      if (ref == NULL || data == NULL) return false;
      double dist_sum = 0.;
      for (int i = 0; i < size; i++ ){
        if (mask[i] != 0){ // if mask is non-zero
          double temp = std::log(1. + data[i]) - std::log(1. +  + ref[i]);
          dist_sum = temp * temp;
        }
      }
      dist.push_back(dist_sum);
      return true;
    }
};

// normalized sum of absolute differences
class AbsoluteDifferenceNorm : public DistanceMeasure {
  public:
    AbsoluteDifferenceNorm() { }
    ~AbsoluteDifferenceNorm() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist) const {
      if(ref == NULL || data == NULL) return false;
      double dist_sum = 0.0;
      double ref_sum = 0.0;
      for(int i = 0; i < size; ++ i) {
        dist_sum += mask[i] * fabs(ref[i] - data[i]);
        ref_sum += mask[i] * ref[i];
      } // for
      dist_sum /= ref_sum;
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class AbsoluteDifferenceNorm


/**
 * normalized/scaled:
 *  NL1X : actual param vector / unitvec norm / L1
 *  NCX  : actual param vector / c norm / L2
 *  NL2X : actual param vector / unitvec norm / L2   [default]
 *  NL1M : mean param vector / unitvec norm / L1
 *  NCM  : mean param vector / c norm / L2
 *  NL2M : mean param vector / unitvec norm / L2
 */


// NL1X
class UnitVectorNormL1Distance : public DistanceMeasure {
  public:
    UnitVectorNormL1Distance() { }
    ~UnitVectorNormL1Distance() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l1(ref, mask, size);
      if(mean == NULL) dat_norm = norm_l1(dat, mask, size);
      else dat_norm = norm_l1(mean, mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = dat[i] / dat_norm,
                    n_ref = ref[i] / ref_norm;
        hig::real_t temp = mask[i] * fabs(n_dat - n_ref);
        dist_sum += temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class UnitVectorNormL1Distance


// NL2X
class UnitVectorNormL2DistanceSquare : public DistanceMeasure {
  public:
    UnitVectorNormL2DistanceSquare() { }
    ~UnitVectorNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l2(ref, mask, size);
      if(mean == NULL) dat_norm = norm_l2(dat, mask, size);
      else dat_norm = norm_l2(mean, mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = dat[i] / dat_norm,
                    n_ref = ref[i] / ref_norm;
        hig::real_t temp = mask[i] * (n_dat - n_ref);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class UnitVectorNormL2DistanceSquare


// NCX
class CNormL2DistanceSquare : public DistanceMeasure {
  public:
    CNormL2DistanceSquare() { }
    ~CNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      hig::real_t dat_norm = 0.0, dot_prod = 0.0;
      if(mean == NULL) {
        dat_norm = norm_l2(dat, mask, size);
        dot_prod = vec_dot(dat, ref, mask, size);
      } else {
        dat_norm = norm_l2(mean, mask, size);
        dot_prod = vec_dot(mean, ref, mask, size);
      } // if-else
      hig::real_t c = dot_prod / (dat_norm * dat_norm);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = c * dat[i],
                    n_ref = ref[i];
        hig::real_t temp = mask[i] * fabs(n_ref - n_dat);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class CNormL2DistanceSquare


/**
 * square-rooted versions
 */

// NL1X
class SqrtUnitVectorNormL1Distance : public DistanceMeasure {
  public:
    SqrtUnitVectorNormL1Distance() { }
    ~SqrtUnitVectorNormL1Distance() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> sqrt_ref, sqrt_dat, sqrt_mean;
      sqrt_ref.resize(size); sqrt_dat.resize(size); sqrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        sqrt_ref[i] = sqrt(ref[i]);
        sqrt_dat[i] = sqrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) sqrt_mean[i] = sqrt(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l1(&sqrt_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l1(&sqrt_dat[0], mask, size);
      else dat_norm = norm_l1(&sqrt_mean[0], mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = sqrt_dat[i] / dat_norm,
                    n_ref = sqrt_ref[i] / ref_norm;
        hig::real_t temp = mask[i] * fabs(n_dat - n_ref);
        dist_sum += temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class SqrtUnitVectorNormL1Distance


// NL2X
class SqrtUnitVectorNormL2DistanceSquare : public DistanceMeasure {
  public:
    SqrtUnitVectorNormL2DistanceSquare() { }
    ~SqrtUnitVectorNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> sqrt_ref, sqrt_dat, sqrt_mean;
      sqrt_ref.resize(size); sqrt_dat.resize(size); sqrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        sqrt_ref[i] = sqrt(ref[i]);
        sqrt_dat[i] = sqrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) sqrt_mean[i] = sqrt(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l2(&sqrt_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l2(&sqrt_dat[0], mask, size);
      else dat_norm = norm_l2(&sqrt_mean[0], mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        double n_dat = sqrt_dat[i] / dat_norm,
               n_ref = sqrt_ref[i] / ref_norm;
        double temp = mask[i] * std::fabs(n_ref - n_dat);
        temp = std::isfinite(temp) ? temp : 0.0;
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class SqrtUnitVectorNormL2DistanceSquare


// NCX
class SqrtCNormL2DistanceSquare : public DistanceMeasure {
  public:
    SqrtCNormL2DistanceSquare() { }
    ~SqrtCNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> sqrt_ref, sqrt_dat, sqrt_mean;
      sqrt_ref.resize(size); sqrt_dat.resize(size);
      if(mean != NULL) sqrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        sqrt_ref[i] = sqrt(ref[i]);
        sqrt_dat[i] = sqrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) sqrt_mean[i] = sqrt(mean[i]);
      hig::real_t dat_norm = 0.0, dot_prod = 0.0;
      if(mean == NULL) {
        dat_norm = norm_l2(&sqrt_dat[0], mask, size);
        dot_prod = vec_dot(&sqrt_dat[0], &sqrt_ref[0], mask, size);
      } else {
        dat_norm = norm_l2(&sqrt_mean[0], mask, size);
        dot_prod = vec_dot(&sqrt_mean[0], &sqrt_ref[0], mask, size);
      } // if-else
      hig::real_t c = dot_prod / (dat_norm * dat_norm);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        double n_dat = c * sqrt_dat[i],
               n_ref = sqrt_ref[i];
        double temp = mask[i] * std::fabs(n_ref - n_dat);
        temp = std::isfinite(temp) ? temp : 0.0;
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class SqrtCNormL2DistanceSquare


/**
 * cube-rooted versions
 */

// NL1X
class CbrtUnitVectorNormL1Distance : public DistanceMeasure {
  public:
    CbrtUnitVectorNormL1Distance() { }
    ~CbrtUnitVectorNormL1Distance() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> cbrt_ref, cbrt_dat, cbrt_mean;
      cbrt_ref.resize(size); cbrt_dat.resize(size); cbrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        cbrt_ref[i] = cbrt(ref[i]);
        cbrt_dat[i] = cbrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) cbrt_mean[i] = cbrt(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l1(&cbrt_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l1(&cbrt_dat[0], mask, size);
      else dat_norm = norm_l1(&cbrt_mean[0], mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = cbrt_dat[i] / dat_norm,
                    n_ref = cbrt_ref[i] / ref_norm;
        hig::real_t temp = mask[i] * fabs(n_dat - n_ref);
        dist_sum += temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class CbrtUnitVectorNormL1Distance


// NL2X [default]
class CbrtUnitVectorNormL2DistanceSquare : public DistanceMeasure {
  public:
    CbrtUnitVectorNormL2DistanceSquare() { }
    ~CbrtUnitVectorNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> cbrt_ref, cbrt_dat, cbrt_mean;
      cbrt_ref.resize(size); cbrt_dat.resize(size); cbrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        cbrt_ref[i] = cbrt(ref[i]);
        cbrt_dat[i] = cbrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) cbrt_mean[i] = cbrt(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l2(&cbrt_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l2(&cbrt_dat[0], mask, size);
      else dat_norm = norm_l2(&cbrt_mean[0], mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = cbrt_dat[i] / dat_norm,
                    n_ref = cbrt_ref[i] / ref_norm;
        hig::real_t temp = mask[i] * (n_dat - n_ref);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class CbrtUnitVectorNormL2DistanceSquare


// NCX
class CbrtCNormL2DistanceSquare : public DistanceMeasure {
  public:
    CbrtCNormL2DistanceSquare() { }
    ~CbrtCNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> cbrt_ref, cbrt_dat, cbrt_mean;
      cbrt_ref.resize(size); cbrt_dat.resize(size); cbrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        cbrt_ref[i] = cbrt(ref[i]);
        cbrt_dat[i] = cbrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) cbrt_mean[i] = cbrt(mean[i]);
      hig::real_t dat_norm = 0.0, dot_prod = 0.0;
      if(mean == NULL) {
        dat_norm = norm_l2(&cbrt_dat[0], mask, size);
        dot_prod = vec_dot(&cbrt_dat[0], &cbrt_ref[0], mask, size);
      } else {
        dat_norm = norm_l2(&cbrt_mean[0], mask, size);
        dot_prod = vec_dot(&cbrt_mean[0], &cbrt_ref[0], mask, size);
      } // if-else
      hig::real_t c = dot_prod / (dat_norm * dat_norm);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = c * cbrt_dat[i],
                    n_ref = cbrt_ref[i];
        hig::real_t temp = mask[i] * fabs(n_ref - n_dat);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class CbrtCNormL2DistanceSquare


/**
 * log versions
 */

// NL1X
class LogUnitVectorNormL1Distance : public DistanceMeasure {
  public:
    LogUnitVectorNormL1Distance() { }
    ~LogUnitVectorNormL1Distance() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> log_ref, log_dat, log_mean;
      log_ref.resize(size); log_dat.resize(size); log_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        log_ref[i] = log(ref[i]);
        log_dat[i] = log(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) log_mean[i] = log(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l1(&log_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l1(&log_dat[0], mask, size);
      else dat_norm = norm_l1(&log_mean[0], mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = log_dat[i] / dat_norm,
                    n_ref = log_ref[i] / ref_norm;
        hig::real_t temp = mask[i] * fabs(n_dat - n_ref);
        dist_sum += temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class LogUnitVectorNormL1Distance


// NL2X
class LogUnitVectorNormL2DistanceSquare : public DistanceMeasure {
  public:
    LogUnitVectorNormL2DistanceSquare() { }
    ~LogUnitVectorNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> log_ref, log_dat, log_mean;
      log_ref.resize(size); log_dat.resize(size); log_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        log_ref[i] = log(ref[i]);
        log_dat[i] = log(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) log_mean[i] = log(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l2(&log_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l2(&log_dat[0], mask, size);
      else dat_norm = norm_l2(&log_mean[0], mask, size);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = log_dat[i] / dat_norm,
                    n_ref = log_ref[i] / ref_norm;
        hig::real_t temp = mask[i] * (n_dat - n_ref);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class LogUnitVectorNormL2DistanceSquare


// NCX
class LogCNormL2DistanceSquare : public DistanceMeasure {
  public:
    LogCNormL2DistanceSquare() { }
    ~LogCNormL2DistanceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      std::vector<hig::real_t> log_ref, log_dat, log_mean;
      log_ref.resize(size); log_dat.resize(size); log_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        log_ref[i] = log(ref[i]);
        log_dat[i] = log(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) log_mean[i] = log(mean[i]);
      hig::real_t dat_norm = 0.0, dot_prod = 0.0;
      if(mean == NULL) {
        dat_norm = norm_l2(&log_dat[0], mask, size);
        dot_prod = vec_dot(&log_dat[0], &log_ref[0], mask, size);
      } else {
        dat_norm = norm_l2(&log_mean[0], mask, size);
        dot_prod = vec_dot(&log_mean[0], &log_ref[0], mask, size);
      } // if-else
      hig::real_t c = dot_prod / (dat_norm * dat_norm);
      double dist_sum = 0.0;
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = c * log_dat[i],
                    n_ref = log_ref[i];
        hig::real_t temp = mask[i] * fabs(n_ref - n_dat);
        dist_sum += temp * temp;
      } // for
      dist.clear();
      dist.push_back((hig::real_t) dist_sum);
      return true;
    } // operator()
}; // class LogCNormL2DistanceSquare


/**
 * with residual vector for pounders
 */

// uses unit-length normalization/scaling -- used in pounders
class UnitLengthNormalizedResidualVector : public DistanceMeasure {
  public:
    UnitLengthNormalizedResidualVector() { }
    ~UnitLengthNormalizedResidualVector() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || data == NULL) return false;
      hig::real_t ref_norm = norm_l2(ref, mask, size);
      hig::real_t dat_norm = norm_l2(data, mask, size);
      dist.clear();
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_data = data[i] / dat_norm,
                    n_ref = ref[i] / ref_norm;
        hig::real_t temp = mask[i] * fabs(n_data - n_ref);
        //temp *= temp;
        dist.push_back(temp);
      } // for
      return true;
    } // operator()
}; // class UnitLengthNormalizedResidualVector


/**
 * sqrt unit vector normalized L2 distance with residual vector for pounders
 */
class SqrtUnitVectorNormL2DistanceSquareResidualVector : public DistanceMeasure {
  public:
    SqrtUnitVectorNormL2DistanceSquareResidualVector() { }
    ~SqrtUnitVectorNormL2DistanceSquareResidualVector() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      dist.clear();
      std::vector<hig::real_t> sqrt_ref, sqrt_dat, sqrt_mean;
      sqrt_ref.resize(size); sqrt_dat.resize(size);
      if(mean != NULL) sqrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        sqrt_ref[i] = sqrt(ref[i]);
        sqrt_dat[i] = sqrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) sqrt_mean[i] = sqrt(mean[i]);
      hig::real_t dat_norm = 0.0,
                  ref_norm = norm_l2(&sqrt_ref[0], mask, size);
      if(mean == NULL) dat_norm = norm_l2(&sqrt_dat[0], mask, size);
      else dat_norm = norm_l2(&sqrt_mean[0], mask, size);
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = sqrt_dat[i] / dat_norm,
                    n_ref = sqrt_ref[i] / ref_norm;
        hig::real_t temp = mask[i] * fabs(n_dat - n_ref);
        temp *= temp;
        dist.push_back(temp);
      } // for
      return true;
    } // operator()
}; // class SqrtUnitVectorNormL2DistanceSquareResidualVector


/**
 * sqrt C-norm normalized L2 distance with residual vector for pounders
 */
class SqrtCNormL2DistanceSquareResidualVector : public DistanceMeasure {
  public:
    SqrtCNormL2DistanceSquareResidualVector() { }
    ~SqrtCNormL2DistanceSquareResidualVector() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      dist.clear();
      std::vector<hig::real_t> sqrt_ref, sqrt_dat, sqrt_mean;
      sqrt_ref.resize(size); sqrt_dat.resize(size);
      if(mean != NULL) sqrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        sqrt_ref[i] = sqrt(ref[i]);
        sqrt_dat[i] = sqrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) sqrt_mean[i] = sqrt(mean[i]);
      hig::real_t dat_norm = 0.0, dot_prod = 0.0;
      if(mean == NULL) {
        dat_norm = norm_l2(&sqrt_dat[0], mask, size);
        dot_prod = vec_dot(&sqrt_dat[0], &sqrt_ref[0], mask, size);
      } else {
        dat_norm = norm_l2(&sqrt_mean[0], mask, size);
        dot_prod = vec_dot(&sqrt_mean[0], &sqrt_ref[0], mask, size);
      } // if-else
      hig::real_t c = dot_prod / (dat_norm * dat_norm);
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = c * sqrt_dat[i],
                    n_ref = sqrt_ref[i];
        hig::real_t temp = mask[i] * fabs(n_ref - n_dat);
        temp *= temp;
        dist.push_back(temp);
      } // for
      return true;
    } // operator()
}; // class SqrtCNormL2DistanceSquareResidualVector


class SqrtCNormL2DistanceResidualVector : public DistanceMeasure {
  public:
    SqrtCNormL2DistanceResidualVector() { }
    ~SqrtCNormL2DistanceResidualVector() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& dat,
                    unsigned int*& mask, unsigned int size,
                    std::vector<hig::real_t>& dist,
                    hig::real_t* mean = NULL) const {
      if(ref == NULL || dat == NULL) return false;
      dist.clear();
      std::vector<hig::real_t> sqrt_ref, sqrt_dat, sqrt_mean;
      sqrt_ref.resize(size); sqrt_dat.resize(size);
      if(mean != NULL) sqrt_mean.resize(size);
      for(unsigned int i = 0; i < size; ++ i) {
        sqrt_ref[i] = sqrt(ref[i]);
        sqrt_dat[i] = sqrt(dat[i]);
      } // for
      if(mean != NULL) for(unsigned int i = 0; i < size; ++ i) sqrt_mean[i] = sqrt(mean[i]);
      hig::real_t dat_norm = 0.0, dot_prod = 0.0;
      if(mean == NULL) {
        dat_norm = norm_l2(&sqrt_dat[0], mask, size);
        dot_prod = vec_dot(&sqrt_dat[0], &sqrt_ref[0], mask, size);
      } else {
        dat_norm = norm_l2(&sqrt_mean[0], mask, size);
        dot_prod = vec_dot(&sqrt_mean[0], &sqrt_ref[0], mask, size);
      } // if-else
      hig::real_t c = dot_prod / (dat_norm * dat_norm);
      for(unsigned int i = 0; i < size; ++ i) {
        hig::real_t n_dat = c * sqrt_dat[i],
                    n_ref = sqrt_ref[i];
        hig::real_t temp = mask[i] * fabs(n_ref - n_dat);
        dist.push_back(temp);
      } // for
      return true;
    } // operator()
}; // class SqrtCNormL2DistanceResidualVector


#endif // __DISTANCE_FUNCTIONS_HPP__
