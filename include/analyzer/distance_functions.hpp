/***
  *  Project:
  *
  *  File: distance_functions.hpp
  *  Created: May 17, 2013
  *  Modified: Wed 08 Oct 2014 12:11:36 PM PDT
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
 * Distance Functors
 */


// The base class used everywhere
class DistanceMeasure {
  public:
    virtual bool operator()(hig::real_t*& ref, hig::real_t*& data,
                            unsigned int*& mask, unsigned int size,
                            std::vector<hig::real_t>& dist) const = 0;
  //  virtual hig::real_t operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int size) const { }
}; // class DistanceMeasure


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


#endif // __DISTANCE_FUNCTIONS_HPP__
