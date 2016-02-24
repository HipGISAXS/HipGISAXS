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

/**
 * Distance Functors
 */


// The base class used everywhere
class DistanceMeasure {
  public:
    virtual bool operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int*& mask, unsigned int size,
                std::vector<hig::real_t>& dist) const = 0;
  //  virtual hig::real_t operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int size) const { }
}; // class DistanceMeasure


// sum of absolute differences
class AbsoluteDifferenceError : public DistanceMeasure {
  public:
    AbsoluteDifferenceError() { }
    ~AbsoluteDifferenceError() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int*& mask, unsigned int size,
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


// vector of differences
class ResidualVector : public DistanceMeasure {
  public:
    ResidualVector() { }
    ~ResidualVector() { }

    //std::vector<hig::real_t> operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int size) const {
    bool operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int*& mask, unsigned int size,
            std::vector<hig::real_t>& dist) const {
      dist.clear();
      for(int i = 0; i < size; ++ i) {
        dist.push_back(mask[i] * (ref[i] - data[i]));
      } // for
      return true;
    } // operator()
}; // class ResidualVector


// sum of squares of absolute differences
class AbsoluteDifferenceSquare : public DistanceMeasure {
  public:
    AbsoluteDifferenceSquare() { }
    ~AbsoluteDifferenceSquare() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int*& mask, unsigned int size,
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


// normalized sum of squares of absolute differences
class AbsoluteDifferenceSquareNorm : public DistanceMeasure {
  public:
    AbsoluteDifferenceSquareNorm() { }
    ~AbsoluteDifferenceSquareNorm() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int*& mask, unsigned int size,
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
        double ref_sum = 0.;
        for (int i = 0; i < size; i++ ){
          if (std::abs(mask[i]) > 1.E-18){ // if mask is non-zero
            double temp = std::log(ref[i]) - std::log(data[i]);
            double logref = std::log(ref[i]);
            dist_sum = temp * temp;
            ref_sum = mask[i] * logref * logref;
        }
      }
      dist_sum /= ref_sum;
      dist.push_back(dist_sum);
      return true;
    }
};

// normalized sum of absolute differences
class AbsoluteDifferenceNorm : public DistanceMeasure {
  public:
    AbsoluteDifferenceNorm() { }
    ~AbsoluteDifferenceNorm() { }

    bool operator()(hig::real_t*& ref, hig::real_t*& data, unsigned int*& mask, unsigned int size,
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
