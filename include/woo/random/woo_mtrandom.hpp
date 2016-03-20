/***
  *  Project: WOO Random Number Generator Library
  *
  *  File: woo_mtrandom.hpp
  *  Created: Aug 25, 2013
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __WOO_RANDOM_MT_HPP__
#define __WOO_RANDOM_MT_HPP__

#include "woorandomnumbers.hpp"
#include <random>

namespace woo {

  enum MTRandomDistribution {
    mt_uniform_dist,    /* uniform distribution */
    mt_normal_dist,     /* normal/gaussian distribution */
    mt_cauchy_dist      /* cauchy distribution */
  };

  // C++ std Mersenne-Twister random number generator
  class MTRandomNumberGenerator : public WooRandomNumberGenerator {
    private:
      // random number generator
      std::mt19937_64 mt_rand_gen_;

      // return a random number in (0,1)
      double mt_rand_01() {
        return ((double) (mt_rand_gen_() - min_) / (max_ - min_));
      } // mt_rand_01()

    public:
      // construct with 0 as seed
      MTRandomNumberGenerator() {
        mt_rand_gen_.seed(0);
        min_ = mt_rand_gen_.min();
        max_ = mt_rand_gen_.max();
        last_ = -1.0;  // nothing
      } // MTRandomNumberGenerator()

      // construct with a given seed
      MTRandomNumberGenerator(unsigned int seed) {
        mt_rand_gen_.seed(seed);
        min_ = mt_rand_gen_.min();
        max_ = mt_rand_gen_.max();
        last_ = -1.0;  // nothing
      } // MTRandomNumberGenerator()

      ~MTRandomNumberGenerator() { }

      void reset() {
        mt_rand_gen_.seed(0);
        last_ = -1.0;
      } // reset()

      void reset(unsigned int seed) {
        mt_rand_gen_.seed(seed);
        last_ = -1.0;
      } // reset()

      //double min() { return min_; }

      //double max() { return max_; }

      // returns the next random number
      double rand() {
        last_ = mt_rand_01();
        return last_;
      } // rand()

      double rand_last() { return last_; }
  }; // class WooRandomNumberGenerator

  // C++ std Mersenne-Twister random number generator with normal distribution
  class MTNormalRandomNumberGenerator : public WooRandomNumberGenerator {
    private:
      // random number generator
      std::mt19937_64 mt_rand_gen_;
      std::normal_distribution<double> dist_;
      double mean_;
      double sd_;

      // return a random number
      double mt_rand() {
        return dist_(mt_rand_gen_);
      } // mt_rand_01()

    public:
      MTNormalRandomNumberGenerator(double mean, double sd): dist_(mean, sd) {
        mean_ = mean;
        sd_ = sd;
        min_ = mt_rand_gen_.min();
        max_ = mt_rand_gen_.max();
        last_ = -1.0;  // nothing
      } // MTRandomNumberGenerator()

      // construct with a given seed
      MTNormalRandomNumberGenerator(double mean, double sd, unsigned int seed): dist_(mean, sd) {
        mt_rand_gen_.seed(seed);
        mean_ = mean;
        sd_ = sd;
        min_ = mt_rand_gen_.min();
        max_ = mt_rand_gen_.max();
        last_ = -1.0;  // nothing
      } // MTRandomNumberGenerator()

      ~MTNormalRandomNumberGenerator() { }

      void reset() {
        mt_rand_gen_.seed(0);
        last_ = -1.0;
      } // reset()

      void reset(unsigned int seed) {
        mt_rand_gen_.seed(seed);
        last_ = -1.0;
      } // reset()

      // returns the next random number
      double rand() {
        last_ = mt_rand();
        return last_;
      } // rand()

      // returns N random numbers
      std::vector <double> rand (int n) {
        std::vector <double> rand(n);
        for(int i = 0; i <  n; i++) rand[i] = mt_rand();
        last_ = rand[n - 1];
        return rand;
      } // rand()

      double rand_last() { return last_; }
  }; // class WooNormalRandomNumberGenerator

  // MT random number generator with Cauchy distribution
  class MTCauchyRandomNumberGenerator : public WooRandomNumberGenerator {
    private:
      // random number generator
      std::mt19937_64 mt_rand_gen_;
      std::cauchy_distribution <double> dist_;
      double location_;
      double scale_;

      // return a random number
      double mt_rand() {
        return dist_(mt_rand_gen_);
      } // mt_rand_01()

    public:
      MTCauchyRandomNumberGenerator(double location, double scale): dist_(location, scale) {
        location_ = location;
        scale_ = scale;
        min_ = mt_rand_gen_.min();
        max_ = mt_rand_gen_.max();
        last_ = -1.0;  // nothing
      } // MTRandomNumberGenerator()

      // construct with a given seed
      MTCauchyRandomNumberGenerator(double location, double scale, unsigned int seed) :
            dist_(location, scale) {
        mt_rand_gen_.seed(seed);
        location_ = location;
        scale_ = scale;
        min_ = mt_rand_gen_.min();
        max_ = mt_rand_gen_.max();
        last_ = -1.0;  // nothing
      } // MTRandomNumberGenerator()

      ~MTCauchyRandomNumberGenerator() { }

      void reset() {
        mt_rand_gen_.seed(0);
        dist_.reset();
        last_ = -1.0;
      } // reset()

      void reset(unsigned int seed) {
        mt_rand_gen_.seed(seed);
        dist_.reset();
        last_ = -1.0;
      } // reset()

      // returns the next random number
      double rand() {
        last_ = mt_rand();
        return last_;
      } // rand()

      // returns N random numbers
      std::vector <double> rand (int n) {
        std::vector <double> rand(n);
        for(int i = 0; i <  n; ++ i) rand[i] = mt_rand();
        last_ = rand[n - 1];
        return rand;
      } // rand()

      double rand_last() { return last_; }
  }; // class WooCauchyRandomNumberGenerator

} // namespace woo

#endif // __WOO_RANDOM_MT_HPP__
