/***
  *  Project:
  *
  *  File: distance_functions.hpp
  *  Created: May 17, 2013
  *  Modified: Wed 05 Feb 2014 10:25:33 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __DISTANCE_FUNCTIONS_HPP__
#define __DISTANCE_FUNCTIONS_HPP__

#include <vector>
#include <cmath>

/**
 * Distance Functors
 */


// The base class used everywhere
class DistanceMeasure {
	public:
		virtual bool operator()(float*& ref, float*& data, unsigned int size,
								std::vector<float>& err) const = 0;
	//	virtual float operator()(float*& ref, float*& data, unsigned int size) const { }
}; // class DistanceMeasure


// sum of absolute differences
class AbsoluteDifferenceError : public DistanceMeasure {
	public:
		AbsoluteDifferenceError() { }
		~AbsoluteDifferenceError() { }

		bool operator()(float*& ref, float*& data, unsigned int size, std::vector<float>& err) const {
			if(ref == NULL || data == NULL) return false;
			double err_sum = 0.0;
			for(int i = 0; i < size; ++ i) {
				err_sum += fabs(ref[i] - data[i]);
			} // for
			err.clear();
			err.push_back((float)err_sum);
			return true;
		} // operator()

		/*float operator()(float*& ref, float*& data, unsigned int size) const {
			if(ref == NULL || data == NULL) return 0;
			double err_sum = 0.0;
			for(int i = 0; i < size; ++ i) {
				err_sum += fabs(ref[i] - data[i]);
			} // for
			return (float) err_sum;
		} // operator()
*/
}; // class AbsoluteDifferenceError


// vector of differences
class ResidualVector : public DistanceMeasure {
	public:
		ResidualVector() { }
		~ResidualVector() { }

		//std::vector<float> operator()(float*& ref, float*& data, unsigned int size) const {
		bool operator()(float*& ref, float*& data, unsigned int size, std::vector<float>& err) const {
			err.clear();
			for(int i = 0; i < size; ++ i) {
				err.push_back(ref[i] - data[i]);
			} // for
			return true;
		} // operator()
		/*
		float operator()(float*& ref, float*& data, unsigned int size) const {
			if(ref == NULL || data == NULL) return 0;
			double err_sum = 0.0;
			for(int i = 0; i < size; ++ i) {
				err_sum += fabs(ref[i] - data[i]);
			} // for
			return (float) err_sum;
		} // operator()
		*/
}; // class ResidualVector


// normalized sum of absolute differences
class AbsoluteDifferenceNorm : public DistanceMeasure {
	public:
		AbsoluteDifferenceNorm() { }
		~AbsoluteDifferenceNorm() { }

		bool operator()(float*& ref, float*& data, unsigned int size, std::vector<float>& err) const {
			if(ref == NULL || data == NULL) return false;
			double err_sum = 0.0;
			double ref_sum = 0.0;
			for(int i = 0; i < size; ++ i) {
				err_sum += fabs(ref[i] - data[i]);
				ref_sum += ref[i];
			} // for
			err_sum /= ref_sum;
			err.clear();
			err.push_back((float) err_sum);
			return true;
		} // operator()
}; // class AbsoluteDifferenceNorm


#endif // __DISTANCE_FUNCTIONS_HPP__
