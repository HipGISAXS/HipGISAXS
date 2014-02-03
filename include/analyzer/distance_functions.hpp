/***
  *  Project:
  *
  *  File: distance_functions.hpp
  *  Created: May 17, 2013
  *  Modified: Sun 02 Feb 2014 06:18:54 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __DISTANCE_FUNCTIONS_HPP__
#define __DISTANCE_FUNCTIONS_HPP__

#include <vector>
#include <cmath>


class DistanceMeasure {
	public:
		virtual bool operator()(float*& ref, float*& data, unsigned int size,
								std::vector<float>& err) const { }
		virtual float operator()(float*& ref, float*& data, unsigned int size) const { }
}; // class DistanceMeasure


class AbsoluteDifferenceError : public DistanceMeasure {

	public:

		AbsoluteDifferenceError() { }
		~AbsoluteDifferenceError() { }

		float operator()(float*& ref, float*& data, unsigned int size) const {
			if(ref == NULL || data == NULL) return 0;
			double err_sum = 0.0;
			for(int i = 0; i < size; ++ i) {
				err_sum += fabs(ref[i] - data[i]);
			} // for
			return (float) err_sum;
		} // operator()

}; // class AbsoluteDifferenceError


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


#endif // __DISTANCE_FUNCTIONS_HPP__
