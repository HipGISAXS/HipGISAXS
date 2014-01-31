/***
  *  Project:
  *
  *  File: error_functions.hpp
  *  Created: May 17, 2013
  *  Modified: Fri 31 Jan 2014 01:27:55 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __ERROR_FUNCTIONS_HPP__
#define __ERROR_FUNCTIONS_HPP__

#include <vector>


class DistanceMeasure {
	public:
		virtual std::vector<float> operator()(float*&, float*&, unsigned int) const = 0;
}; // class DistanceMeasure


class AbsoluteDifferenceError : public DistanceMeasure {

	public:

	AbsoluteDifferenceError() { }
	~AbsoluteDifferenceError() { }

	std::vector<float> operator()(float*& ref, float*& data, unsigned int size) const {
		//if(ref == NULL || data == NULL) return 0;
		std::vector<float> err;
		for(int i = 0; i < size; ++ i) {
			err.push_back(fabs(ref[i] - data[i]));
		} // for
		return err;
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
}; // class AbsoluteDifferenceError


#endif // __ERROR_FUNCTIONS_HPP__
