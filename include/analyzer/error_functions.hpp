/***
  *  Project:
  *
  *  File: error_functions.hpp
  *  Created: May 17, 2013
  *  Modified: Wed 29 Jan 2014 04:12:22 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __ERROR_FUNCTIONS_HPP__
#define __ERROR_FUNCTIONS_HPP__


class DistanceMeasure {
	public:
		virtual float operator()(float*&, float*&, unsigned int) const = 0;
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


#endif // __ERROR_FUNCTIONS_HPP__
