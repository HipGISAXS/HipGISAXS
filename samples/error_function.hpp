/***
  *  Project:
  *
  *  File: error_function.hpp
  *  Created: May 17, 2013
  *  Modified: Fri 10 Jan 2014 10:06:57 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __ERROR_FUNCTION_HPP__
#define __ERROR_FUNCTION_HPP__

class ErrorFunction {

	public:

	ErrorFunction() { }
	~ErrorFunction() { }

	double operator()(float* &ref, float* &data,
						unsigned int nqx, unsigned int nqy, unsigned int nqz) const {
		if(ref == NULL || data == NULL) return 0;
		double err_sum = 0.0;
		for(int z = 0; z < nqz; ++ z) {
			for(int y = 0; y < nqy; ++ y) {
				for(int x = 0; x < nqx; ++ x) {
					unsigned int index = nqx * nqy * z + nqx * y + x;
					err_sum += fabs(ref[index] - data[index]);
				} // for x
			} // for
		} // for z
		return err_sum;
	} // operator()

}; // class ErrorFunction

#endif // __ERROR_FUNCTION_HPP__
