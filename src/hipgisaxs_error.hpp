/***
  *  Project:
  *
  *  File: hipgisaxs_error.hpp
  *  Created: May 17, 2013
  *  Modified: Sat 18 May 2013 10:21:12 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __HIPGISAXS_ERROR_HPP__
#define __HIPGISAXS_ERROR_HPP__

#include "typedefs.hpp"
#include "hipgisaxs_main.hpp"

namespace hig {

	class HipGISAXSError {

		public:

		HipGISAXSError() { }
		~HipGISAXSError() { }

		double operator()(const float_t* ref, const float_t* data,
							unsigned int nqx, unsigned int nqy, unsigned int nqz) {
			// there is no error checking for memory allocated for ref or data ... FIXME ...
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

	}; // class HipGISAXSError

} // namespace

#endif // __HIPGISAXS_ERROR_HPP__
