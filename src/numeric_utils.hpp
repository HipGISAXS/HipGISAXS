/***
  *  $Id$
  *
  *  Project:
  *
  *  File: numeric_utils.hpp
  *  Created: Oct 08, 2012
  *  Modified: Tue 09 Oct 2012 12:17:13 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <complex>
#include <cmath>

#include "typedefs.hpp"

namespace hig {

	const int MAXK = 20;

	double gamma(double x) {
		double coef[6];
		coef[0] = 76.18009173;
		coef[1] = -86.50532033;
		coef[2] = 24.01409822;
		coef[3] = -1.231739516;
		coef[4] = 0.120858003e-2;
		coef[5] = -0.536382e-5;
		double stp = 2.50662827465;
		double temp = x + 5.5;
		temp = (x + 0.5) * log(temp) - temp;
		double ser = 1.0;
		for(int j = 0; j < 6; ++ j) {
			x += 1.0;
			ser += coef[j] / x;
		} // for
		return exp(temp + log(stp * ser));
	} // gamma()

	/*--------------------------------------------------
	                        inf.     (-z^2/4)^k
    	Jnu(z) = (z/2)^nu x Sum  ------------------
	                        k=0  k! x Gamma(nu+k+1)
		(nu must be >= 0). Here k=15.
	---------------------------------------------------*/
	complex_t cbessj(complex_t zz, int order) {
		std::complex<long double> z = zz;
		std::complex<long double> temp1 = pow(z / (long double) 2.0, order);
		std::complex<long double> z2 = z * z / (long double) 4.0;
		std::complex<long double> sum(0.0, 0.0);
		long double factorial_k = 0.0;
		std::complex<long double> pow_z2_k = 1.0;
		for(int k = 0; k <= MAXK; ++ k, pow_z2_k *= z2) {
			if(k == 1) factorial_k = 1.0;	// base case
			factorial_k *= k;				// compute k!
			std::complex<long double> temp2 =
						(pow_z2_k / factorial_k) * (long double) gamma((double) order + k);
			sum += temp2;
		} // for
		temp1 *= sum;
		return complex_t((float_t) temp1.real(), (float_t) temp1.imag());
	} // cbessj()

} // namespace hig
