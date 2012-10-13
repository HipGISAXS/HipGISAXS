/***
  *  $Id$
  *
  *  Project:
  *
  *  File: numeric_utils.hpp
  *  Created: Oct 08, 2012
  *  Modified: Sat 13 Oct 2012 02:32:42 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <complex>
#include <cmath>

#include "typedefs.hpp"

namespace hig {

	const int MAXK = 20;

	extern double gamma(double x);
	extern complex_t cbessj(complex_t zz, int order);

} // namespace hig
