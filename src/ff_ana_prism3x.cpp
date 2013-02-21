/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_prism3x.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 05:08:50 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/math/special_functions/fpclassify.hpp>
//#include <boost/timer/timer.hpp>

#include "woo/timer/woo_boostchronotimers.hpp"

#include "ff_ana.hpp"
#include "shape.hpp"
#include "enums.hpp"
#include "qgrid.hpp"
#include "utilities.hpp"
#include "numeric_utils.hpp"

namespace hig {

	/**
	 * triangular grating in the x-direction
	 */
	bool AnalyticFormFactor::compute_prism3x() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_prism3x"
					<< std::endl;
		return false;
		// ...
		/*for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
				case param_height:
				case param_radius:
				case param_baseangle:
				default:
			} // switch
		} // for */

/*            Lx = 2* dims(1);
            H = dims(2);
            L = 2* dims(3);
            d =0.85; %;  0.9512;
            gamma =0.0 * pi/180;
            FF = Triangular_Grating_FF_Matrix(QX,QY,QZ,H,L,d,gamma,Lx) .* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3))) ; */

	} // AnalyticFormFactor::compute_prism3x()

} // namespace hig
