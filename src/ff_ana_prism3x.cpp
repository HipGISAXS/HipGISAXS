/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_prism3x.cpp
  *  Created: Jul 12, 2012
  *  Modified: Thu 21 Feb 2013 10:57:24 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/math/special_functions/fpclassify.hpp>

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
	bool AnalyticFormFactor::compute_prism3x(shape_param_list_t& params, std::vector<complex_t>& ff,
												float_t tau, float_t eta, vector3_t transvec) {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_prism3x"
					<< std::endl;
		return false;

		std::vector <float_t> l, distr_l;
		std::vector <float_t> h, distr_h;
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			switch((*i).second.type()) {
				case param_edge:
					param_distribution((*i).second, l, distr_l);
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_xsize:
				case param_ysize:
				case param_radius:
				case param_baseangle:
					std::cerr << "warning: ignoring unwanted values for shape type 'box'" << std::endl;
					break;
				default:
					std::cerr << "warning: unknown parameters for shape type 'box'. ignoring" << std::endl;
			} // switch
		} // for

		// TODO ...
/*            Lx = 2* dims(1);		// what is are L and Lx, why are they multiplied by 2
            H = dims(2);
            L = 2* dims(3);
            d =0.85; %;  0.9512;
            gamma =0.0 * pi/180;
            FF = Triangular_Grating_FF_Matrix(QX,QY,QZ,H,L,d,gamma,Lx) .* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3))) ; */

		return true;
	} // AnalyticFormFactor::compute_prism3x()

} // namespace hig
