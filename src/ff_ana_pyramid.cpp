/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_pyramid.cpp
  *  Created: Jul 12, 2012
  *  Modified: Tue 19 Feb 2013 11:42:54 AM PST
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
	 * pyramid
	 */
	bool AnalyticFormFactor::compute_pyramid() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_pyramid"
					<< std::endl;
		return false;
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
	} // AnalyticFormFactor::compute_pyramid()

} // namespace hig

