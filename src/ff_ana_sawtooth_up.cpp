/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_sawtooth_up.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 05:09:27 PM PST
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
	 * upwards sawtooth
	 */
	bool AnalyticFormFactor::compute_sawtooth_up() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_sawtooth_up"
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

/*      RD = Dimension_Distr( dims(1,:) );[rows, NRR] = size(RD); RR = RD(1,:); RRD = RD(2,:);
        HD = Dimension_Distr( dims(2,:) );[rows, NHH] = size(HD); HH = HD(1,:); HHD = HD(2,:);
        WD = Dimension_Distr( dims(3,:) );[rows, NWW] = size(WD); WW = WD(1,:); WWD = WD(2,:);
        AD = Dimension_Distr( dims(4,:) );[rows, AWW] = size(AD); AA = AD(1,:); AAD = AD(2,:);
        Lx = 2* RR(1);
        H = HH(1);
        L = 2* WW(1);
        d =AA(1);

            gamma =0.0;
            FF = Sawtooth_Fsup_Matrix(QX,QY,QZ,H,L,d,gamma,Lx) .* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3))) ; */

	} // AnalyticFormFactor::compute_sawtooth_up()

} // namespace hig

