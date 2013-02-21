/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_rand_cylinder.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 05:05:17 PM PST
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
	 * random cylinders
	 */
	bool AnalyticFormFactor::compute_random_cylinders() { // for saxs
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_random_cylinders"
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

/*            RD = Dimension_Distr( dims(1,:) );[rows, NRR] = size(RD); RR = RD(1,:); RRD = RD(2,:);
            HD = Dimension_Distr( dims(2,:) );[rows, NHH] = size(HD); HH = HD(1,:); HHD = HD(2,:);

            dx= 0.001;
            X = 0:dx:1-dx;
            nx = length(X);

            FFX = zeros(NQZ, 1);

            for id1 =1: NRR
                R = RR(id1);
                for id2=1:NHH
                    H = HH(id2);

                    FFx = zeros(NQZ, nx);
                    for ix=1:nx
                        x = X(ix);
                        FFx(:, ix) = SINC_Matrix( QZ * H * x /2 ) .* besselj(1, QZ * R * (1-x^2)^(1/2) )./( QZ * R * (1-x^2)^(1/2) ) ;
                    end

                    FFX = FFX + RRD(id1) * HHD(id2) * 4 * sum(FFx, 2).^2;
                end
            end


            FFPAR = ones(NQX, NQY);
            for iqz=1:NQZ
                FF(:,:,iqz) = FFX(iqz) * FFPAR;
            end */

	} // AnalyticFormFactor::compute_random_cylinders()

} // namespace hig
