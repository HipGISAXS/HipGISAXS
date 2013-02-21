/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_trunc_cone.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 05:06:27 PM PST
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
	 * truncated cone
	 */
	bool AnalyticFormFactor::compute_truncated_cone() {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_cone"
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
            AD = Dimension_Distr( dims(4,:) );[rows, NAA] = size(AD); AA = AD(1,:) * pi/180; AAD = AD(2,:);

            NZ = 40;

            qpar = sqrt(qx.^2 + qy.^2);

            for id4=1:NAA
                ang = AA(id4);
                tanang = tan(ang);
                for id2=1:NHH
                    H = HH(id2);
                    dz = H/(NZ-1);
                    Z = 0: dz: H;
                    for id1 =1: NRR
                        R = RR(id1);

                        if tanang < 1e14
                            Rz = R - Z/tanang;
                        else
                            Rz = R * ones(1, NZ);
                        end

                        FF_z = zeros(NQX, NQY, NQZ, NZ);

                        for n=1:NZ %NZ:-1:1
                            FF_z(:,:,:,n) = FF_z(:,:,:,n) +  Rz(n)^2 * exp(1i* qz * Rz(n)) .* (besselj(1, qpar * Rz(n)) ./ (qpar * Rz(n)) )  ;
                        end

                        FF = FF  + RRD(id1) * HHD(id2) * AAD(id4) *  2*pi * sum(FF_z , 4);
                    end
                end
            end

            FF = FF.* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3))); */

	} // AnalyticFormFactor::compute_pyramid()

} // namespace hig

