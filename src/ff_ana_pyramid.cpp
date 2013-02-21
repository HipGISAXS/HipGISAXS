/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_pyramid.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 20 Feb 2013 05:07:33 PM PST
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

/*            RD = Dimension_Distr( dims(1,:) );[rows, NRR] = size(RD); RR = RD(1,:); RRD = RD(2,:);
            HD = Dimension_Distr( dims(2,:) );[rows, NHH] = size(HD); HH = HD(1,:); HHD = HD(2,:);
            WD = Dimension_Distr( dims(3,:) );[rows, NWW] = size(WD); WW = WD(1,:); WWD = WD(2,:);
            AD = Dimension_Distr( dims(4,:) );[rows, NAA] = size(AD); AA = AD(1,:) * pi/180 ; AAD = AD(2,:);

            for id4=1:NAA
                ang = AA(id4);
                tanang = tan(ang);
                qxy_s = (qx + qy)/tanang;
                qxy_d = (qx - qy)/tanang;
                q1 = ( qz + qxy_d )/2;
                q2 = (-qz + qxy_d )/2;
                q3 = ( qz + qxy_s )/2;
                q4 = (-qz + qxy_s )/2;

                for id2=1:NHH
                    H = HH(id2);

                    sinceiq1 = SINC_Matrix(q1 * H) .* exp( 1i*q1*H) ;%./ (q1 *H);
                    sinceiq2 = SINC_Matrix(q2 * H) .* exp(-1i*q2*H) ;%./ (q2 *H);
                    sinceiq3 = SINC_Matrix(q3 * H) .* exp( 1i*q3*H) ;%./ (q3 *H);
                    sinceiq4 = SINC_Matrix(q4 * H) .* exp(-1i*q4*H) ;%./ (q4 *H);

                    K1 =      sinceiq1 + sinceiq2;
                    K2 = -1i*(sinceiq1 - sinceiq2);
                    K3 =      sinceiq3 + sinceiq4;
                    K4 = -1i*(sinceiq3 - sinceiq4);

                    for id1 =1: NRR
                        R = RR(id1);

                        for id3=1:NWW
                            W = WW(id3);

                            FF = FF +   RRD(id1) * HHD(id2) * WWD(id3) * AAD(id4) * (H ./ (qx .* qy ) ) .* ( ...
                                cos(qx*R - qy*W ) .* K1 + ...
                                sin(qx*R - qy*W ) .* K2 - ...
                                cos(qx*R + qy*W ) .* K3 - ...
                                sin(qx*R + qy*W ) .* K4 );
                        end
                    end
                end
            end


            %%special points
            dx= abs(qx);
            Zx = find(dx < 1e-14 );
            [Ix, Jx] = ind2sub(size(dx),Zx);

            dy= abs(qy);
            Zy = find(dy < 1e-14 );
            [Iy, Jy , Ky] = ind2sub(size(dy),Zy);

            dz= abs(qz);
            Zz = find(dz < 1e-14 );
            [Iz, Jz, Kz] = ind2sub(size(dz),Zz);

            for n=1:length(Zy)
                i0 = Iy(n); j0=Jy(n);  k0=Ky(n);
                if abs(qy(i0,j0,k0)) < 1e-14  %% qy=0
                    qy0  = 1e-6;
                    qx0 =  qx(i0,j0,k0);
                    qz0 =  qz(i0,j0,k0);

                    qxy_s0 = (qx0 + qy0)/tanang;
                    qxy_d0 = (qx0 - qy0)/tanang;
                    q1 = ( qz0 + qxy_d0 )/2;
                    q2 = (-qz0 + qxy_d0 )/2;
                    q3 = ( qz0 + qxy_s0 )/2;
                    q4 = (-qz0 + qxy_s0)/2;

                    sinceiq1 = SINC_Matrix(q1 * H) .* exp( 1i*q1*H) ;%./ (q1 *H);
                    sinceiq2 = SINC_Matrix(q2 * H) .* exp(-1i*q2*H) ;%./ (q2 *H);
                    sinceiq3 = SINC_Matrix(q3 * H) .* exp( 1i*q3*H) ;%./ (q3 *H);
                    sinceiq4 = SINC_Matrix(q4 * H) .* exp(-1i*q4*H) ;%./ (q4 *H);

                    K1 =      sinceiq1 + sinceiq2;
                    K2 = -1i*(sinceiq1 - sinceiq2);
                    K3 =      sinceiq3 + sinceiq4;
                    K4 = -1i*(sinceiq3 - sinceiq4);

                    FF(i0,j0,k0) = (H ./ (qx0 .* qy0 ) ) .* ( ...
                        cos(qx0*R - qy0*W ) .* K1 + ...
                        sin(qx0*R - qy0*W ) .* K2 - ...
                        cos(qx0*R + qy0*W ) .* K3 - ...
                        sin(qx0*R + qy0*W ) .* K4 );
                end
            end
            %%%%%%
            FF = FF  .* exp(1i* (qx*T(1) + qy*T(2) + qz*T(3))); */

	} // AnalyticFormFactor::compute_pyramid()

} // namespace hig

