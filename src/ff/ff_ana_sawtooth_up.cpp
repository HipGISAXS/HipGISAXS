/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_ana_sawtooth_up.cpp
 *  Created: Jul 12, 2012
 *  Modified: Sun 26 Jan 2014 10:24:33 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <boost/math/special_functions/fpclassify.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

#include <ff/ff_ana.hpp>
#include <common/enums.hpp>
#include <model/shape.hpp>
#include <model/qgrid.hpp>
#include <utils/utilities.hpp>
#include <numerics/numeric_utils.hpp>

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

