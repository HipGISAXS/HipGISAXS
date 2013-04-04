/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana_trunc_cone.cpp
  *  Created: Jul 12, 2012
  *  Modified: Wed 03 Apr 2013 07:43:06 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/math/special_functions/fpclassify.hpp>

#include "woo/timer/woo_boostchronotimers.hpp"

#include "constants.hpp"
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
	bool AnalyticFormFactor::compute_truncated_cone(shape_param_list_t& params, float_t tau, float_t eta,
													std::vector<complex_t>& ff, vector3_t transvec) {
		std::cerr << "uh-oh: you reach an unimplemented part of the code, compute_truncated_cone"
			                    << std::endl;
		return false;

		// this is not working ... something is not right in matlab code ...

		std::vector <float_t> h, distr_h;	// for h dimension: param_height
		std::vector <float_t> r, distr_r;	// for r dimension: param_radius
		std::vector <float_t> a, distr_a;	// for a angle: param_baseangle
		for(shape_param_iterator_t i = params.begin(); i != params.end(); ++ i) {
			if(!(*i).second.isvalid()) {
				std::cerr << "warning: ignoring invalid shape parameter" << std::endl;
				continue;
			} // if
			switch((*i).second.type()) {
				case param_edge:
				case param_xsize:
				case param_ysize:
					std::cerr << "warning: ignoring unwanted input parameters for 'cylinder'" << std::endl;
					break;
				case param_baseangle:
					param_distribution((*i).second, a, distr_a);
					break;
				case param_height:
					param_distribution((*i).second, h, distr_h);
					break;
				case param_radius:
					param_distribution((*i).second, r, distr_r);
					break;
				default:
					std::cerr << "error: unknown parameter type given" << std::endl;
					return false;
			} // switch
		} // for

		if(h.size() < 1 || r.size() < 1 || a.size() < 1) {
			std::cerr << "error: missing parameters for truncated cone" << std::endl;
			return false;
		} // if

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
		// on cpu
		std::cout << "-- Computing truncated cone FF on CPU ..." << std::endl;

		ff.clear();
		ff.reserve(nqx_ * nqy_ * nqz_);
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) ff.push_back(complex_t(0.0, 0.0));

		unsigned int nz = 40;	// FIXME: hard-coded ... what is this???
		for(unsigned z = 0; z < nqz_; ++ z) {
			for(unsigned y = 0; y < nqy_; ++ y) {
				for(unsigned x = 0; x < nqx_; ++ x) {
					complex_t mqx, mqy, mqz;
					compute_meshpoints(QGrid::instance().qx(x), QGrid::instance().qy(y),
										QGrid::instance().qz_extended(z), rot_, mqx, mqy, mqz);
					complex_t qpar = sqrt(mqx * mqx + mqy * mqy);
					complex_t temp_ff(0.0, 0.0);
					for(unsigned int i_a = 0; i_a < a.size(); ++ i_a) {
						float_t temp1 = tan(a[i_a]);
						for(unsigned int i_h = 0; i_h < h.size(); ++ i_h) {
							float_t dz = h[i_h] / (float_t)(nz - 1);
							for(unsigned int i_r = 0; i_r < r.size(); ++ i_r) {
								float_t z_val = 0.0;
								complex_t temp_ffz(0.0, 0.0);
								for(unsigned int i_z = 0; i_z < nz; ++ i_z, z_val += dz) {
									float_t rz = r[i_r] - z_val / temp1;
									complex_t temp2 = exp(complex_t(-(mqz * rz).imag(), (mqz * rz).real()));
									complex_t temp3 = cbessj(qpar * rz, 1) / (qpar * rz);
									temp_ffz += rz * rz * temp2 * temp3;
									if((boost::math::fpclassify(temp_ffz.real()) == FP_ZERO &&
										boost::math::fpclassify(temp_ffz.imag()) == FP_ZERO) ||
										!boost::math::isfinite(temp_ffz.real()) ||
										!boost::math::isfinite(temp_ffz.imag())) {
										std::cout << "--------------------here it is: " << a[i_a] << ", "
													<< r[i_r] << ", " << h[i_h] << ", " << x << ", " << y
													<< ", " << z << ", " << i_z << std::endl;
										exit(1);
									} // if
								} // for
								temp_ff += 2 * PI_ * distr_r[i_r] * distr_h[i_h] * distr_a[i_a] * temp_ffz;
							} // for r
						} // for h
					} // for a
					complex_t temp1 = mqx * transvec[0] + mqy * transvec[1] + mqz * transvec[2];
					complex_t temp2 = exp(complex_t(-temp1.imag(), temp1.real()));
					unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
					ff[curr_index] = temp_ff * temp2;
				} // for x
			} // for y
		} // for z

	} // AnalyticFormFactor::compute_pyramid()

} // namespace hig

