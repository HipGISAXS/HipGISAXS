/***
  *  $Id: ff_ana.cpp 37 2012-08-09 22:59:59Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: ff_ana.cpp
  *  Created: Jul 12, 2012
  *  Modified: Thu 21 Feb 2013 12:46:25 PM PST
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

	// TODO: decompose into two init functions:
	// 	one for overall (sets qgrid in gff_),
	// 	other for each run/invocation (sets rotation matrices)
	bool AnalyticFormFactor::init(vector3_t &rot1, vector3_t &rot2, vector3_t &rot3,
				std::vector<complex_t> &ff) {
		nqx_ = QGrid::instance().nqx();
		nqy_ = QGrid::instance().nqy();
		nqz_ = QGrid::instance().nqz_extended();

		// first make sure there is no residue from any previous computations
		ff.clear();

#ifdef FF_ANA_GPU
		gff_.init(nqx_, nqy_, nqz_);
#endif // FF_ANA_GPU

		// rotation matrices are new for each ff calculation
		rot_ = new (std::nothrow) float_t[9];
		rot_[0] = rot1[0]; rot_[1] = rot1[1]; rot_[2] = rot1[2];
		rot_[3] = rot2[0]; rot_[4] = rot2[1]; rot_[5] = rot2[2];
		rot_[6] = rot3[0]; rot_[7] = rot3[1]; rot_[8] = rot3[2];

		return true;
	} // AnalyticFormFactor::init()


	void AnalyticFormFactor::clear() {
		nqx_ = nqy_ = nqz_ = 0;
	} // AnalyticFormFactor::clear()


	bool AnalyticFormFactor::compute(ShapeName shape, float_t tau, float_t eta, vector3_t transvec,
									std::vector<complex_t>& ff,
									shape_param_list_t& params, float_t single_layer_thickness,
									vector3_t rot1, vector3_t rot2, vector3_t rot3,
									MPI::Intracomm& world_comm) {

		std::cout << "-- Computing form factor analytically ... " << std::endl;
#ifdef TIME_DETAIL_1
		woo::BoostChronoTimer compute_timer;
		compute_timer.start();
#endif // TIME_DETAIL_1

		switch(shape) {
			case shape_box:						// cube or box
				if(!compute_box(nqx_, nqy_, nqz_, ff, shape, params, tau, eta, transvec, rot1, rot2, rot3)) {
					std::cerr << "error: something went wrong while computing FF for a box"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_cylinder:				// standing cylinder
				if(!compute_cylinder(params, tau, eta, ff, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a cylinder"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_random_cylinders:		// randomly oriented cylinders, for saxs
				if(!compute_random_cylinders(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for random cylinders"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_horizontal_cylinder:		// horizontal cylinder
				if(!compute_horizontal_cylinder(tau, eta, params, transvec, ff)) {
					std::cerr << "error: something went wrong while computing FF for a horizontal cylinder"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sphere:					// simple sphere
				if(!compute_sphere(params, ff, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a sphere"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism3:					// triangular prism (prism with 3 sides)
				if(!compute_prism(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a prism3"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism6:					// hexagonal prism (prism with 6 sides)
				if(!compute_prism6(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a prism6"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sawtooth_down:			// downwards sawtooth
				if(!compute_sawtooth_down()) {
					std::cerr << "error: something went wrong while computing FF for a sawtooth down"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_sawtooth_up:				// upwards sawtooth
				if(!compute_sawtooth_up()) {
					std::cerr << "error: something went wrong while computing FF for a sawtooth up"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_prism3x:					// triangular grating in x direction
				if(!compute_prism3x(params, ff, tau, eta, transvec)) {
					std::cerr << "error: something went wrong while computing FF for a prism3x"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_pyramid:					// pyramid
				if(!compute_pyramid()) {
					std::cerr << "error: something went wrong while computing FF for a pyramid"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_trunccone:				// truncated cone
				if(!compute_truncated_cone()) {
					std::cerr << "error: something went wrong while computing FF for a truncated cone"
								<< std::endl;
					return false;
				} // if
				break;
			case shape_truncpyr:				// truncated pyramid
				if(!compute_truncated_pyramid()) {
					std::cerr << "error: something went wrong while computing FF for a truncated pyramid"
								<< std::endl;
					return false;
				} // if
				break;
			default:
				std::cerr << "error: invalid shape. given shape is not supported" << std::endl;
				return false;
		} // switch

#ifdef TIME_DETAIL_1
		compute_timer.stop();
		std::cout << "**               FF compute time: " << compute_timer.elapsed_msec() << " ms."
					<< std::endl;
#endif // TIME_DETAIL_1

		/*int rank;
		MPI_Comm_rank(world_comm, &rank);
		if(rank == 0) {
			int naninfs = count_naninfs(nqx_, nqy_, nqz_, ff);
			std::cout << " ------ " << naninfs << " / " << nqx_ * nqy_ * nqz_ << " nans or infs" << std::endl;
		} // if
*/
		return true;
	} // AnalyticFormFactor::compute()


	/**
	 * matrix computation helpers
	 */

	bool AnalyticFormFactor::mat_fq_inv_in(unsigned int x_size,
											unsigned int y_size,
											unsigned int z_size,
											complex_vec_t& matrix, float_t y) {
		for(complex_vec_t::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			*i = fq_inv(*i, y);
		} // for
		
		return true;
	} // AnalyticFormFactor::mat_fq_inv()

	bool AnalyticFormFactor::mat_fq_inv(unsigned int x_size, unsigned int y_size, unsigned int z_size,
										const complex_vec_t& matrix, float_t y, complex_vec_t& result) {
		result.clear();
		for(complex_vec_t::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			result.push_back(fq_inv(*i, y));
		} // for
		return true;
	} // AnalyticFormFactor::mat_fq_inv()

	complex_t AnalyticFormFactor::fq_inv(complex_t value, float_t y) {
		complex_t unitc(0, 1.0);
		complex_t temp = 2.0 * exp(unitc * value * y / (float_t) 2.0) *
							sin(value * y / (float_t) 2.0) / value;
		if(fabs(temp) <= 1e-14) temp = y;
		return temp;
	} // AnalyticFormFactor::fq_inv()


	bool AnalyticFormFactor::mat_sinc(unsigned int x_size, unsigned int y_size, unsigned int z_size,
										const complex_vec_t& matrix, complex_vec_t& result) {
		result.clear();
		for(std::vector<complex_t>::const_iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			result.push_back(sinc(*i));
		} // for
		return true;
	} // AnalyticFormFactor::mat_sinc()

	bool AnalyticFormFactor::mat_sinc_in(unsigned int x_size,
											unsigned int y_size,
											unsigned int z_size,
											std::vector<complex_t>& matrix) {
		for(std::vector<complex_t>::iterator i = matrix.begin(); i != matrix.end(); ++ i) {
			*i = sinc(*i);
		} // for
		return true;
	} // AnalyticFormFactor::mat_sinc()


	/**
	 * Other helpers
	 */

	complex_t AnalyticFormFactor::sinc(complex_t value) {
		complex_t temp;
		if(fabs(value.real()) <= 1e-14 && fabs(value.imag()) <= 1e-14) temp = complex_t(1.0, 0.0);
		else temp = sin(value) / value;
		return temp;
	} // AnalyticFormFactor::sinc()


	/**
	 * generate parameter's distribution
	 */
	bool AnalyticFormFactor::param_distribution(ShapeParam& param, std::vector<float_t>& dim,
												std::vector<float_t>& dim_vals) {
		if(!param.isvalid()) {
			std::cerr << "error: invalid shape parameter encountered" << std::endl;
			return false;
		} // if

		if(param.nvalues() < 1) {
			std::cerr << "error: empty parameter found (nvalues = 0)" << std::endl;
			return false;
		} // if
		float_t pmax = param.max(), pmin = param.min();
		if(pmax < pmin) pmax = pmin;
		if(param.nvalues() > 1) {
			float_t step = fabs(pmax - pmin) / (param.nvalues() - 1);
			float_t curr = pmin;
			do {
				dim.push_back(curr);
				curr += step;
			} while(curr < param.max());	// assumes min < max ...
		} else {
			dim.push_back(pmin);
		} // if-else

		if(param.stat() == stat_none || param.stat() == stat_null) {	// just one value
			dim_vals.push_back(1.0);
		} else if(param.stat() == stat_uniform) {
			for(unsigned int i = 0; i < dim.size(); ++ i) {
				dim_vals.push_back(1.0);
			} // for
		} else if(param.stat() == stat_gaussian) {
			float_t mean = param.mean();
			if(!boost::math::isfinite(mean)) {
				mean = (dim[0] + dim[dim.size() - 1]) / 2;
			} // if
			for(unsigned int i = 0; i < dim.size(); ++ i) {
				dim_vals.push_back(exp(-1.0 * pow((dim[i] - mean), 2) / (2 * pow(param.deviation(), 2)))
									/ (sqrt(2 * PI_) * param.deviation()));
			} // for
		} else if(param.stat() == stat_random) {
			std::cerr << "uh-oh: random statistic has not been implemented yet" << std::endl;
			return false;
			// ...
		} else {
			std::cerr << "error: an invalid statistic value given for shape parameter" << std::endl;
			return false;
		} // if-else

		return true;
	} // AnalyticFormFactor::param_distribution()


	void AnalyticFormFactor::compute_meshpoints(const float_t qx, const float_t qy, const complex_t qz,
												const float_t* rot,
	    		                                complex_t& mqx, complex_t& mqy, complex_t& mqz) {
		mqx = qy * rot[0] + qx * rot[1] + qz * rot[2];
		mqy = qy * rot[3] + qx * rot[4] + qz * rot[5];
		mqz = qy * rot[6] + qx * rot[7] + qz * rot[8];
	} // AnalyticFormFactor::compute_meshpoints()

} // namespace hig

