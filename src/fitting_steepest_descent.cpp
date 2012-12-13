/***
  *  Project:
  *
  *  File: fitting_steepest_descent.cpp
  *  Created: Dec 06, 2012
  *  Modified: Thu 13 Dec 2012 05:19:45 AM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <gsl/gsl_linalg.h>

#include "hipgisaxs_main.hpp"


namespace hig {

	// for now m = n = 2
	bool mldivide(int m, float_t* a_mat, float_t* &c_mat) {
	/*	double *a_data = new (std::nothrow) double[m * m];
		for(int i = 0; i < m; ++ i) {
			for(int j = 0; j < m; ++ j) {
				a_data[m * i + j] = (double) a_mat[m * i + j];
			} // for
		} // for
		gsl_matrix_view a = gsl_matrix_view_array(a_data, m, m);
							      
		int s;
		gsl_matrix *x = gsl_matrix_alloc(m, m);
		gsl_permutation* p = gsl_permutation_alloc(m);
		gsl_linalg_LU_decomp(&a.matrix, p, &s);
		gsl_linalg_LU_invert(&a.matrix, p, x);

		c_mat = new (std::nothrow) float_t[m * m];
		for(int i = 0; i < m; ++ i) {
			for(int j = 0; j < m; ++ j) {
				c_mat[m * i + j] = (float_t) gsl_matrix_get(x, i, j);
			} // for
		} // for
																     
		gsl_permutation_free(p);
		gsl_matrix_free(x);
		*/

		// shortcut for 2x2 matrix a_mat
		if(m != 2) {
			std::cerr << "error: only 2x2 matrix case is implemented" << std::endl;
			return false;
		} // if
		float_t det_a = a_mat[0] * a_mat[3] - a_mat[1] * a_mat[2];
		c_mat = new (std::nothrow) float_t[m * m];
		c_mat[0] = a_mat[3] / det_a;
		c_mat[1] = - a_mat[1] / det_a;
		c_mat[2] = - a_mat[2] / det_a;
		c_mat[3] = a_mat[0] / det_a;
	} // mldivide()

	// currently only hacked for spheres, with radius and sd as two parameters
	bool HipGISAXS::fit_steepest_descent(unsigned int zcut,
					float_t radius_min, float_t radius_max, float_t radius_num,
					float_t sd_min, float_t sd_max, float_t sd_num,
					unsigned int dim, MPI::Intracomm& world_comm,
					int x_min, int x_max, int x_step) {
		int mpi_rank = world_comm.Get_rank();

		if(!init(world_comm)) return false;

		int num_alphai = 0, num_phi = 0, num_tilt = 0;;

		float_t alphai_min, alphai_max, alphai_step;
		HiGInput::instance().scattering_alphai(alphai_min, alphai_max, alphai_step);
		if(alphai_max < alphai_min) alphai_max = alphai_min;
		if(alphai_min == alphai_max || alphai_step == 0) num_alphai = 1;
		else num_alphai = (alphai_max - alphai_min) / alphai_step + 1;

		float_t phi_min, phi_max, phi_step;
		HiGInput::instance().scattering_inplanerot(phi_min, phi_max, phi_step);
		if(phi_step == 0) num_phi = 1;
		else num_phi = (phi_max - phi_min) / phi_step + 1;

		float_t tilt_min, tilt_max, tilt_step;
		HiGInput::instance().scattering_tilt(tilt_min, tilt_max, tilt_step);
		if(tilt_step == 0) num_tilt = 1;
		else num_tilt = (tilt_max - tilt_min) / tilt_step + 1;

		std::cout << "**                    Num alphai: " << num_alphai << std::endl
					<< "**                       Num phi: " << num_phi << std::endl
					<< "**                      Num tilt: " << num_tilt << std::endl;

		// prepare parameters

		std::vector<std::vector<float_t> > params;
		int num_params = 2;
		//for(int i = 0; i < num_params; ++ i) {
		//	std::vector<float_t> temp;
		//} // for
		std::vector<float_t> temp;
		float_t deltay = 0.0;
		if(radius_num <= 1)
			temp.push_back(radius_min);
		else {
			deltay = fabs(radius_max - radius_min) / (radius_num - 1);
			for(int i = 0; i < radius_num; ++ i) {
				temp.push_back(radius_min + i * deltay);
			} // for
		} // if-else
		params.push_back(temp);
		temp.clear();
		if(sd_num <= 1)
			temp.push_back(sd_min);
		else {
			float_t delta = fabs(sd_max - sd_min) / (sd_num - 1);
			for(int i = 0; i < sd_num; ++ i) {
				temp.push_back(sd_min + i * delta);
			} // for
		} // if-else
		params.push_back(temp);
		temp.clear();

		// this will work only on one shape and one structure

		const float_t err_threshold = 1e-3;
		const unsigned int max_iter = 200;

		std::vector<float_t> param_vals;
		param_vals.push_back(23.0);
		param_vals.push_back(2.0);
		std::vector<float_t> param_deltas;
		param_deltas.push_back(0.05);
		param_deltas.push_back(0.05);
		float_t gamma_const = 0.05;

		float_t alpha_i = alphai_min;
		// high level of parallelism here (alphai, phi, tilt) for dynamicity ...
		for(int i = 0; i < num_alphai; i ++, alpha_i += alphai_step) {
			float_t alphai = alpha_i * PI_ / 180;
			float_t phi = phi_min;
			for(int j = 0; j < num_phi; j ++, phi += phi_step) {
				float_t tilt = tilt_min;
				for(int k = 0; k < num_tilt; k ++, tilt += tilt_step) {

					std::cout << "-- Computing reference GISAXS "
								<< i * num_phi * num_tilt + j * num_tilt + k + 1 << " / "
								<< num_alphai * num_phi * num_tilt
								<< " [alphai = " << alpha_i << ", phi = " << phi
								<< ", tilt = " << tilt << "] ..." << std::endl;

					/* run the reference gisaxs simulation using input params */
					float_t* ref_data = NULL;
					if(!run_gisaxs(alpha_i, alphai, phi, tilt, ref_data, world_comm)) {
						if(mpi_rank == 0) std::cerr << "error: could not finish successfully" << std::endl;
						return false;
					} // if

					if(dim != 1) {
						std::cerr << "uh-oh: only 1D is supported for now" << std::endl;
						return false;
					} // if

					float_t* ref_z_cut = new (std::nothrow) float_t[nqy_];
					for(unsigned int iy = 0; iy < nqy_; ++ iy) {
						ref_z_cut[iy] = ref_data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
					} // for

					delete[] ref_data;

					// this will store z cut values for each iteration for plotting later
					float_t* z_cuts = new (std::nothrow) float_t[nqy_ * max_iter * 13];

					// do some preprocessing
					// start the main loop, bound by max_iter and err_threshold
					// 	compute gisaxs for current parameter values
					// 	compute the neighbors parameter values
					// 	for 12 combinations of current and neighbors, compute gisaxs and error
					// 	compute the derivatives (gradient) and error stuff
					// 	update parameter values
					// compute the error surface

					float_t err = 10.0;
					std::vector<float_t> param1_list;
					std::vector<float_t> param2_list;
					structure_iterator_t structure_iter = HiGInput::instance().structure_begin();
					Structure* structure = &((*structure_iter).second);
					Shape* shape = HiGInput::instance().shape(*structure);
					shape_param_iterator_t shape_param = (*shape).param_begin();
					float_t* data = NULL;
					for(unsigned int iter = 0; iter < max_iter; ++ iter) {
						param1_list.clear();
						param1_list.push_back(param_vals[0] - 2 * param_deltas[0]);	// p1mm
						param1_list.push_back(param_vals[0] - param_deltas[0]);		// p1m
						param1_list.push_back(param_vals[0]);							// p1
						param1_list.push_back(param_vals[0] + param_deltas[0]);		// p1p
						param1_list.push_back(param_vals[0] + 2 * param_deltas[0]);	// p1pp
						param2_list.clear();
						param2_list.push_back(param_vals[1] - 2 * param_deltas[1]);	// p2mm
						param2_list.push_back(param_vals[1] - param_deltas[1]);		// p2m
						param2_list.push_back(param_vals[1]);							// p2
						param2_list.push_back(param_vals[1] + param_deltas[1]);		// p2p
						param2_list.push_back(param_vals[1] + 2 * param_deltas[1]);	// p2pp

						// current point
						(*shape_param).second.mean(param1_list[2]);
						(*shape_param).second.deviation(param2_list[2]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 0 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err22 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 0 * nqy_,
																ref_z_cut, deltay);

						// 12 neighbors
						(*shape_param).second.mean(param1_list[0]);
						(*shape_param).second.deviation(param2_list[2]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 1 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err02 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 1 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[1]);
						(*shape_param).second.deviation(param2_list[1]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 2 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err11 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 2 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[1]);
						(*shape_param).second.deviation(param2_list[2]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 3 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err12 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 3 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[1]);
						(*shape_param).second.deviation(param2_list[3]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 4 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err13 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 4 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[2]);
						(*shape_param).second.deviation(param2_list[0]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 5 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err20 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 5 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[2]);
						(*shape_param).second.deviation(param2_list[1]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 6 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err21 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 6 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[2]);
						(*shape_param).second.deviation(param2_list[3]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 7 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err23 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 7 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[2]);
						(*shape_param).second.deviation(param2_list[4]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 8 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err24 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 8 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[3]);
						(*shape_param).second.deviation(param2_list[1]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 9 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err31 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 9 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[3]);
						(*shape_param).second.deviation(param2_list[2]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 10 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err32 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 10 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[3]);
						(*shape_param).second.deviation(param2_list[3]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 11 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err33 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 11 * nqy_,
																ref_z_cut, deltay);

						(*shape_param).second.mean(param1_list[4]);
						(*shape_param).second.deviation(param2_list[2]);
						if(!run_gisaxs(alpha_i, alphai, phi, tilt, data, world_comm)) {
							if(mpi_rank == 0)
								std::cerr << "error: could not finish successfully" << std::endl;
							return false;
						} // if
						for(unsigned int iy = 0; iy < nqy_; ++ iy)
							z_cuts[13 * iter * nqy_ + 12 * nqy_ + iy] =
													data[nqx_ * nqy_ * zcut + nqx_ * iy + 0];
						delete[] data; data = NULL;
						float_t err42 = compute_cut_fit_error(z_cuts + 13 * iter * nqy_ + 12 * nqy_,
																ref_z_cut, deltay);

						// 22	0
						// 02	1mm
						// 11	1m2m
						// 12	1m
						// 13	1m2p
						// 20	2mm
						// 21	2m
						// 23	2p
						// 24	2pp
						// 31	1p2m
						// 32	1p
						// 33	1p2p
						// 42	1pp

						float_t derr1 = (err32 - err12) / (2 * param_deltas[0]);
						float_t derr2 = (err23 - err21) / (2 * param_deltas[1]);
						err = sqrt(derr1 * derr1 + derr2 * derr2);
						if(err < err_threshold) break;

						float_t herr11 = (err42 - err02 - 2 * err22) /
											(4 * param_deltas[0] * param_deltas[0]);
						float_t herr12 = (err33 - err13 - (err31 - err11)) /
											(4 * param_deltas[0] * param_deltas[1]);
						float_t herr21 = (err33 - err13 - (err31 - err11)) /
											(4 * param_deltas[0] * param_deltas[1]);
						float_t herr22 = (err24 - err20 - 2 * err22) /
											(4 * param_deltas[1] * param_deltas[1]);
						float_t* herr = new (std::nothrow) float_t[2 * 2];
						herr[0] = herr11;
						herr[1] = herr12;
						herr[2] = herr21;
						herr[3] = herr22;
						float_t* herrinv;
						mldivide(2, herr, herrinv);

						param_vals[0] =
								param_vals[0] - gamma_const * (herrinv[0] * derr1 + herrinv[1] * derr2);
						param_vals[1] =
								param_vals[1] - gamma_const * (herrinv[2] * derr1 + herrinv[3] * derr2);

						delete[] herrinv;
						delete[] herr;
					} // for

					(*shape_param).second.mean(22.0);
					(*shape_param).second.deviation(7.0);

					delete[] z_cuts;
					delete[] ref_z_cut;

					std::cout << "parameter values: " << param_vals[0] << ", " << param_vals[1]
								<< " (err: " << err << ")" << std::endl;

					// synchronize all procs after each run
					world_comm.Barrier();
				} // for tilt
			} // for phi
		} // for alphai

		return true;
	} // HipGISAXS::fit_all_gisaxs()


	float_t HipGISAXS::compute_cut_fit_error(float_t* cut, float_t* ref, float_t dy) {
		float_t err_sum = 0.0;
		for(unsigned int i = 0; i < nqy_; ++ i) {
			err_sum += fabs(cut[i] - ref[i]);
		} // for
		return err_sum * dy;
	} // HipGISAXS::compute_cut_fit_error()



} // namespace hig
