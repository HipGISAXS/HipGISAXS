/***
  *  Project:
  *
  *  File: fitting_steepest_descent.cpp
  *  Created: Dec 06, 2012
  *  Modified: Thu 06 Dec 2012 10:04:41 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include "hipgisaxs_main.hpp"


namespace hig {

	// currently only hacked for spheres, with radius and sd as two parameters
	bool HipGISAXS::fit_steepest_descent(float_t zcut,
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

		float_t alpha_i = alphai_min;
		// high level of parallelism here (alphai, phi, tilt) for dynamicity ...
		for(int i = 0; i < num_alphai; i ++, alpha_i += alphai_step) {
			float_t alphai = alpha_i * PI_ / 180;
			float_t phi = phi_min;
			for(int j = 0; j < num_phi; j ++, phi += phi_step) {
				float_t tilt = tilt_min;
				for(int k = 0; k < num_tilt; k ++, tilt += tilt_step) {

					std::cout << "-- Computing GISAXS "
								<< i * num_phi * num_tilt + j * num_tilt + k + 1 << " / "
								<< num_alphai * num_phi * num_tilt
								<< " [alphai = " << alpha_i << ", phi = " << phi
								<< ", tilt = " << tilt << "] ..." << std::endl;

					/* run a gisaxs simulation */

					float_t* final_data = NULL;
					if(!run_gisaxs(alpha_i, alphai, phi, tilt, final_data, world_comm)) {
						if(mpi_rank == 0) std::cerr << "error: could not finish successfully" << std::endl;
						return false;
					} // if

					if(mpi_rank == 0) {
						// note that final_data stores 3d info
						// for 2d, just take a slice of the data
						std::cout << "-- Constructing GISAXS image ... " << std::flush;
						Image img(nqx_, nqy_, nqz_);
						img.construct_image(final_data, 0); // merge this into the contructor ...
						std::cout << "done." << std::endl;

						if(x_max < x_min) x_max = x_min;
						// define output filename
						std::stringstream alphai_b, phi_b, tilt_b;
						std::string alphai_s, phi_s, tilt_s;
						alphai_b << alpha_i; alphai_s = alphai_b.str();
						phi_b << phi; phi_s = phi_b.str();
						tilt_b << tilt; tilt_s = tilt_b.str();
						std::string output(HiGInput::instance().param_pathprefix() +
											"/" + HiGInput::instance().runname() +
											"/img_ai=" + alphai_s + "_rot=" + phi_s +
											"_tilt=" + tilt_s + ".tif");

						std::cout << "**                    Image size: " << nqy_  << " x " << nqz_
									<< std::endl;
						std::cout << "-- Saving image in " << output << " ... " << std::flush;
						img.save(output);
						std::cout << "done." << std::endl;

						// save the actual data into a file also
						std::string data_file(HiGInput::instance().param_pathprefix() +
										"/" + HiGInput::instance().runname() +
										"/gisaxs_ai=" + alphai_s + "_rot=" + phi_s +
										"_tilt=" + tilt_s + ".out");
						std::cout << "-- Saving raw data in " << data_file << " ... "
								<< std::flush;
						save_gisaxs(final_data, data_file);
						std::cout << "done." << std::endl;
					} // if

					delete[] final_data;

					// synchronize all procs after each run
					world_comm.Barrier();
				} // for tilt
			} // for phi
		} // for alphai

		return true;
	} // HipGISAXS::fit_all_gisaxs()

} // namespace hig
