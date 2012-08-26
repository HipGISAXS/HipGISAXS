/***
  *  $Id: hipgisaxs_main.cpp 47 2012-08-23 21:05:16Z asarje $
  *
  *  Project: HipGISAXS
  *
  *  File: hipgisaxs_main.cpp
  *  Created: Jun 14, 2012
  *  Modified: Thu 23 Aug 2012 01:58:32 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
//#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <mpi.h>

#include "hipgisaxs_main.hpp"
#include "typedefs.hpp"
#include "utilities.hpp"

namespace hig {

	bool HipGISAXS::init(MPI::Intracomm& world_comm) {
						// is called at the beginning of the runs (after input is read)
		// first check if the input has been constructed ...

		int mpi_rank = world_comm.Get_rank();

		//photon conversion
		float_t photon = 0.0;
		std::string unit;
		freq_ = 0; k0_ = 0;
		HiGInput::instance().photon_energy(photon, unit);
		if(unit == "ev") {
		//	freq_ = photon * 1.60217646e9 / 6.626068;
			photon = photon / 1000;		// in keV
			freq_ = 1e-9 * photon * 1.60217646e-19 * 1000 / 6.626068e-34;
		} else { /* do something else ? ... */
			if(mpi_rank == 0) std::cerr << "error: photon energy is not given in 'ev'" << std::endl;
			return false;
		} // if-else

		k0_ = 2 * PI_ * freq_ / LIGHT_SPEED_;

		// create output directory
		if(mpi_rank == 0) {		// this is not quite good for mpi ... improve ...
			const std::string p = HiGInput::instance().path() + "/" + HiGInput::instance().runname();
			if(!boost::filesystem::create_directory(p)) {
				std::cerr << "error: could not create output directory " << p << std::endl;
				return false;
			} // if
		} // if

		world_comm.Barrier();

		// create Q-grid
		float_t min_alphai = HiGInput::instance().scattering_min_alpha_i() * PI_ / 180;
		if(!QGrid::instance().create(freq_, min_alphai, k0_)) {
			if(mpi_rank == 0) std::cerr << "error: could not create Q-grid" << std::endl;
			return false;
		} // if

		nqx_ = QGrid::instance().nqx();
		nqy_ = QGrid::instance().nqy();
		nqz_ = QGrid::instance().nqz();

		return true;
	} // HipGISAXS::init()


	bool HipGISAXS::run_init(float_t alphai, float_t phi, float_t tilt, MPI::Intracomm& world_comm) {
					// this is called for each config-run during the main run

		int mpi_rank = world_comm.Get_rank();
		// get all the variable values from the input structures	(can't we do some of this in init()?)

		/* get initialization data from layers */

		if(HiGInput::instance().has_substrate_layer()) {	//
			substrate_refindex_ = HiGInput::instance().substrate_layer().refindex();
		} else {
			substrate_refindex_.delta(0); substrate_refindex_.beta(0);
		} // if-else

		num_layers_ = HiGInput::instance().num_layers();	// this excludes substrate layer
		if(num_layers_ == 1) {		// is this really needed? ...
			single_layer_refindex_ = HiGInput::instance().single_layer().refindex();
			single_layer_thickness_ = HiGInput::instance().single_layer().thickness();
		} else {
			// initialize
			single_layer_refindex_.delta(0); single_layer_refindex_.beta(0);
			single_layer_thickness_ = 0.0;
		} // if-else

		complex_t temp(substrate_refindex_.delta(), substrate_refindex_.beta());
		dns2_ = ((float_t) 2.0) * temp - (complex_t) pow(temp, 2);

		/* get initialization data from structures */
		
		num_structures_ = HiGInput::instance().num_structures();
		/* construct lattice vectors in each structure, if needed */
		if(!HiGInput::instance().construct_lattice_vectors()) {	// this can also be done at input reading ...
			if(mpi_rank == 0) std::cerr << "error: could not construct lattice vectors" << std::endl;
			return false;
		} // if
		if(!illuminated_volume(alphai, HiGInput::instance().scattering_spot_area(),
				HiGInput::instance().min_layer_order(), HiGInput::instance().substrate_refindex())) {
			if(mpi_rank == 0) std::cerr << "error: something went wrong in illuminatedvolume()" << std::endl;
			return false;
		} // if
		if(!HiGInput::instance().construct_layer_profile()) {	// also can be done at input reading ...
			if(mpi_rank == 0) std::cerr << "error: could not construct layer profile" << std::endl;
			return false;
		} // if

		/* domain cell size */
		vector3_t min_vec(0.0, 0.0, 0.0), max_vec(0.0, 0.0, 0.0);
		float_t z_min_0 = 0.0, z_max_0 = 0.0;
		if(!HiGInput::instance().compute_domain_size(min_vec, max_vec, z_min_0, z_max_0)) {
			if(mpi_rank == 0) std::cerr << "error: could not construct domain sizes" << std::endl;
			return false;
		} // if
		// min_vec and max_vec are the domain
		cell_[0] = fabs(max_vec[0] - min_vec[0]);
		cell_[1] = fabs(max_vec[1] - min_vec[1]);
		cell_[2] = fabs(z_max_0 - z_min_0);

		/* rotation matrices */
		// compute r_phi = sample rotation by phi
		vector3_t rp1(0.0, 0.0, 0.0), rp2(0.0, 0.0, 0.0), rp3(0.0, 0.0, 0.0);
		vector3_t rt1(0.0, 0.0, 0.0), rt2(0.0, 0.0, 0.0), rt3(0.0, 0.0, 0.0);
		compute_rotation_matrix_z(phi, rp1, rp2, rp3);
		compute_rotation_matrix_x(tilt, rt1, rt2, rt3);
		mat_mul_3x3(rp1, rp2, rp3, rt1, rt2, rt3,
				rotation_matrix_.r1_, rotation_matrix_.r2_, rotation_matrix_.r3_);

		return true;
	} // HipGISAXS::run_init()


	bool HipGISAXS::run_all_gisaxs(MPI::Intracomm& world_comm,
									int x_min = 0, int x_max = 0, int x_step = 0) {
		int mpi_rank = world_comm.Get_rank();

		if(!init(world_comm)) return false;

		int num_alphai = 0, num_phi = 0, num_tilt = 0;;

		float_t alphai_min, alphai_max, alphai_step;
		HiGInput::instance().scattering_alphai(alphai_min, alphai_max, alphai_step);
		if(alphai_step == 0) num_alphai = 1;
		else num_alphai = (alphai_max - alphai_min) / alphai_step + 1;

		float_t phi_min, phi_max, phi_step;
		HiGInput::instance().scattering_inplanerot(phi_min, phi_max, phi_step);
		if(phi_step == 0) num_phi = 1;
		else num_phi = (phi_max - phi_min) / phi_step + 1;

		float_t tilt_min, tilt_max, tilt_step;
		HiGInput::instance().scattering_tilt(tilt_min, tilt_max, tilt_step);
		if(tilt_step == 0) num_tilt = 1;
		else num_tilt = (tilt_max - tilt_min) / tilt_step + 1;

		//std::cout << "num alphai: " << num_alphai << ", num phi: " << num_phi
		//			<< ", num tilt: " << num_tilt << std::endl;

		float_t alpha_i = alphai_min;
		for(int i = 0; i < num_alphai; i ++, alpha_i += alphai_step) {
			float_t alphai = alpha_i * PI_ / 180;
			float_t phi = phi_min;
			for(int j = 0; j < num_phi; j ++, phi += phi_step) {
				float_t tilt = tilt_min;
				for(int k = 0; k < num_tilt; k ++, tilt += tilt_step) {

					// run one gisaxs simulation
					float_t* final_data = NULL;
					if(!run_gisaxs(alpha_i, alphai, phi, tilt, final_data, world_comm)) {
						if(mpi_rank == 0) std::cerr << "error: could not finish successfully" << std::endl;
						return false;
					} // if

					// for future ...
					//Image img3d(nqx_, nqy_, nqz_);
					//if(!run_gisaxs(alphai, phi, tilt, img3d)) {
					
					//printfr("final_data", final_data, nqx_ * nqy_ * nqz_);

					if(mpi_rank == 0) {
						/*for(int i = 0; i < nqz_; ++ i) {
							for(int j = 0; j < nqy_; ++ j) {
								for(int k = 0; k < nqx_; ++ k) {
									std::cout << final_data[nqx_ * nqy_ * i + nqx_ * j + k] << "\t";
								} // for
								std::cout << std::endl;
							} // for
							std::cout << std::endl;
						} // for */

						//Image img(nqx_, nqy_, nqz_);
						// note that final_data stores 3d info
						// for 2d, just take a slice of the data
						std::cout << "-- Constructing image ..." << std::endl;
						Image img(nqx_, nqy_, nqz_);
						img.construct_image(final_data, 0); // merge this into the contructor ...

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

						//img.save(output, x_min, x_max);
						std::cout << "-- Saving image in " << output << " ..." << std::endl;
						img.save(output);

						// for future ...
						/*for(int x = x_min; x <= x_max; x += x_step) {
							Image *img2d = NULL;
							img3d.slice(x, img2d);
	
							// define output names ...
							std::ostream temp;
							std::string x_s, alphai_s, phi_s, tilt_s;
							temp << std::setw(4) << x; temp >> x_s;
							temp << alphai; temp >> alphai_s;
							temp << phi; temp >> phi_s;
							temp << tilt; temp >> tilt_s;
							std::string output(HiGInput::instance().outputdir() +
											"/img" + x_s + "_ai=" + alphai_s + "_rot=" + phi_s +
											"_tilt=" + tilt_s + ".img");
							img2d.save(output);
								// save in buffer too ... ? not for now
						 //} for x */
					} // if

					delete[] final_data;

					// synchronize all procs after each run
					world_comm.Barrier();
				} // for tilt
			} // for phi
		} // for alphai

		return true;
	} // HipGISAXS::run_all_gisaxs()
	
	
	void HipGISAXS::printfr(const char* name, float_t* arr, unsigned int size) {
		std::cout << name << ":" << std::endl;
		if(arr == NULL) { std::cout << "NULL" << std::endl; return; }
		for(unsigned int i = 0; i < size; ++ i) {
			std::cout << arr[i] << "\t";
		} // for
		std::cout << std::endl;
	} // HipGISAXS::printfc()
 

	void HipGISAXS::printfc(const char* name, complex_t* arr, unsigned int size) {
		std::cout << name << ":" << std::endl;
		if(arr == NULL) { std::cout << "NULL" << std::endl; return; }
		for(unsigned int i = 0; i < size; ++ i) {
			std::cout << arr[i].real() << "," << arr[i].imag() << "\t";
		} // for
		std::cout << std::endl;
	} // HipGISAXS::printfc()
 

	/* all the real juice is here */
	bool HipGISAXS::run_gisaxs(float_t alpha_i, float_t alphai, float_t phi, float_t tilt, float_t* &img3d,
								MPI::Intracomm& world_comm, int corr_doms) {

		if(!run_init(alphai, phi, tilt, world_comm)) return false;

		//std::cout << "nqx: " << nqx_ << ", nqy: " << nqy_ << ", nqz: "
		//		<< nqz_ << ", nqz_e: " << nqz_extended_ << std::endl;

		int mpi_rank = world_comm.Get_rank();

		// compute propagation coefficients/fresnel coefficients
		// this can also go into run_init() ...
		complex_t *amm = NULL, *apm = NULL, *amp = NULL, *app = NULL;
		complex_t *rk1 = NULL, *rk2 = NULL, *rk1rk2 = NULL, *tk1tk2 = NULL, *h0 = NULL;
		if(!compute_propagation_coefficients(alphai, amm, apm, amp, app,
											rk1, rk2, rk1rk2, tk1tk2, h0)) return false;
		//printfc("amm", amm, nqx_ * nqy_ * nqz_);
		//printfc("apm", apm, nqx_ * nqy_ * nqz_);
		//printfc("amp", amp, nqx_ * nqy_ * nqz_);
		//printfc("app", app, nqx_ * nqy_ * nqz_);
		//printfc("rk1", rk1, nqx_ * nqy_ * nqz_);
		//printfc("rk2", rk2, nqx_ * nqy_ * nqz_);
		//printfc("rk1rk2", rk1rk2, nqx_ * nqy_ * nqz_);
		//printfc("tk1tk2", tk1tk2, nqx_ * nqy_ * nqz_);
		//printfc("h0", h0, nqx_ * nqy_ * nqz_);

		// initialize memory for struct_intensity, ff and sf
		float_t* struct_intensity = new (std::nothrow) float_t[num_structures_ * nqx_ * nqy_ * nqz_];
		//complex_t *ff, *sf;
		// memory usage can be reduced here ...
		//if(HiGInput::instance().experiment() == "saxs") {
		//	ff = new (std::nothrow) complex_t[num_structures_ * nqx_ * nqy_ * nqz_];
		//	sf = new (std::nothrow) complex_t[num_structures_ * nqx_ * nqy_ * nqz_];
		//} else if(HiGInput::instance().experiment() == "gisaxs") {
		//	ff = new (std::nothrow) complex_t[num_structures_ * nqx_ * nqy_ * nqz_extended_];
		//	sf = new (std::nothrow) complex_t[num_structures_ * nqx_ * nqy_ * nqz_extended_];
		//} else {
		//	if(mpi_rank == 0) std::cerr << "error: experiment type '"
		//								<< HiGInput::instance().experiment()
		//								<< "' is either unknown or has not been implemented."
		//								<< std::endl;
		//	return false;
		//} // if-else

		/* loop over all structures and domains/distributions */
		int s_num = 0;
		for(structure_iterator_t s = HiGInput::instance().structure_begin();
				s != HiGInput::instance().structure_end(); ++ s, ++ s_num) {
			// get all repetitions of the structure in the volume
			// with position dr and orientation (tau, eta)
			// i.e. compute matrix that defines all num_domains domains inside the volume
			// distr = [ drx_i dry_i drz_i  tau_i eta_i ] (NDISTR x 5 )
			// get params of the shape repeated in this structure: [dn2, id, dims, invar, t, seta, stau]
			// allocate id(nqy, nqz, num_domains)

			// get the shape
			// compute t, lattice, ndoms, dd, nn, id etc.

			//std::string shape_k = (*s).grain_shape_key();
			//Shape *curr_shape = shape(shape_k);
			//Lattice *curr_lattice = (*s).lattice();

			Structure *curr_struct = &((*s).second);
			Shape *curr_shape = HiGInput::instance().shape(*curr_struct);
			Lattice *curr_lattice = (Lattice*) HiGInput::instance().lattice(*curr_struct);

			vector3_t grain_repeats = (*s).second.grain_repetition();

			int num_domains = 0, num_coords = 0;
			int num_nn = 0;
			float_t *dd = NULL, *nn = NULL;		// come back to this ...
												// these structures can be improved ...
			float_t tz = 0;
			int num_dimen = 3;
			int ndx = 0, ndy = 0;
			// compute dd and nn
			spatial_distribution(s, tz, num_dimen, ndx, ndy, dd);
			orientation_distribution(s, dd, ndx, ndy, nn); // num_nn = dim(nn.x)

			num_nn = ndx;
			num_domains = num_nn;

			complex_t *id = NULL;		// check this ... come back and improve ...
										// this may be reduced on the fly to reduce mem usage ...
			id = new (std::nothrow) complex_t[num_domains * nqx_ * nqy_ * nqz_];
			if(id == NULL) {
				if(mpi_rank == 0) std::cerr << "error: could not allocate memory for 'id'" << std::endl;
				return false;
			} // if

			vector3_t curr_transvec = (*s).second.grain_transvec();
			curr_transvec = curr_transvec - vector3_t(0, 0, single_layer_thickness_);
			ShapeName shape_name = HiGInput::instance().shape_name((*s).second);
			float_t shape_tau = HiGInput::instance().shape_tau((*s).second);
			float_t shape_eta = HiGInput::instance().shape_eta((*s).second);
			vector3_t shape_originvec = HiGInput::instance().shape_originvec((*s).second);
			std::string shape_file = HiGInput::instance().shape_filename((*s).second);
			shape_param_list_t shape_params = HiGInput::instance().shape_params((*s).second);

			/* computing dwba ff for each domain in structure (*s) with
			 * ensemble containing num_domains grains */
			for(int j = 0; j < num_domains; j ++) {	// or distributions
				// define r_norm (domain orientation by tau and eta)
				// define full domain rotation matrix r_total = r_phi * r_norm
				// ... i think these tau eta zeta can be computed on the fly to save memory ...
				float_t tau = nn[0 * num_nn + j];
				float_t eta = nn[1 * num_nn + j];
				float_t zeta = nn[2 * num_nn + j];
				vector3_t z1, z2, z3, e1, e2, e3, t1, t2, t3;
				compute_rotation_matrix_z(zeta, z1, z2, z3);
				compute_rotation_matrix_y(eta, e1, e2, e3);
				compute_rotation_matrix_x(tau, t1, t2, t3);

				vector3_t temp1, temp2, temp3;
				vector3_t r_norm1, r_norm2, r_norm3;
				mat_mul_3x3(z1, z2, z3, e1, e2, e3, temp1, temp2, temp3);
				mat_mul_3x3(temp1, temp2, temp3, t1, t2, t3, r_norm1, r_norm2, r_norm3);

				vector3_t r_tot1, r_tot2, r_tot3;
				mat_mul_3x3(rotation_matrix_.r1_, rotation_matrix_.r2_, rotation_matrix_.r3_,
							r_norm1, r_norm2, r_norm3, r_tot1, r_tot2, r_tot3);

				/* center of unit cell replica */
				vector3_t curr_dd_vec(dd[3 * j + 0], dd[3 * j + 1], dd[3 * j + 2]);
				vector3_t result(0.0, 0.0, 0.0);
				mat_mul_3x1(rotation_matrix_.r1_, rotation_matrix_.r2_, rotation_matrix_.r3_,
							curr_dd_vec, result);
				vector3_t center = result + curr_transvec;

				//std::cout << "nqx: " << nqx_ << ", nqy: " << nqy_ << ", nqz: "
				//		<< nqz_ << ", nqz_e: " << nqz_extended_ << std::endl;

				/* compute the structure factor and form factor */
//				if(HiGInput::instance().experiment() == "saxs") {
//					structure_factor(HiGInput::instance().experiment(), center,
//								curr_lattice, grain_repeats, r_tot1, r_tot2, r_tot3, world_comm);
//					form_factor(shape_name, shape_file, shape_params, curr_transvec,
//								shape_tau, shape_eta, r_tot1, r_tot2, r_tot3, world_comm);
//				} else if(HiGInput::instance().experiment() == "gisaxs") {
				structure_factor(HiGInput::instance().experiment(), center,
							curr_lattice, grain_repeats, r_tot1, r_tot2, r_tot3, world_comm);
				//sf_.printsf();

				// write q grid
				//write_qgrid("current_qgrid");
				//exit(0);

				//read_form_factor("curr_ff.out");
				form_factor(shape_name, shape_file, shape_params, curr_transvec,
							shape_tau, shape_eta, r_tot1, r_tot2, r_tot3, world_comm);
				//ff_.print_ff(nqx_, nqy_, nqz_extended_);
				//ff_.printff(nqx_, nqy_, nqz_extended_);
				std::stringstream alphai_b, phi_b, tilt_b;
				std::string alphai_s, phi_s, tilt_s;
				alphai_b << alpha_i; alphai_s = alphai_b.str();
				phi_b << phi; phi_s = phi_b.str();
				tilt_b << tilt; tilt_s = tilt_b.str();
				std::string ff_output(HiGInput::instance().param_pathprefix() +
								"/" + HiGInput::instance().runname() +
								"/ff_ai=" + alphai_s + "_rot=" + phi_s +
								"_tilt=" + tilt_s + ".out");
				std::cout << "-- Saving form factor in " << ff_output << " ..." << std::endl;
				ff_.save_ff(nqx_, nqy_, nqz_extended_, ff_output.c_str());

//				} else {
//					if(mpi_rank == 0)
//						std::cerr << "error: experiment type '" << HiGInput::instance().experiment()
//									<< "' is either unknown or has not been implemented." << std::endl;
//					return false;
//				} // if-else

				/*std::cout << "SF: nqx = " << nqx_ << ", nqy = " << nqy_ << ", nqz = " << nqz_ << std::endl;
				for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i)
					std::cout << sf_[i].real() << "+" << sf_[i].imag() << "i ";
				std::cout << std::endl;

				std::cout << "FF: nqx = " << nqx_ << ", nqy = " << nqy_ << ", nqz = " << nqz_ << std::endl;
				for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i)
					std::cout << ff_[i].x << "+" << ff_[i].y << "i ";
				std::cout << std::endl;*/

				complex_t* base_id = id + j * nqx_ * nqy_ * nqz_;

				unsigned int nslices = HiGInput::instance().param_nslices();
				if(nslices <= 1) {
					/* without slicing */
					if(HiGInput::instance().structure_layer_order((*s).second) == 1) {
													// structure is inside a layer, on top of substrate
						if(single_layer_refindex_.delta() < 0 || single_layer_refindex_.beta() < 0) {
							// this should never happen
							if(mpi_rank == 0)
								std::cerr << "error: single layer information not correctly set" << std::endl;
							return false;
						} // if
						complex_t dn2(-2.0 * ((*s).second.grain_refindex().delta() -
										single_layer_refindex_.delta()),
										-2.0 * ((*s).second.grain_refindex().beta() -
										single_layer_refindex_.beta()));

						//std::cout << "DN2 = " << dn2 << std::endl;

						// base_id = dn2 * (amm .* sf() .* ff() + amp .* sf() .* ff() +
						//				apm .* sf() .* ff() + app .* sf() .* ff());
						for(int z = 0; z < nqz_; ++ z) {
							for(int y = 0; y < nqy_; ++ y) {
								for(int x = 0; x < nqx_; ++ x) {
									unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
									unsigned int curr_index_0 = curr_index;
									unsigned int curr_index_1 = nqx_ * nqy_ * nqz_ + curr_index;
									unsigned int curr_index_2 = 2 * nqx_ * nqy_ * nqz_ + curr_index;
									unsigned int curr_index_3 = 3 * nqx_ * nqy_ * nqz_ + curr_index;

									// what happends in case of "saxs" ? ... when nqz_ext == nqz
									base_id[curr_index] = dn2 *
										(amm[curr_index] * sf_[curr_index_0] * ff_[curr_index_0] +
										amp[curr_index] * sf_[curr_index_1] * ff_[curr_index_1] +
										apm[curr_index] * sf_[curr_index_2] * ff_[curr_index_2] +
										app[curr_index] * sf_[curr_index_3] * ff_[curr_index_3]);
								} // for x
							} // for y
						} // for z
					} else if(HiGInput::instance().structure_layer_order((*s).second) == -1) {
														// structure burried in substrate
						complex_t dn2(-2.0 * (substrate_refindex_.delta() -
										(*s).second.grain_refindex().delta()),
										-2.0 * (substrate_refindex_.beta() -
										(*s).second.grain_refindex().beta()));

						// what happens in the case when sf and ff have nqz_extended = 4 * nqz  ... ?
						for(int z = 0; z < nqz_; ++ z) {
							for(int y = 0; y < nqy_; ++ y) {
								for(int x = 0; x < nqx_; ++ x) {
									unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
									base_id[curr_index] = dn2 * sf_[curr_index] * ff_[curr_index];
								} // for x
							} // for y
						} // for z
					} else if(HiGInput::instance().structure_layer_order((*s).second) == 0) {
														// structure on top of substrate
						complex_t dn2(-2.0 * (*s).second.grain_refindex().delta(),
										-2.0 * (*s).second.grain_refindex().beta());

						for(int z = 0; z < nqz_; ++ z) {
							for(int y = 0; y < nqy_; ++ y) {
								for(int x = 0; x < nqx_; ++ x) {
									unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
									unsigned int curr_index_0 = curr_index;
									unsigned int curr_index_1 = nqx_ * nqy_ * nqz_ + curr_index;
									unsigned int curr_index_2 = 2 * nqx_ * nqy_ * nqz_ + curr_index;
									unsigned int curr_index_3 = 3 * nqx_ * nqy_ * nqz_ + curr_index;
									base_id[curr_index] = dn2 *
										(h0[curr_index] * sf_[curr_index_0] * ff_[curr_index_0] +
										rk2[curr_index] * sf_[curr_index_1] * ff_[curr_index_1] +
										rk1[curr_index] * sf_[curr_index_2] * ff_[curr_index_2] +
										rk1rk2[curr_index] * sf_[curr_index_3] * ff_[curr_index_3]);
								} // for x
							} // for y
						} // for z
					} else {
						if(mpi_rank == 0)
							std::cerr << "error: unable to determine sample structure. "
										<< "make sure the layer order is correct" << std::endl;
						return false;
					} // if-else
				} else {
					/* perform slicing */
					// not yet implemented ...
					if(mpi_rank == 0)
						std::cout << "uh-oh: ever thought about implementing the slicing scheme?"
									<< std::endl;
					return false;
				} // if-else
			} // for num_domains

			delete[] nn;
			delete[] dd;

			/*std::cout << "ID: nqx = " << nqx_ << ", nqy = " << nqy_ << ", nqz = " << nqz_ << std::endl;
			for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
				std::cout << id[i] << " ";
			} // for
			std::cout << std::endl; */

			//printfc("id", id, nqx_ * nqy_ * nqz_);

			// in slim's code, this is inside num_domains loop ...
			// i think it should be outside because it is summing
			// over the 4th dimension which is the set of domains ...
			// nqz_extended ... ?? ...
			if(corr_doms == 1) {		// note: currently this is hardcoded as 0
				unsigned int soffset = s_num * nqx_ * nqy_ * nqz_;
				for(int z = 0; z < nqz_; ++ z) {
					for(int y = 0; y < nqy_; ++ y) {
						for(int x = 0; x < nqx_; ++ x) {
							unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
							complex_t sum(0.0, 0.0);
							for(int d = 0; d < num_domains; ++ d) {
								unsigned int id_index = d * nqx_ * nqy_ * nqz_ + curr_index;
								sum += id[id_index];
							} // for d
							struct_intensity[soffset + curr_index] = sum.real() * sum.real() +
																		sum.imag() * sum.imag();
						} // for x
					} // for y
				} // for z
			} else {
				unsigned int soffset = s_num * nqx_ * nqy_ * nqz_;
				for(unsigned int z = 0; z < nqz_; ++ z) {
					for(unsigned int y = 0; y < nqy_; ++ y) {
						for(unsigned int x = 0; x < nqx_; ++ x) {
							unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
							float_t sum = 0.0;
							for(int d = 0; d < num_domains; ++ d) {
								unsigned int id_index = d * nqx_ * nqy_ * nqz_ + curr_index;
								sum += id[id_index].real() * id[id_index].real() +
										id[id_index].imag() * id[id_index].imag();
							} // for d
							struct_intensity[soffset + curr_index] = sum;
						} // for x
					} // for y
				} // for z
			} // if-elseo

			//printfr("struct_intensity", struct_intensity, nqx_ * nqy_ * nqz_);

			delete[] id;
		} // for num_structs
		int num_structs = s_num;

		// arrays struct_intensity and id etc can be eliminated/reduced in size to nqx*nqy*nqz only
		// ...

/*		std::cout << "STRUCT_INTENSITY: num_structs = " << num_structs
				<< ", nqx = " << nqx_ << ", nqy = " << nqy_ << ", nqz = " << nqz_ << std::endl;
		for(int i = 0; i < num_structs * nqx_ * nqy_ * nqz_; ++ i) {
			std::cout << struct_intensity[i] << " ";
		} // for
		std::cout << std::endl;
*/
		img3d = new (std::nothrow) float_t[nqx_ * nqy_ * nqz_];
		// sum of struct_intensity into intensity
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				for(unsigned int x = 0; x < nqx_; ++ x) {
					unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
					float_t sum = 0.0;
					for(int s = 0; s < num_structs; ++ s) {
						unsigned int index = s * nqx_ * nqy_ * nqz_ + curr_index;
						sum += struct_intensity[index];
					} // for d
					img3d[curr_index] = sum;
				} // for x
			} // for y
		} // for z

/*		std::cout << "IMG3D: nqx = " << nqx_ << ", nqy = " << nqy_ << ", nqz = " << nqz_ << std::endl;
		for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
			std::cout << img3d[i] << " ";
		} // for
		std::cout << std::endl;
*/
		delete[] struct_intensity;
		delete[] fc_;

		return true;
	} // HipGISAXS::run_gisaxs()


	// merge the two ... ?
	//
	bool HipGISAXS::structure_factor(std::string expt, vector3_t& center, Lattice* &curr_lattice,
									vector3_t& grain_repeats, vector3_t& r_tot1,
									vector3_t& r_tot2, vector3_t& r_tot3,
									MPI::Intracomm& world_comm) {
		return sf_.compute_structure_factor(expt, center, curr_lattice, grain_repeats,
											r_tot1, r_tot2, r_tot3, world_comm);
	} // HipGISAXS::structure_factor()


	//bool HipGISAXS::structure_factor(float_t*& qz_extended, vector3_t& center, Lattice* &curr_lattice,
	//								vector3_t& grain_repeats, vector3_t& r_tot1,
	//								vector3_t& r_tot2, vector3_t& r_tot3) {
	//	return sf_.compute_structure_factor(center, curr_lattice, grain_repeats,
	//										r_tot1, r_tot2, r_tot3);
	//} // HipGISAXS::structure_factor()


	//template <typename float_t, typename complex_t>
	bool HipGISAXS::form_factor(ShapeName shape_name, std::string shape_file,
								shape_param_list_t& shape_params, vector3_t &curr_transvec,
								float_t shp_tau, float_t shp_eta,
								vector3_t &r_tot1, vector3_t &r_tot2, vector3_t &r_tot3,
								MPI::Intracomm& world_comm) {
		return ff_.compute_form_factor(shape_name, shape_file, shape_params, single_layer_thickness_,
										curr_transvec, shp_tau, shp_eta, r_tot1, r_tot2, r_tot3,
										world_comm);
	} // HipGISAXS::form_factor()


	//bool HipGISAXS::form_factor(float_t* &qz_extended, ShapeName shape_name, std::string shape_file,
	//							const shape_param_list_t& shape_params, const vector3_t &originvec,
	//							const vector3_t &curr_transvec, float_t shp_tau, float_t shp_eta,
	//							const vector3_t &r_tot1, const vector3_t &r_tot2, const vector3_t &r_tot3) {
	//	return ff_.compute_form_factor(shape_name, shape_file, shape_params, single_layer_thickness_,
	//									curr_transvec, shp_tau, shp_eta, r_tot1, r_tot2, r_tot3);
	//} // HipGISAXS::form_factor()


	bool HipGISAXS::compute_rotation_matrix_z(float_t angle, vector3_t& r1, vector3_t& r2, vector3_t& r3) {
		float_t s = sin(angle);
		float_t c = cos(angle);

		r1[0] = c; r1[1] = -s; r1[2] = 0.0;
		r2[0] = s; r2[1] = c; r2[2]  = 0.0;
		r3[0] = 0.0; r3[1] = 0.0; r3[2]  = 1.0;

		return true;
	} // HipGISAXS::comptue_rotation_matrix_z()


	bool HipGISAXS::compute_rotation_matrix_y(float_t angle, vector3_t& r1, vector3_t& r2, vector3_t& r3) {
		float_t s = sin(angle);
		float_t c = cos(angle);

		r1[0] = c; r1[1] = 0.0; r1[2] = -s;
		r2[0] = 0.0; r2[1] = 1.0; r2[2]  = 0.0;
		r3[0] = s; r3[1] = 0.0; r3[2]  = c;

		return true;
	} // HipGISAXS::compute_rotation_matrix_x()


	bool HipGISAXS::compute_rotation_matrix_x(float_t angle, vector3_t& r1, vector3_t& r2, vector3_t& r3) {
		float_t s = sin(angle);
		float_t c = cos(angle);

		r1[0] = 1.0; r1[1] = 0.0; r1[2] = 0.0;
		r2[0] = 0.0; r2[1] = c; r2[2]  = -s;
		r3[0] = 0.0; r3[1] = s; r3[2]  = c;

		return true;
	} // HipGISAXS::compute_rotation_matrix_x()


	bool HipGISAXS::illuminated_volume(float_t alpha_i, float_t spot_area, int min_layer_order,
										RefractiveIndex substrate_refindex) {
		float_t spot_diameter = 2.0 * sqrt(spot_area / PI_) * 1e6;	// in nm
		float_t substr_delta = substrate_refindex.delta();
		float_t substr_beta = substrate_refindex.beta();

		if(HiGInput::instance().experiment() == "saxs") {	// SAXS
			vol_[0] = vol_[1] = vol_[2] = spot_diameter;
		} else if(HiGInput::instance().experiment() == "gisaxs") {		// GISAXS
			complex_t c_max_depth = complex_t(MAX_DEPTH_, 0);
			complex_t penetration_depth_layer;
			if(HiGInput::instance().num_layers() == 1 && min_layer_order > 0) {
				RefractiveIndex r1 = HiGInput::instance().single_layer().refindex();
				float_t alpha_c = sqrt(2.0 * r1.delta());
				penetration_depth_layer = -1.0 /
											(2.0 * k0_ * (sqrt(alpha_i * alpha_i - alpha_c * alpha_c -
											complex_t(0, 2.0 * r1.beta()))).imag());
			} else if(HiGInput::instance().num_layers() == 0) {
				float_t alpha_c = sqrt(2.0 * substr_delta);
				penetration_depth_layer = -1.0 /
											(2.0 * k0_ * (sqrt(alpha_i * alpha_i - alpha_c * alpha_c -
											complex_t(0, 2.0 * substr_beta))).imag());
			} else {
				// the sample is described by 2 or more layers, slicing scheme will be applied
				// NOTE: the case where a structure is implicitly in 2 different layers is
				// not currently handled.
				std::cerr << "uh-oh: this case (num_layers > 1) has not yet been implemented yet"
							<< std::endl;
				std::cerr << "go get yourself a nice cup of yummmy hot chocolate instead!" << std::endl;
				return false;
			} // if-else
			vol_[0] = vol_[1] = spot_diameter;
			vol_[2] = std::real((penetration_depth_layer < c_max_depth) ? penetration_depth_layer : c_max_depth);
		} else {
			std::cerr << "error: experiment type '" << HiGInput::instance().experiment()
						<< "' is either unknown or has not been implemented." << std::endl;
			return false;
		} // if-else

		return true;
	} // HipGISAXS::illuminated_volume()


	bool HipGISAXS::compute_propagation_coefficients(float_t alpha_i,
					complex_t* &amm, complex_t* &apm, complex_t* &amp, complex_t* &app,
					complex_t* &rk1, complex_t* &rk2, complex_t* &rk1rk2, complex_t* &tk1tk2,
					complex_t* &h0) {
		amm = apm = amp = app = NULL;
		rk1 = rk2 = rk1rk2 = tk1tk2 = NULL;
		h0 = NULL;

		// doesnt this also depend on the type of polarization of the light? ... ?

		if(HiGInput::instance().param_nslices() <= 1) {	/* computing without sample slicing */
				// what is "single layer" when there are > 1 layers ? ... ?
			complex_t dnl2 = 2.0 * complex_t(single_layer_refindex_.delta(),
											single_layer_refindex_.beta());
			if(!layer_qgrid_qz(alpha_i, dnl2)) {
				std::cerr << "error: could not compute extended qz" << std::endl;
				return false;
			} // if

			if(HiGInput::instance().num_layers() == 1) {
				/* compute fresnel coefficients for a structure
				 * embedded inside a layer on a substrate */
				if(!compute_fresnel_coefficients_embedded(alpha_i)) {
					std::cerr << "error: could not compute fresnel coefficients" << std::endl;
					return false;
				} // if

				// set aliases
				amm = fc_;
				apm = fc_ + nqx_ * nqy_ * nqz_;
				amp = fc_ + 2 * nqx_ * nqy_ * nqz_;
				app = fc_ + 3 * nqx_ * nqy_ * nqz_;
				h0 = fc_ + 4 * nqx_ * nqy_ * nqz_;
			} else if(HiGInput::instance().num_layers() == 0) {
				/* compute fresnel coefficients for a structure
				 * on top of or buried inside a substrate */
				if(!compute_fresnel_coefficients_top_buried(alpha_i)) {
					std::cerr << "error: could not compute fresnel coefficients" << std::endl;
					return false;
				} // if
				// set aliases
				rk1 = fc_;							// re-check these ...
				rk2 = fc_ + nqx_ * nqy_ * nqz_;
				rk1rk2 = fc_ + 2 * nqx_ * nqy_ * nqz_;
				tk1tk2 = fc_ + 3 * nqx_ * nqy_ * nqz_;
				h0 = fc_ + 4 * nqx_ * nqy_ * nqz_;
			} else {
				std::cerr << "error: invalid number of layers in non-slicing scheme" << std::endl;
				return false;
			} // if-else
		} else {
			/* computing using the sample slicing scheme */

			// not yet implemented ...
			std::cout << "uh-oh: you hit a slicing part that has not been implemented yet" << std::endl;
			return false;
		} // if

/*		std::cout << "FC: nqx = " << nqx_ << ", nqy = " << nqy_
				<< ", nqz extended = " << nqz_extended_ << std::endl;
		for(int ex = 0; ex < 5; ++ ex) {
			for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
				std::cout << fc_[ex * nqx_ * nqy_ * nqz_ + i] << " ";
			} // for
			std::cout << std::endl;
		} // for
		std::cout << std::endl;
*/
		//exit(1);

		return true;
	} // HipGISAXS::compute_propagation_coefficients()


	bool HipGISAXS::layer_qgrid_qz(float_t alpha_i, complex_t dnl_j) {
		float_t kzi_0 = -1.0 * k0_ * sin(alpha_i);
		complex_t kzi_j = -1.0 * k0_ * sqrt((float_t)pow(sin(alpha_i), 2) - dnl_j);

		if(!QGrid::instance().create_qz_extended(k0_, kzi_0, kzi_j, dnl_j)) {
			std::cerr << "error: something went wrong while creating qz_extended" << std::endl;
			return false;
		} // if
		nqz_extended_ = QGrid::instance().nqz_extended();

		return true;
	} // HipGISAXS::layer_qgrid_qz()


	// optimize this later ...
	bool HipGISAXS::compute_fresnel_coefficients_embedded(float_t alpha_i) {
		RefractiveIndex nl = single_layer_refindex_;
		RefractiveIndex ns = substrate_refindex_;
		float_t lt = single_layer_thickness_;

		complex_t dnl2(2.0 * nl.delta(), 2.0 * nl.beta());
		complex_t dns2(2.0 * ns.delta(), 2.0 * ns.beta());

		float_t sinai = sin(alpha_i);
		float_t kiz0 = -1.0 * k0_ * sinai;
		complex_t kiz1 = -1.0 * k0_ * sqrt(sinai * sinai - dnl2);
		complex_t kiz2 = -1.0 * k0_ * sqrt(sinai * sinai - dns2);

		complex_t r01_kiz1 = (kiz0 - kiz1) / (kiz0 + kiz1);
		complex_t r12_kiz1 = (kiz1 - kiz2) / (kiz1 + kiz2);
		complex_t t01_kiz1 = 2.0 * (kiz0 / (kiz0 + kiz1));

		complex_t a1m_kiz1 = t01_kiz1 /
							((float_t) 1.0 + r01_kiz1 * r12_kiz1 * exp(complex_t(0, 2) * kiz1 * lt));
		complex_t a1p_kiz1 = a1m_kiz1 * r12_kiz1 * exp(complex_t(0, 2) * kiz1 * lt);

		complex_t *a1mi = NULL, *a1pi = NULL;
		a1mi = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];
		a1pi = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];
		if(a1mi == NULL || a1pi == NULL) {
			std::cerr << "error: failed to allocate memory for a1mi, a1pi" << std::endl;
			return false;
		} // if
//		std::cout << "A1MI+A1PI: " << std::endl;
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
			a1mi[i] = a1m_kiz1; a1pi[i] = a1p_kiz1;
//			std::cout << a1mi[i] << "+" << a1pi[i] << " ";
		} // for
//		std::cout << std::endl;

		complex_t *a1mf = NULL, *a1pf = NULL;
		a1mf = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];
		a1pf = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];
		if(a1mf == NULL || a1pf == NULL) {
			std::cerr << "error: failed to allocate memory for a1mf, a1pf" << std::endl;
			return false;
		} // if
		for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
			a1mf[i] = a1pf[i] = complex_t(0.0, 0.0);
		} // for

		// allocate fc memory
		fc_ = new (std::nothrow) complex_t[5 * nqx_ * nqy_ * nqz_];	// 5 sets
		if(fc_ == NULL) {
			std::cerr << "error: failed to allocate memory for fc" << std::endl;
			return false;
		} // if

		// test
/*		std::cout << "NQZ_EXTENDED: " << nqz_extended_ << std::endl;
		for(int i = 0; i < nqz_; ++ i) {
			std::cout << QGrid::instance().qz(i) << " ";
		} // for
		std::cout << std::endl;
*/
		//for(int z = nqz_ - 1; z >= 0; z --) {
		for(int z = 0; z < nqz_; ++ z) {
			complex_t a1m_nkfz1, a1p_nkfz1;
			float_t kfz0 = QGrid::instance().qz(z) + kiz0;

			if(kfz0 < 0) {
				a1m_nkfz1 = complex_t(0.0, 0.0);
				a1p_nkfz1 = complex_t(0.0, 0.0);

				for(int y = 0; y < nqy_; y ++) {
					for(int x = 0; x < nqx_; x ++) {
						fc_[4 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(0.0, 0.0);
					} // for x
				} // for y
			} else {	// kfz0 >= 0
				float_t nkfz0 = -1.0 * kfz0;
				complex_t nkfz1 = -1.0 * sqrt(kfz0 * kfz0 - k0_ * k0_ * dnl2);
				complex_t nkfz2 = -1.0 * sqrt(kfz0 * kfz0 - k0_ * k0_ * dns2);

				complex_t r01_nkfz1 = (nkfz0 - nkfz1) / (nkfz0 + nkfz1);
				complex_t r12_nkfz1 = (nkfz1 - nkfz2) / (nkfz1 + nkfz2);
				complex_t t01_nkfz1 = 2.0 * (nkfz0 / (nkfz0 + nkfz1));

				complex_t uniti = complex_t(0, 1);
				complex_t temp0 = 2.0 * nkfz1 * lt;
				float_t temp1 = exp(-1.0 * temp0.imag());
				complex_t temp2 = exp(uniti * temp0.real());
									//((float_t)cos(temp0.real()) + uniti * sin(temp0.real()));
				complex_t temp = temp1 * temp2;
				a1m_nkfz1 = t01_nkfz1 /
							((float_t) 1.0 + r01_nkfz1 * r12_nkfz1 * temp);
				a1p_nkfz1 = a1m_nkfz1 * r12_nkfz1 * temp;

				for(int y = 0; y < nqy_; y ++) {
					for(int x = 0; x < nqx_; x ++) {
						fc_[4 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(1.0, 0.0);
					} // for x
				} // for y
			} // if-else

			for(int y = 0; y < nqy_; y ++) {							// these can be aliminated ...
				for(int x = 0; x < nqx_; x ++) {
					a1mf[z * nqx_ * nqy_ + y * nqx_ + x] = a1m_nkfz1;
					a1pf[z * nqx_ * nqy_ + y * nqx_ + x] = a1p_nkfz1;
				} // for x
			} // for y
		} // for z

/*		std::cout << "A1MF+A1PF: " << std::endl;
		for(int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
			std::cout << a1mf[i] << "+" << a1pf[i] << " ";
		} // for
		std::cout << std::endl;
*/
		// the element-element products
		for(int z = 0; z < nqz_; z ++) {
			for(int y = 0; y < nqy_; y ++) {
				for(int x = 0; x < nqx_; x ++) {
					fc_[0 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1mi[z * nqx_ * nqy_ + y * nqx_ + x] * a1mf[z * nqx_ * nqy_ + y * nqx_ + x];
					fc_[1 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1pi[z * nqx_ * nqy_ + y * nqx_ + x] * a1mf[z * nqx_ * nqy_ + y * nqx_ + x];
					fc_[2 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1mi[z * nqx_ * nqy_ + y * nqx_ + x] * a1pf[z * nqx_ * nqy_ + y * nqx_ + x];
					fc_[3 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1pi[z * nqx_ * nqy_ + y * nqx_ + x] * a1pf[z * nqx_ * nqy_ + y * nqx_ + x];
				} // for x
			} // for y
		} // for z

		delete[] a1pf;
		delete[] a1mf;
		delete[] a1pi;
		delete[] a1mi;

		return true;
	} // HipGISAXS::compute_fresnel_coefficients_embedded()


	// optimize ...
	bool HipGISAXS::compute_fresnel_coefficients_top_buried(float_t alpha_i) {
		//complex_t tk1 = ((float_t) 2.0 * sin(alpha_i)) / ((float_t)sin(alpha_i) +
		//										sqrt((float_t)pow(sin(alpha_i), 2) - dns2_));
		complex_t tk1 = ((float_t) (2.0 * sin(alpha_i))) / ((float_t) sin(alpha_i) + sqrt((float_t)pow(sin(alpha_i), 2) - dns2_));
		complex_t rk1 = tk1 - complex_t(1.0, 0.0);

		fc_ = new (std::nothrow) complex_t[5 * nqx_ * nqy_ * nqz_];
		if(fc_ == NULL) {
			std::cerr << "error: failed to allocate memory for fc" << std::endl;
			return false;
		} // if
		for(unsigned int z = 0; z < nqz_; z ++) {
			for(unsigned int y = 0; y < nqy_; y ++) {
				for(unsigned int x = 0; x < nqx_; x ++) {
					fc_[z * nqy_ * nqx_ + y * nqx_ + x] = rk1;
				} // for x
			} // for y
		} // for z

		float_t k1z = -1.0 * k0_ * sin(alpha_i);
		complex_t tk2, rk2;
		for(unsigned int z = nqz_ - 1; z >= 0; -- z) {
			complex_t k2z = QGrid::instance().qz(z) + k1z;
			if(k2z < 0) {
				tk2 = complex_t(0.0, 0.0);
				rk2 = complex_t(0.0, 0.0);
				for(unsigned int y = 0; y < nqy_; ++ y) {
					for(unsigned int x = 0; x < nqx_; ++ x) {
						fc_[0 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(0.0, 0.0);
						fc_[4 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(0.0, 0.0);
					} // for x
				} // for y
			} else {
				complex_t sf = k2z / k0_;
				tk2 = 2.0 * sf / (sf + sqrt(sf * sf - dns2_));
				rk2 = tk2 - (float_t) 1.0;
				for(unsigned int y = 0; y < nqy_; y ++) {
					for(unsigned int x = 0; x < nqx_; x ++) {
						fc_[4 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(1.0, 0.0);
					} // for x
				} // for y
			} // if-else

			for(unsigned int y = 0; y < nqy_; y ++) {
				for(unsigned int x = 0; x < nqx_; x ++) {
					fc_[1 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = rk2;
					fc_[2 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = rk1 * rk2;
					fc_[3 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] = tk1 * tk2;
				} // for x
			} // for y
		} // for z

		return true;
	} // HipGISAXS::compute_fresnel_coefficients_top_buried()


	bool HipGISAXS::spatial_distribution(structure_iterator_t s, float_t tz, int dim,
										int& rand_dim_x, int& rand_dim_y, float_t* &d) {
		vector3_t spacing = (*s).second.ensemble_spacing();
		vector3_t maxgrains = (*s).second.ensemble_maxgrains();
		std::string distribution = (*s).second.ensemble_distribution();
		// vol_, cell_

		vector3_t spaced_cell = cell_ + spacing;

		if(distribution == "random") {
			srand(time(NULL));
			if(dim == 3) {
				// find max density - number of domains in vol
				vector3_t max_density = min(floor(vol_ / cell_) + 1, maxgrains);
				rand_dim_x = max(1, (int)std::floor(max_density[0] * max_density[1] * max_density[2] / 4));
				rand_dim_y = 3;

				// construct random matrix
				float_t *d_rand = new (std::nothrow) float_t[rand_dim_x * rand_dim_y];
				srand(time(NULL));
				for(int i = 0; i < rand_dim_x * rand_dim_y; ++ i) d_rand[i] = ((float_t)rand() / RAND_MAX);

				d = new (std::nothrow) float_t[rand_dim_x * rand_dim_y * 4];

				int base_index = 0;
				float_t mul_val1 = vol_[0] / 2;
				float_t mul_val2 = vol_[1] / 2;
				float_t mul_val3 = vol_[2];
				int y = 0;
				for(int x = 0; x < rand_dim_x; ++ x) {	// d1
					d[4 * rand_dim_x * y + x] = d_rand[rand_dim_x * y + x] * mul_val1;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d2
					d[4 * rand_dim_x * y + rand_dim_x + x] =
							d_rand[rand_dim_x * y + rand_dim_x + x] * mul_val1 * -1.0;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d3
					d[4 * rand_dim_x * y + 3 * rand_dim_x + x] =
							d_rand[rand_dim_x * y + 2 * rand_dim_x + x] * mul_val1;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d4
					d[4 * rand_dim_x * y + 3 * rand_dim_x + x] =
							d_rand[rand_dim_x * y + 3 * rand_dim_x + x] * mul_val1 * -1.0;
				} // for x
				y = 1;
				for(int x = 0; x < rand_dim_x; ++ x) {	// d1
					d[4 * rand_dim_x * y + x] = d_rand[rand_dim_x * y + x] * mul_val2;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d2
					d[4 * rand_dim_x * y + rand_dim_x + x] =
							d_rand[rand_dim_x * y + rand_dim_x + x] * mul_val2;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d3
					d[4 * rand_dim_x * y + 3 * rand_dim_x + x] =
							d_rand[rand_dim_x * y + 2 * rand_dim_x + x] * mul_val2 * -1.0;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d4
					d[4 * rand_dim_x * y + 3 * rand_dim_x + x] =
							d_rand[rand_dim_x * y + 3 * rand_dim_x + x] * mul_val2 * -1.0;
				} // for x
				y = 2;
				for(int x = 0; x < rand_dim_x; ++ x) {	// d1
					d[4 * rand_dim_x * y + x] = d_rand[rand_dim_x * y + x] * mul_val3 + tz;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d2
					d[4 * rand_dim_x * y + rand_dim_x + x] =
							d_rand[rand_dim_x * y + rand_dim_x + x] * mul_val3 + tz;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d3
					d[4 * rand_dim_x * y + 3 * rand_dim_x + x] =
							d_rand[rand_dim_x * y + 2 * rand_dim_x + x] * mul_val3 + tz;
				} // for x
				for(int x = 0; x < rand_dim_x; ++ x) {	// d4
					d[4 * rand_dim_x * y + 3 * rand_dim_x + x] =
							d_rand[rand_dim_x * y + 3 * rand_dim_x + x] * mul_val3 + tz;
				} // for x
			} else if(dim == 2) {
				std::cerr << "error: dim == 2 case not implemented" << std::endl;
				return false;
			} else if(dim == 1) {
				std::cerr << "error: dim == 1 case not implemented" << std::endl;
				return false;
			} else {
				std::cerr << "error: invalid dim size " << dim << std::endl;
				return false;
			} // if-else

		} else if(distribution == "regular") {
			if(dim == 3) {
				vector3_t nd = min(floor(vol_ / spaced_cell) + 1, maxgrains);

				int size = nd[0] * nd[1] * nd[2];
				d = new (std::nothrow) float_t[3 * size];
				float_t* d1 = d;
				float_t* d2 = d + size;
				float_t* d3 = d + 2 * size;

				rand_dim_x = size;		// total number of rows in d
				rand_dim_y = 3;			// 3 dims, number of cols in d

				// this is more correct I think:
				float_t* dx = new (std::nothrow) float_t[(int)nd[0]];
				float_t* dy = new (std::nothrow) float_t[(int)nd[1]];
				float_t* dz = new (std::nothrow) float_t[(int)nd[2]];
				float_t val_dx = 0, step_dx = spaced_cell[0], max_dx = spaced_cell[0] * (nd[0] - 1);
				for(int i = 0; i < (int)nd[0]; ++ i, val_dx += step_dx) {
					dx[i] = val_dx - max_dx / 2;
				} // for
				float_t val_dy = 0, step_dy = spaced_cell[1], max_dy = spaced_cell[1] * (nd[1] - 1);
				for(int i = 0; i < (int)nd[1]; ++ i, val_dy += step_dy) {
					dy[i] = val_dy - max_dy / 2;
				} // for
				float_t val_dz = 0, step_dz = spaced_cell[2];
				for(int i = 0; i < (int)nd[2]; ++ i, val_dz += step_dz) {
					dz[i] = val_dz + tz;
				} // for

				for(int x = 0; x < (int)nd[0]; ++ x) {
					for(int yz = 0; yz < (int)nd[1] * (int)nd[2]; ++ yz) {
						d1[(int)nd[1] * (int)nd[2] * x + yz] = dx[x];
					} // for
				} // for
				for(int x = 0; x < (int)nd[0]; ++ x) {
					for(int y = 0; y < (int)nd[1]; ++ y) {
						for(int z = 0; z < (int)nd[2]; ++ z) {
							d2[(int)nd[1] * (int)nd[2] * x + (int)nd[2] * y + z] = dy[y];
						} // for
					} // for
				} // for
				for(int xy = 0; xy < (int)nd[0] * (int)nd[1]; ++ xy) {
					for(int z = 0; z < nd[2]; ++ z) {
						d3[(int)nd[2] * xy + z] = dz[z];
					} // for z
				} // for xy
				delete[] dz;
				delete[] dy;
				delete[] dx;
			} else if(dim == 2) {
				std::cerr << "error: dim == 2 case not implemented" << std::endl;
				return false;
			} else if(dim == 1) {
				std::cerr << "error: dim == 1 case not implemented" << std::endl;
				return false;
			} else {
				std::cerr << "error: invalid dim size " << dim << std::endl;
				return false;
			} // if-else
		} else {
			// read .spa file ...
			std::cerr << "uh-oh: seems like you wanted to read distribution from a file" << std::endl;
			std::cerr << "sorry dear, this has not been implemented yet" << std::endl;
			return false;
		} // if-else

		return true;
	} // HipGISAXS::spatial_distribution()


	bool HipGISAXS::orientation_distribution(structure_iterator_t s, float_t* dd,
												int ndx, int ndy, float_t* &nn) {
		std::string distribution = (*s).second.grain_orientation();
		vector2_t tau = (*s).second.rotation_tau();
		vector2_t eta = (*s).second.rotation_eta();
		vector2_t zeta = (*s).second.rotation_zeta();

		nn = new (std::nothrow) float_t[ndx * ndy];
												// i believe constructing nn may not be needed ...
		if(distribution == "single") {			// single
			for(int x = 0; x < ndx; ++ x) {
				nn[x] = tau[0];
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				nn[ndx + x] = eta[0];
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				nn[2 * ndx + x] = zeta[0];
			} // for x
		} else if(distribution == "random") {	// random
			for(int x = 0; x < 3 * ndx; ++ x) {
				nn[x] = (float_t(rand()) / RAND_MAX) * 2 * PI_;
			} // for x
		} else if(distribution == "range") {	// range
			float_t dtau = fabs(tau[1] - tau[0]);
			float_t deta = fabs(eta[1] - eta[0]);
			float_t dzeta = fabs(zeta[1] - zeta[0]);
			for(int x = 0; x < ndx; ++ x) {
				nn[x] = tau[0] + (float_t(rand()) / RAND_MAX) * dtau;
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				nn[ndx + x] = eta[0] + (float_t(rand()) / RAND_MAX) * deta;
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				nn[2 * ndx + x] = zeta[0] + (float_t(rand()) / RAND_MAX) * dzeta;
			} // for x
		} else {
			// read .ori file ...
			std::cerr << "uh-oh: I guess you wanted to read orientations from a file" << std::endl;
			std::cerr << "too bad, its not implemented yet" << std::endl;
			return false;
		} // if-else

		return true;
	} // HipGISAXS::orientation_distribution()


	bool HipGISAXS::write_qgrid(char* filename) {
		std::ofstream qout(filename);

		qout << nqx_ << " " << nqy_ << " " << nqz_extended_ << std::endl;
		for(int i = 0; i < nqx_; ++ i) {
			qout << QGrid::instance().qx(i) << " ";
		} // for
		qout << std::endl;
		for(int i = 0; i < nqy_; ++ i) {
			qout << QGrid::instance().qy(i) << " ";
		} // for
		qout << std::endl;
		for(int i = 0; i < nqz_extended_; ++ i) {
			qout << QGrid::instance().qz_extended(i).real() << " "
					<< QGrid::instance().qz_extended(i).imag() << " ";
		} // for
		qout << std::endl;

		qout.close();
		return true;
	} // HipGISAXS::write_qgrid()


	bool HipGISAXS::read_form_factor(const char* filename) {
		return ff_.read_form_factor(filename, nqx_, nqy_, nqz_extended_);
	} // HipGISAXS::read_form_factor()

} // namespace hig


/* The main for HipGISAXS
 */
int main(int narg, char** args) {

	if(narg != 2) {
		std::cout << "usage: hipgisaxs <input_config>" << std::endl;
		return 1;
	} // if

	/* initialize MPI */
	MPI::Init(narg, args);
	int mpi_rank = MPI::COMM_WORLD.Get_rank();
	int mpi_num_procs = MPI::COMM_WORLD.Get_size();

	/* read input file and construct input structures */
	hig::HipGISAXS my_gisaxs;
	if(!my_gisaxs.construct_input(args[1])) {
		if(mpi_rank == 0) std::cerr << "error: failed to construct input containers" << std::endl;
		MPI::Finalize();
		return 1;
	} // if
	//hig::HiGInput::instance().print_all();	// for testing

	/* run the simulation */
	if(!my_gisaxs.run_all_gisaxs(MPI::COMM_WORLD)) {
		std::cerr << "error: could not run the simulation - some error occured" << std::endl;
		MPI::Finalize();
		return 1;
	} // if

	MPI::Finalize();
	return 0;
} // main()
