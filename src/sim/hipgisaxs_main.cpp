/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_main.cpp
 *  Created: Jun 14, 2012
 *  Modified: Fri 10 Jan 2014 10:05:30 AM PST
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

#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "../woo/timer/woo_boostchronotimers.hpp"

#include "hipgisaxs_main.hpp"
#include "../common/typedefs.hpp"
#include "../utils/utilities.hpp"

#if defined USE_GPU || defined FF_ANA_GPU || defined FF_NUM_GPU
	#include "../init/gpu/init_gpu.cuh"
#elif defined USE_MIC
	#include "../init/mic/init_mic.hpp"
#endif


namespace hig {

	HipGISAXS::HipGISAXS(int narg, char** args): freq_(0.0), k0_(0.0),
				num_layers_(0), num_structures_(0),
				nqx_(0), nqy_(0), nqz_(0), nqz_extended_(0)
				#ifdef USE_MPI
					, multi_node_(narg, args)
				#endif
				//#ifdef FF_NUM_GPU   // use GPU
				//	#ifdef FF_NUM_GPU_FUSED
				//		ff_(64, 8)
				//	#elif defined KERNEL2
				//		ff_(2, 4, 4)
				//	#else
				//		ff_(64)
				//	#endif
				//#else   // use CPU
				//	ff_()
				//#endif
					{
		single_layer_refindex_.delta(0.0);
		single_layer_refindex_.beta(0.0);
		single_layer_thickness_ = 0.0;
		HiGInput::instance();
		QGrid::instance();
	} // HipGISAXS::HipGISAXS()


	HipGISAXS::~HipGISAXS() {
		// nothing to do here yet ...
	} // HipGISAXS::~HipGISAXS()


	bool HipGISAXS::init() {
						// is called at the beginning of the runs (after input is read)
						// it does the following:
						// 	+ set detector/system stuff
						// 	+ initialize output dir
						// 	+ create q-grid
						// 	+ construct layer profile
						// 	+ get layer profile data
						// 	+ get structure info
						// 	+ construct lattice vectors
						// 	+ compute cell size
		// TODO first check if the input has been constructed ...

		#ifdef USE_MPI
			int mpi_rank = multi_node_.rank();
			bool master = multi_node_.is_master();
		#else
			int mpi_rank = 0;
			bool master = true;
		#endif

		if(master) {
			std::cout << std::endl
					<< "*******************************************************************" << std::endl
					<< "*********************** HipGISAXS v0.9b ***************************" << std::endl
					<< "*******************************************************************" << std::endl
					<< std::endl;
		} // if

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
			if(master) std::cerr << "error: photon energy is not given in 'ev'" << std::endl;
			return false;
		} // if-else

		k0_ = 2 * PI_ * freq_ / LIGHT_SPEED_;

		// create output directory
		if(master) {		// this is not quite good for mpi ... improve ...
			const std::string p = HiGInput::instance().path() + "/" + HiGInput::instance().runname();
			if(!boost::filesystem::create_directory(p)) {
				std::cerr << "error: could not create output directory " << p << std::endl;
				return false;
			} // if
		} // if

		#ifdef USE_MPI
			multi_node_.barrier();
		#endif

		// create Q-grid
		float_t min_alphai = HiGInput::instance().scattering_min_alpha_i() * PI_ / 180;
		if(!QGrid::instance().create(freq_, min_alphai, k0_, mpi_rank)) {
			if(master) std::cerr << "error: could not create Q-grid" << std::endl;
			return false;
		} // if

		nqx_ = QGrid::instance().nqx();
		nqy_ = QGrid::instance().nqy();
		nqz_ = QGrid::instance().nqz();

		/* construct layer profile */
		if(!HiGInput::instance().construct_layer_profile()) {	// also can be done at input reading ...
			if(master) std::cerr << "error: could not construct layer profile" << std::endl;
			return false;
		} // if

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
		// TODO: this can also be done at input reading ...
		if(!HiGInput::instance().construct_lattice_vectors()) {
			if(master) std::cerr << "error: could not construct lattice vectors" << std::endl;
			return false;
		} // if

		/* domain size */
		vector3_t min_vec(0.0, 0.0, 0.0), max_vec(0.0, 0.0, 0.0);
		float_t z_min_0 = 0.0, z_max_0 = 0.0;
		if(!HiGInput::instance().compute_domain_size(min_vec, max_vec, z_min_0, z_max_0)) {
			if(master) std::cerr << "error: could not compute domain size" << std::endl;
			return false;
		} // if
//		std::cout << "++ Domain min point: " << min_vec[0] << ", " << min_vec[1]
//					<< ", " << min_vec[2] << std::endl;
//		std::cout << "++ Domain max point: " << max_vec[0] << ", " << max_vec[1]
//					<< ", " << max_vec[2] << std::endl;
//		std::cout << "++ Domain z max and min: " << z_max_0 << ", " << z_min_0 << std::endl;
//		std::cout << "++ Domain dimensions: " << max_vec[0] - min_vec[0] << " x "
//					<< max_vec[1] - min_vec[1] << " x "
//					<< max_vec[2] - min_vec[2] << " ("
//					<< z_max_0 - z_min_0 << ")" << std::endl;

		// cell size ... check ... FIXME
		cell_[0] = fabs(max_vec[0] - min_vec[0]);
		cell_[1] = fabs(max_vec[1] - min_vec[1]);
		cell_[2] = fabs(z_max_0 - z_min_0);


		#ifdef _OPENMP
			if(master)
				std::cout << "++ Number of host OpenMP threads: "
							<< omp_get_max_threads() << std::endl;
		#endif

		#if defined USE_GPU || defined FF_ANA_GPU || defined FF_NUM_GPU
			if(master) std::cout << "-- Waking up GPU(s) ..." << std::flush;
			init_gpu();
			if(master) std::cout << " it woke up!" << std::endl;
		#elif defined USE_MIC
			if(master) std::cout << "-- Waking up MIC(s) ..." << std::flush;
			init_mic();
			if(master) std::cout << " done." << std::endl;
		#else
			if(master) std::cout << "-- Not set up to use any accelerator!" << std::endl;
		#endif

		return true;
	} // HipGISAXS::init()


	/*
	// this is temporary, for newton's fit method
	bool HipGISAXS::init_steepest_fit(float_t qzcut) {
						// is called at the beginning of the runs (after input is read)
		// first check if the input has been constructed ...

		#ifdef USE_MPI
			int mpi_rank = multi_node_.rank();
			bool master = multi_node_.is_master();
		#else
			bool master = true;
		#endif

		//photon conversion
		float_t photon = 0.0;
		std::string unit;
		freq_ = 0; k0_ = 0;
		HiGInput::instance().photon_energy(photon, unit);
		if(unit == "ev") {
			photon = photon / 1000;		// in keV
			freq_ = 1e-9 * photon * 1.60217646e-19 * 1000 / 6.626068e-34;
		} else { // do something else ? ...
			if(master) std::cerr << "error: photon energy is not given in 'ev'" << std::endl;
			return false;
		} // if-else

		k0_ = 2 * PI_ * freq_ / LIGHT_SPEED_;

		// create output directory
		if(master) {		// this is not quite good for mpi ... improve ...
			const std::string p = HiGInput::instance().path() + "/" + HiGInput::instance().runname();
			if(!boost::filesystem::create_directory(p)) {
				std::cerr << "error: could not create output directory " << p << std::endl;
				return false;
			} // if
		} // if

		#ifdef USE_MPI
			multi_node_.barrier();
		#endif

		// create Q-grid
		float_t min_alphai = HiGInput::instance().scattering_min_alpha_i() * PI_ / 180;
		if(!QGrid::instance().create_z_cut(freq_, min_alphai, k0_, qzcut)) {
			if(master) std::cerr << "error: could not create Q-grid" << std::endl;
			return false;
		} // if

		nqx_ = QGrid::instance().nqx();
		nqy_ = QGrid::instance().nqy();
		nqz_ = QGrid::instance().nqz();

		// construct layer profile
		if(!HiGInput::instance().construct_layer_profile()) {	// also can be done at input reading ...
			if(master) std::cerr << "error: could not construct layer profile" << std::endl;
			return false;
		} // if

		return true;
	} // HipGISAXS::init_steepest_fit()
	*/


	bool HipGISAXS::run_init(float_t alphai, float_t phi, float_t tilt, SampleRotation& rot_matrix) {
					// this is called for each config-run during the main run
					// it does the following:
					// 	+ construct the illuminated volume
					//	+ construct rotation matrices

		#ifdef USE_MPI
			bool master = multi_node_.is_master();
		#else
			bool master = true;
		#endif

		if(!illuminated_volume(alphai, HiGInput::instance().scattering_spot_area(),
				HiGInput::instance().min_layer_order(), HiGInput::instance().substrate_refindex())) {
			if(master) std::cerr << "error: something went wrong in illuminated_volume()" << std::endl;
			return false;
		} // if

		/* rotation matrices */
		// compute r_phi = sample rotation by phi
		vector3_t rp1(0.0, 0.0, 0.0), rp2(0.0, 0.0, 0.0), rp3(0.0, 0.0, 0.0);
		vector3_t rt1(0.0, 0.0, 0.0), rt2(0.0, 0.0, 0.0), rt3(0.0, 0.0, 0.0);
		compute_rotation_matrix_z(phi, rp1, rp2, rp3);
		compute_rotation_matrix_x(tilt, rt1, rt2, rt3);
		mat_mul_3x3(rp1, rp2, rp3, rt1, rt2, rt3, rot_matrix.r1_, rot_matrix.r2_, rot_matrix.r3_);

		return true;
	} // HipGISAXS::run_init()


	/**
	 * This is the main function called from outside
	 * It loops over all configurations and calls the simulation routine
	 */
	bool HipGISAXS::run_all_gisaxs(int x_min, int x_max, int x_step) {
		#ifdef USE_MPI
			// this is for the whole comm world
			const char* world_comm = "world";
			bool master = multi_node_.is_master(world_comm);
		#else
			bool master = true;
		#endif

		if(!init()) return false;

		woo::BoostChronoTimer sim_timer;
		sim_timer.start();

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

		if(master) {
			std::cout << "**                  Num alphai: " << num_alphai << std::endl
						<< "**                     Num phi: " << num_phi << std::endl
						<< "**                    Num tilt: " << num_tilt << std::endl;
		} // if

		// loop over all alphai, phi, and tilt

		#ifdef USE_MPI
			// divide among processors
			int num_procs = multi_node_.size(world_comm);
			int rank = multi_node_.rank(world_comm);
			int alphai_color = 0;
			if(num_procs > num_alphai) {
				alphai_color = rank % num_alphai;
				alphai_min = alphai_min + alphai_step * alphai_color;
				num_alphai = 1;
			} else {
				alphai_color = rank;
				alphai_min = alphai_min + alphai_step *
							((num_alphai / num_procs) * rank + min(rank, num_alphai % num_procs));
				num_alphai = (num_alphai / num_procs) + (rank < num_alphai % num_procs);
			} // if-else
			const char* alphai_comm = "alphai";
			multi_node_.split(alphai_comm, world_comm, alphai_color);

			bool amaster = multi_node_.is_master(alphai_comm);
			int temp_amaster = amaster;
			int *amasters = new (std::nothrow) int[multi_node_.size(world_comm)];
			// all alphai masters tell the world master about who they are
			multi_node_.allgather(world_comm, &temp_amaster, 1, amasters, 1);
		#else
			bool amaster = true;
		#endif // USE_MPI

		float_t alpha_i = alphai_min;
		for(int i = 0; i < num_alphai; i ++, alpha_i += alphai_step) {
			float_t alphai = alpha_i * PI_ / 180;

			float_t* averaged_data = NULL;		// to hold summation of all (if needed)

			#ifdef USE_MPI
				// divide among processors
				int num_procs = multi_node_.size(alphai_comm);
				int rank = multi_node_.rank(alphai_comm);
				int phi_color = 0;
				if(num_procs > num_phi) {
					phi_color = rank % num_phi;
					phi_min = phi_min + phi_step * phi_color;
					num_phi = 1;
				} else {
					phi_color = rank;
					phi_min = phi_min + phi_step *
								((num_phi / num_procs) * rank + min(rank, num_phi % num_procs));
					num_phi = (num_phi / num_procs) + (rank < num_phi % num_procs);
				} // if-else
				const char* phi_comm = "phi";
				multi_node_.split(phi_comm, alphai_comm, phi_color);

				bool pmaster = multi_node_.is_master(phi_comm);
				int temp_pmaster = pmaster;
				int *pmasters = new (std::nothrow) int[multi_node_.size(alphai_comm)];
				// all phi masters tell the alphai master about who they are
				multi_node_.allgather(alphai_comm, &temp_pmaster, 1, pmasters, 1);
			#else
				bool pmaster = true;
			#endif // USE_MPI

			float_t phi = phi_min;
			for(int j = 0; j < num_phi; j ++, phi += phi_step) {
				float_t phi_rad = phi * PI_ / 180;

				#ifdef USE_MPI
					// divide among processors
					int num_procs = multi_node_.size(phi_comm);
					int rank = multi_node_.rank(phi_comm);
					int tilt_color = 0;
					if(num_procs > num_tilt) {
						tilt_color = rank % num_tilt;
						tilt_min = tilt_min + tilt_step * tilt_color;
						num_tilt = 1;
					} else {
						tilt_color = rank;
						tilt_min = tilt_min + tilt_step *
								((num_tilt / num_procs) * rank + min(rank, num_tilt % num_procs));
						num_tilt = (num_tilt / num_procs) + (rank < num_tilt % num_procs);
					} // if-else
					const char* tilt_comm = "tilt";
					multi_node_.split(tilt_comm, phi_comm, tilt_color);

					bool tmaster = multi_node_.is_master(tilt_comm);
					int temp_tmaster = tmaster;
					int *tmasters = new (std::nothrow) int[multi_node_.size(phi_comm)];
					// all tilt masters tell the phi master about who they are
					multi_node_.allgather(phi_comm, &temp_tmaster, 1, tmasters, 1);
				#else
					bool tmaster = true;
				#endif // USE_MPI

				float_t tilt = tilt_min;
				for(int k = 0; k < num_tilt; k ++, tilt += tilt_step) {
					float_t tilt_rad = tilt * PI_ / 180;

					if(tmaster) {
						std::cout << "-- Computing GISAXS "
									<< i * num_phi * num_tilt + j * num_tilt + k + 1 << " / "
									<< num_alphai * num_phi * num_tilt
									<< " [alphai = " << alpha_i << ", phi = " << phi
									<< ", tilt = " << tilt << "] ..." << std::endl << std::flush;
					} // if

					/* run a gisaxs simulation */

					float_t* final_data = NULL;
					if(!run_gisaxs(alpha_i, alphai, phi_rad, tilt_rad, final_data,
								#ifdef USE_MPI
									tilt_comm,
								#endif
								0)) {
						if(tmaster)
							std::cerr << "error: could not finish successfully" << std::endl;
						return false;
					} // if

					// for future - 3d image ...
					//Image img3d(nqx_, nqy_, nqz_);
					//if(!run_gisaxs(alphai, phi, tilt, img3d))
					
					if(tmaster) {
						// note that final_data stores 3d info
						// for 2d, just take a slice of the data
						std::cout << "-- Constructing GISAXS image ... " << std::flush;
						//Image img(nqx_, nqy_, nqz_);
						// testing ...
						//Image img(nqx_, nqy_, nqz_, 37, 36, 27);
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

						//img.save(output, x_min, x_max);
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

					// also compute averaged values over phi and tilt
					if(num_phi > 1 || num_tilt > 1) {
						if(tmaster) {
							if(averaged_data == NULL) {
								averaged_data = new (std::nothrow) float_t[nqx_ * nqy_ * nqz_];
								memset(averaged_data, 0, nqx_ * nqy_ * nqz_ * sizeof(float_t));
							} // if
							int i = 0;
							add_data_elements(averaged_data, final_data, averaged_data,
												nqx_ * nqy_ * nqz_);
						} // if
						delete[] final_data;
					} else {
						averaged_data = final_data;
					} // if-else

					#ifdef USE_MPI
						multi_node_.barrier(tilt_comm);
					#endif

				} // for tilt
				#ifdef USE_MPI
					multi_node_.free(tilt_comm);

//					if(num_phi > 1 || num_tilt > 1) {
						// get data from all other phi_comm processors which were tmasters
						int psize = multi_node_.size(phi_comm);
						float_t* temp_data = NULL;
						if(psize > 1) {
							int *proc_sizes = new (std::nothrow) int[psize];
							int *proc_displacements = new (std::nothrow) int[psize];
							if(pmaster) {
								temp_data = new (std::nothrow) float_t[psize * nqx_ * nqy_ * nqz_];
							} // if
							for(int i = 0; i < psize; ++ i) {
								proc_sizes[i] = tmasters[i] * nqx_ * nqy_ * nqz_;
							} // for
							proc_displacements[0] = 0;
							for(int i = 1; i < psize; ++ i) {
								proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
							} // for
							int prank = multi_node_.rank(phi_comm);
							multi_node_.gatherv(phi_comm, averaged_data, proc_sizes[prank],
												temp_data, proc_sizes, proc_displacements);
							if(pmaster) {
								for(int i = 1; i < psize; ++ i) {
									add_data_elements(averaged_data, temp_data + proc_displacements[i],
														averaged_data, proc_sizes[i]);
								} // for
								delete[] temp_data;
							} // if
							delete[] proc_displacements;
							delete[] proc_sizes;
						} // if
//					} // if

					multi_node_.barrier(phi_comm);

				#endif
			} // for phi
			#ifdef USE_MPI
				multi_node_.free(phi_comm);

//				if(num_phi > 1 || num_tilt > 1) {
					// get data from all other phi_comm processors which were tmasters
					int asize = multi_node_.size(alphai_comm);
					float_t* temp_data = NULL;
					if(asize > 1) {
						int *proc_sizes = new (std::nothrow) int[asize];
						int *proc_displacements = new (std::nothrow) int[asize];
						if(amaster) {
							temp_data = new (std::nothrow) float_t[asize * nqx_ * nqy_ * nqz_];
						} // if
						for(int i = 0; i < asize; ++ i) {
							proc_sizes[i] = pmasters[i] * nqx_ * nqy_ * nqz_;
						} // for
						proc_displacements[0] = 0;
						for(int i = 1; i < asize; ++ i) {
							proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
						} // for
						int arank = multi_node_.rank(alphai_comm);
						multi_node_.gatherv(alphai_comm, averaged_data, proc_sizes[arank],
											temp_data, proc_sizes, proc_displacements);
						if(amaster) {
							for(int i = 1; i < asize; ++ i) {
								add_data_elements(averaged_data, temp_data + proc_displacements[i],
													averaged_data, nqx_ * nqy_ * nqz_);
							} // for
							delete[] temp_data;
						} // if
						delete[] proc_displacements;
						delete[] proc_sizes;
					} // if
//				} // if

				multi_node_.barrier(alphai_comm);

			#endif

			if(amaster && (num_phi > 1 || num_tilt > 1)) {
				if(averaged_data != NULL) {
					Image img(nqx_, nqy_, nqz_);
					img.construct_image(averaged_data, 0); // slice x = 0

					// define output filename
					std::stringstream alphai_b;
					std::string alphai_s;
					alphai_b << alpha_i; alphai_s = alphai_b.str();
					std::string output(HiGInput::instance().param_pathprefix() +
										"/" + HiGInput::instance().runname() +
										"/img_ai=" + alphai_s + "_averaged.tif");
					std::cout << "-- Saving averaged image in " << output << " ... " << std::flush;
					img.save(output);
					std::cout << "done." << std::endl;

					// save the actual data into a file also
					std::string data_file(HiGInput::instance().param_pathprefix() +
									"/" + HiGInput::instance().runname() +
									"/gisaxs_ai=" + alphai_s + "_averaged.out");
					std::cout << "-- Saving averaged raw data in " << data_file << " ... " << std::flush;
					save_gisaxs(averaged_data, data_file);
					std::cout << "done." << std::endl;

					delete[] averaged_data;
				} // if
			} // if

		} // for alphai
		#ifdef USE_MPI
			multi_node_.free(alphai_comm);
		#endif

		sim_timer.stop();
		if(master) {
			std::cout << "**         Total simulation time: " << sim_timer.elapsed_msec() << " ms."
						<< std::endl;
		} // if

		return true;
	} // HipGISAXS::run_all_gisaxs()


	/**
	 * used for fitting
	 */
	bool HipGISAXS::fit_init() { return init(); }
	bool HipGISAXS::compute_gisaxs(float_t* &final_data) {
		#ifdef USE_MPI
			// this is for the whole comm world
			const char* world_comm = "world";
			bool master = multi_node_.is_master(world_comm);
		#else
			bool master = true;
		#endif

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
		if(num_alphai > 1 || num_phi > 1 || num_tilt > 1) {
			if(master)
				std::cerr << "error: currently you can simulate only for single "
							<< "alpha_i, phi and tilt angles"
							<< std::endl;
			return -1.0;
		} // if

		woo::BoostChronoTimer sim_timer;
		sim_timer.start();
		float_t alpha_i = alphai_min;
		float_t alphai = alpha_i * PI_ / 180;
		float_t phi_rad = phi_min * PI_ / 180;
		float_t tilt_rad = tilt_min * PI_ / 180;
		if(master) std::cout << "-- Computing GISAXS ... " << std::endl << std::flush;
		/* run a gisaxs simulation */
		if(!run_gisaxs(alpha_i, alphai, phi_rad, tilt_rad, final_data,
					#ifdef USE_MPI
						world_comm,
					#endif
					0)) {
			if(master) std::cerr << "error: could not finish successfully" << std::endl;
			return -1.0;
		} // if
		sim_timer.stop();
		if(master)
			std::cout << "**        Total Simulation time: " << sim_timer.elapsed_msec()
						<< " ms." << std::endl;

		return true;
	} // HipGISAXS::compute_gisaxs()

	
	/**
	 * run an experiment configuration
	 * called for each configuration
	 */
	/* all the real juice is here */
	bool HipGISAXS::run_gisaxs(float_t alpha_i, float_t alphai, float_t phi, float_t tilt,
								float_t* &img3d, const char* comm_key, int corr_doms) {

		SampleRotation rotation_matrix;
		if(!run_init(alphai, phi, tilt, rotation_matrix)) return false;

		#ifdef USE_MPI
			bool master = multi_node_.is_master(comm_key);
			int ss = multi_node_.size(comm_key);
		#else
			bool master = true;
		#endif

		// compute propagation coefficients/fresnel coefficients
		// this can also go into run_init() ...
		complex_t *amm = NULL, *apm = NULL, *amp = NULL, *app = NULL;
		complex_t *rk1 = NULL, *rk2 = NULL, *rk1rk2 = NULL, *tk1tk2 = NULL, *h0 = NULL;
		// TODO : where is fc used ?????????? ....................
		complex_t* fc = NULL;
		if(!compute_propagation_coefficients(alphai, amm, apm, amp, app,
												rk1, rk2, rk1rk2, tk1tk2, h0, fc)) {
			if(master) std::cerr << "error: failed to compute propogation coefficients" << std::endl;
			return false;
		} // if

		/* loop over all structures and grains/grains */
		structure_iterator_t s = HiGInput::instance().structure_begin();
		int num_structs = num_structures_;
		#ifdef USE_MPI
			// divide among processors
			int num_procs = multi_node_.size(comm_key);
			int rank = multi_node_.rank(comm_key);
			int struct_color = 0;
			int soffset = 0;
			if(num_procs > num_structs) {
				struct_color = rank % num_structs;
				soffset = struct_color;
				num_structs = 1;
			} else {
				struct_color = rank;
				soffset = ((num_structs / num_procs) * rank + min(rank, num_structs % num_procs));
				num_structs = (num_structs / num_procs) + (rank < num_structs % num_procs);
			} // if-else
			const char* struct_comm = "structure";
			multi_node_.split(struct_comm, comm_key, struct_color);
			for(int i = 0; i < soffset; ++ i) ++ s;

			bool smaster = multi_node_.is_master(struct_comm);
			int temp_smaster = smaster;
			int *smasters = new (std::nothrow) int[multi_node_.size(comm_key)];
			// all grain masters tell the structure master about who they are
			multi_node_.gather(comm_key, &temp_smaster, 1, smasters, 1);
		#else
			bool smaster = true;
		#endif // USE_MPI

		// initialize memory for struct_intensity
		unsigned int size = nqx_ * nqy_ * nqz_;
		float_t* struct_intensity = NULL;
		complex_t* c_struct_intensity = NULL;
		if(smaster) {
			struct_intensity = new (std::nothrow) float_t[num_structs * size];
			c_struct_intensity = new (std::nothrow) complex_t[num_structs * size];
		} // master

		// loop over structures
		for(int s_num = 0; s_num < num_structs && s != HiGInput::instance().structure_end();
				++ s, ++ s_num) {
			// get all repetitions of the structure in the volume
			// with position dr and orientation (tau, eta)
			// i.e. compute matrix that defines all num_grains grains inside the volume
			// distr = [ drx_i dry_i drz_i  tau_i eta_i ] (NDISTR x 5 )
			// get params of the shape repeated in this structure: [dn2, id, dims, invar, t, seta, stau]
			// allocate id(nqy, nqz, num_grains)

			// get the shape
			// compute t, lattice, ndoms, dd, nn, id etc.

			if(smaster) {
				std::cout << "-- Processing structure " << s_num + 1 << " ..." << std::endl;
			} // if

			Structure *curr_struct = &((*s).second);
			Shape *curr_shape = HiGInput::instance().shape(*curr_struct);
			Lattice *curr_lattice = (Lattice*) HiGInput::instance().lattice(*curr_struct);
			bool struct_in_layer = (*s).second.grain_in_layer();

			vector3_t grain_repeats = (*s).second.grain_repetition();

			float_t *dd = NULL, *nn = NULL;		// come back to this ...
												// these structures can be improved ...
			float_t tz = 0;
			int num_dimen = 3;
			int ndx = 0, ndy = 0;
			// compute dd and nn
			spatial_distribution(s, tz, num_dimen, ndx, ndy, dd);
			orientation_distribution(s, dd, ndx, ndy, nn);
			int num_grains = ndx;

			int r1axis = (int) (*s).second.rotation_rot1()[0];
			int r2axis = (int) (*s).second.rotation_rot2()[0];
			int r3axis = (int) (*s).second.rotation_rot3()[0];

			if(smaster) {
				std::cout << "-- Grains: " << num_grains << std::endl;
			} // if

			//std::cout << "DD: " << std::endl;
			//for(unsigned int d = 0; d < num_grains; ++ d) {
			//	std::cerr << dd[d] << "\t" << dd[num_grains + d] << "\t" << dd[num_grains * 2 + d]
			//				<< std::endl;
			//} // for*/

			vector3_t curr_transvec = (*s).second.grain_transvec();
			if(HiGInput::instance().param_nslices() <= 1) {
				curr_transvec[2] = curr_transvec[2] - single_layer_thickness_;
			} else {
				curr_transvec[2] = curr_transvec[2]; // TODO/FIXME... for more than 1 layers ... 
													 // ... incomplete
			} // if-else
			ShapeName shape_name = HiGInput::instance().shape_name((*s).second);
			float_t shape_tau = HiGInput::instance().shape_tau((*s).second);
			float_t shape_eta = HiGInput::instance().shape_eta((*s).second);
			std::string shape_file = HiGInput::instance().shape_filename((*s).second);
			shape_param_list_t shape_params = HiGInput::instance().shape_params((*s).second);

			/* computing dwba ff for each grain in structure (*s) with
			 * ensemble containing num_grains grains */

			int grain_min = 0;
			int num_gr = num_grains;
			int grain_max = num_gr;

			#ifdef USE_MPI
				// divide among processors
				int num_procs = multi_node_.size(struct_comm);
				int rank = multi_node_.rank(struct_comm);
				int grain_color = 0;
				if(num_procs > num_gr) {
					grain_color = rank % num_gr;
					grain_min = grain_color;
					num_gr = 1;
				} else {
					grain_color = rank;
					grain_min = ((num_gr / num_procs) * rank + min(rank, num_gr % num_procs));
					num_gr = (num_gr / num_procs) + (rank < num_gr % num_procs);
				} // if-else
				grain_max = grain_min + num_gr;
				const char* grain_comm = "grain";
				multi_node_.split(grain_comm, struct_comm, grain_color);

				bool gmaster = multi_node_.is_master(grain_comm);
				int temp_gmaster = gmaster;
				int *gmasters = new (std::nothrow) int[multi_node_.size(struct_comm)];
				// all grain masters tell the structure master about who they are
				multi_node_.gather(struct_comm, &temp_gmaster, 1, gmasters, 1);
			#else
				bool gmaster = true;
			#endif // USE_MPI

			complex_t *grain_ids = NULL;
			if(gmaster) {		// only the structure master needs this
				grain_ids = new (std::nothrow) complex_t[num_gr * nqx_ * nqy_ * nqz_];
				if(grain_ids == NULL) {
					std::cerr << "error: could not allocate memory for 'id'" << std::endl;
					return false;
				} // if
				// initialize to 0
				memset(grain_ids, 0 , num_gr * nqx_ * nqy_ * nqz_ * sizeof(complex_t));
			} // if


			// loop over grains - each process processes num_gr grains
			for(int grain_i = grain_min; grain_i < grain_max; grain_i ++) {	// or distributions

				if(gmaster) {
					std::cout << "-- Processing grain " << grain_i + 1 << " / " << num_grains << " ..."
								<< std::endl;
				} // if

				// define r_norm (grain orientation by tau and eta)
				// define full grain rotation matrix r_total = r_phi * r_norm
				// TODO: ... i think these tau eta zeta can be computed on the fly to save memory ...
				float_t rot1 = nn[0 * num_gr + grain_i];
				float_t rot2 = nn[1 * num_gr + grain_i];
				float_t rot3 = nn[2 * num_gr + grain_i];
				vector3_t z1, z2, z3, e1, e2, e3, t1, t2, t3;
				switch(r1axis) {
					case 0:
						compute_rotation_matrix_x(rot1, z1, z2, z3);
						break;
					case 1:
						compute_rotation_matrix_y(rot1, z1, z2, z3);
						break;
					case 2:
						compute_rotation_matrix_z(rot1, z1, z2, z3);
						break;
					default:
						std::cerr << "error: unknown axis: " << r1axis << std::endl;
						return false;
				} // switch
				switch(r2axis) {
					case 0:
						compute_rotation_matrix_x(rot2, e1, e2, e3);
						break;
					case 1:
						compute_rotation_matrix_y(rot2, e1, e2, e3);
						break;
					case 2:
						compute_rotation_matrix_z(rot2, e1, e2, e3);
						break;
					default:
						std::cerr << "error: unknown axis: " << r2axis << std::endl;
						return false;
				} // switch
				switch(r3axis) {
					case 0:
						compute_rotation_matrix_x(rot3, t1, t2, t3);
						break;
					case 1:
						compute_rotation_matrix_y(rot3, t1, t2, t3);
						break;
					case 2:
						compute_rotation_matrix_z(rot3, t1, t2, t3);
						break;
					default:
						std::cerr << "error: unknown axis: " << r3axis << std::endl;
						return false;
				} // switch

				vector3_t temp1, temp2, temp3;
				vector3_t r_norm1, r_norm2, r_norm3;
				mat_mul_3x3(z1, z2, z3, e1, e2, e3, temp1, temp2, temp3);
				mat_mul_3x3(temp1, temp2, temp3, t1, t2, t3, r_norm1, r_norm2, r_norm3);

				vector3_t r_tot1, r_tot2, r_tot3;
				mat_mul_3x3(rotation_matrix.r1_, rotation_matrix.r2_, rotation_matrix.r3_,
							r_norm1, r_norm2, r_norm3, r_tot1, r_tot2, r_tot3);

				/* center of unit cell replica */
				vector3_t curr_dd_vec(dd[grain_i + 0], dd[grain_i + num_gr], dd[grain_i + 2 * num_gr]);
				vector3_t result(0.0, 0.0, 0.0);

				mat_mul_3x1(rotation_matrix.r1_, rotation_matrix.r2_, rotation_matrix.r3_,
							curr_dd_vec, result);
				vector3_t center = result + curr_transvec;

				/* compute structure factor and form factor */

				StructureFactor sf;
				#ifdef FF_NUM_GPU   // use GPU
					#ifdef FF_NUM_GPU_FUSED
						FormFactor ff(64, 8);
					#elif defined KERNEL2
						FormFactor ff(2, 4, 4);
					#else
						FormFactor ff(64);
					#endif
				#else   // use CPU or MIC
					FormFactor ff;
				#endif

				// TODO: parallel tasks ...

				structure_factor(sf, HiGInput::instance().experiment(), center, curr_lattice,
									grain_repeats, r_tot1, r_tot2, r_tot3
									#ifdef USE_MPI
										, grain_comm
									#endif
									);
				//sf.printsf();

				/*if(master) {
					std::stringstream alphai_b, phi_b, tilt_b;
					std::string alphai_s, phi_s, tilt_s;
					alphai_b << alpha_i; alphai_s = alphai_b.str();
					phi_b << phi; phi_s = phi_b.str();
					tilt_b << tilt; tilt_s = tilt_b.str();
					std::string sf_output(HiGInput::instance().param_pathprefix() +
									"/" + HiGInput::instance().runname() +
									"/sf_ai=" + alphai_s + "_rot=" + phi_s +
									"_tilt=" + tilt_s + ".out");
					std::cout << "-- Saving structure factor in " << sf_output << " ... " << std::flush;
					sf.save_sf(nqx_, nqy_, nqz_extended_, sf_output.c_str());
					std::cout << "done." << std::endl;
				} // if*/

				//read_form_factor("curr_ff.out");
				form_factor(ff, shape_name, shape_file, shape_params, curr_transvec,
							shape_tau, shape_eta, r_tot1, r_tot2, r_tot3
							#ifdef USE_MPI
								, grain_comm
							#endif
							);
				//ff.print_ff(nqx_, nqy_, nqz_extended_);
				//ff.printff(nqx_, nqy_, nqz_extended_);

				// save form factor data on disk
				/*if(master) {
					std::stringstream alphai_b, phi_b, tilt_b;
					std::string alphai_s, phi_s, tilt_s;
					alphai_b << alpha_i; alphai_s = alphai_b.str();
					phi_b << phi; phi_s = phi_b.str();
					tilt_b << tilt; tilt_s = tilt_b.str();
					std::string ff_output(HiGInput::instance().param_pathprefix() +
									"/" + HiGInput::instance().runname() +
									"/ff_ai=" + alphai_s + "_rot=" + phi_s +
									"_tilt=" + tilt_s + ".out");
					std::cout << "-- Saving form factor in " << ff_output << " ... " << std::flush;
					ff.save_ff(nqx_, nqy_, nqz_extended_, ff_output.c_str());
					std::cout << "done." << std::endl;
				} // if*/

				/* compute intensities using sf and ff */
				// TODO: parallelize ...
				if(gmaster) {	// grain master
					complex_t* base_id = grain_ids + grain_i * nqx_ * nqy_ * nqz_;

					unsigned int nslices = HiGInput::instance().param_nslices();
					if(nslices <= 1) {
						/* without slicing */
						if(struct_in_layer &&
								HiGInput::instance().structure_layer_order((*s).second) == 1) {
													// structure is inside a layer, on top of substrate
							if(single_layer_refindex_.delta() < 0 || single_layer_refindex_.beta() < 0) {
								// this should never happen
								std::cerr << "error: single layer information not correctly set"
										<< std::endl;
								return false;
							} // if
							complex_t dn2(-2.0 * ((*s).second.grain_refindex().delta() -
											single_layer_refindex_.delta()),
											-2.0 * ((*s).second.grain_refindex().beta() -
											single_layer_refindex_.beta()));

							for(unsigned int z = 0; z < nqz_; ++ z) {
								for(unsigned int y = 0; y < nqy_; ++ y) {
									for(unsigned int x = 0; x < nqx_; ++ x) {
										unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
										unsigned int curr_index_0 = curr_index;
										unsigned int curr_index_1 = nqx_ * nqy_ * nqz_ + curr_index;
										unsigned int curr_index_2 = 2 * nqx_ * nqy_ * nqz_ + curr_index;
										unsigned int curr_index_3 = 3 * nqx_ * nqy_ * nqz_ + curr_index;

										base_id[curr_index] = dn2 *
											(amm[curr_index] * sf[curr_index_0] * ff[curr_index_0] +
											amp[curr_index] * sf[curr_index_1] * ff[curr_index_1] +
											apm[curr_index] * sf[curr_index_2] * ff[curr_index_2] +
											app[curr_index] * sf[curr_index_3] * ff[curr_index_3]);
									} // for x
								} // for y
							} // for z
						} else if(HiGInput::instance().structure_layer_order((*s).second) == -1) {
															// structure burried in substrate
							complex_t dn2(-2.0 * (substrate_refindex_.delta() -
											(*s).second.grain_refindex().delta()),
											-2.0 * (substrate_refindex_.beta() -
											(*s).second.grain_refindex().beta()));

							for(unsigned int z = 0; z < nqz_; ++ z) {
								for(unsigned int y = 0; y < nqy_; ++ y) {
									for(unsigned int x = 0; x < nqx_; ++ x) {
										unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
										base_id[curr_index] = dn2 * sf[curr_index] * ff[curr_index];
									} // for x
								} // for y
							} // for z
						} else if(HiGInput::instance().structure_layer_order((*s).second) == 0) {
															// structure on top of substrate
							complex_t dn2(-2.0 * (*s).second.grain_refindex().delta(),
											-2.0 * (*s).second.grain_refindex().beta());
	
							for(unsigned int z = 0; z < nqz_; ++ z) {
								for(unsigned int y = 0; y < nqy_; ++ y) {
									for(unsigned int x = 0; x < nqx_; ++ x) {
										unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
										unsigned int curr_index_0 = curr_index;
										unsigned int curr_index_1 = nqx_ * nqy_ * nqz_ + curr_index;
										unsigned int curr_index_2 = 2 * nqx_ * nqy_ * nqz_ + curr_index;
										unsigned int curr_index_3 = 3 * nqx_ * nqy_ * nqz_ + curr_index;
										base_id[curr_index] = dn2 *
											(h0[curr_index] * sf[curr_index_0] * ff[curr_index_0] +
											rk2[curr_index] * sf[curr_index_1] * ff[curr_index_1] +
											rk1[curr_index] * sf[curr_index_2] * ff[curr_index_2] +
											rk1rk2[curr_index] * sf[curr_index_3] * ff[curr_index_3]);
									} // for x
								} // for y
							} // for z
						} else {
							if(gmaster)
								std::cerr << "error: unable to determine sample structure. "
											<< "make sure the layer order is correct" << std::endl;
							return false;
						} // if-else
					} else {
						/* perform slicing */
						// not yet implemented ...
						std::cout << "uh-oh: ever thought about implementing the slicing scheme?"
									<< std::endl;
						return false;
					} // if-else
				} // if gmaster

				// clean everything before going to next
				//ff_.clear();
				//sf_.clear();
				ff.clear();
				sf.clear();

			} // for num_gr

			complex_t* id = NULL;
			#ifdef USE_MPI
				if(multi_node_.size(struct_comm) > 1) {
					// collect grain_ids from all procs in struct_comm
					if(smaster) {
						id = new (std::nothrow) complex_t[num_grains * nqx_ * nqy_ * nqz_];
					} // if
					int *proc_sizes = new (std::nothrow) int[multi_node_.size(struct_comm)];
					int *proc_displacements = new (std::nothrow) int[multi_node_.size(struct_comm)];
					multi_node_.gather(struct_comm, &num_gr, 1, proc_sizes, 1);
					if(smaster) {
						// make sure you get data only from the gmasters
						for(int i = 0; i < multi_node_.size(struct_comm); ++ i)
							proc_sizes[i] *= (gmasters[i] * nqx_ * nqy_ * nqz_);
						proc_displacements[0] = 0;
						for(int i = 1; i < multi_node_.size(struct_comm); ++ i)
							proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
					} // if
					multi_node_.gatherv(struct_comm, grain_ids, gmaster * num_gr * nqx_ * nqy_ * nqz_,
										id, proc_sizes, proc_displacements);
					delete[] proc_displacements;
					delete[] proc_sizes;
				} else {
					id = grain_ids;
				} // if-else

				delete[] gmasters;
				multi_node_.free(grain_comm);
			#else
				id = grain_ids;
			#endif

			delete[] nn;
			delete[] dd;

			if(smaster) {
				// FIXME: double check the following correlation stuff ...
				// new stuff for grain/ensemble correlation
				unsigned int soffset = 0;
				switch(HiGInput::instance().param_structcorrelation()) {
					case structcorr_null:	// default
					case structcorr_nGnE:	// no correlation
						// struct_intensity = sum_grain(abs(grain_intensity)^2)
						// intensity = sum_struct(struct_intensity)
						soffset = s_num * nqx_ * nqy_ * nqz_;
						for(unsigned int z = 0; z < nqz_; ++ z) {
							for(unsigned int y = 0; y < nqy_; ++ y) {
								for(unsigned int x = 0; x < nqx_; ++ x) {
									unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
									float_t sum = 0.0;
									for(int d = 0; d < num_grains; ++ d) {
										unsigned int id_index = d * nqx_ * nqy_ * nqz_ + curr_index;
										sum += id[id_index].real() * id[id_index].real() +
												id[id_index].imag() * id[id_index].imag();
									} // for d
									struct_intensity[soffset + curr_index] = sum;
								} // for x
							} // for y
						} // for z
						break;

					case structcorr_nGE:	// non corr grains, corr ensemble
						// grain_intensity = abs(sum_struct(struct_intensity))^2
						// intensity = sum_grain(grain_itensity)
						// TODO ...
						std::cerr << "uh-oh: this nGE correlation is not yet implemented" << std::endl;
						return false;
						break;

					case structcorr_GnE:	// corr grains, non corr ensemble
						// struct_intensity = abs(sum_grain(grain_intensity))^2
						// intensty = sum_struct(struct_intensity)
						soffset = s_num * nqx_ * nqy_ * nqz_;
						for(unsigned int z = 0; z < nqz_; ++ z) {
							for(unsigned int y = 0; y < nqy_; ++ y) {
								for(unsigned int x = 0; x < nqx_; ++ x) {
									unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
									complex_t sum(0.0, 0.0);
									for(int d = 0; d < num_grains; ++ d) {
										unsigned int id_index = d * nqx_ * nqy_ * nqz_ + curr_index;
										sum += id[id_index];
									} // for d
									struct_intensity[soffset + curr_index] = sum.real() * sum.real() +
																				sum.imag() * sum.imag();
								} // for x
							} // for y
						} // for z
						break;

					case structcorr_GE:		// both correlated
						// struct_intensity = sum_grain(grain_intensity)
						// intensity = abs(sum_struct(struct_intensity))^2
						soffset = s_num * nqx_ * nqy_ * nqz_;
						for(unsigned int z = 0; z < nqz_; ++ z) {
							for(unsigned int y = 0; y < nqy_; ++ y) {
								for(unsigned int x = 0; x < nqx_; ++ x) {
									unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
									complex_t sum(0.0, 0.0);
									for(int d = 0; d < num_grains; ++ d) {
										unsigned int id_index = d * nqx_ * nqy_ * nqz_ + curr_index;
										sum += id[id_index];
									} // for d
									c_struct_intensity[soffset + curr_index] = sum;
								} // for x
							} // for y
						} // for z
						break;

					default:				// error
						std::cerr << "error: unknown correlation type." << std::endl;
						return false;
				} // switch

				delete[] id;
			} // if master

		} // for num_structs

		float_t* all_struct_intensity = NULL;
		complex_t* all_c_struct_intensity = NULL;
		#ifdef USE_MPI
			if(multi_node_.size(comm_key) > 1) {
				// collect struct_intensity from all procs in comm_key
				if(master) {
					if(HiGInput::instance().param_structcorrelation() == structcorr_GE) {
						all_c_struct_intensity =
								new (std::nothrow) complex_t[num_structures_ * nqx_ * nqy_ * nqz_];
					} else {
						all_struct_intensity =
								new (std::nothrow) float_t[num_structures_ * nqx_ * nqy_ * nqz_];
					} // if-else
				} // if
				int *proc_sizes = new (std::nothrow) int[multi_node_.size(comm_key)];
				int *proc_displacements = new (std::nothrow) int[multi_node_.size(comm_key)];
				multi_node_.gather(comm_key, &num_structs, 1, proc_sizes, 1);
				if(master) {
					// make sure you get data only from the smasters
					for(int i = 0; i < multi_node_.size(comm_key); ++ i)
						proc_sizes[i] *= (smasters[i] * nqx_ * nqy_ * nqz_);
					proc_displacements[0] = 0;
					for(int i = 1; i < multi_node_.size(comm_key); ++ i)
						proc_displacements[i] = proc_displacements[i - 1] + proc_sizes[i - 1];
				} // if
				if(HiGInput::instance().param_structcorrelation() == structcorr_GE) {
					multi_node_.gatherv(comm_key, c_struct_intensity,
										smaster * num_structs * nqx_ * nqy_ * nqz_,
										all_c_struct_intensity, proc_sizes, proc_displacements);
				} else {
					multi_node_.gatherv(comm_key, struct_intensity,
										smaster * num_structs * nqx_ * nqy_ * nqz_,
										all_struct_intensity, proc_sizes, proc_displacements);
				} // if-else
				delete[] proc_displacements;
				delete[] proc_sizes;
			} else {
				all_struct_intensity = struct_intensity;
				all_c_struct_intensity = c_struct_intensity;
			} // if-else
		#else
			all_struct_intensity = struct_intensity;
			all_c_struct_intensity = c_struct_intensity;
		#endif

		#ifdef USE_MPI
			multi_node_.free(struct_comm);
			multi_node_.barrier(comm_key);
		#endif

		if(master) {
			img3d = new (std::nothrow) float_t[nqx_ * nqy_ * nqz_];
			// sum of all struct_intensity into intensity

			// new stuff for correlation
			switch(HiGInput::instance().param_structcorrelation()) {
				case structcorr_null:	// default
				case structcorr_nGnE:	// no correlation
					// struct_intensity = sum_grain(abs(grain_intensity)^2)
					// intensity = sum_struct(struct_intensity)
					for(unsigned int z = 0; z < nqz_; ++ z) {
						for(unsigned int y = 0; y < nqy_; ++ y) {
							for(unsigned int x = 0; x < nqx_; ++ x) {
								unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
								float_t sum = 0.0;
								for(int s = 0; s < num_structures_; ++ s) {
									unsigned int index = s * nqx_ * nqy_ * nqz_ + curr_index;
									sum += all_struct_intensity[index];
								} // for d
								img3d[curr_index] = sum;
							} // for x
						} // for y
					} // for z
					break;

				case structcorr_nGE:	// non corr grains, corr ensemble
					// grain_intensity = abs(sum_struct(struct_intensity))^2
					// intensity = sum_grain(grain_itensity)
					// TODO ...
					std::cerr << "uh-oh: this nGE correlation is not yet implemented" << std::endl;
					return false;
					break;

				case structcorr_GnE:	// corr grains, non corr ensemble
					// struct_intensity = abs(sum_grain(grain_intensity))^2
					// intensty = sum_struct(struct_intensity)
					for(unsigned int z = 0; z < nqz_; ++ z) {
						for(unsigned int y = 0; y < nqy_; ++ y) {
							for(unsigned int x = 0; x < nqx_; ++ x) {
								unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
								float_t sum = 0.0;
								for(int s = 0; s < num_structures_; ++ s) {
									unsigned int index = s * nqx_ * nqy_ * nqz_ + curr_index;
									sum += all_struct_intensity[index];
								} // for d
								img3d[curr_index] = sum;
							} // for x
						} // for y
					} // for z
					break;

				case structcorr_GE:		// both correlated
					// struct_intensity = sum_grain(grain_intensity)
					// intensity = abs(sum_struct(struct_intensity))^2
					for(unsigned int z = 0; z < nqz_; ++ z) {
						for(unsigned int y = 0; y < nqy_; ++ y) {
							for(unsigned int x = 0; x < nqx_; ++ x) {
								unsigned int curr_index = nqx_ * nqy_ * z + nqx_ * y + x;
								complex_t sum = 0.0;
								for(int s = 0; s < num_structures_; ++ s) {
									unsigned int index = s * nqx_ * nqy_ * nqz_ + curr_index;
									sum += all_c_struct_intensity[index];
								} // for d
								img3d[curr_index] = sum.real() * sum.real() + sum.imag() * sum.imag();
							} // for x
						} // for y
					} // for z
					break;

				default:				// error
					std::cerr << "error: unknown correlation type." << std::endl;
					return false;
			} // switch

			delete[] all_c_struct_intensity;
			delete[] all_struct_intensity;
		} // if master

		delete[] fc;

		return true;
	} // HipGISAXS::run_gisaxs()


	bool HipGISAXS::structure_factor(StructureFactor& sf,
									std::string expt, vector3_t& center, Lattice* &curr_lattice,
									vector3_t& grain_repeats, vector3_t& r_tot1,
									vector3_t& r_tot2, vector3_t& r_tot3
									#ifdef USE_MPI
										, const char* comm_key
									#endif
									) {
		#ifndef GPUSF
			return sf.compute_structure_factor(expt, center, curr_lattice, grain_repeats,
											r_tot1, r_tot2, r_tot3
											#ifdef USE_MPI
												, multi_node_, comm_key
											#endif
											);
		#else
			return sf.compute_structure_factor_gpu(expt, center, curr_lattice, grain_repeats,
											r_tot1, r_tot2, r_tot3
											#ifdef USE_MPI
												, multi_node_, comm_key
											#endif
											);
		#endif
	} // HipGISAXS::structure_factor()


	bool HipGISAXS::form_factor(FormFactor& ff,
								ShapeName shape_name, std::string shape_file,
								shape_param_list_t& shape_params, vector3_t &curr_transvec,
								float_t shp_tau, float_t shp_eta,
								vector3_t &r_tot1, vector3_t &r_tot2, vector3_t &r_tot3
								#ifdef USE_MPI
									, const char* comm_key
								#endif
								) {
		#ifndef GPUSF
			return ff.compute_form_factor(shape_name, shape_file, shape_params,
											single_layer_thickness_,
											curr_transvec, shp_tau, shp_eta, r_tot1, r_tot2, r_tot3
											#ifdef USE_MPI
												, multi_node_, comm_key
											#endif
											);
		#else
			return ff.compute_form_factor_gpu(shape_name, shape_file, shape_params,
											single_layer_thickness_,
											curr_transvec, shp_tau, shp_eta, r_tot1, r_tot2, r_tot3
											#ifdef USE_MPI
												, multi_node_, comm_key
											#endif
											);
		#endif
	} // HipGISAXS::form_factor()


	bool HipGISAXS::compute_rotation_matrix_z(float_t angle,
												vector3_t& r1, vector3_t& r2, vector3_t& r3) {
		float_t s = sin(angle);
		float_t c = cos(angle);

		r1[0] = c; r1[1] = -s; r1[2] = 0.0;
		r2[0] = s; r2[1] = c; r2[2]  = 0.0;
		r3[0] = 0.0; r3[1] = 0.0; r3[2]  = 1.0;

		return true;
	} // HipGISAXS::comptue_rotation_matrix_z()


	bool HipGISAXS::compute_rotation_matrix_y(float_t angle,
												vector3_t& r1, vector3_t& r2, vector3_t& r3) {
		float_t s = sin(angle);
		float_t c = cos(angle);

		r1[0] = c; r1[1] = 0.0; r1[2] = -s;
		r2[0] = 0.0; r2[1] = 1.0; r2[2]  = 0.0;
		r3[0] = s; r3[1] = 0.0; r3[2]  = c;

		return true;
	} // HipGISAXS::compute_rotation_matrix_y()


	bool HipGISAXS::compute_rotation_matrix_x(float_t angle,
												vector3_t& r1, vector3_t& r2, vector3_t& r3) {
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
			if(HiGInput::instance().is_single_layer()) {
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
				// TODO: not currently handled ...
				std::cerr << "uh-oh: this case (num_layers > 1) has not yet been implemented yet"
							<< std::endl;
				std::cerr << "go get yourself a nice cup of yummmy hot chocolate instead!" << std::endl;
				return false;
			} // if-else
			vol_[0] = vol_[1] = spot_diameter;
			vol_[2] = (penetration_depth_layer.real() < c_max_depth.real()) ?
								penetration_depth_layer.real() : c_max_depth.real();
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
					complex_t* &h0, complex_t* &fc) {
		amm = apm = amp = app = NULL;
		rk1 = rk2 = rk1rk2 = tk1tk2 = NULL;
		h0 = NULL;

		// doesnt this also depend on the type of polarization of the light? ... ?

		if(HiGInput::instance().param_nslices() <= 1) {	/* computing without sample slicing */
			complex_t dnl2 = 2.0 * complex_t(single_layer_refindex_.delta(),
											single_layer_refindex_.beta());
			if(!layer_qgrid_qz(alpha_i, dnl2)) {
				std::cerr << "error: could not compute extended qz" << std::endl;
				return false;
			} // if

			if(HiGInput::instance().struct_in_layer() && HiGInput::instance().num_layers() == 1) {
				/* compute fresnel coefficients for a structure
				 * embedded inside a layer on a substrate */
				if(!compute_fresnel_coefficients_embedded(alpha_i, fc)) {
					std::cerr << "error: could not compute fresnel coefficients" << std::endl;
					return false;
				} // if
				// set aliases
				amm = fc;
				apm = fc + nqx_ * nqy_ * nqz_;
				amp = fc + 2 * nqx_ * nqy_ * nqz_;
				app = fc + 3 * nqx_ * nqy_ * nqz_;
				h0 = fc + 4 * nqx_ * nqy_ * nqz_;
			} else if(HiGInput::instance().num_layers() == 1 || HiGInput::instance().num_layers() == 0) {
				/* compute fresnel coefficients for a structure
				 * on top of or buried inside a substrate */
				if(!compute_fresnel_coefficients_top_buried(alpha_i, fc)) {
					std::cerr << "error: could not compute fresnel coefficients" << std::endl;
					return false;
				} // if
				// set aliases
				rk1 = fc;							// re-check these ...
				rk2 = fc + nqx_ * nqy_ * nqz_;
				rk1rk2 = fc + 2 * nqx_ * nqy_ * nqz_;
				tk1tk2 = fc + 3 * nqx_ * nqy_ * nqz_;
				h0 = fc + 4 * nqx_ * nqy_ * nqz_;
			} else {
				std::cerr << "error: invalid number of layers in non-slicing scheme" << std::endl;
				return false;
			} // if-else
		} else {
			/* computing using the sample slicing scheme */
			// TODO: not yet implemented ...
			std::cerr << "uh-oh: you hit a slicing part that has not been implemented yet" << std::endl;
			return false;
		} // if

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


	// TODO optimize this later ...
	bool HipGISAXS::compute_fresnel_coefficients_embedded(float_t alpha_i, complex_t* &fc) {
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
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
			a1mi[i] = a1m_kiz1; a1pi[i] = a1p_kiz1;
		} // for

		complex_t *a1mf = NULL, *a1pf = NULL;
		a1mf = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];
		a1pf = new (std::nothrow) complex_t[nqx_ * nqy_ * nqz_];
		if(a1mf == NULL || a1pf == NULL) {
			std::cerr << "error: failed to allocate memory for a1mf, a1pf" << std::endl;
			return false;
		} // if
		for(unsigned int i = 0; i < nqx_ * nqy_ * nqz_; ++ i) {
			a1mf[i] = a1pf[i] = complex_t(0.0, 0.0);
		} // for

		// allocate fc memory
		fc = new (std::nothrow) complex_t[5 * nqx_ * nqy_ * nqz_];	// 5 sets
		if(fc == NULL) {
			std::cerr << "error: failed to allocate memory for fc" << std::endl;
			return false;
		} // if

		unsigned int nqxyz = nqx_ * nqy_ * nqz_;
		unsigned int nqxy = nqx_ * nqy_;

		for(unsigned int z = 0; z < nqz_; ++ z) {
			complex_t a1m_nkfz1, a1p_nkfz1;
			float_t kfz0 = QGrid::instance().qz(z) + kiz0;
			unsigned int temp0 = nqxy * z;
			unsigned int temp1 = 4 * nqxyz + temp0;

			if(kfz0 < 0) {
				a1m_nkfz1 = complex_t(0.0, 0.0);
				a1p_nkfz1 = complex_t(0.0, 0.0);

				for(unsigned int y = 0; y < nqy_; y ++) {
					unsigned int temp2 = temp1 + nqx_ * y;
					for(unsigned int x = 0; x < nqx_; x ++) {
						fc[temp2 + x] = complex_t(0.0, 0.0);
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

				for(unsigned int y = 0; y < nqy_; y ++) {
					unsigned int temp2 = temp1 + nqx_ * y;
					for(unsigned int x = 0; x < nqx_; x ++) {
						fc[temp2 + x] = complex_t(1.0, 0.0);
					} // for x
				} // for y
			} // if-else

			for(unsigned int y = 0; y < nqy_; y ++) {			// these can be aliminated ...
				unsigned int temp2 = temp0 + nqx_ * y;
				for(unsigned int x = 0; x < nqx_; x ++) {
					a1mf[temp2 + x] = a1m_nkfz1;
					a1pf[temp2 + x] = a1p_nkfz1;
				} // for x
			} // for y
		} // for z

		// the element-element products
		for(unsigned int z = 0; z < nqz_; z ++) {
			for(unsigned int y = 0; y < nqy_; y ++) {
				for(unsigned int x = 0; x < nqx_; x ++) {
					fc[0 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1mi[z * nqx_ * nqy_ + y * nqx_ + x] * a1mf[z * nqx_ * nqy_ + y * nqx_ + x];
					fc[1 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1pi[z * nqx_ * nqy_ + y * nqx_ + x] * a1mf[z * nqx_ * nqy_ + y * nqx_ + x];
					fc[2 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
						a1mi[z * nqx_ * nqy_ + y * nqx_ + x] * a1pf[z * nqx_ * nqy_ + y * nqx_ + x];
					fc[3 * nqx_ * nqy_ * nqz_ + z * nqx_ * nqy_ + y * nqx_ + x] =
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
	bool HipGISAXS::compute_fresnel_coefficients_top_buried(float_t alpha_i, complex_t* &fc) {
		//complex_t tk1 = ((float_t) 2.0 * sin(alpha_i)) / ((float_t)sin(alpha_i) +
		//										sqrt((float_t)pow(sin(alpha_i), 2) - dns2_));
		complex_t tk1 = ((float_t) (2.0 * sin(alpha_i))) /
						((float_t) sin(alpha_i) + sqrt((float_t)pow(sin(alpha_i), 2) - dns2_));
		complex_t rk1 = tk1 - complex_t(1.0, 0.0);

		fc = new (std::nothrow) complex_t[5 * nqx_ * nqy_ * nqz_];
		if(fc == NULL) {
			std::cerr << "error: failed to allocate memory for fc" << std::endl;
			return false;
		} // if
		for(unsigned int z = 0; z < nqz_; z ++) {
			for(unsigned int y = 0; y < nqy_; y ++) {
				for(unsigned int x = 0; x < nqx_; x ++) {
					fc[z * nqy_ * nqx_ + y * nqx_ + x] = rk1;
				} // for x
			} // for y
		} // for z

		float_t k1z = -1.0 * k0_ * sin(alpha_i);
		complex_t tk2, rk2;
		unsigned long int nqxyz = nqx_ * nqy_ * nqz_;
		for(unsigned int z = 0; z < nqz_; ++ z) {
			complex_t k2z = QGrid::instance().qz(z) + k1z;
			if(k2z < 0) {
				tk2 = complex_t(0.0, 0.0);
				rk2 = complex_t(0.0, 0.0);
				for(unsigned int y = 0; y < nqy_; ++ y) {
					for(unsigned int x = 0; x < nqx_; ++ x) {
						fc[0 * nqxyz + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(0.0, 0.0);
						fc[4 * nqxyz + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(0.0, 0.0);
					} // for x
				} // for y
			} else {
				complex_t sf = k2z / k0_;
				tk2 = 2.0 * sf / (sf + sqrt(sf * sf - dns2_));
				rk2 = tk2 - (float_t) 1.0;
				for(unsigned int y = 0; y < nqy_; y ++) {
					for(unsigned int x = 0; x < nqx_; x ++) {
						fc[4 * nqxyz + z * nqx_ * nqy_ + y * nqx_ + x] = complex_t(1.0, 0.0);
					} // for x
				} // for y
			} // if-else

			for(unsigned int y = 0; y < nqy_; y ++) {
				for(unsigned int x = 0; x < nqx_; x ++) {
					fc[1 * nqxyz + z * nqx_ * nqy_ + y * nqx_ + x] = rk2;
					fc[2 * nqxyz + z * nqx_ * nqy_ + y * nqx_ + x] = rk1 * rk2;
					fc[3 * nqxyz + z * nqx_ * nqy_ + y * nqx_ + x] = tk1 * tk2;
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
				// find max density - number of grains in vol
				vector3_t max_density = min(floor(vol_ / spaced_cell) + 1, maxgrains);
				rand_dim_x = max(1,
							(int)std::floor(max_density[0] * max_density[1] * max_density[2] / 4.0));
				rand_dim_y = 3;

				// construct random matrix
				float_t *d_rand = new (std::nothrow) float_t[rand_dim_x * rand_dim_y];
				srand(time(NULL));
				for(int i = 0; i < rand_dim_x * rand_dim_y; ++ i)
					d_rand[i] = ((float_t)rand() / RAND_MAX);

				d = new (std::nothrow) float_t[rand_dim_x * rand_dim_y * 4];

				//int base_index = 0;
				float_t mul_val1 = vol_[0] / 2;
				float_t mul_val2 = vol_[1] / 2;
				float_t mul_val3 = vol_[2];

				// D0
				if(maxgrains[0] == 1) {
					for(int x = 0; x < rand_dim_x; ++ x) d[x] = 0;
				} else {
					for(int x = 0; x < rand_dim_x; ++ x)
						d[x] = d_rand[x] * mul_val1;
				} // if-else
				if(maxgrains[1] == 1) {
					for(int x = 0; x < rand_dim_x; ++ x) d[4 * rand_dim_x + x] = 0;
				} else {
					for(int x = 0; x < rand_dim_x; ++ x)
						d[4 * rand_dim_x + x] = d_rand[rand_dim_x + x] * mul_val2;
				} // if-else
				if(maxgrains[2] == 1) {
					for(int x = 0; x < rand_dim_x; ++ x) d[8 * rand_dim_x + x] = tz;
				} else {
					for(int x = 0; x < rand_dim_x; ++ x)
						d[8 * rand_dim_x + x] = d_rand[2 * rand_dim_x + x] * mul_val3 + tz;
				} // if-else
				// D1
				for(int x = 0; x < rand_dim_x; ++ x)
					d[rand_dim_x + x] = - d[x];
				for(int x = 0; x < rand_dim_x; ++ x)
					d[4 * rand_dim_x + rand_dim_x + x] = d[4 * rand_dim_x + x];
				for(int x = 0; x < rand_dim_x; ++ x)
					d[8 * rand_dim_x + rand_dim_x + x] = d[8 * rand_dim_x + x];
				// D2
				for(int x = 0; x < rand_dim_x; ++ x)
					d[2 * rand_dim_x + x] = d[x];
				for(int x = 0; x < rand_dim_x; ++ x)
					d[4 * rand_dim_x + 2 * rand_dim_x + x] = - d[4 * rand_dim_x + x];
				for(int x = 0; x < rand_dim_x; ++ x)
					d[8 * rand_dim_x + 2 * rand_dim_x + x] = d[8 * rand_dim_x + x];
				// D3
				for(int x = 0; x < rand_dim_x; ++ x)
					d[3 * rand_dim_x + x] = - d[x];
				for(int x = 0; x < rand_dim_x; ++ x)
					d[4 * rand_dim_x + 3 * rand_dim_x + x] = - d[4 * rand_dim_x + x];
				for(int x = 0; x < rand_dim_x; ++ x)
					d[8 * rand_dim_x + 3 * rand_dim_x + x] = d[8 * rand_dim_x + x];

				// fix the rand_dim_x
				rand_dim_x *= 4;

			} else if(dim == 2) {
				std::cerr << "error: dim == 2 case not implemented" << std::endl;
				// ...
				return false;
			} else if(dim == 1) {
				std::cerr << "error: dim == 1 case not implemented" << std::endl;
				// ...
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
				// ...
				return false;
			} else if(dim == 1) {
				std::cerr << "error: dim == 1 case not implemented" << std::endl;
				// ...
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
		//vector2_t tau = (*s).second.rotation_tau();
		//vector2_t eta = (*s).second.rotation_eta();
		//vector2_t zeta = (*s).second.rotation_zeta();
		vector3_t rot1 = (*s).second.rotation_rot1();
		vector3_t rot2 = (*s).second.rotation_rot2();
		vector3_t rot3 = (*s).second.rotation_rot3();

		nn = new (std::nothrow) float_t[ndx * ndy];
											// TODO i believe constructing nn may not be needed ...
		if(distribution == "single") {		// single
			for(int x = 0; x < ndx; ++ x) {
				//nn[x] = tau[0] * PI_ / 180;
				nn[x] = rot1[1] * PI_ / 180;
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				//nn[ndx + x] = eta[0] * PI_ / 180;
				nn[ndx + x] = rot2[1] * PI_ / 180;
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				//nn[2 * ndx + x] = zeta[0] * PI_ / 180;
				nn[2 * ndx + x] = rot3[1] * PI_ / 180;
			} // for x
		} else if(distribution == "random") {	// random
			for(int x = 0; x < 3 * ndx; ++ x) {
				nn[x] = (float_t(rand()) / RAND_MAX) * 2 * PI_;
			} // for x
		} else if(distribution == "range") {	// range
			//float_t dtau = fabs(tau[1] - tau[0]);
			//float_t deta = fabs(eta[1] - eta[0]);
			//float_t dzeta = fabs(zeta[1] - zeta[0]);
			float_t drot1 = fabs(rot1[2] - rot1[1]);
			float_t drot2 = fabs(rot2[2] - rot2[1]);
			float_t drot3 = fabs(rot3[2] - rot3[1]);
			for(int x = 0; x < ndx; ++ x) {
				//nn[x] = (tau[0] + (float_t(rand()) / RAND_MAX) * dtau) * PI_ / 180;
				nn[x] = (rot1[1] + (float_t(rand()) / RAND_MAX) * drot1) * PI_ / 180;
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				//nn[ndx + x] = (eta[0] + (float_t(rand()) / RAND_MAX) * deta) * PI_ / 180;
				nn[ndx + x] = (rot2[1] + (float_t(rand()) / RAND_MAX) * drot2) * PI_ / 180;
			} // for x
			for(int x = 0; x < ndx; ++ x) {
				//nn[2 * ndx + x] = (zeta[0] + (float_t(rand()) / RAND_MAX) * dzeta) * PI_ / 180;
				nn[2 * ndx + x] = (rot3[1] + (float_t(rand()) / RAND_MAX) * drot3) * PI_ / 180;
			} // for x
			/*float_t da = PI_ / (2 * (ndx - 1));
			for(int x = 0; x < ndx; ++ x) {
				nn[x] = x * da;
			} // for x
			for(int x = 0; x < 2 * ndx; ++ x) nn[ndx + x] = 0;*/
		} else {
			// TODO read .ori file ...
			std::cerr << "uh-oh: I guess you wanted to read orientations from a file" << std::endl;
			std::cerr << "too bad, its not implemented yet" << std::endl;
			return false;
		} // if-else

		return true;
	} // HipGISAXS::orientation_distribution()


	/**
	 * fitting related functions
	 */

	bool HipGISAXS::update_params(const map_t& params) {
		return HiGInput::instance().update_params(params);
	} // HipGISAXS::update_params()


	/**
	 * miscellaneous functions
	 */

	void HipGISAXS::save_gisaxs(float_t *final_data, std::string output) {
		std::ofstream f(output.c_str());
		for(unsigned int z = 0; z < nqz_; ++ z) {
			for(unsigned int y = 0; y < nqy_; ++ y) {
				unsigned int index = nqy_ * z + y;
				f << final_data[index] << "\t";
			} // for
			f << std::endl;
		} // for
		f.close();
	} // HipGISAXS::save_gisaxs()
	
	
	void HipGISAXS::printfr(const char* name, float_t* arr, unsigned int size) {
		std::cout << name << ":" << std::endl;
		if(arr == NULL) { std::cout << "NULL" << std::endl; return; }
		for(unsigned int i = 0; i < size; ++ i) {
			std::cout << arr[i] << "\t";
		} // for
		std::cout << std::endl;
	} // HipGISAXS::printfr()
 

	void HipGISAXS::printfc(const char* name, complex_t* arr, unsigned int size) {
		std::cout << name << ":" << std::endl;
		if(arr == NULL) { std::cout << "NULL" << std::endl; return; }
		for(unsigned int i = 0; i < size; ++ i) {
			std::cout << arr[i].real() << "," << arr[i].imag() << "\t";
		} // for
		std::cout << std::endl;
	} // HipGISAXS::printfc()
 

	// for testing
	bool HipGISAXS::write_qgrid(char* filename) {
		std::ofstream qout(filename);

		qout << nqx_ << " " << nqy_ << " " << nqz_extended_ << std::endl;
		for(unsigned int i = 0; i < nqx_; ++ i) {
			qout << QGrid::instance().qx(i) << " ";
		} // for
		qout << std::endl;
		for(unsigned int i = 0; i < nqy_; ++ i) {
			qout << QGrid::instance().qy(i) << " ";
		} // for
		qout << std::endl;
		for(unsigned int i = 0; i < nqz_extended_; ++ i) {
			qout << QGrid::instance().qz_extended(i).real() << " "
					<< QGrid::instance().qz_extended(i).imag() << " ";
		} // for
		qout << std::endl;

		qout.close();
		return true;
	} // HipGISAXS::write_qgrid()


	bool HipGISAXS::read_form_factor(FormFactor& ff, const char* filename) {
		return ff.read_form_factor(filename, nqx_, nqy_, nqz_extended_);
	} // HipGISAXS::read_form_factor()


} // namespace hig
