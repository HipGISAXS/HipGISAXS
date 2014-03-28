/**
 *  Project:
 *
 *  File: fit_pso.cpp
 *  Created: Jan 13, 2014
 *  Modified: Sun 23 Mar 2014 12:42:07 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>
#include <limits>
#include <ctime>

#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <hipgisaxs.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

namespace hig {

	//ParticleSwarmOptimization::ParticleSwarmOptimization(int narg, char** args, ObjectiveFunction* obj) :
	ParticleSwarmOptimization::ParticleSwarmOptimization(int narg, char** args, ObjectiveFunction* obj,
			float_t omega, float_t phi1, float_t phi2, int npart, int ngen,
			bool tune_omega = false, bool foresee = false) :
   			rand_(time(NULL))	{
		name_ = algo_pso;
		max_hist_ = 100;			// not used ...
		tol_ = 1e-8;

		foresee_num_ = 5;			// TODO: make it modifiable

		pso_omega_ = omega;
		pso_phi1_ = phi1;
		pso_phi2_ = phi2;
		obj_func_ = obj;

		tune_omega_ = tune_omega;
		foresee_ = foresee;

		std::cout << "** PSO Parameters: " << pso_omega_ << ", " << pso_phi1_ << ", " << pso_phi2_
					<< std::endl;

		params_ = (*obj_func_).fit_param_keys();
		num_params_ = (*obj_func_).num_fit_params();
		x0_ = (*obj_func_).fit_param_init_values();
		#ifdef USE_MPI
			multi_node_ = (*obj_func_).multi_node_comm();
			root_comm_ = (*multi_node_).universe_key();
			std::cout << "** MPI Size: " << (*multi_node_).size(root_comm_);
			std::cout << ", my rank: " << (*multi_node_).rank(root_comm_) << std::endl;
			rand_.reset(time(NULL) + (*multi_node_).rank(root_comm_));
		#else
			root_comm_ = "world";
		#endif

		best_values_.resize(num_params_, 0.0);
		best_fitness_ = std::numeric_limits<float_t>::max();

		std::vector <std::pair <float_t, float_t> > limits = obj_func_->fit_param_limits();
		parameter_data_list_t max_limits, min_limits;
		for(std::vector <std::pair <float_t, float_t> >::iterator iter = limits.begin();
				iter != limits.end(); ++ iter) {
			float_t min_val = iter->first;
			float_t max_val = iter->second;
			if(max_val < min_val) std::swap(min_val, max_val);
			min_limits.push_back(min_val);
			max_limits.push_back(max_val);
		} // for
		if(min_limits.size() < num_params_) min_limits.resize(num_params_, 0.0);
		if(max_limits.size() < num_params_) max_limits.resize(num_params_,
																std::numeric_limits<float_t>::max());
		constraints_.param_values_min_ = min_limits;
		constraints_.param_values_max_ = max_limits;

		// set constraints
		// for now let min max velocity be -max to max param (assuming all params are +ve) ...
		for(int i = 0; i < num_params_; ++ i) min_limits[i] = - max_limits[i];
		constraints_.velocity_min_ = min_limits;
		constraints_.velocity_max_ = max_limits;

		// get number of particles
		num_particles_ = npart;
		// max number of generations
		max_iter_ = ngen;

		#ifdef USE_MPI
			int num_procs = (*multi_node_).size(root_comm_);
			int rank = (*multi_node_).rank(root_comm_);
			// two cases:	num_procs <= num_part (easy case)
			// 				num_procs > num_part
			particle_comm_ = "pso_particle";
			if(num_procs <= num_particles_) {
				// assign multiple particles per proc
				num_particles_ = num_particles_ / num_procs +
										(rank < (num_particles_ % num_procs) ? 1 : 0);
				// each proc will now be in its own new world
				(*multi_node_).split(particle_comm_, root_comm_, rank);
			} else {
				// assign multiple procs per particle
				// procs for a particle will be in their common new world
				// TODO ...
				std::cout << "error: this case has not been implemented yet" << std::endl;
				exit(-1);
			} // if-else
		#else
			particle_comm_ = root_comm_;
		#endif

		// initialize particles
		particles_.clear();
		for(int i = 0; i < num_particles_; ++ i) {
			particles_.push_back(PSOParticle(num_params_, PSO_UNIFORM, rand_, constraints_));
		} // for
	} // ParticleSwarmOptimization::ParticleSwarmOptimization()


	ParticleSwarmOptimization::~ParticleSwarmOptimization() {

	} // ParticleSwarmOptimization::~ParticleSwarmOptimization()


	parameter_map_t ParticleSwarmOptimization::get_best_values() const {
		parameter_map_t values;
		for(int i = 0; i < num_params_; ++ i)
			values[params_[i]] = best_values_[i];
		return values;
	} // ParticleSwarmOptimization::get_best_values()


	bool ParticleSwarmOptimization::run(int narg, char** args, int img_num) {
		std::cout << "Running Particle Swarm Optimization ..." << std::endl;

		(*obj_func_).set_reference_data(img_num);

		woo::BoostChronoTimer gen_timer;
		double total_time = 0;

		for(int gen = 0; gen < max_iter_; ++ gen) {
			std::cout << (*multi_node_).rank(root_comm_) << ": Generation " << gen << std::endl;
			gen_timer.start();
			if(foresee_) {
				if(!simulate_soothsayer_generation()) {
					std::cerr << "error: failed in generation " << gen << std::endl;
					return false;
				} // if
			} else {
				if(!simulate_generation()) {
					std::cerr << "error: failed in generation " << gen << std::endl;
					return false;
				} // if
			} // if-else
			gen_timer.stop();
			double elapsed = gen_timer.elapsed_msec();
			total_time += elapsed;
			#ifdef USE_MPI
			if((*multi_node_).is_master(root_comm_))
			#endif
				std::cout << "@@@@@@ Generation time: " << elapsed << " ms. [total: "
							<< total_time << " ms.]" << std::endl;
		} // for
		
		// set the final values
		xn_ = best_values_;

		return true;
	} // ParticleSwarmOptimization::run()


	bool ParticleSwarmOptimization::simulate_generation() {
		// for each particle, simulate
		for(int i = 0; i < num_particles_; ++ i) {
			std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
			// construct param map
			//parameter_map_t curr_particle;
			float_vec_t curr_particle;
			for(int j = 0; j < num_params_; ++ j)
				curr_particle.push_back(particles_[i].param_values_[j]);
				//curr_particle[params_[j]] = particles_[i].param_values_[j];
			// compute the fitness
			(*obj_func_).update_sim_comm(particle_comm_);
			float_vec_t curr_fitness = (*obj_func_)(curr_particle);

			// update particle fitness
			if(particles_[i].best_fitness_ > curr_fitness[0]) {
				particles_[i].best_fitness_ = curr_fitness[0];
				particles_[i].best_values_ = particles_[i].param_values_;
			} // if
			// update global best fitness locally
			if(best_fitness_ > curr_fitness[0]) {
				best_fitness_ = curr_fitness[0];
				best_values_ = particles_[i].param_values_;
			} // if
			//std::cout << "@@@@ " << i + 1 << ": " << curr_fitness[0] << " "
			//	<< particles_[i].best_fitness_ << " [ ";
			//for(int j = 0; j < num_params_; ++ j) std::cout << particles_[i].param_values_[j] << " ";
			//std::cout << "]" << std::endl;
		} // for

		#ifdef USE_MPI
			// communicate and find globally best fitness
			float_t new_best_fitness; int best_rank;
			if((*multi_node_).size("pso_particle") > 1 && num_particles_ == 1) {
				// in the case when multiple procs work on one particle,
				// only the master needs to communicate with other masters (for other particles)
				// TODO ...
				bool pmaster = (*multi_node_).is_master(particle_comm_);
				std::cout << "error: this case has not been implemented yet" << std::endl;
				return false;
			} else {
				if(!(*multi_node_).allreduce(root_comm_, best_fitness_, new_best_fitness, best_rank,
												woo::comm::minloc))
					return false;
				best_fitness_ = new_best_fitness;
			} // if-else
			// the best rank proc broadcasts its best values to all others
			if(!(*multi_node_).broadcast(root_comm_, best_values_, best_rank)) return false;
		#endif

		// update each particle using new global best
		for(int i = 0; i < num_particles_; ++ i) {
			if(!particles_[i].update_particle(pso_omega_, pso_phi1_, pso_phi2_, best_values_,
											constraints_, rand_)) {
				std::cerr << "error: failed to update particle " << i << std::endl;
				return false;
			} // if
		} // for

		if(tune_omega_) pso_omega_ /= 1.5;
		std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

		/*std::cout << "~~~~ Generation best: ";
		for(int i = 0; i < num_particles_; ++ i) {
			std::cout << "[ ";
			for(int j = 0; j < num_params_; ++ j)
				std::cout << particles_[i].param_values_[j] << " ";
			std::cout << "]\t";
		} // for
		std::cout << std::endl;*/

		#ifdef USE_MPI
		if((*multi_node_).is_master(root_comm_)) {
		#endif
			std::cout << "@@@@@@ Global best: ";
			std::cout << best_fitness_ << " [ ";
			for(int j = 0; j < num_params_; ++ j)
				std::cout << best_values_[j] << " ";
			std::cout << "]\t";
			std::cout << std::endl;
		#ifdef USE_MPI
		} // if
		#endif

		// write out stuff to output files
		/*std::string prefix(HiGInput::instance().param_pathprefix()+"/"+HiGInput::instance().runname());
		for(int i = 0; i < num_particles_; ++ i) {
			std::stringstream outss; outss << "/pso_params_particle_" << i << ".dat";
			std::ofstream out(prefix + outss.str(), std::ios::app);
			for(int j = 0; j < num_params_; ++ j)
				out << particles_[i].param_values_[j] << "\t";
			out << std::endl;
			out.close();
			std::stringstream out2ss; out2ss << "/pso_fitness_particle_" << i << ".dat";
			std::ofstream out2(prefix + out2ss.str(), std::ios::app);
			out2 << particles_[i].best_fitness_ << std::endl;
			out2.close();
		} // for
		std::ofstream out(prefix + "/pso_params_best.dat", std::ios::app);
		for(int j = 0; j < num_params_; ++ j)
			out << best_values_[j] << "\t";
		out << std::endl;
		out.close();
		std::ofstream out2(prefix + "/pso_fitness_best.dat", std::ios::app);
		out2 << best_fitness_ << std::endl;
		out2.close();*/

		return true;
	} // ParticleSwarmOptimization::simulate_generation()


	bool ParticleSwarmOptimization::simulate_soothsayer_generation() {
		// for each particle, simulate
		for(int i = 0; i < num_particles_; ++ i) {
			std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
			// first save the current state of the particle
			parameter_data_list_t orig_param_vals = particles_[i].param_values_;
			parameter_data_list_t orig_velocity = particles_[i].velocity_;
			float_t foresee_best_fitness = std::numeric_limits<float_t>::max();
			parameter_data_list_t foresee_best_params, foresee_best_velocity;
			// for each forseeing step
			for(int k = 0; k < foresee_num_; ++ k) {
				// create a new position from the actual current position
				particles_[i].compute_and_set_values(orig_param_vals, orig_velocity,
														pso_omega_, pso_phi1_, pso_phi2_,
														best_values_, constraints_, rand_);
				// compute the fitness
				float_vec_t curr_particle = particles_[i].param_values_;
				(*obj_func_).update_sim_comm(particle_comm_);
				float_vec_t curr_fitness = (*obj_func_)(curr_particle);

				if(foresee_best_fitness > curr_fitness[0]) {
					foresee_best_fitness = curr_fitness[0];
					foresee_best_params = curr_particle;
					foresee_best_velocity = particles_[i].velocity_;
				} // if
			} // if
			// update particle
			particles_[i].param_values_ = foresee_best_params;
			particles_[i].velocity_ = foresee_best_velocity;
			if(particles_[i].best_fitness_ > foresee_best_fitness) {
				particles_[i].best_fitness_ = foresee_best_fitness;
				particles_[i].best_values_ = foresee_best_params;
			} // if
			// update global best fitness locally
			if(best_fitness_ > foresee_best_fitness) {
				best_fitness_ = foresee_best_fitness;
				best_values_ = foresee_best_params;
			} // if
		} // for

		#ifdef USE_MPI
			// communicate and find globally best fitness
			float_t new_best_fitness; int best_rank;
			if((*multi_node_).size("pso_particle") > 1 && num_particles_ == 1) {
				// in the case when multiple procs work on one particle,
				// only the master needs to communicate with other masters (for other particles)
				// TODO ...
				bool pmaster = (*multi_node_).is_master(particle_comm_);
				std::cout << "error: this case has not been implemented yet" << std::endl;
				return false;
			} else {
				if(!(*multi_node_).allreduce(root_comm_, best_fitness_, new_best_fitness, best_rank,
												woo::comm::minloc))
					return false;
				best_fitness_ = new_best_fitness;
			} // if-else
			// the best rank proc broadcasts its best values to all others
			if(!(*multi_node_).broadcast(root_comm_, best_values_, best_rank)) return false;
		#endif

		if(tune_omega_) pso_omega_ /= 2.0;
		std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

		#ifdef USE_MPI
		if((*multi_node_).is_master(root_comm_)) {
		#endif
			std::cout << "@@@@@@ Global best:\t";
			std::cout << best_fitness_ << "\t[ ";
			for(int j = 0; j < num_params_; ++ j) std::cout << best_values_[j] << " ";
			std::cout << "]" << std::endl;
		#ifdef USE_MPI
		} // if
		#endif

		return true;
	} // ParticleSwarmOptimization::simulate_soothsayer_generation()


} // namespace hig
