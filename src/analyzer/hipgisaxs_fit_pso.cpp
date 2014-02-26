/**
 *  Project:
 *
 *  File: fit_pso.cpp
 *  Created: Jan 13, 2014
 *  Modified: Wed 26 Feb 2014 01:39:49 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>
#include <limits>
#include <ctime>

#include <analyzer/hipgisaxs_fit_pso.hpp>

namespace hig {

	//ParticleSwarmOptimization::ParticleSwarmOptimization(int narg, char** args, ObjectiveFunction* obj) :
	ParticleSwarmOptimization::ParticleSwarmOptimization(int narg, char** args, ObjectiveFunction* obj,
			float_t omega, float_t phi1, float_t phi2, int npart, int ngen) :
   			rand_(time(NULL))	{
		name_ = algo_pso;
		max_hist_ = 100;			// not used ...
		tol_ = 1e-8;

		//if(narg < 7) {
		//	std::cerr << "error: for PSO please specify <num_particles> <num_gen> "
		//				<< "<omega> <phi1> <phi2>" << std::endl;
		//	exit(1);
		//} // if

		//pso_omega_ = atof(args[4]);
		//pso_phi1_ = atof(args[5]);
		//pso_phi2_ = atof(args[6]);
		pso_omega_ = omega;
		pso_phi1_ = phi1;
		pso_phi2_ = phi2;
		obj_func_ = obj;

		std::cout << "** PSO Parameters: " << pso_omega_ << ", " << pso_phi1_ << ", " << pso_phi2_
					<< std::endl;

		params_ = (*obj_func_).fit_param_keys();
		num_params_ = (*obj_func_).num_fit_params();
		x0_ = (*obj_func_).fit_param_init_values();

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
		//num_particles_ = atoi(args[2]);
		num_particles_ = npart;
		// max number of generations
		//max_iter_ = atoi(args[3]);
		max_iter_ = ngen;

		// initialize particles
		particles_.clear();
		for(int i = 0; i < num_particles_; ++ i) {
			particles_.push_back(PSOParticle(num_params_, PSO_UNIFORM, rand_, constraints_));
		} // for
	} // ParticleSwarmOptimization::ParticleSwarmOptimization()


	ParticleSwarmOptimization::~ParticleSwarmOptimization() {

	} // ParticleSwarmOptimization::~ParticleSwarmOptimization()


	bool ParticleSwarmOptimization::run(int narg, char** args, int img_num) {
		std::cout << "Running Particle Swarm Optimization ..." << std::endl;

		(*obj_func_).set_reference_data(img_num);

		for(int gen = 0; gen < max_iter_; ++ gen) {
			if(!simulate_generation()) {
				std::cerr << "error: failed in generation " << gen << std::endl;
				return false;
			} // if
		} // for
		
		// set the final values
		xn_ = best_values_;

		return true;
	} // ParticleSwarmOptimization::run()


	parameter_map_t ParticleSwarmOptimization::get_best_values() const {
		parameter_map_t values;
		for(int i = 0; i < num_params_; ++ i)
			values[params_[i]] = best_values_[i];
		return values;
	} // ParticleSwarmOptimization::get_best_values()


	bool ParticleSwarmOptimization::simulate_generation() {
		// for each particle, simulate
		for(int i = 0; i < num_particles_; ++ i) {
			// construct param map
			//parameter_map_t curr_particle;
			float_vec_t curr_particle;
			for(int j = 0; j < num_params_; ++ j)
				curr_particle.push_back(particles_[i].param_values_[j]);
				//curr_particle[params_[j]] = particles_[i].param_values_[j];
			// compute the fitness
			float_vec_t curr_fitness = (*obj_func_)(curr_particle);
			// update particle fitness
			if(particles_[i].best_fitness_ > curr_fitness[0]) {
				particles_[i].best_fitness_ = curr_fitness[0];
				particles_[i].best_values_ = particles_[i].param_values_;
			} // if
			// update global best fitness
			if(best_fitness_ > curr_fitness[0]) {
				best_fitness_ = curr_fitness[0];
				best_values_ = particles_[i].param_values_;
			} // if
			std::cout << "@@@@ " << i + 1 << ": " << curr_fitness[0] << " "
				<< particles_[i].best_fitness_ << " [ ";
			for(int j = 0; j < num_params_; ++ j) std::cout << particles_[i].param_values_[j] << " ";
			std::cout << "]" << std::endl;

			if(!particles_[i].update_particle(pso_omega_, pso_phi1_, pso_phi2_, best_values_,
											constraints_, rand_)) {
				std::cerr << "error: failed to update particle " << i << std::endl;
				return false;
			} // if
		} // for
		std::cout << "~~~~ Generation best: ";
		for(int i = 0; i < num_particles_; ++ i) {
			std::cout << "[ ";
			for(int j = 0; j < num_params_; ++ j)
				std::cout << particles_[i].param_values_[j] << " ";
			std::cout << "]\t";
		} // for
		std::cout << std::endl;
		std::cout << "~~~~ Global best: ";
		std::cout << best_fitness_ << " [ ";
		for(int j = 0; j < num_params_; ++ j)
			std::cout << best_values_[j] << " ";
		std::cout << "]\t";
		std::cout << std::endl;

		// write out stuff to output files
		std::string prefix(HiGInput::instance().param_pathprefix()+"/"+HiGInput::instance().runname());
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
		out2.close();

		return true;
	} // ParticleSwarmOptimization::simulate_generation()


} // namespace hig
