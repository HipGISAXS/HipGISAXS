/**
 *  Project:
 *
 *  File: hipgisaxs_fit_pso_particle.cpp
 *  Created: Jan 13, 2014
 *  Modified: Mon 14 Apr 2014 09:58:17 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */


#include <analyzer/hipgisaxs_fit_pso.hpp>

namespace hig {

	PSOParticle::PSOParticle(unsigned int id, unsigned int num_param, pso_parameter_dist_t dist,
							woo::MTRandomNumberGenerator& rng,
							const PSOParticleConstraints& constraints) {
		index_ = id;
		num_parameters_ = num_param;
		init(dist, rng, constraints);
	} // PSOParticle::PSOParticle()


	PSOParticle::~PSOParticle() {

	} // PSOParticle::~PSOParticle()


	bool PSOParticle::init(pso_parameter_dist_t dist, woo::MTRandomNumberGenerator& rng,
							const PSOParticleConstraints& constraints) {
		param_values_.clear();
		velocity_.clear();
		best_values_.clear();
		best_fitness_ = std::numeric_limits<float_t>::max();
		switch(dist) {
			case PSO_DEFAULT:
			case PSO_UNIFORM:
				return init_random_uniform(rng, constraints);

			case PSO_GAUSSIAN:
				return init_random_gaussian(rng, constraints);

			case PSO_SINGLE:
				return init_single(constraints);

			default:
				std::cerr << "error: invalid parameter distribution" << std::endl;
				return false;
		} // switch
		return true;
	} // PSOParticle::init()


	bool PSOParticle::init_random_uniform(woo::MTRandomNumberGenerator& rng,
											const PSOParticleConstraints& constraints) {
		for(int i = 0; i < num_parameters_; ++ i) {
			float_t val = constraints.param_values_min_[i] +
							rng.rand() *
							(constraints.param_values_max_[i] - constraints.param_values_min_[i]);
			float_t vel = constraints.velocity_min_[i] +
							rng.rand() *
							(constraints.velocity_max_[i] - constraints.velocity_min_[i]);
			param_values_.push_back(val);
			velocity_.push_back(vel);
			best_values_.push_back(val);
		} // for
		return true;
	} // PSOParticle::init_random_uniform()


	bool PSOParticle::init_random_gaussian(woo::MTRandomNumberGenerator& rng,
											const PSOParticleConstraints& constraints) {
		std::cerr << "error: init_random_gaussian() is not currently implemented" << std::endl;
		return false;
	} // PSOParticle::init_random_gaussian()


	bool PSOParticle::init_single(const PSOParticleConstraints& constraints) {
		std::cerr << "error: init_single() is not currently implemented" << std::endl;
		return false;
	} // PSOParticle::init_single()


	// basic
	bool PSOParticle::update_particle(float_t omega, float_t phi1, float_t phi2,
										const parameter_data_list_t& global_best,
										const PSOParticleConstraints& constraints,
										woo::MTRandomNumberGenerator& rng) {
		for(int i = 0; i < num_parameters_; ++ i) {
			float_t r1 = rng.rand();
			float_t r2 = rng.rand();
			float_t new_vel = omega * velocity_[i] +
								phi1 * r1 * (best_values_[i] - param_values_[i]) +
								phi2 * r2 * (global_best[i] - param_values_[i]);
			velocity_[i] = std::min(std::max(new_vel,
										constraints.velocity_min_[i]),
										constraints.velocity_max_[i]);
			param_values_[i] = std::min(std::max(param_values_[i] + velocity_[i],
										constraints.param_values_min_[i]),
										constraints.param_values_max_[i]);
		} // for
		return true;
	} // PSOParticle::update_particle()


	// FIPS
	bool PSOParticle::update_fips_particle(float_t omega, float_t phi1, float_t phi2,
										const float_vec_t& best_values,
										const PSOParticleConstraints& constraints,
										woo::MTRandomNumberGenerator& rng) {
		int tot_num_parts = best_values.size() / num_parameters_;
		for(int i = 0; i < num_parameters_; ++ i) {
			float_t new_vel = omega * velocity_[i];
			float_t sum = 0.0;
			for(int j = 0; j < tot_num_parts; ++ j) {
				float_t r1 = rng.rand();
				sum += phi1 * r1 * (best_values[num_parameters_ * j + i] - param_values_[i]);
			} // for
			new_vel += sum / tot_num_parts;
			velocity_[i] = std::min(std::max(new_vel,
										constraints.velocity_min_[i]),
										constraints.velocity_max_[i]);
			param_values_[i] = std::min(std::max(param_values_[i] + velocity_[i],
										constraints.param_values_min_[i]),
										constraints.param_values_max_[i]);
		} // for
		return true;
	} // PSOParticle::update_fips_particle()


	// FDR
	bool PSOParticle::update_fdr_particle(float_t omega, float_t phi1, float_t phi2,
											const parameter_data_list_t& global_best,
											const float_vec_t& all_values,
											const float_vec_t& all_fitness,
											const PSOParticleConstraints& constraints,
											woo::MTRandomNumberGenerator& rng) {
		// for each dimension i, find among all other particles p the one with largest FDR value
		// fdr_p_i = (fitness - fitness_p) / (param_values_i - param_values_p_i)
		int tot_num_parts = all_values.size() / num_parameters_;	// total number of particles
		for(int i = 0; i < num_parameters_; ++ i) {
			// for each particle compute fdr and find best one
			float_t best_fdr = 0, best_val = 0;
			for(int p = 0; p < tot_num_parts; ++ p) {
				// need to ignore self
				if(index_ == p) continue;	// make sure they are ordered in same way
				float_t val = all_values[num_parameters_ * p + i];
				float_t fdr = fabs(fitness_ - all_fitness[i]) / fabs(param_values_[i] - val);
				if(fdr > best_fdr) { best_fdr = fdr; best_val = val; }
			} // for
			// use best particle parameter to compute new velocity
			float_t r1 = rng.rand(), r2 = rng.rand();
			float_t new_vel = omega * velocity_[i] +
								phi1 * r1 * (best_val - param_values_[i]) +
								phi2 * r2 * (global_best[i] - param_values_[i]);
			velocity_[i] = std::min(std::max(new_vel,
										constraints.velocity_min_[i]),
										constraints.velocity_max_[i]);
			param_values_[i] = std::min(std::max(param_values_[i] + velocity_[i],
										constraints.param_values_min_[i]),
										constraints.param_values_max_[i]);
		} // for
		return true;
	} // PSOParticle::update_fdr_particle()


	bool PSOParticle::compute_and_set_values(const parameter_data_list_t& start_pos,
										const parameter_data_list_t& start_vel,
										float_t omega, float_t phi1, float_t phi2,
										const parameter_data_list_t& global_best,
										const PSOParticleConstraints& constraints,
										woo::MTRandomNumberGenerator& rng) {
		for(int i = 0; i < num_parameters_; ++ i) {
			float_t r1 = rng.rand();
			float_t r2 = rng.rand();
			float_t new_vel = omega * start_vel[i] +
								phi1 * r1 * (best_values_[i] - start_pos[i]) +
								phi2 * r2 * (global_best[i] - start_pos[i]);
			velocity_[i] = std::min(std::max(new_vel,
										constraints.velocity_min_[i]),
										constraints.velocity_max_[i]);
			param_values_[i] = std::min(std::max(start_pos[i] + new_vel,
										constraints.param_values_min_[i]),
										constraints.param_values_max_[i]);
		} // for
		return true;
	} // PSOParticle::compute_and_set_values()

} // namespace hig
