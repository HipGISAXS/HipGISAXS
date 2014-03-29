/**
 *  Project:
 *
 *  File: hipgisaxs_fit_pso.hpp
 *  Created: Jan 13, 2014
 *  Modified: Fri 28 Mar 2014 04:53:21 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __HIPGISAXS_FIT_PSO_HPP__
#define __HIPGISAXS_FIT_PSO_HPP__

#include <analyzer/analysis_algorithm.hpp>
#include <analyzer/hipgisaxs_fit_pso_typedefs.hpp>
#include <woo/random/woo_mtrandom.hpp>

#ifdef USE_MPI
#include <woo/comm/multi_node_comm.hpp>
#endif

namespace hig {

	// class defining constrainsts for particles
	class PSOParticleConstraints {
		private:
			parameter_data_list_t param_values_min_;	// minimum parameter values (inclusive)
			parameter_data_list_t param_values_max_;	// maximum parameter values (inclusive)
			parameter_data_list_t velocity_min_;		// minimum particle velocity components
			parameter_data_list_t velocity_max_;		// maximum particle velocity components

		protected:
			PSOParticleConstraints() { }
			~PSOParticleConstraints() { }

		friend class PSOParticle;
		friend class ParticleSwarmOptimization;
	}; // class PSOParticleConstraints


	// class defining one particle
	class PSOParticle {
		private:
			unsigned int num_parameters_;				// number of parameters
			parameter_data_list_t param_values_;		// list of all current parameter values
			parameter_data_list_t velocity_;			// current particle velocity componenets
			parameter_data_list_t best_values_;			// particle's best known parameter values
			float_t best_fitness_;						// the fitness value for best parameter values

			bool init(pso_parameter_dist_t, woo::MTRandomNumberGenerator&,
						const PSOParticleConstraints&);
			bool init_random_uniform(woo::MTRandomNumberGenerator&,
						const PSOParticleConstraints&);
			bool init_random_gaussian(woo::MTRandomNumberGenerator&,
						const PSOParticleConstraints&);
			bool init_single(const PSOParticleConstraints&);

		protected:
			PSOParticle(unsigned int, pso_parameter_dist_t, woo::MTRandomNumberGenerator&,
						const PSOParticleConstraints&);
			bool update_particle(float_t, float_t, float_t, const parameter_data_list_t&,
						const PSOParticleConstraints&, woo::MTRandomNumberGenerator&);
			bool update_fips_particle(float_t, float_t, float_t, const parameter_data_list_t&,
						const PSOParticleConstraints&, woo::MTRandomNumberGenerator&);
			bool compute_and_set_values(const parameter_data_list_t&, const parameter_data_list_t&,
						float_t, float_t, float_t, const parameter_data_list_t&,
						const PSOParticleConstraints&, woo::MTRandomNumberGenerator&);
		public:
			~PSOParticle();

		friend class ParticleSwarmOptimization;
	}; // class PSOParticle

	typedef std::vector <PSOParticle> pso_particle_list_t;


	class ParticleSwarmOptimization : public AnalysisAlgorithm {
		private:
			unsigned int num_parameters_;				// number of parameters (particle dimension)
			parameter_name_list_t params_;				// list of parameter names
			parameter_data_list_t best_values_;			// current overall best parameter values
			float_t best_fitness_;						// fitness value for best parameter values

			unsigned int num_particles_;				// total number of particles (population size)
			pso_particle_list_t particles_;				// list of particles

			PSOParticleConstraints constraints_;		// particle constraints

			// PSO related parameters
			float_t pso_omega_;
			float_t pso_phi1_;
			float_t pso_phi2_;

			bool tune_omega_;							// whether to vary omega or keep constant
			bool foresee_;
			unsigned int foresee_num_;					// number of samples to choose from

			// helpers
			woo::MTRandomNumberGenerator rand_;			// random number generator

			// for multiple node usage
			#ifdef USE_MPI
				// multinode communicator pointer
				woo::MultiNode* multi_node_;
			#endif
			std::string root_comm_;
			std::string particle_comm_;

			bool simulate_generation();					// simulate single generation
			bool simulate_fips_generation();			// simulate single generation
			bool simulate_soothsayer_generation();		// simulate single generation with foresee

		public:
			ParticleSwarmOptimization(int, char**, ObjectiveFunction*,
										float_t, float_t, float_t, int, int, bool, bool);
			~ParticleSwarmOptimization();

			bool run(int, char**, int);							// simulate
			parameter_map_t get_best_values() const;
			#ifdef USE_MPI
				bool is_master() const { return (*multi_node_).is_master(root_comm_); }
			#else
				bool is_master() const { return true; }
			#endif
	}; // class ParticleSwarmOptimization

} // namespace hig

#endif // __HIPGISAXS_FIT_PSO_HPP__
