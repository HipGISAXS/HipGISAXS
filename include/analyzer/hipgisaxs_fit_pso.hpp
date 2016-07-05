/**
 *  Project:
 *
 *  File: hipgisaxs_fit_pso.hpp
 *  Created: Jan 13, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __HIPGISAXS_FIT_PSO_HPP__
#define __HIPGISAXS_FIT_PSO_HPP__

#include <set>

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
      parameter_data_list_t param_values_min_;  // minimum parameter values (inclusive)
      parameter_data_list_t param_values_max_;  // maximum parameter values (inclusive)
      parameter_data_list_t velocity_min_;    // minimum particle velocity components
      parameter_data_list_t velocity_max_;    // maximum particle velocity components
      // list of stepping values ...

    protected:
      PSOParticleConstraints() { }
      ~PSOParticleConstraints() { }

    friend class PSOParticle;
    friend class ParticleSwarmOptimization;
  }; // class PSOParticleConstraints


  // class defining one particle
  class PSOParticle {
    private:
      unsigned int index_;            // global index/id of this particle
      unsigned int num_parameters_;        // number of parameters
      parameter_data_list_t param_values_;    // list of all current parameter values
      parameter_data_list_t velocity_;      // current particle velocity componenets
      real_t fitness_;              // fitness for current parameter values
      parameter_data_list_t best_values_;      // particle's best known parameter values
      real_t best_fitness_;            // the fitness value for best parameter values

      // these are used with custom topologies only
      real_t best_fitness_global_;        // best fitness globally wrt this particle
      parameter_data_list_t best_values_global_;  // particle's global best values wrt this particle

      // for non-fully-connected topologies

      struct ParticleNeighbor {
        int proc_;      // processor rank where the neighbor lies
        int index_;      // the index (global) of the neighbor particle

        bool operator<(const ParticleNeighbor b) const {
          return (proc_ == b.proc_) ? index_ < b.index_ : proc_ < b.proc_;
        } // operator<()
      };
      typedef std::set <ParticleNeighbor> particle_neighbor_set_t;
      particle_neighbor_set_t neighbors_;

      // functions

      bool init(pso_parameter_dist_t, woo::MTRandomNumberGenerator&,
            const PSOParticleConstraints&);
      bool init_random_uniform(woo::MTRandomNumberGenerator&,
            const PSOParticleConstraints&);
      bool init_random_gaussian(woo::MTRandomNumberGenerator&,
            const PSOParticleConstraints&);
      bool init_single(const PSOParticleConstraints&);

    protected:
      PSOParticle(unsigned int, unsigned int, pso_parameter_dist_t,
            woo::MTRandomNumberGenerator&, const PSOParticleConstraints&);
      bool insert_neighbor(unsigned int, int);
      bool update_particle(real_t, real_t, real_t, const parameter_data_list_t&,
            const PSOParticleConstraints&, woo::MTRandomNumberGenerator&);
      bool update_fips_particle(real_t, real_t, real_t, const parameter_data_list_t&,
            const PSOParticleConstraints&, woo::MTRandomNumberGenerator&);
      bool update_fdr_particle(real_t, real_t, real_t, const parameter_data_list_t&,
            const real_vec_t&, const real_vec_t&, const PSOParticleConstraints&,
            woo::MTRandomNumberGenerator&);
      bool update_barebones_particle(const parameter_data_list_t&, const PSOParticleConstraints&,
            woo::MTRandomNumberGenerator&);
      bool compute_and_set_values(const parameter_data_list_t&, const parameter_data_list_t&,
            real_t, real_t, real_t, const parameter_data_list_t&,
            const PSOParticleConstraints&, woo::MTRandomNumberGenerator&);
      parameter_data_list_t get_param_values() { return param_values_; }

      void print_neighbors() {
        for(particle_neighbor_set_t::iterator i = neighbors_.begin(); i != neighbors_.end(); ++ i)
          std::cout << "(" << (*i).proc_ << "," << (*i).index_ << ")";
        std::cout << std::endl;
      } // print_neighbors()
    public:
      ~PSOParticle();

    friend class ParticleSwarmOptimization;
  }; // class PSOParticle

  typedef std::vector <PSOParticle> pso_particle_list_t;


  class ParticleSwarmOptimization : public AnalysisAlgorithm {
    private:
      //unsigned int num_parameters_;        // number of parameters (particle dimension)
      parameter_name_list_t params_;        // list of parameter names
      parameter_data_list_t best_values_;      // current overall best parameter values
      real_t best_fitness_;            // fitness value for best parameter values

      unsigned int num_particles_global_;      // total number of particles globally
      unsigned int num_particles_;        // total number of local particles
      unsigned int start_particle_index_;      // the first particle index i have
      pso_particle_list_t particles_;        // list of particles

      PSOParticleConstraints constraints_;    // particle constraints

      // PSO related parameters
      real_t pso_omega_;
      real_t pso_phi1_;
      real_t pso_phi2_;

      bool tune_omega_;              // whether to vary omega or keep constant
      unsigned int type_;              // pso type
      unsigned int foresee_num_;          // number of samples to choose from

      // helpers
      woo::MTRandomNumberGenerator rand_;      // random number generator

      // for multiple node usage
      typedef std::map <int, std::set <unsigned int> > comm_list_t;
      typedef comm_list_t::iterator comm_list_iter_t;
      typedef woo::comm_t comm_t;
      #ifdef USE_MPI
        // multinode communicator pointer
        woo::MultiNode *multi_node_;
        comm_list_t neighbor_send_lists_;
        comm_list_t neighbor_recv_lists_;
        std::vector <unsigned int> start_indices_;  // starting index for each proc
      #endif
      comm_t root_comm_;
      comm_t particle_comm_;
      comm_t pmasters_comm_;

      bool construct_neighbor_lists();
      bool neighbor_data_exchange();

      bool simulate_generation();          // simulate single generation
      bool simulate_fips_generation();      // simulate single generation
      bool simulate_fdr_generation();        // simulate single generation
      bool simulate_barebones_generation();    // simulate single generation
      bool simulate_soothsayer_generation();    // simulate single generation with foresee

      // other static topologies
      bool simulate_lbest_generation();      // lbest version with K = 2
      bool simulate_vonnewmann_generation();    // lbest grid/torus with K = 4
      bool simulate_random_generation();      // random connectivity with degree k

      unsigned int decode_pso_algo_type(real_t val) { return (unsigned int) val; }

    public:
      ParticleSwarmOptimization(int, char**, ObjectiveFunction*,
                                real_t, real_t, real_t, int, int, bool, int);
      ParticleSwarmOptimization(int, char**, ObjectiveFunction*, unsigned int, bool, int);
      ~ParticleSwarmOptimization();

      bool init();

      bool run(int, char**, int, int);              // simulate
      parameter_map_t get_best_values() const;
      #ifdef USE_MPI
        bool is_master() const { return (*multi_node_).is_master(root_comm_); }
      #else
        bool is_master() const { return true; }
      #endif
  }; // class ParticleSwarmOptimization

} // namespace hig

#endif // __HIPGISAXS_FIT_PSO_HPP__
