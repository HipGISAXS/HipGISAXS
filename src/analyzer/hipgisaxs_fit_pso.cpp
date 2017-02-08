/**
 *  Project:
 *
 *  File: fit_pso.cpp
 *  Created: Jan 13, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <ctime>

#include <analyzer/hipgisaxs_fit_pso.hpp>
#include <hipgisaxs.hpp>

#include <woo/timer/woo_boostchronotimers.hpp>

namespace hig {

  ParticleSwarmOptimization::ParticleSwarmOptimization(int narg, char** args, ObjectiveFunction* obj,
      real_t omega, real_t phi1, real_t phi2, int npart, int ngen,
      bool tune_omega = false, int type = 0) :
        rand_(time(NULL)), type_(type) {
    name_ = algo_pso;
    max_hist_ = 100;      // not used in pso
    tol_ = 1e-6;          // default?

    obj_func_ = obj;

    //if(narg < 7) {
    //  std::cerr << "error: for PSO please specify <num_particles> <num_gen> "
    //        << "<omega> <phi1> <phi2>" << std::endl;
    //  exit(1);
    //} // if
    foresee_num_ = 5;      // TODO: make it modifiable

    pso_omega_ = omega;
    pso_phi1_ = phi1;
    pso_phi2_ = phi2;
    tune_omega_ = tune_omega;

    // get number of particles
    num_particles_ = npart;
    num_particles_global_ = npart;
    // max number of generations
    max_iter_ = ngen;

    init();
  } // ParticleSwarmOptimization::ParticleSwarmOptimization()


  ParticleSwarmOptimization::ParticleSwarmOptimization(int narg, char** args,
                                                       ObjectiveFunction* obj,
                                                       unsigned int algo_num,
                                                       bool tune_omega = false,
                                                       int type = 0) :
      rand_(time(NULL)), type_(type) {
    name_ = algo_pso;
    max_hist_ = 200;
    obj_func_ = obj;

    tol_ = (*obj_func_).analysis_tolerance(algo_num);

    foresee_num_ = 5;      // TODO: make it modifiable

    real_t temp_val = 0.0;

    // mandatory PSO parameters
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_omega", pso_omega_)) {
      std::cerr << "error: mandatory PSO parameter 'pso_omega' missing" << std::endl;
      exit(-1);
    } // if
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_phi1", pso_phi1_)) {
      std::cerr << "error: mandatory PSO parameter 'pso_phi1' missing" << std::endl;
      exit(-1);
    } // if
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_phi2", pso_phi2_)) {
      std::cerr << "error: mandatory PSO parameter 'pso_phi2' missing" << std::endl;
      exit(-1);
    } // if
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_num_particles", temp_val)) {
      std::cerr << "error: mandatory PSO parameter 'pso_num_particles' missing" << std::endl;
      exit(-1);
    } else {
      num_particles_ = (unsigned int) temp_val;
    } // if-else
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_num_generations", temp_val)) {
      std::cerr << "error: mandatory PSO parameter 'pso_num_generations' missing" << std::endl;
      exit(-1);
    } else {
      max_iter_ = (unsigned int) temp_val;
    } // if
    num_particles_global_ = num_particles_;

    // optional PSO parameters
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_tune_omega", temp_val)) {
      std::cerr << "warning: optional PSO parameter 'pso_tune_omega' defaulted" << std::endl;
      tune_omega_ = false;
    } else {
      if(temp_val > 0.5) tune_omega_ = true;
      else tune_omega_ = false;
    } // if-else
    if(!(*obj_func_).analysis_algo_param(algo_num, "pso_type", temp_val)) {
      std::cerr << "warning: optional PSO parameter 'pso_type' defaulted" << std::endl;
      type_ = 0;  // base algorithm
    } else {
      type_ = decode_pso_algo_type(temp_val);
    } // if-else

    init();
  } // ParticleSwarmOptimization::ParticleSwarmOptimization()


  ParticleSwarmOptimization::~ParticleSwarmOptimization() {
    // nothing to do here yet ...
  } // ParticleSwarmOptimization::~ParticleSwarmOptimization()


  bool ParticleSwarmOptimization::init() {
    params_ = (*obj_func_).fit_param_keys();
    num_params_ = (*obj_func_).num_fit_params();
    x0_ = (*obj_func_).fit_param_init_values();
    #ifdef USE_MPI
      // obtain the multi node object from hipgisaxs
      multi_node_ = (*obj_func_).multi_node_comm();
      // the root communicator
      root_comm_ = (*multi_node_).universe_key();
      if((*multi_node_).is_master(root_comm_))
        std::cout << "** MPI Size: " << (*multi_node_).size(root_comm_)
                  << ", my rank: " << (*multi_node_).rank(root_comm_) << std::endl;
      rand_.reset(time(NULL) + (*multi_node_).rank(root_comm_));
    #else
      root_comm_ = "world";
    #endif

    best_values_.resize(num_params_, 0.0);
    best_fitness_ = std::numeric_limits<real_t>::max();

    std::vector <std::pair <real_t, real_t> > limits = (*obj_func_).fit_param_limits();
    parameter_data_list_t max_limits, min_limits;
    for(std::vector <std::pair <real_t, real_t> >::iterator iter = limits.begin();
        iter != limits.end(); ++ iter) {
      real_t min_val = iter->first;
      real_t max_val = iter->second;
      if(max_val < min_val) std::swap(min_val, max_val);
      min_limits.push_back(min_val);
      max_limits.push_back(max_val);
    } // for
    if(min_limits.size() < num_params_) min_limits.resize(num_params_, 0.0);
    if(max_limits.size() < num_params_) max_limits.resize(num_params_,
                                                          std::numeric_limits<real_t>::max());
    constraints_.param_values_min_ = min_limits;
    constraints_.param_values_max_ = max_limits;

    // set constraints
    // for now let min max velocity be -max to max param (assuming all params are +ve) ...
    for(int i = 0; i < num_params_; ++ i) min_limits[i] = - max_limits[i];
    constraints_.velocity_min_ = min_limits;
    constraints_.velocity_max_ = max_limits;

    unsigned int index_offset = 0;
    #ifdef USE_MPI
      int num_procs = (*multi_node_).size(root_comm_);
      int rank = (*multi_node_).rank(root_comm_);
      // two cases:  num_procs <= num_part (easy case)
      //         num_procs > num_part
      particle_comm_ = "pso_particle";
      int particle_master = 1;  // am I the master proc for a particle
      if(num_procs <= num_particles_global_) {
        // assign multiple particles per proc
        num_particles_ = num_particles_global_ / num_procs +
                    (rank < (num_particles_global_ % num_procs) ? 1 : 0);
        // each proc will now be in its own new world
        (*multi_node_).split(particle_comm_, root_comm_, rank);
        (*multi_node_).scan_sum(root_comm_, num_particles_, index_offset);
        index_offset -= num_particles_;
        start_particle_index_ = index_offset;
      } else {
        // assign multiple procs per particle
        // procs for a particle will be in their common new world
        // split the world
        num_particles_ = 1;    // no fractional particles :P
        int color = rank % num_particles_global_;  // round-robin distribution for ease
        (*multi_node_).split(particle_comm_, root_comm_, color);
        particle_master = (*multi_node_).is_master(particle_comm_);
        start_particle_index_ = color;  // the index of particle i am responsible for
        index_offset = color;      // it is same as the index offset
      } // if-else
      // create communicator for just the particle masters (which may be all the procs)
      pmasters_comm_ = "pso_pmasters";
      (*multi_node_).split(pmasters_comm_, root_comm_, particle_master);
      // allgather to collect the start indices for each proc
      unsigned int * buff = new (std::nothrow) unsigned int[(*multi_node_).size(root_comm_)];
      (*multi_node_).allgather(root_comm_, &start_particle_index_, 1, buff, 1);
      for(int i = 0; i < (*multi_node_).size(root_comm_); ++ i)
        start_indices_.push_back(buff[i]);
      delete[] buff;
    #else
      particle_comm_ = root_comm_;
      pmasters_comm_ = root_comm_;
      start_particle_index_ = 0;
    #endif

    // initialize particles
    particles_.clear();
    for(int i = 0; i < num_particles_; ++ i) {
      particles_.push_back(PSOParticle(index_offset + i, num_params_,
                        PSO_UNIFORM, rand_, constraints_));
    } // for

    if(type_ > 4) {    // need to construct neighbor list (non-fully-connected topology)
      if(!construct_neighbor_lists()) exit(2);
    } // if

    if(is_master()) {
      std::cout << "** Fitting Parameters: Tolerance = " << tol_ << std::endl
                << "**     PSO Parameters:     Omega = " << pso_omega_ << std::endl
                << "                     Tune Omega? = " << tune_omega_ << std::endl
                << "                            Phi1 = " << pso_phi1_ << std::endl
                << "                            Phi2 = " << pso_phi2_ << std::endl
                << "                Number of Agents = " << num_particles_ << std::endl
                << "           Number of Generations = " << max_iter_ << std::endl
                << "            PSO Algorithm Flavor = " << type_ << std::endl;
    } // if

    return true;
  } // ParticleSwarmOptimization::init()


  bool ParticleSwarmOptimization::construct_neighbor_lists() {
    if((*multi_node_).is_master(root_comm_))
      std::cout << "-- Constructing neighbor lists for each particle ..." << std::endl;
    #ifdef USE_MPI
      // easy way:
      // find location of all particles through one allgatherv operation
      unsigned int *start_indices = new unsigned int[(*multi_node_).size(root_comm_) + 1];
      (*multi_node_).allgather(root_comm_, &start_particle_index_, 1, start_indices, 1);
      start_indices[(*multi_node_).size(root_comm_)] = num_particles_global_;
      std::map<unsigned int, int> proc_loc;
      unsigned int i = 0;
      for(int curr_proc = 0; curr_proc < (*multi_node_).size(root_comm_); ++ curr_proc) {
        while(i < start_indices[curr_proc + 1]) proc_loc[i ++] = curr_proc;
      } // while
      delete [] start_indices;
    #endif

    int x = 0, y = 0, k = 0, size = 0, parts = 0, overall_size = 0, my_size = 0, idx = 0, nidx = 0;
    int * counts, * offsets, * my_neighbors, * all_neighbors, * sizes, * tot_counts, * num_neighbors;
    std::vector<std::vector<int> > neighbors(num_particles_global_);
    switch(type_) {
      case 5:      // lbest, K = 2
        for(int i = 0; i < num_particles_; ++ i) {
          unsigned int myid = particles_[i].index_;
          unsigned int left;
          if(myid == 0) left = num_particles_global_ - 1;    // wrap around
          else left = myid - 1;
          unsigned int right = (myid + 1) % num_particles_global_;
          #ifdef USE_MPI
            particles_[i].insert_neighbor(left, proc_loc[left]);
            particles_[i].insert_neighbor(right, proc_loc[right]);
            neighbor_send_lists_[proc_loc[left]].insert(myid);
            neighbor_send_lists_[proc_loc[right]].insert(myid);
            neighbor_recv_lists_[proc_loc[left]].insert(left);
            neighbor_recv_lists_[proc_loc[right]].insert(right);
          #else
            particles_[i].insert_neighbor(left, 0);
            particles_[i].insert_neighbor(right, 0);
          #endif
        } // for
        break;

      case 6:      // von Newmann. K = 4
        // arrange particles in a 2d grid
        x = sqrt(num_particles_global_);
        y = num_particles_global_ / x;
        std::cout << "======================= " << x << " x " << y << std::endl;
        if(x * y != num_particles_global_) {
          std::cerr << "error: number of particles is not compatible with specified topology"
                << std::endl;
          return false;
        } // if
        for(int i = 0; i < num_particles_; ++ i) {
          unsigned int myid = particles_[i].index_;
          int myx = myid % x;
          int myy = myid / x;
          int up = (myy - 1 < 0 ? y - 1 : myy - 1) * x + myx;
          int down = ((myy + 1) % y) * x + myx;
          int left = myy * x + (myx - 1 < 0 ? x - 1 : myx - 1);
          int right = myy * x + (myx + 1) % x;
          #ifdef USE_MPI
            particles_[i].insert_neighbor(up, proc_loc[up]);
            particles_[i].insert_neighbor(down, proc_loc[down]);
            particles_[i].insert_neighbor(left, proc_loc[left]);
            particles_[i].insert_neighbor(right, proc_loc[right]);
            // send lists
            neighbor_send_lists_[proc_loc[up]].insert(myid);
            neighbor_send_lists_[proc_loc[down]].insert(myid);
            neighbor_send_lists_[proc_loc[left]].insert(myid);
            neighbor_send_lists_[proc_loc[right]].insert(myid);
            // recv lists
            neighbor_recv_lists_[proc_loc[up]].insert(up);
            neighbor_recv_lists_[proc_loc[down]].insert(down);
            neighbor_recv_lists_[proc_loc[left]].insert(left);
            neighbor_recv_lists_[proc_loc[right]].insert(right);
          #else
            particles_[i].insert_neighbor(up, 0);
            particles_[i].insert_neighbor(down, 0);
            particles_[i].insert_neighbor(left, 0);
            particles_[i].insert_neighbor(right, 0);
          #endif
        } // for
        break;

      case 7:      // random, K = k
        k = 0.2 * (num_particles_global_ - 1);  // 20% connectivity
        // let just one proc do this and scatter to others
        if(k < 2 || k > num_particles_global_ - 1) {
          std::cerr << "error: you cannot have value of k greater than "
                << "a fully connected topology. use the base case instead"
                << std::endl;
          return false;
        } // if
        for(int i = 0; i < num_particles_global_; ++ i) {
          // for each particle, pick k other particles and connect those
          for(int j = 0; j < k; ++ j) {
            int idx = rand_.rand() * (num_particles_global_ - 2);
            if(idx > i) ++ idx;
            if(idx > num_particles_global_ - 1) {
              std::cerr << "error: something went very wrong in picking randoms"
                    << std::endl;
              return false;
            } // if
            neighbors[i].push_back(idx);
            neighbors[idx].push_back(i);
          } // for
        } // for
        counts = new (std::nothrow) int[(*multi_node_).size(root_comm_)];
        tot_counts = new (std::nothrow) int[(*multi_node_).size(root_comm_)];
        sizes = new (std::nothrow) int[num_particles_global_];
        (*multi_node_).allgather(root_comm_, (int*) &num_particles_, 1, counts, 1);
        parts = 0; overall_size = 0;
        for(int i = 0; i < (*multi_node_).size(root_comm_); ++ i) {
          tot_counts[i] = 0;
          for(int j = 0; j < counts[i]; ++ j, ++ parts) {
            sizes[parts] = neighbors[parts].size();
            tot_counts[i] += neighbors[parts].size();
          } // for
          overall_size += tot_counts[i];
        } // for
        offsets = new (std::nothrow) int[(*multi_node_).size(root_comm_)];
        offsets[0] = 0;
        for(int i = 1; i < (*multi_node_).size(root_comm_); ++ i)
          offsets[i] = offsets[i - 1] + counts[i - 1];
        num_neighbors = new (std::nothrow) int[num_particles_];
        (*multi_node_).scatterv(root_comm_, sizes, counts, offsets, num_neighbors,
                    num_particles_);
        my_size = 0;
        for(int i = 0; i < num_particles_; ++ i) my_size += num_neighbors[i];
        all_neighbors = new (std::nothrow) int[overall_size];
        idx = 0;
        for(std::vector<std::vector<int> >::iterator i = neighbors.begin();
            i != neighbors.end(); ++ i)
          for(std::vector<int>::iterator j = (*i).begin(); j != (*i).end(); ++ j)
            all_neighbors[idx ++] = *j;
        my_neighbors = new (std::nothrow) int[my_size];
        (*multi_node_).scatterv(root_comm_, all_neighbors, tot_counts, offsets, my_neighbors,
                    my_size);
        // now we have received the neighbor information, save them
        nidx = 0;
        for(int i = 0; i < num_particles_; ++ i) {
          unsigned int myid = particles_[i].index_;
          for(int j = 0; j < num_neighbors[i]; ++ j) {
            int val = my_neighbors[nidx ++];
            particles_[i].insert_neighbor(val, proc_loc[val]);
            neighbor_send_lists_[proc_loc[val]].insert(myid);
            neighbor_recv_lists_[proc_loc[val]].insert(val);
          } // for
        } // for
        delete[] my_neighbors;
        delete[] all_neighbors;
        delete[] offsets;
        delete[] sizes;
        delete[] tot_counts;
        delete[] counts;
        break;

      default:
        std::cerr << "error: given type is not custom topology" << std::endl;
        return false;
    } // switch

    // once the particle neighbors are known, construct the neighbor processor communication lists
    // TODO ...

    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << " -> ";
      particles_[i].print_neighbors();
    } // for
    for(comm_list_iter_t i = neighbor_send_lists_.begin(); i != neighbor_send_lists_.end(); ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": send list " << (*i).first << " -> ";
      for(std::set<unsigned int>::iterator j = (*i).second.begin(); j != (*i).second.end(); ++ j)
        std::cout << (*multi_node_).rank(root_comm_) << ":" << *j << " ";
      std::cout << std::endl;
    } // for
    for(comm_list_iter_t i = neighbor_recv_lists_.begin(); i != neighbor_recv_lists_.end(); ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": recv list " << (*i).first << " -> ";
      for(std::set<unsigned int>::iterator j = (*i).second.begin(); j != (*i).second.end(); ++ j)
        std::cout << (*multi_node_).rank(root_comm_) << ":" << *j << " ";
      std::cout << std::endl;
    } // for

    return true;
  } // ParticleSwarmOptimization::construct_neighbor_lists()


  parameter_map_t ParticleSwarmOptimization::get_best_values() const {
    parameter_map_t values;
    for(int i = 0; i < num_params_; ++ i)
      values[params_[i]] = best_values_[i];
    return values;
  } // ParticleSwarmOptimization::get_best_values()


  bool ParticleSwarmOptimization::run(int narg, char** args, int algo_num, int img_num) {
    #ifdef USE_MPI
    if((*multi_node_).is_master(root_comm_))
    #endif
      std::cout << "Running Particle Swarm Optimization ..." << std::endl;

    if(!(*obj_func_).set_reference_data(img_num)) return false;

    woo::BoostChronoTimer gen_timer;
    double total_time = 0;

    for(int gen = 0; gen < max_iter_; ++ gen) {

      #ifdef USE_MPI
      if((*multi_node_).is_master(root_comm_))
      #endif
        std::cout << "** Generation " << gen << std::endl;

      gen_timer.start();
      if(type_ == 2) {        // soothsayer
        if(!simulate_soothsayer_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 1) {    // fips
        if(!simulate_fips_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 3) {    // fdr
        if(!simulate_fdr_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 4) {    // barebones
        if(!simulate_barebones_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 5) {    // lbest
        if(!simulate_lbest_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 6) {    // von newmann
        if(!simulate_vonnewmann_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 7) {    // random
        if(!simulate_random_generation()) {
          std::cerr << "error: failed in generation " << gen << std::endl;
          return false;
        } // if
      } else if(type_ == 0) {    // base (inertia weight)
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
        std::cout << "@@@ Generation time: " << elapsed << " ms. [total: "
                  << total_time << " ms.]" << std::endl;
    } // for
    
    // set the final values
    xn_ = best_values_;

    // print the global best - only root master does this
    #ifdef USE_MPI
    if((*multi_node_).is_master(root_comm_)) {
    #endif
      std::cout << "@@@ Global best: ";
      std::cout << best_fitness_ << " [ ";
      for(int j = 0; j < num_params_; ++ j) std::cout << params_[j] << ": " << best_values_[j] << " ";
      std::cout << "]" << std::endl;
      std::ofstream outfile("hipgisaxs_fit_params.txt");
      for(int j = 0; j < num_params_; ++ j)
        outfile << params_[j] << ":\t" << best_values_[j] << std::endl;
      outfile.close();
    #ifdef USE_MPI
    } // if
    #endif

    return true;
  } // ParticleSwarmOptimization::run()


  // Base case and with constriction coeff, tuned omega
  bool ParticleSwarmOptimization::simulate_generation() {
    // for each particle, simulate
    int myrank = 0;
    #ifdef USE_MPI
    myrank = (*multi_node_).rank(root_comm_);
    #endif

    for(int i = 0; i < num_particles_; ++ i) {

      #ifdef USE_MPI
      if((*multi_node_).is_master(root_comm_))
      #endif
        std::cout << "** Particle " << i << std::endl;

      // construct param map
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j) curr_particle.push_back(particles_[i].param_values_[j]);

      // tell hipgisaxs about the communicator to work with
      (*obj_func_).update_sim_comm(particle_comm_);

      // compute the fitness
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);

      // update particle fitness
      // this is meaningful only at the particle masters
      #ifdef USE_MPI
      (*multi_node_).barrier(particle_comm_);
      if((*multi_node_).is_master(particle_comm_)) {
      #endif
        if(particles_[i].best_fitness_ > curr_fitness[0]) {
          particles_[i].best_fitness_ = curr_fitness[0];
          particles_[i].best_values_ = particles_[i].param_values_;
        } // if
        // update global best fitness locally
        if(best_fitness_ > curr_fitness[0]) {
          best_fitness_ = curr_fitness[0];
          best_values_ = particles_[i].param_values_;
        } // if
      #ifdef USE_MPI
      } // if
      #endif

      // write out the current status
      #ifdef USE_MPI
      if((*multi_node_).is_master(particle_comm_)) {  // only particle masters do it
      #endif
        std::stringstream cfilename_s;
        cfilename_s << "convergence." << myrank << "." << i << ".dat";
        std::string prefix((*obj_func_).param_pathprefix() + "/" + (*obj_func_).runname());
        std::ofstream out(prefix + "/" + cfilename_s.str(), std::ios::app);
        out.precision(10);
        out << myrank << "\t" << i << "\t" << curr_fitness[0] << "\t";
        for(int j = 0; j < num_params_; ++ j) out << particles_[i].param_values_[j] << "\t";
        out << particles_[i].best_fitness_ << "\t";
        for(int j = 0; j < num_params_; ++ j) out << particles_[i].best_values_[j] << "\t";
        out << std::endl;
        out.close();
      #ifdef USE_MPI
      } // if
      #endif
    } // for

    #ifdef USE_MPI
      // communicate and find globally best fitness
      real_t new_best_fitness = 0.0; int best_rank = -1;
      // all particle masters find the proc with global min
      if(!(*multi_node_).allreduce(pmasters_comm_, best_fitness_, new_best_fitness, best_rank,
                                   woo::comm::minloc))
        return false;
      best_fitness_ = new_best_fitness;   // this makes sense only on particle masters, but doesnt hurt
      // the best rank proc broadcasts its best values to all other particle masters
      if(!(*multi_node_).broadcast(pmasters_comm_, best_values_, best_rank)) return false;
      // also broadcast to other procs responsible for same particle, if is the case
      if((*multi_node_).size(particle_comm_) > 1) {
        if(!(*multi_node_).broadcast(particle_comm_, best_fitness_,
                                     (*multi_node_).master(particle_comm_))) return false;
        if(!(*multi_node_).broadcast(particle_comm_, best_values_,
                                     (*multi_node_).master(particle_comm_))) return false;
      } // if
    #endif

    // update each particle using new global best
    // only the particle masters need to do this
    #ifdef USE_MPI
    if((*multi_node_).is_master(particle_comm_)) {
    #endif
      for(int i = 0; i < num_particles_; ++ i) {
        if(!particles_[i].update_particle(pso_omega_, pso_phi1_, pso_phi2_, best_values_,
                        constraints_, rand_)) {
          std::cerr << "error: failed to update particle " << i << std::endl;
          return false;
        } // if
      } // for
    #ifdef USE_MPI
    } // if
    // the particle master then broadcasts this updated info to all
    // other procs responsible for the same particle, if any
    if((*multi_node_).size(particle_comm_) > 1 && num_particles_ == 1) {
      parameter_data_list_t pvals = particles_[0].get_param_values();
      if(!(*multi_node_).broadcast(particle_comm_, pvals,
                    (*multi_node_).master(particle_comm_))) return false;
    } // if
    #endif

    if(tune_omega_) {
      pso_omega_ /= 2.0;
      #ifdef USE_MPI
      if((*multi_node_).is_master(root_comm_))
      #endif
        std::cout << myrank << ": Omega = " << pso_omega_ << std::endl;
    } // if

    // print the global best - only root master does this
    #ifdef USE_MPI
    if((*multi_node_).is_master(root_comm_)) {
    #endif
      std::cout << "@@@ Global best: " << best_fitness_ << " [ ";
      for(int j = 0; j < num_params_; ++ j) std::cout << params_[j] << ": " << best_values_[j] << " ";
      std::cout << "]" << std::endl;
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


  // Bare-bones PSO
  bool ParticleSwarmOptimization::simulate_barebones_generation() {
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // construct param map
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j)
        curr_particle.push_back(particles_[i].param_values_[j]);
      // compute the fitness
      (*obj_func_).update_sim_comm(particle_comm_);
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);

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
    } // for

    #ifdef USE_MPI
      // communicate and find globally best fitness
      real_t new_best_fitness; int best_rank;
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
      if(!particles_[i].update_barebones_particle(best_values_, constraints_, rand_)) {
        std::cerr << "error: failed to update particle " << i << std::endl;
        return false;
      } // if
    } // for

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

    return true;
  } // ParticleSwarmOptimization::simulate_barebones_generation()


  // FDR: Fitness Distance Ratio
  // [ratio of difference in fitness of nearest neighbor to distance, in each dimension]
  bool ParticleSwarmOptimization::simulate_fdr_generation() {
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // construct param map
      //parameter_map_t curr_particle;
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j)
        curr_particle.push_back(particles_[i].param_values_[j]);
        //curr_particle[params_[j]] = particles_[i].param_values_[j];
      // compute the fitness
      (*obj_func_).update_sim_comm(particle_comm_);
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);
      particles_[i].fitness_ = curr_fitness[0];    // in PSO we get only one fitness value

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
    } // for

    real_vec_t all_values, all_fitness;
    #ifdef USE_MPI
      // communicate with everyone to get their current fitness and params (alltoall)
      real_t new_best_fitness; int best_rank;
      if((*multi_node_).size("pso_particle") > 1 && num_particles_ == 1) {
        // TODO ...
        bool pmaster = (*multi_node_).is_master(particle_comm_);
        std::cout << "error: this case has not been implemented yet" << std::endl;
        return false;
      } else {
        // communicate all the values to all
        // and the fitnesses
        real_vec_t all_local_values(num_particles_ * num_params_);
        real_vec_t all_local_fitness(num_particles_);
        for(int i = 0; i < num_particles_; ++ i) {
          for(int j = 0; j < num_params_; ++ j)
            all_local_values[i * num_params_ + j] = particles_[i].param_values_[j];
          all_local_fitness[i] = particles_[i].fitness_;
        } // for
        if(!(*multi_node_).allgatherv(root_comm_, all_local_values, all_values)) return false;
        if(!(*multi_node_).allgatherv(root_comm_, all_local_fitness, all_fitness)) return false;

        // also keep track of global best
        if(!(*multi_node_).allreduce(root_comm_, best_fitness_, new_best_fitness, best_rank,
                        woo::comm::minloc)) return false;
        best_fitness_ = new_best_fitness;
      } // if-else
      // the best rank proc broadcasts its best values to all others
      if(!(*multi_node_).broadcast(root_comm_, best_values_, best_rank)) return false;
    #else
      for(int i = 0; i < num_particles_; ++ i) {
        for(int j = 0; j < num_params_; ++ j)
          all_values.push_back(particles_[i].param_values_[j]);
        all_fitness.push_back(particles_[i].fitness_);
      } // for
    #endif

    // update each particle using all best values
    for(int i = 0; i < num_particles_; ++ i) {
      if(!particles_[i].update_fdr_particle(pso_omega_, pso_phi1_, pso_phi2_, best_values_,
                          all_values, all_fitness,
                          constraints_, rand_)) {
        std::cerr << "error: failed to update particle " << i << std::endl;
        return false;
      } // if
    } // for

    if(tune_omega_) pso_omega_ /= 1.5;
    std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

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

    return true;
  } // ParticleSwarmOptimization::simulate_fdr_generation()


  // FIPS: Fully Informed Particle Swarm
  bool ParticleSwarmOptimization::simulate_fips_generation() {
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // construct param map
      //parameter_map_t curr_particle;
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j)
        curr_particle.push_back(particles_[i].param_values_[j]);
        //curr_particle[params_[j]] = particles_[i].param_values_[j];
      // compute the fitness
      (*obj_func_).update_sim_comm(particle_comm_);
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);

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
    } // for

    real_vec_t all_best_values;
    #ifdef USE_MPI
      // communicate with everyone to get their best (alltoall)
      real_t new_best_fitness; int best_rank;
      if((*multi_node_).size("pso_particle") > 1 && num_particles_ == 1) {
        // TODO ...
        bool pmaster = (*multi_node_).is_master(particle_comm_);
        std::cout << "error: this case has not been implemented yet" << std::endl;
        return false;
      } else {
        // communicate all the values to all
        real_vec_t all_local_best_values(num_particles_ * num_params_);
        for(int i = 0; i < num_particles_; ++ i)
          for(int j = 0; j < num_params_; ++ j)
            all_local_best_values[i * num_params_ + j] = particles_[i].best_values_[j];
        if(!(*multi_node_).allgatherv(root_comm_, all_local_best_values, all_best_values))
          return false;

        // also keep track of global best
        // first find out who has the best one
        if(!(*multi_node_).allreduce(root_comm_, best_fitness_, new_best_fitness, best_rank,
                        woo::comm::minloc))
          return false;
        best_fitness_ = new_best_fitness;
      } // if-else
      // the best rank proc broadcasts its best values to all others
      if(!(*multi_node_).broadcast(root_comm_, best_values_, best_rank)) return false;
    #else
      for(int i = 0; i < num_particles_; ++ i)
        for(int j = 0; j < num_params_; ++ j)
          all_best_values.push_back(particles_[i].best_values_[j]);
    #endif

    // update each particle using all best values (from all particles)
    for(int i = 0; i < num_particles_; ++ i) {
      if(!particles_[i].update_fips_particle(pso_omega_, pso_phi1_, pso_phi2_, all_best_values,
                      constraints_, rand_)) {
        std::cerr << "error: failed to update particle " << i << std::endl;
        return false;
      } // if
    } // for

    if(tune_omega_) pso_omega_ /= 1.5;
    std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

    // print current global best
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

    return true;
  } // ParticleSwarmOptimization::simulate_fips_generation()


  // Soothsayer: look ahead few possibilities and pick the best
  bool ParticleSwarmOptimization::simulate_soothsayer_generation() {
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // first save the current state of the particle
      parameter_data_list_t orig_param_vals = particles_[i].param_values_;
      parameter_data_list_t orig_velocity = particles_[i].velocity_;
      real_t foresee_best_fitness = std::numeric_limits<real_t>::max();
      parameter_data_list_t foresee_best_params, foresee_best_velocity;
      // for each forseeing step
      for(int k = 0; k < foresee_num_; ++ k) {
        // create a new position from the actual current position
        particles_[i].compute_and_set_values(orig_param_vals, orig_velocity,
                            pso_omega_, pso_phi1_, pso_phi2_,
                            best_values_, constraints_, rand_);
        // compute the fitness
        real_vec_t curr_particle = particles_[i].param_values_;
        (*obj_func_).update_sim_comm(particle_comm_);
        real_vec_t curr_fitness = (*obj_func_)(curr_particle);

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
      real_t new_best_fitness; int best_rank;
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


  // other topologies:

  // TODO: optimize communication later ...
  // for each proc in the send list, collect/pack data and send
  // for each proc in the recv list, receive data and unpack
  bool ParticleSwarmOptimization::neighbor_data_exchange() {
    // Irecv
    // first receive fitness values and parameter values
    std::vector<MPI_Request> rreq;
    std::vector<real_t*> recv_bufs;
    int myrank = (*multi_node_).rank(root_comm_);
    for(comm_list_iter_t r = neighbor_recv_lists_.begin(); r != neighbor_recv_lists_.end(); ++ r) {
//      if(myrank == (*r).first || (*r).second.size() == 0) continue;
      real_t * recv_fitness = new (std::nothrow) real_t[(*r).second.size()];
      real_t * recv_values = new (std::nothrow) real_t[(*r).second.size() * num_params_];
      recv_bufs.push_back(recv_fitness); recv_bufs.push_back(recv_values);
      MPI_Request r1, r2;
      (*multi_node_).irecv(root_comm_, recv_fitness, (*r).second.size(), (*r).first, r1);
      (*multi_node_).irecv(root_comm_, recv_values, (*r).second.size() * num_params_,
                  (*r).first, r2);
      rreq.push_back(r1); rreq.push_back(r2);
    } // for

    // Isend
    // first send fitness values, then send parameter values
    std::vector <MPI_Request> sreq;
    std::vector<real_t*> send_bufs;
    typedef std::map <int, std::set <unsigned int> >::const_iterator comm_list_iter;
    for(comm_list_iter s = neighbor_send_lists_.begin(); s != neighbor_send_lists_.end(); ++ s) {
//      if(myrank == (*s).first || (*s).second.size() == 0) continue;
      real_t * send_fitness = new (std::nothrow) real_t[(*s).second.size()];
      real_t * send_values = new (std::nothrow) real_t[(*s).second.size() * num_params_];
      send_bufs.push_back(send_fitness); send_bufs.push_back(send_values);
      std::set<unsigned int>::const_iterator iter = (*s).second.begin();
      for(int i = 0; iter != (*s).second.end(); ++ i, ++ iter) {
        unsigned int index = *iter - start_particle_index_;
        send_fitness[i] = particles_[index].best_fitness_global_;
        for(int j = 0; j < num_params_; ++ j) {
          send_values[i * num_params_ + j] = particles_[index].best_values_global_[j];
        } // for
      } // for
      MPI_Request s1, s2;
      (*multi_node_).isend(root_comm_, send_fitness, (*s).second.size(), (*s).first, s1);
      (*multi_node_).isend(root_comm_, send_values, (*s).second.size() * num_params_,
                  (*s).first, s2);
      sreq.push_back(s1); sreq.push_back(s2);
    } // for

    // wait for all recvs to complete
    (*multi_node_).waitall(root_comm_, rreq.size(), &rreq[0]);
    // then wait for all sends to complete
    (*multi_node_).waitall(root_comm_, sreq.size(), &sreq[0]);

    // for each particle update data using received neighbor information
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << "--------------- Particle " << i << "(" << particles_[i].index_ << "): ";
      // for each neighbor
      //for(particle_neighbor_set_t::iterator n = particles_[i].neighbors_.begin();
      for(std::set<PSOParticle::ParticleNeighbor>::iterator n = particles_[i].neighbors_.begin();
          n != particles_[i].neighbors_.end(); ++ n) {
        std::cout << "nbr " << (*n).index_ << " ";
        bool found = false;
        comm_list_iter_t r = neighbor_recv_lists_.begin();
        for(int j = 0; r != neighbor_recv_lists_.end(); ++ r, ++ j) {
          if((*r).first != (*n).proc_) continue;
          std::set<unsigned int>::iterator f = (*r).second.begin();
          for(int k = 0; f != (*r).second.end(); ++ k, ++ f) {
            if((*n).index_ == *f) {    // this is what we are looking for
              found = true;
              std::cout << (*multi_node_).rank(root_comm_) << ": comparing "
                    << particles_[i].best_fitness_global_ << " with "
                    << recv_bufs[2 * j][k] << std::endl;
              if(particles_[i].best_fitness_global_ > recv_bufs[2 * j][k]) {
                particles_[i].best_fitness_global_ = recv_bufs[2 * j][k];
                for(int l = 0; l < num_params_; ++ l)
                  particles_[i].best_values_global_[l] =
                              recv_bufs[2 * j + 1][k * num_params_ + l];
              } // if
              break;
            } // if
          } // for k
          if(found) break;
        } // for recv list
      } // for neighbors
    } // for particle

    // free all buffers
    for(std::vector<real_t*>::iterator i = recv_bufs.begin(); i != recv_bufs.end(); ++ i)
      delete[] *i;
    for(std::vector<real_t*>::iterator i = send_bufs.begin(); i != send_bufs.end(); ++ i)
      delete[] *i;

    return true;
  } // ParticleSwarmOptimization::neighbor_data_exchange()


  // lbest, with K = 2 in a ring
  bool ParticleSwarmOptimization::simulate_lbest_generation() {
    int k = 2;
    if(k != 2) {
      std::cout << "error: case k = " << k << " not implemented yet!" << std::endl;
      return false;
    } // if
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // construct param map
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j)
        curr_particle.push_back(particles_[i].param_values_[j]);
      // compute the fitness
      (*obj_func_).update_sim_comm(particle_comm_);
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);
      // update particle fitness
      if(particles_[i].best_fitness_ > curr_fitness[0]) {
        particles_[i].best_fitness_ = curr_fitness[0];
        particles_[i].best_values_ = particles_[i].param_values_;
      } // if
      if(particles_[i].best_fitness_global_ > curr_fitness[0]) {
        particles_[i].best_fitness_global_ = curr_fitness[0];
        particles_[i].best_values_global_ = particles_[i].param_values_;
      } // if
      std::cout << (*multi_node_).rank(root_comm_) << ": " << curr_fitness[0] << " "
            << particles_[i].best_fitness_ << " " << particles_[i].best_fitness_global_
            << std::endl;
    } // for

    #ifdef USE_MPI
      // exchange data between neighbors
      if(!neighbor_data_exchange()) return false;
    #else
      // TODO ...
      std::cerr << "error: this has not been implemented!!!" << std::endl;
      return false;
    #endif

    // update each particle using new best values
    for(int i = 0; i < num_particles_; ++ i) {
      if(!particles_[i].update_particle(pso_omega_, pso_phi1_, pso_phi2_,
                        particles_[i].best_values_global_,
                        constraints_, rand_)) {
        std::cerr << "error: failed to update particle " << i << std::endl;
        return false;
      } // if
      // also keep track of global best :P -- not really needed ...
      if(best_fitness_ > particles_[i].best_fitness_global_) {
        best_fitness_ = particles_[i].best_fitness_global_;
        best_values_ = particles_[i].best_values_global_;
      } // if
    } // for

    if(tune_omega_) pso_omega_ /= 2.0;
    std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

    #ifdef USE_MPI
    std::cout << (*multi_node_).rank(root_comm_) << " :: ";
    #endif
    std::cout << "@@@@@@ Global best: ";
    std::cout << best_fitness_ << " [ ";
    for(int j = 0; j < num_params_; ++ j) std::cout << best_values_[j] << " ";
    std::cout << "]\t";
    std::cout << std::endl;

    return true;
  } // ParticleSwarmOptimization::simulate_lbest_generation()


  // vonnewmann, with K = 4 in a grid/torus
  bool ParticleSwarmOptimization::simulate_vonnewmann_generation() {
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // construct param map
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j)
        curr_particle.push_back(particles_[i].param_values_[j]);
      // compute the fitness
      (*obj_func_).update_sim_comm(particle_comm_);
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);
      // update particle fitness
      if(particles_[i].best_fitness_ > curr_fitness[0]) {
        particles_[i].best_fitness_ = curr_fitness[0];
        particles_[i].best_values_ = particles_[i].param_values_;
      } // if
      if(particles_[i].best_fitness_global_ > curr_fitness[0]) {
        particles_[i].best_fitness_global_ = curr_fitness[0];
        particles_[i].best_values_global_ = particles_[i].param_values_;
      } // if
    } // for

    #ifdef USE_MPI
      // exchange data between neighbors
      if(!neighbor_data_exchange()) return false;
    #else
      // TODO ...
      std::cerr << "error: this has not been implemented!!!" << std::endl;
      return false;
    #endif

    // update each particle using new best values
    for(int i = 0; i < num_particles_; ++ i) {
      if(!particles_[i].update_particle(pso_omega_, pso_phi1_, pso_phi2_,
                        particles_[i].best_values_global_,
                        constraints_, rand_)) {
        std::cerr << "error: failed to update particle " << i << std::endl;
        return false;
      } // if
      // also keep track of global best :P -- not really needed ...
      if(best_fitness_ > particles_[i].best_fitness_global_) {
        best_fitness_ = particles_[i].best_fitness_global_;
        best_values_ = particles_[i].best_values_global_;
      } // if
    } // for

    if(tune_omega_) pso_omega_ /= 2.0;
    std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

    #ifdef USE_MPI
    std::cout << (*multi_node_).rank(root_comm_) << " :: ";
    #endif
    std::cout << "@@@@@@ Global best: ";
    std::cout << best_fitness_ << " [ ";
    for(int j = 0; j < num_params_; ++ j) std::cout << best_values_[j] << " ";
    std::cout << "]\t";
    std::cout << std::endl;

    return true;
  } // ParticleSwarmOptimization::simulate_vonnewmann_generation()


  // random topology
  bool ParticleSwarmOptimization::simulate_random_generation() {
    // for each particle, simulate
    for(int i = 0; i < num_particles_; ++ i) {
      std::cout << (*multi_node_).rank(root_comm_) << ": Particle " << i << std::endl;
      // construct param map
      real_vec_t curr_particle;
      for(int j = 0; j < num_params_; ++ j)
        curr_particle.push_back(particles_[i].param_values_[j]);
      // compute the fitness
      (*obj_func_).update_sim_comm(particle_comm_);
      real_vec_t curr_fitness = (*obj_func_)(curr_particle);
      // update particle fitness
      if(particles_[i].best_fitness_ > curr_fitness[0]) {
        particles_[i].best_fitness_ = curr_fitness[0];
        particles_[i].best_values_ = particles_[i].param_values_;
      } // if
      if(particles_[i].best_fitness_global_ > curr_fitness[0]) {
        particles_[i].best_fitness_global_ = curr_fitness[0];
        particles_[i].best_values_global_ = particles_[i].param_values_;
      } // if
    } // for

    #ifdef USE_MPI
      // exchange data between neighbors
      if(!neighbor_data_exchange()) return false;
    #else
      // TODO ...
      std::cerr << "error: this has not been implemented!!!" << std::endl;
      return false;
    #endif

    // update each particle using new best values
    for(int i = 0; i < num_particles_; ++ i) {
      if(!particles_[i].update_particle(pso_omega_, pso_phi1_, pso_phi2_,
                        particles_[i].best_values_global_,
                        constraints_, rand_)) {
        std::cerr << "error: failed to update particle " << i << std::endl;
        return false;
      } // if
      // also keep track of global best :P -- not really needed ...
      if(best_fitness_ > particles_[i].best_fitness_global_) {
        best_fitness_ = particles_[i].best_fitness_global_;
        best_values_ = particles_[i].best_values_global_;
      } // if
    } // for

    if(tune_omega_) pso_omega_ /= 2.0;
    std::cout << (*multi_node_).rank(root_comm_) << ": Omega = " << pso_omega_ << std::endl;

    #ifdef USE_MPI
    std::cout << (*multi_node_).rank(root_comm_) << " :: ";
    #endif
    std::cout << "@@@@@@ Global best: ";
    std::cout << best_fitness_ << " [ ";
    for(int j = 0; j < num_params_; ++ j) std::cout << best_values_[j] << " ";
    std::cout << "]\t";
    std::cout << std::endl;

    return true;
  } // ParticleSwarmOptimization::simulate_random_generation()



} // namespace hig
