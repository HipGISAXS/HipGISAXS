/***
  *  Project:
  *
  *  File: multi_node_comm.hpp
  *  Created: Mar 18, 2013
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __MULTI_NODE_COMM_HPP__
#define __MULTI_NODE_COMM_HPP__

#ifdef USE_MPI

#include <iostream>
#include <mpi.h>
#include <complex>
#include <map>
#include <vector>
#include <string>
#include <iostream>

namespace woo {

	static const int MASTER_RANK = 0;

	/**
	 * Some enums and types
	 */
	namespace comm {

		enum ReduceOp {
			min,
			max,
			sum,
			prod,
			land,
			band,
			lor,
			bor,
			lxor,
			bxor,
			minloc,
			maxloc
		}; // enum CommOp

	} // namespace comm
	typedef std::string comm_t;


	/**
	 * Communicator
	 */
	class MultiNodeComm {

		public:

			MultiNodeComm() :
				valid_(false), num_procs_(0), idle_(false), master_rank_(MASTER_RANK), rank_(0) {
				// create communicator ...
			} // MultiNodeComm()

			MultiNodeComm(const MPI_Comm& comm) {
				world_ = comm;
				num_procs_ = size();
				valid_ = true;
				idle_ = false;
				master_rank_ = MASTER_RANK;
				MPI_Comm_rank(world_, &rank_);
			} // operator=()

			~MultiNodeComm() {
				// destroy communicator ...
			} // ~MultiNodeComm()

		protected:

			inline int size() const { int size; MPI_Comm_size(world_, &size); return size; }
			inline int rank() const { return rank_; }
			inline bool is_idle() const { return idle_; }
			inline bool is_valid() const { return valid_; }
			inline int master() const { return master_rank_; }
			inline int is_master() const { return (master_rank_ == rank_); }

			MultiNodeComm& operator=(const MPI_Comm& comm) {
				if(valid_) {
					std::cerr << "warning: assigning communicator to an already valid communicator. "
								<< "proceeding as is." << std::endl;
				} // if
				world_ = comm;
				num_procs_ = size();
				valid_ = true;
				idle_ = false;
				master_rank_ = MASTER_RANK;
				MPI_Comm_rank(world_, &rank_);
        return *this;
			} // operator=()

			MultiNodeComm& operator=(const MultiNodeComm& comm) {
				if(valid_) {
					std::cerr << "warning: assigning communicator to an already valid communicator. "
								<< "proceeding as is." << std::endl;
				} // if
				world_ = comm.world_;
				num_procs_ = size();
				valid_ = comm.valid_;
				idle_ = comm.idle_;
				master_rank_ = comm.master_rank_;
				rank_ = comm.rank_;
        return *this;
			} // operator=()

			inline void set_idle() { idle_ = true; }

			inline MPI_Op reduce_op_map(comm::ReduceOp op) {
				switch(op) {
					case comm::min: return MPI_MIN;
					case comm::max: return MPI_MAX;
					case comm::sum: return MPI_SUM;
					case comm::prod: return MPI_PROD;
					case comm::land: return MPI_LAND;
					case comm::band: return MPI_BAND;
					case comm::lor: return MPI_LOR;
					case comm::bor: return MPI_BOR;
					case comm::lxor: return MPI_LXOR;
					case comm::bxor: return MPI_BXOR;
					case comm::minloc: return MPI_MINLOC;
					case comm::maxloc: return MPI_MAXLOC;
					default: return NULL;
				} // switch
			} // reduce_op_map()

			inline MultiNodeComm split(int color) {
				MPI_Comm new_comm;
				MPI_Comm_split(world_, color, rank_, &new_comm);
				MultiNodeComm new_mncomm(new_comm);
				return new_mncomm;
			} // split()

			inline MultiNodeComm dup() {
				MPI_Comm new_comm;
				MPI_Comm_dup(world_, &new_comm);
				MultiNodeComm new_mncomm(new_comm);
				int rankie;
				MPI_Comm_rank(new_comm, &rankie);
				return new_mncomm;
			} // dup()

			inline bool free() {
				MPI_Comm_free(&world_);
				valid_ = false;
				return true;
			} // free()

			/*inline MultiNodeComm* dup() {
				MPI_Comm new_comm;
				MPI_Comm_dup(world_, &new_comm);
				MultiNodeComm *new_mncomm = new MultiNodeComm(new_comm);
				int rankie;
				MPI_Comm_rank(new_comm, &rankie);
				return new_mncomm;
			} // dup()*/

			inline bool broadcast(float& data) {
				MPI_Bcast(&data, 1, MPI_FLOAT, master_rank_, world_);
				return true;
			} // broadcast()

			inline bool broadcast(double& data) {
				MPI_Bcast(&data, 1, MPI_DOUBLE, master_rank_, world_);
				return true;
			} // broadcast()

			inline bool broadcast(float* data, int size) {
				MPI_Bcast(&(*data), size, MPI_FLOAT, master_rank_, world_);
				return true;
			} // broadcast()

			inline bool broadcast(double* data, int size) {
				MPI_Bcast(&(*data), size, MPI_DOUBLE, master_rank_, world_);
				return true;
			} // broadcast()

			inline bool broadcast(unsigned int* data, int size) {
				MPI_Bcast(&(*data), size, MPI_UNSIGNED, master_rank_, world_);
				return true;
			} // broadcast()

			inline bool broadcast(float& data, int root) {
				if(MPI_Bcast(&data, 1, MPI_FLOAT, root, world_) != MPI_SUCCESS) return false;
				return true;
			} // broadcast()

			inline bool broadcast(double& data, int root) {
				if(MPI_Bcast(&data, 1, MPI_DOUBLE, root, world_) != MPI_SUCCESS) return false;
				return true;
			} // broadcast()

			inline bool broadcast(float* data, int size, int root) {
				if(MPI_Bcast(&(*data), size, MPI_FLOAT, root, world_) != MPI_SUCCESS) return false;
				return true;
			} // broadcast()

			inline bool broadcast(double* data, int size, int root) {
				if(MPI_Bcast(&(*data), size, MPI_DOUBLE, root, world_) != MPI_SUCCESS) return false;
				return true;
			} // broadcast()

			inline bool allreduce(float sendval, float& recvval, int& proc_rank, comm::ReduceOp op) {
				if(op != comm::minloc && op != comm::maxloc) {
					std::cerr << "error: invalid reduction operation" << std::endl;
					return false;
				} // if
				MPI_Op mpi_op = reduce_op_map(op);
				struct {
					float val;
					int rank;
				} send, recv;
				send.val = sendval;
				send.rank = rank_;
				if(MPI_Allreduce(&send, &recv, 1, MPI_FLOAT_INT, mpi_op, world_) != MPI_SUCCESS)
					return false;
				recvval = recv.val;
				proc_rank = recv.rank;
				return true;
			} // allreduce()

			inline bool allreduce(double sendval, double & recvval, int& proc_rank, comm::ReduceOp op) {
				if(op != comm::minloc && op != comm::maxloc) {
					std::cerr << "error: invalid reduction operation" << std::endl;
					return false;
				} // if
				MPI_Op mpi_op = reduce_op_map(op);
				struct {
					double val;
					int rank;
				} send, recv;
				send.val = sendval;
				send.rank = rank_;
				if(MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE_INT, mpi_op, world_) != MPI_SUCCESS)
					return false;
				recvval = recv.val;
				proc_rank = recv.rank;
				return true;
			} // allreduce()

			inline bool scan_sum(unsigned int in, unsigned int& out) {
				if(MPI_Scan(&in, &out, 1, MPI_UNSIGNED, MPI_SUM, world_) != MPI_SUCCESS)
					return false;
				return true;
			} // scan_sum()

			inline bool gather(int* sbuf, int scount, int* rbuf, int rcount) {
				MPI_Gather(sbuf, scount, MPI_INT, rbuf, rcount, MPI_INT, master_rank_, world_);
				return true;
			} // gather()

			inline bool gather(float* sbuf, int scount, float* rbuf, int rcount) {
				MPI_Gather(sbuf, scount, MPI_FLOAT, rbuf, rcount, MPI_FLOAT, master_rank_, world_);
				return true;
			} // gather()

			inline bool gather(double* sbuf, int scount, double* rbuf, int rcount) {
				MPI_Gather(sbuf, scount, MPI_DOUBLE, rbuf, rcount, MPI_DOUBLE, master_rank_, world_);
				return true;
			} // gather()

			inline bool allgather(int* sbuf, int scount, int* rbuf, int rcount) {
				MPI_Allgather(sbuf, scount, MPI_INT, rbuf, rcount, MPI_INT, world_);
				return true;
			} // allgather()

			inline bool allgather(unsigned int* sbuf, int scount, unsigned int* rbuf, int rcount) {
				MPI_Allgather(sbuf, scount, MPI_UNSIGNED, rbuf, rcount, MPI_UNSIGNED, world_);
				return true;
			} // allgather()

			inline bool allgatherv(float* sbuf, int scount, float* rbuf, int* rcount) {
				int* displs = new (std::nothrow) int[num_procs_];
				displs[0] = 0;
				for(int i = 1; i < num_procs_; ++ i) displs[i] = displs[i - 1] + rcount[i - 1];
				MPI_Allgatherv(sbuf, scount, MPI_FLOAT, rbuf, rcount, displs, MPI_FLOAT, world_);
				delete[] displs;
				return true;
			} // allgatherv()

			inline bool allgatherv(double* sbuf, int scount, double* rbuf, int* rcount) {
				int* displs = new (std::nothrow) int[num_procs_];
				displs[0] = 0;
				for(int i = 1; i < num_procs_; ++ i) displs[i] = displs[i - 1] + rcount[i - 1];
				MPI_Allgatherv(sbuf, scount, MPI_DOUBLE, rbuf, rcount, displs, MPI_DOUBLE, world_);
				delete[] displs;
				return true;
			} // allgatherv()

			inline bool gatherv(float* sbuf, int scount, float* rbuf, int* rcount, int* displs) {
				MPI_Gatherv(sbuf, scount, MPI_FLOAT, rbuf, rcount, displs, MPI_FLOAT,
							master_rank_, world_);
				return true;
			} // gatherv()

			inline bool gatherv(double* sbuf, int scount, double* rbuf, int* rcount, int* displs) {
				MPI_Gatherv(sbuf, scount, MPI_DOUBLE, rbuf, rcount, displs, MPI_DOUBLE,
							master_rank_, world_);
				return true;
			} // gatherv()

			inline bool gatherv(std::complex<float>* sbuf, int scount,
									std::complex<float>* rbuf, int* rcount, int* displs) {
				MPI_Gatherv(sbuf, scount, MPI_COMPLEX, rbuf, rcount, displs, MPI_COMPLEX,
							master_rank_, world_);
				return true;
			} // gatherv()

			inline bool gatherv(std::complex<double>* sbuf, int scount,
									std::complex<double>* rbuf, int* rcount, int* displs) {
				MPI_Gatherv(sbuf, scount, MPI_DOUBLE_COMPLEX,
								rbuf, rcount, displs, MPI_DOUBLE_COMPLEX, master_rank_, world_);
				return true;
			} // gatherv()

			inline bool scatter(int* sbuf, int scount, int* rbuf, int rcount) {
				MPI_Scatter(sbuf, scount, MPI_INT, rbuf, rcount, MPI_INT, master_rank_, world_);
				return true;
			} // scatter()

			inline bool scatterv(int* sbuf, int* scounts, int* displs, int* rbuf, int rcount) {
				MPI_Scatterv(sbuf, scounts, displs, MPI_INT, rbuf, rcount, MPI_INT, master_rank_, world_);
				return true;
			} // scatterv()

			inline bool barrier() {
				MPI_Barrier(world_);
				return true;
			} // barrier()

			/**
			 * Point-to-point
			 */

			bool isend(float* sbuf, int scount, int to, MPI_Request& req) {
				if(MPI_Isend(sbuf, scount, MPI_FLOAT, to, rank_, world_, &req) != MPI_SUCCESS)
					return false;
				return true;
			} // send()

			bool isend(double * sbuf, int scount, int to, MPI_Request& req) {
				if(MPI_Isend(sbuf, scount, MPI_DOUBLE, to, rank_, world_, &req) != MPI_SUCCESS)
					return false;
				return true;
			} // send()

			bool irecv(float* rbuf, int rcount, int from, MPI_Request& req) {
				if(MPI_Irecv(rbuf, rcount, MPI_FLOAT, from, from, world_, &req) != MPI_SUCCESS)
					return false;
				return true;
			} // send()

			bool irecv(double * rbuf, int rcount, int from, MPI_Request& req) {
				if(MPI_Irecv(rbuf, rcount, MPI_DOUBLE, from, from, world_, &req) != MPI_SUCCESS)
					return false;
				return true;
			} // send()

			bool waitall(int count, MPI_Request* req) {
				MPI_Status* stats = new (std::nothrow) MPI_Status[count];
				MPI_Waitall(count, req, stats);
				delete[] stats;
				return true;
			} // waitall()

		private:
			unsigned int num_procs_;	// number of processes in this communicator
			MPI_Comm world_;			// MPI communicator for this communicator
			bool valid_;
			bool idle_;					// status of a node
			int rank_;					// my rank
			int master_rank_;			// who is the master here

			friend class MultiNode;

	}; // class MultiNodeComm


	//typedef std::map <const char*, MultiNodeComm> multi_node_comm_map_t;
	typedef std::map <comm_t, MultiNodeComm> multi_node_comm_map_t;

	/**
	 * For communication - TODO: make it singleton
	 */
	class MultiNode {

		public:
			MultiNode(int narg, char** args): universe_key_("world") {
				MPI_Init(&narg, &args);
				comms_.clear();
				MultiNodeComm universe(MPI_COMM_WORLD);
				universe_num_procs_ = universe.size();
				comms_[universe_key_] = universe;
			} // MultiNodeComm()

			MultiNode(int narg, char** args, const comm_t& univ_key): universe_key_(univ_key) {
				MPI_Init(&narg, &args);
				comms_.clear();
				MultiNodeComm universe(MPI_COMM_WORLD);
				universe_num_procs_ = universe.size();
				comms_[universe_key_] = universe;
			} // MultiNodeComm()

			~MultiNode() {
				MPI_Finalize();
			} // ~MultiNodeComm()


			bool init() {
				return true;
			} // init()

			// number of procs in the world
			/*inline int size() { return comms_.at("world").size(); }
			inline int rank() { return comms_.at("world").rank(); }
			inline bool is_master() { return (comms_.at("world").rank() == comms_.at("world").master()); }
			inline bool is_idle() { return comms_.at("world").is_idle(); }
			inline int master() { return comms_.at("world").master(); }*/

			inline comm_t universe_key() const { return universe_key_; }

			inline int size(comm_t key) { return comms_.at(key).size(); }
			inline int rank(comm_t key) { return comms_.at(key).rank(); }
			inline bool is_master(comm_t key) { return comms_.at(key).is_master(); }
			inline bool is_idle(comm_t key) { return comms_.at(key).is_idle(); }
			inline int master(comm_t key) { return comms_.at(key).master(); }

			inline void set_idle(comm_t key) { comms_.at(key).set_idle(); }

			// number of communicators including the world
			inline int num_comms() const { return comms_.size(); }

			/**
			 * Create new communicators
			 */

			inline bool split(comm_t new_key, comm_t key, int color) {
				comms_[new_key] = comms_[key].split(color);
				return true;
			} // split()

			inline bool dup(comm_t new_key, comm_t key) {
				comms_[new_key] = comms_[key].dup();
				return true;
			} // dup()

			inline bool free(comm_t key) {
				comms_[key].free();
				comms_.erase(key);
				return true;
			} // free

			/**
			 * Broadcasts
			 */

			inline bool broadcast(comm_t key, float& data) {
				return comms_[key].broadcast(data);
			} // send_broadcast()

			inline bool broadcast(comm_t key, double& data) {
				return comms_[key].broadcast(data);
			} // send_broadcast()

			inline bool broadcast(comm_t key, float* data, int size) {
				return comms_[key].broadcast(data, size);
			} // send_broadcast()

			inline bool broadcast(comm_t key, double* data, int size) {
				return comms_[key].broadcast(data, size);
			} // send_broadcast()

			inline bool broadcast(comm_t key, unsigned int* data, int size) {
				return comms_[key].broadcast(data, size);
			} // send_broadcast()

			inline bool broadcast(comm_t key, float& data, int rank) {
				return comms_[key].broadcast(data, rank);
      } // broadcast()

			inline bool broadcast(comm_t key, double& data, int rank) {
				return comms_[key].broadcast(data, rank);
      } // broadcast()

			inline bool broadcast(comm_t key, std::vector<float>& data, int rank) {
				float* temp_data = new (std::nothrow) float[data.size()];
				for(int i = 0; i < data.size(); ++ i) temp_data[i] = data[i];
				bool success = comms_[key].broadcast(temp_data, data.size(), rank);
				for(int i = 0; i < data.size(); ++ i) data[i] = temp_data[i];
				delete[] temp_data;
				return success;
			} // broadcast()

			inline bool broadcast(comm_t key, std::vector<double>& data, int rank) {
				double* temp_data = new (std::nothrow) double[data.size()];
				for(int i = 0; i < data.size(); ++ i) temp_data[i] = data[i];
				bool success = comms_[key].broadcast(temp_data, data.size(), rank);
				for(int i = 0; i < data.size(); ++ i) data[i] = temp_data[i];
				delete[] temp_data;
				return success;
			} // broadcast()

			inline bool allreduce(comm_t key, float sendval, float& recvval,
									int& rank, comm::ReduceOp op) {
				return comms_[key].allreduce(sendval, recvval, rank, op);
			} // allreduce()

			inline bool allreduce(comm_t key, double sendval, double& recvval,
									int& rank, comm::ReduceOp op) {
				return comms_[key].allreduce(sendval, recvval, rank, op);
			} // allreduce()

			/**
			 * Scans
			 */

			inline bool scan_sum(comm_t key, unsigned int in, unsigned int& out) {
				return comms_[key].scan_sum(in, out);
			} // scan_sum()

			/**
			 * Gathers
			 */

			inline bool allgather(comm_t key, int* sbuf, int scount, int* rbuf, int rcount) {
				return comms_[key].allgather(sbuf, scount, rbuf, rcount);
			} // allgather()

			inline bool allgather(comm_t key, unsigned int* sbuf, int scount,
									unsigned int* rbuf, int rcount) {
				return comms_[key].allgather(sbuf, scount, rbuf, rcount);
			} // allgather()

			inline bool allgatherv(comm_t key, std::vector<float>& sbuf,
									std::vector<float>& rbuf) {
				int size = sbuf.size();
				int* rsize = new (std::nothrow) int[comms_[key].size()];
				comms_[key].allgather(&size, 1, rsize, 1);
				int recv_tot = 0; for(int i = 0; i < comms_[key].size(); ++ i) recv_tot += rsize[i];
				float* temp_rbuf = new (std::nothrow) float[recv_tot];
				rbuf.clear();
				if(!comms_[key].allgatherv(&sbuf[0], size, temp_rbuf, rsize)) return false;
				for(int i = 0; i < recv_tot; ++ i) rbuf.push_back(temp_rbuf[i]);
				delete[] temp_rbuf;
				delete[] rsize;
				return true;
			} // allgather()

			inline bool allgatherv(comm_t key, std::vector<double>& sbuf,
									std::vector<double>& rbuf) {
				int size = sbuf.size();
				int* rsize = new (std::nothrow) int[comms_[key].size()];
				comms_[key].allgather(&size, 1, rsize, 1);
				int recv_tot = 0; for(int i = 0; i < comms_[key].size(); ++ i) recv_tot += rsize[i];
				double* temp_rbuf = new (std::nothrow) double[recv_tot];
				rbuf.clear();
				if(!comms_[key].allgatherv(&sbuf[0], size, temp_rbuf, rsize)) return false;
				for(int i = 0; i < recv_tot; ++ i) rbuf.push_back(temp_rbuf[i]);
				delete[] temp_rbuf;
				delete[] rsize;
				return true;
			} // allgather()

			inline bool gather(comm_t key, int* sbuf, int scount, int* rbuf, int rcount) {
				return comms_[key].gather(sbuf, scount, rbuf, rcount);
			} // gather()

			inline bool gather(comm_t key, float* sbuf, int scount, float* rbuf, int rcount) {
				return comms_[key].gather(sbuf, scount, rbuf, rcount);
			} // gather()

			inline bool gather(comm_t key, double* sbuf, int scount, double* rbuf, int rcount) {
				return comms_[key].gather(sbuf, scount, rbuf, rcount);
			} // gather()

			inline bool gatherv(comm_t key, float* sbuf, int scount,
									float* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			inline bool gatherv(comm_t key, double* sbuf, int scount,
									double* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			inline bool gatherv(comm_t key, std::complex<float>* sbuf, int scount,
									std::complex<float>* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			inline bool gatherv(comm_t key, std::complex<double>* sbuf, int scount,
									std::complex<double>* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			/**
			 * Scatters
			 */

			inline bool scatter(comm_t key, int* sbuf, int scount, int* rbuf, int rcount) {
				return comms_[key].scatter(sbuf, scount, rbuf, rcount);
			} // scatter()

			inline bool scatterv(comm_t key, int* sbuf, int* scounts, int* displs,
									int* rbuf, int rcount) {
				return comms_[key].scatterv(sbuf, scounts, displs, rbuf, rcount);
			} // scatterv()

			/**
			 * Barrier
			 */

			/*bool barrier() {
				return comms_["world"].barrier();
			} // barrier()*/

			inline bool barrier(comm_t key) {
				return comms_[key].barrier();
			} // barrier()

			/**
			 * Point-to-point
			 */

			inline bool isend(comm_t key, float * sbuf, int scount, int to, MPI_Request& req) {
				return comms_[key].isend(sbuf, scount, to, req);
			} // send()

			inline bool isend(comm_t key, double * sbuf, int scount, int to, MPI_Request& req) {
				return comms_[key].isend(sbuf, scount, to, req);
			} // send()

			inline bool irecv(comm_t key, float * rbuf, int rcount, int from, MPI_Request& req) {
				return comms_[key].irecv(rbuf, rcount, from, req);
			} // send()

			inline bool irecv(comm_t key, double * rbuf, int rcount, int from, MPI_Request& req) {
				return comms_[key].irecv(rbuf, rcount, from, req);
			} // send()

			inline bool waitall(comm_t key, int count, MPI_Request* req) {
				return comms_[key].waitall(count, req);
			} // waitall()

		private:
			unsigned int universe_num_procs_;		// total number of processes
			comm_t universe_key_;
			multi_node_comm_map_t comms_;		// list of all communicators in the world

	}; // class MultiNode

} // namespace woo

#else

namespace woo {

	typedef std::string comm_t;

} // namespace woo

#endif // USE_MPI

#endif // __MULTI_NODE_COMM_HPP__
