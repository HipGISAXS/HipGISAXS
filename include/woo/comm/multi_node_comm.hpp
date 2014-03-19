/***
  *  Project:
  *
  *  File: multi_node_comm.hpp
  *  Created: Mar 18, 2013
  *  Modified: Mon 17 Mar 2014 06:50:30 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __MULTI_NODE_COMM_HPP__
#define __MULTI_NODE_COMM_HPP__

#ifdef USE_MPI

#include <mpi.h>
#include <complex>
#include <map>
#include <vector>
#include <string>

namespace woo {

	static const int MASTER_RANK = 0;

	/**
	 * Some enums
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

			inline bool broadcast(float* data, int size, int root) {
				if(MPI_Bcast(&(*data), size, MPI_FLOAT, root, world_) != MPI_SUCCESS) return false;
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

			inline bool barrier() {
				MPI_Barrier(world_);
				return true;
			} // barrier()

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
	typedef std::map <std::string, MultiNodeComm> multi_node_comm_map_t;

	/**
	 * For communication - TODO: make it singleton
	 */
	class MultiNode {

		public:
			MultiNode(int narg, char** args) {
				MPI_Init(&narg, &args);
				comms_.clear();
				MultiNodeComm universe(MPI_COMM_WORLD);
				world_num_procs_ = universe.size();
				comms_["world"] = universe;
			} // MultiNodeComm()

			~MultiNode() {
				MPI_Finalize();
			} // ~MultiNodeComm()


			bool init() {
				return true;
			} // init()

			// number of procs in the world
			inline int size() { return comms_.at("world").size(); }
			inline int rank() { return comms_.at("world").rank(); }
			inline bool is_master() { return (comms_.at("world").rank() == comms_.at("world").master()); }
			inline bool is_idle() { return comms_.at("world").is_idle(); }
			inline int master() { return comms_.at("world").master(); }

			inline int size(const char* key) { return comms_.at(key).size(); }
			inline int rank(const char* key) { return comms_.at(key).rank(); }
			inline bool is_master(const char* key) { return comms_.at(key).is_master(); }
			inline bool is_idle(const char* key) { return comms_.at(key).is_idle(); }
			inline int master(const char* key) { return comms_.at(key).master(); }

			inline void set_idle(const char* key) { comms_.at(key).set_idle(); }

			// number of communicators including the world
			inline int num_comms() const { return comms_.size(); }

			/**
			 * Create new communicators
			 */

			bool split(const char* new_key, const char* key, int color) {
				comms_[new_key] = comms_[key].split(color);
				return true;
			} // split()

			bool dup(const char* new_key, const char* key) {
				comms_[new_key] = comms_[key].dup();
				return true;
			} // dup()

			bool free(const char* key) {
				comms_[key].free();
				comms_.erase(key);
				return true;
			} // free

			/**
			 * Broadcasts
			 */

			bool broadcast(const char* key, float* data, int size) {
				return comms_[key].broadcast(data, size);
			} // send_broadcast()

			bool broadcast(const char* key, double* data, int size) {
				return comms_[key].broadcast(data, size);
			} // send_broadcast()

			bool broadcast(const char* key, unsigned int* data, int size) {
				return comms_[key].broadcast(data, size);
			} // send_broadcast()

			bool broadcast(const char* key, std::vector<float>& data, int rank) {
				float* temp_data = new (std::nothrow) float[data.size()];
				for(int i = 0; i < data.size(); ++ i) temp_data[i] = data[i];
				bool success = comms_[key].broadcast(temp_data, data.size(), rank);
				for(int i = 0; i < data.size(); ++ i) data[i] = temp_data[i];
				delete[] temp_data;
				return success;
			} // broadcast()

			bool allreduce(const char* key, float sendval, float& recvval, int& rank, comm::ReduceOp op) {
				return comms_[key].allreduce(sendval, recvval, rank, op);
			} // allreduce()

			/**
			 * Scans
			 */

			bool scan_sum(const char* key, unsigned int in, unsigned int& out) {
				return comms_[key].scan_sum(in, out);
			} // scan_sum()

			/**
			 * Gathers
			 */

			bool allgather(const char* key, int* sbuf, int scount, int* rbuf, int rcount) {
				return comms_[key].allgather(sbuf, scount, rbuf, rcount);
			} // allgather()

			inline bool gather(const char* key, int* sbuf, int scount, int* rbuf, int rcount) {
				return comms_[key].gather(sbuf, scount, rbuf, rcount);
			} // gather()

			inline bool gather(const char* key, float* sbuf, int scount, float* rbuf, int rcount) {
				return comms_[key].gather(sbuf, scount, rbuf, rcount);
			} // gather()

			inline bool gather(const char* key, double* sbuf, int scount, double* rbuf, int rcount) {
				return comms_[key].gather(sbuf, scount, rbuf, rcount);
			} // gather()

			inline bool gatherv(const char* key, float* sbuf, int scount,
									float* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			inline bool gatherv(const char* key, double* sbuf, int scount,
									double* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			inline bool gatherv(const char* key, std::complex<float>* sbuf, int scount,
									std::complex<float>* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			inline bool gatherv(const char* key, std::complex<double>* sbuf, int scount,
									std::complex<double>* rbuf, int* rcount, int* displs) {
				return comms_[key].gatherv(sbuf, scount, rbuf, rcount, displs);
			} // gatherv()

			/**
			 * Barrier
			 */

			bool barrier() {
				return comms_["world"].barrier();
			} // barrier()

			bool barrier(const char* key) {
				return comms_[key].barrier();
			} // barrier()

		private:
			unsigned int world_num_procs_;		// total number of processes
			multi_node_comm_map_t comms_;		// list of all communicators in the world

	}; // class MultiNode

} // namespace woo

#endif // USE_MPI

#endif // __MULTI_NODE_COMM_HPP__
