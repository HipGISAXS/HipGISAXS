/**
 *  $Id: ff_num.hpp 47 2012-08-23 21:05:16Z asarje $
 *
 *  Project: HipGISAXS
 *
 *  File: ff_num.hpp
 *  Created: Nov 05, 2011
 *  Modified: Thu 23 Aug 2012 02:02:59 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef _NUMERIC_FF_HPP_
#define _NUMERIC_FF_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <mpi.h>
#include <cmath>
#include <complex>
#include <cstring>

#include "object2hdf5.h"
#include "qgrid.hpp"

namespace hig {

	/**
	 * Some declarations.
	 */
	template<typename float_t, typename complex_t>
	extern unsigned int compute_form_factor_wrap(int rank,
							std::vector<float_t> &shape_def, unsigned int num_triangles,
							std::vector<short int> &axes,
							float_t* &qx_h, int nqx,
							float_t* &qy_h, int nqy,
							complex_t* &qz_h, int nqz,
							complex_t* &ff,
							float_t& kernel_time, float_t& red_time, float_t& mem_time
	#ifdef FINDBLOCK
							, const int block_x, const int block_y, const int block_z, const int block_t
	#endif
	#ifdef KERNEL2
							, unsigned int cuda_t, unsigned int cuda_y, unsigned int cuda_z
	#else // default kernel
							, unsigned int cuda_t
	#endif
							);
	
	template<typename complex_t>
	void write_slice_to_file(complex_t *ff, int nqx, int nqy, int nqz, char* filename, int axis, int slice);
	
	
	/**
	 * Templatized class for computing Form Factor across multiple nodes using MPI.
	 */
	template<typename float_t, typename complex_t>
	class NumericFormFactor {

		public:
		NumericFormFactor(int block_cuda): block_cuda_(block_cuda) { }
		NumericFormFactor(int block_cuda_t, int block_cuda_y, int block_cuda_z):
			block_cuda_t_(block_cuda_t), block_cuda_y_(block_cuda_y), block_cuda_z_(block_cuda_z) { }
		~NumericFormFactor() { }

		bool init() { return true; } 	// ...

		// old
		//bool compute(char* filename, complex_t* &ff, float_t* &qx_h, unsigned int nqx,
		//			float_t* &qy_h, unsigned int nqy, complex_t* qz_h, unsigned int nqz,
		//			MPI::Intracomm& world_comm);
		// new
		bool compute(const char* filename, complex_t* &ff, MPI::Intracomm& world_comm);
	
		private:
		unsigned int block_cuda_;
	
		/* cuda block size for main kernel */
		unsigned int block_cuda_t_;
		unsigned int block_cuda_y_;
		unsigned int block_cuda_z_;
	
		unsigned int read_shape_surface_file(const char* filename, std::vector<float_t>& shape_def);
		unsigned int read_shapes_hdf5(const char* filename, std::vector<float_t>& shape_def,
										MPI::Intracomm& comm);
		void find_axes_orientation(std::vector<float_t> &shape_def, std::vector<short int> &axes);
		void construct_ff(int rank, int num_procs, MPI::Comm &comm,
							MPI::Comm &col_comm, MPI::Comm &row_comm,
							int p_nqx, int p_nqy, int p_nqz,
							int nqx, int nqy, int nqz,
							int p_y, int p_z,
							complex_t* p_ff, complex_t* &ff,
							float_t&, float_t&);
		void gather_all(std::complex<float> *cast_p_ff, unsigned long int local_qpoints,
							std::complex<float> *cast_ff,
							int *recv_counts, int *displacements, MPI::Comm &comm);
		void gather_all(std::complex<double> *cast_p_ff, unsigned long int local_qpoints,
							std::complex<double> *cast_ff,
							int *recv_counts, int *displacements, MPI::Comm &comm);
	}; // class NumericFormFactor
	
	
	/**
	 * main host function
	 */
	template<typename float_t, typename complex_t>
	bool NumericFormFactor<float_t, complex_t>::compute(
			const char* filename, complex_t* &ff,
			/*float_t* &qx_h, unsigned int nqx,
			float_t* &qy_h, unsigned int nqy,
			complex_t* qz_h, unsigned int nqz,*/
			MPI::Intracomm& world_comm) {
		float_t comp_start = 0.0, comp_end = 0.0, comm_start = 0.0, comm_end = 0.0;
		float_t mem_start = 0.0, mem_end = 0.0;
		float_t comp_time = 0.0, comm_time = 0.0, mem_time = 0.0, kernel_time = 0.0, red_time = 0.0;
		float_t total_start = 0.0, total_end = 0.0, total_time = 0.0;

		unsigned int nqx = QGrid::instance().nqx();
		unsigned int nqy = QGrid::instance().nqy();
		unsigned int nqz = QGrid::instance().nqz_extended();
		
		comm_start = MPI::Wtime();
	
		world_comm.Barrier();	
		int rank = world_comm.Get_rank();
		int num_procs = world_comm.Get_size();
		
		comm_end = MPI::Wtime();
		comm_time += comm_end - comm_start;
	
		mem_start = MPI::Wtime();
	
		// all procs read the shape file
		// improve to parallel IO, or one proc reading and sending to all ...
		std::vector<float_t> shape_def;
		// use the new file reader instead ...
//		unsigned int num_triangles = read_shape_surface_file(filename, shape_def);
		unsigned int num_triangles = read_shapes_hdf5(filename, shape_def, world_comm);
						// ... <--- all procs read this? IMPROVE!!!
	
		std::vector<short int> axes(4);			// axes[i] = j
												// i: x=0 y=1 z=2
												// j: 0=a 1=b 2=c
		find_axes_orientation(shape_def, axes);
	
		if(rank == 0) {
			std::cout << std::endl
				<< "**        Using input shape file: " << filename << std::endl
				<< "**     Number of input triangles: " << num_triangles << std::endl
				<< "**  Q-grid resolution (q-points): " << nqx * nqy * nqz << std::endl
	            << "**               NQX x NQY x NQZ: " << nqx << " x " << nqy << " x " << nqz << std::endl
				<< "** Number of processes requested: " << num_procs << std::endl << std::flush;
		} // if
		if(num_triangles < 1) {
			//MPI::Finalize();
			std::cerr << "error: no triangles found in specified definition file" << std::endl;
			return false;
		} // if
	
		// decompose along y and z directions into blocks
		int p_y = floor(sqrt(num_procs));	// some procs may be idle ...
		int p_z = num_procs / p_y;
		
		int p_nqx = nqx;
		int p_nqy = nqy / p_y + (((rank / p_z) < nqy % p_y) ? 1 : 0);
		int p_nqz = nqz / p_z + (((rank % p_z) < nqz % p_z) ? 1 : 0);
	
		mem_end = MPI::Wtime();
		mem_time += mem_end - mem_start;
		
		// create row-wise and column-wise communicators
		comm_start = MPI::Wtime();
		MPI::Intracomm row_comm, col_comm;
		int row = rank / p_z, col = rank % p_z;
		if(row >= p_y) col = p_z;	// idle processes
		row_comm = world_comm.Split(row, rank);
		col_comm = world_comm.Split(col, rank);
	
		// construct a communicator for procs with rank < p_y * p_z
		// so that the idle procs do not participate at all
		MPI::Intracomm main_comm; int idle = 0;
		if(rank >= p_y * p_z) idle = 1;
		main_comm = world_comm.Split(idle, rank);
		comm_end = MPI::Wtime();
		comm_time += comm_end - comm_start;
	
	#ifdef FINDBLOCK
		int block_x = 0, block_y = 0, block_z = 0, block_t = 0;
		int block_x_max = 0, block_y_max = 0, block_z_max = 0, block_t_max = 0;
		block_x_max = (nqx < 100) ? nqx : 100;
		block_y_max = (nqy < 100) ? nqy : 100;
		block_z_max = (nqz < 100) ? nqz : 100;
		block_t_max = (num_triangles < 2000) ? num_triangles : 2000;
		block_t = block_t_max;
		for(block_x = block_x_max; block_x > 0; block_x -= 5) {
		for(block_y = 5; block_y < block_y_max; block_y += 5) {
		for(block_z = 5; block_z < block_z_max; block_z += 5) {
	#endif
		
		if(rank < p_y * p_z) {
			total_start = MPI::Wtime();
	
			int value = 1;
			int nidle_num_procs = 0;
	
			comm_start = MPI::Wtime();
			main_comm.Allreduce(&value, &nidle_num_procs, 1, MPI::INT, MPI::SUM);
			if(rank == 0) {
				std::cout << "++  Number of MPI processes used: " << nidle_num_procs << std::endl
						  << "++                 MPI grid size: 1 x " << p_y << " x " << p_z
						  << std::endl << std::flush;
			} // if
	
			// perform MPI scan operation to compute y_offset and z_offset
			unsigned int y_offset = 0, z_offset = 0;
			col_comm.Scan(&p_nqy, &y_offset, 1, MPI::INT, MPI::SUM);
			row_comm.Scan(&p_nqz, &z_offset, 1, MPI::INT, MPI::SUM);
			comm_end = MPI::Wtime();
			comm_time += comm_end - comm_start;
	
			mem_start = MPI::Wtime();
			y_offset -= p_nqy;
			z_offset -= p_nqz;

			// this is a temporary fix ... fix properly ...
			float_t* qx_h = new (std::nothrow) float_t[nqx]();
			float_t* qy_h = new (std::nothrow) float_t[nqy]();
			complex_t* qz_h = new (std::nothrow) complex_t[nqz]();
			// create qy_h and qz_h using qgrid instance
			for(int i = 0; i < nqx; ++ i) {
				qx_h[i] = QGrid::instance().qx(i);
			} // for
			for(int i = 0; i < nqy; ++ i) {
				qy_h[i] = QGrid::instance().qy(i);
			} // for
			for(int i = 0; i < nqz; ++ i) {
				qz_h[i].x = QGrid::instance().qz_extended(i).real();
				qz_h[i].y = QGrid::instance().qz_extended(i).imag();
			} // for
			
			// create p_ff buffers	<----- IMPROVE for all procs!!!
			float_t *p_qy_h = NULL;
			p_qy_h = new (std::nothrow) float_t[p_nqy]();
			if(p_qy_h == NULL) { return 0; }
			memcpy(p_qy_h, (void*) (qy_h + y_offset), p_nqy * sizeof(float_t));
			complex_t *p_qz_h = NULL;
			p_qz_h = new (std::nothrow) complex_t[p_nqz]();
			if(p_qz_h == NULL) { delete[] p_qy_h; return 0; }
			memcpy(p_qz_h, (void*) (qz_h + z_offset), p_nqz * sizeof(complex_t));
	
			mem_end = MPI::Wtime();
			mem_time += mem_end - mem_start;
		
			// compute local
			comp_start = MPI::Wtime();
			complex_t *p_ff = NULL;
			float_t temp_mem_time = 0.0, temp_comm_time = 0.0;
	
			unsigned int ret_nt = 0;
			ret_nt = compute_form_factor_wrap(rank, shape_def, num_triangles, axes,
						qx_h, p_nqx, p_qy_h, p_nqy,	p_qz_h, p_nqz, p_ff,
						kernel_time, red_time, temp_mem_time
	#ifdef FINDBLOCK
						, block_x, block_y, block_z, block_t
	#endif
	#ifdef KERNEL2
							, block_cuda_t_, block_cuda_y_, block_cuda_z_
	#else // default kernel
							, block_cuda_
	#endif
						);
			//p_ff = new (std::nothrow) complex_t[p_nqx * p_nqy * p_nqz];
			//ret_nt = 1;
			mem_time += (temp_mem_time / 1000);
			comp_end = MPI::Wtime();
			comp_time += comp_end - comp_start;
	
			if(ret_nt > 0) {
				temp_mem_time = 0.0; temp_comm_time = 0.0;
				construct_ff(rank, nidle_num_procs, main_comm, col_comm, row_comm,
						p_nqx, p_nqy, p_nqz, nqx, nqy, nqz, p_y, p_z, p_ff, ff,
						temp_mem_time, temp_comm_time);
				mem_time += temp_mem_time;
				comm_time += temp_comm_time;
			} // if
	
			/*if(rank == 0) {
				write_slice_to_file(ff, nqx, nqy, nqz, filename, 0, 0);	// x = 0, y = 1, z = 2
														// only slice along x implemented for now
			} // if*/
	
			comm_start = MPI::Wtime();
			main_comm.Barrier();
			comm_end = MPI::Wtime();
			comm_time += comm_end - comm_start;
	
			mem_start = MPI::Wtime();
	#ifdef FINDBLOCK
			if(ff != NULL) delete[] ff;
	#endif
			if(p_ff != NULL) delete[] p_ff;
			delete[] p_qz_h;
			delete[] p_qy_h;
			delete[] qz_h;
			delete[] qy_h;
			delete[] qx_h;
	
			total_end = MPI::Wtime();
			total_time = total_end - total_start;
			mem_end = MPI::Wtime();
			mem_time += mem_end - mem_start;
	
			if(rank == 0) {
				std::cout << "**                   Kernel time: " << kernel_time / 1000 << "s." << std::endl;
				std::cout << "**                Reduction time: " << red_time / 1000 << "s." << std::endl;
				std::cout << "**            Memory and IO time: " << mem_time << "s." << std::endl;
				std::cout << "**            Communication time: " << comm_time << "s." << std::endl;
				//std::cout << "**                    Total time: " << ((kernel_time + red_time) / 1000 +
				//												mem_time + comm_time) << "s." << std::endl;
				std::cout << "**                    Total time: " << total_time << "s." << std::endl;
				std::cout << std::flush;
				//std::cout << std::endl << std::flush;
			} // if
		} // if
	
		world_comm.Barrier();
	#ifdef FINDBLOCK
		} // block_z
		} // block_y
		} // block_x
	#endif
		//MPI::Finalize();	
		return true;
	} // NumericFormFactor::compute_form_factor()
	

	/**
	 * Function to gather partial FF arrays from all processes to construct the final FF.
	 * This is a bottleneck for large num procs ...
	 */
	template<typename float_t, typename complex_t>
	void NumericFormFactor<float_t, complex_t>::construct_ff(int rank, int num_procs,
			MPI::Comm &comm, MPI::Comm &col_comm, MPI::Comm &row_comm,
			int p_nqx, int p_nqy, int p_nqz,
			int nqx, int nqy, int nqz,
			int p_y, int p_z,
			complex_t* p_ff, complex_t* &ff,
			float_t& mem_time, float_t& comm_time) {
		float_t mem_start = 0, mem_end = 0, comm_start = 0, comm_end = 0;
		mem_time = 0; comm_time = 0;
	
		mem_start = MPI::Wtime();
	
		unsigned long int local_qpoints = p_nqx * p_nqy * p_nqz;
		unsigned long int total_qpoints = nqx * nqy * nqz;
	
		// process 0 creates the main ff, and collects computed p_ff from all others (just use gather)
		if(rank == 0) ff = new (std::nothrow) complex_t[total_qpoints];
	
		// construct stuff for gatherv
		int *recv_counts = new (std::nothrow) int[num_procs]();
		int *displacements = new (std::nothrow) int[num_procs]();
		mem_end = MPI::Wtime();
		mem_time += mem_end - mem_start;
		comm_start = MPI::Wtime();
		comm.Allgather(&local_qpoints, 1, MPI::INT, recv_counts, 1, MPI::INT);
		comm_end = MPI::Wtime();
		comm_time += comm_end - comm_start;
	
		mem_start = MPI::Wtime();
		displacements[0] = 0;
		for(int i = 1; i < num_procs; ++ i) {
			displacements[i] = displacements[i - 1] + recv_counts[i - 1];
		} // for
		complex_t *ff_buffer = new (std::nothrow) complex_t[total_qpoints];
		if(ff_buffer == NULL) {
			std::cerr << "Error allocating memory for ff buffer" << std::endl;
			return;
		} // if
		std::complex<float_t> *cast_p_ff, *cast_ff;
		cast_p_ff = reinterpret_cast<std::complex<float_t>*>(p_ff);
		cast_ff = reinterpret_cast<std::complex<float_t>*>(ff_buffer);
		mem_end = MPI::Wtime();
		mem_time += mem_end - mem_start;
	
		comm_start = MPI::Wtime();
		gather_all(cast_p_ff, local_qpoints,
					cast_ff, recv_counts, displacements,
					comm);
	
		int *recv_p_nqy = new (std::nothrow) int[p_y]();
		col_comm.Gather(&p_nqy, 1, MPI::INT, recv_p_nqy, 1, MPI::INT, 0);
		comm_end = MPI::Wtime();
		comm_time += comm_end - comm_start;
	
		mem_start = MPI::Wtime();
		int *off_p_nqy = new (std::nothrow) int[p_y]();
		off_p_nqy[0] = 0;
		for(int i = 1; i < p_y; ++ i) off_p_nqy[i] = off_p_nqy[i - 1] + recv_p_nqy[i - 1];
	
		// move all the data to correct places
		if(rank == 0) {
			unsigned long int ff_index = 0;
			for(int i_nqz = 0; i_nqz < nqz; ++ i_nqz) {
				for(int i_py = 0; i_py < p_y; ++ i_py) {
					unsigned long int ffb_index = nqx * (i_nqz * recv_p_nqy[i_py] + nqz * off_p_nqy[i_py]);
					memcpy(&ff[ff_index], &ff_buffer[ffb_index], nqx * recv_p_nqy[i_py] * sizeof(complex_t));
					ff_index += nqx * recv_p_nqy[i_py];
				} // for i_py
			} // for i_nqz
		} // if
	
		delete[] off_p_nqy;
		delete[] recv_p_nqy;
		delete[] ff_buffer;
		delete[] displacements;
		delete[] recv_counts;
		mem_end = MPI::Wtime();
		mem_time += mem_end - mem_start;
	} // NumericFormFactor::construct_ff()
	
	
	template<typename float_t, typename complex_t>
	void NumericFormFactor<float_t, complex_t>::gather_all(std::complex<float> *cast_p_ff,
			unsigned long int local_qpoints,
			std::complex<float> *cast_ff, int *recv_counts, int *displacements, MPI::Comm &comm) {
		comm.Gatherv(cast_p_ff, local_qpoints, MPI::COMPLEX, cast_ff, recv_counts,
					displacements, MPI::COMPLEX, 0);
	} // NumericFormFactor::gather_all()
	
	
	template<typename float_t, typename complex_t>
	void NumericFormFactor<float_t, complex_t>::gather_all(std::complex<double> *cast_p_ff,
			unsigned long int local_qpoints,
			std::complex<double> *cast_ff, int *recv_counts, int *displacements, MPI::Comm &comm) {
		comm.Gatherv(cast_p_ff, local_qpoints, MPI::DOUBLE_COMPLEX, cast_ff, recv_counts,
					displacements, MPI::DOUBLE_COMPLEX, 0);
	} // NumericFormFactor::gather_all()
	
	
	/**
	 * Function to read the input shape file.
	 */
	template<typename float_t, typename complex_t>
	unsigned int NumericFormFactor<float_t, complex_t>::read_shape_surface_file(const char* filename,
			std::vector<float_t> &shape_def) {
		std::ifstream f(filename);
		if(!f.is_open()) {
			std::cout << "Cannot open file " << filename << std::endl;
			return 1;
		} // if
		float_t s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;
	
		while(true) {
			f >> s;
			if(f.eof() || !f.good()) break;
			f >> nx; f >> ny; f >> nz;
			f >> cx; f >> cy; f >> cz;
			shape_def.push_back(s);
			shape_def.push_back(nx);
			shape_def.push_back(ny);
			shape_def.push_back(nz);
			shape_def.push_back(cx);
			shape_def.push_back(cy);
			shape_def.push_back(cz);
		} // while
	
		f.close();
		return shape_def.size() / 7;
	} // NumericFormFactor::read_shape_surface_file()
	
	
	template<typename float_t, typename complex_t>
	void NumericFormFactor<float_t, complex_t>::find_axes_orientation(
			std::vector<float_t> &shape_def, std::vector<short int> &axes) {
		float_t min_a = shape_def[4], max_a = shape_def[4];
		float_t min_b = shape_def[5], max_b = shape_def[5];
		float_t min_c = shape_def[6], max_c = shape_def[6];
	
		for(int i = 0; i + 6 < shape_def.size(); i += 7) {
			min_a = (min_a > shape_def[i + 4]) ? shape_def[i + 4] : min_a ;
			max_a = (max_a < shape_def[i + 4]) ? shape_def[i + 4] : max_a ;
			min_b = (min_b > shape_def[i + 5]) ? shape_def[i + 5] : min_b ;
			max_b = (max_b < shape_def[i + 5]) ? shape_def[i + 5] : max_b ;
			min_c = (min_c > shape_def[i + 6]) ? shape_def[i + 6] : min_c ;
			max_c = (max_c < shape_def[i + 6]) ? shape_def[i + 6] : max_c ;
		} // for
	
		float_t diff_a = max_a - min_a;
		float_t diff_b = max_b - min_b;
		float_t diff_c = max_c - min_c;
	
		// axes[i] = j
		// i: x=0 y=1 z=2
		// j: 0=a 1=b 2=c
	
		// the smallest one is x, other two are y and z
		if(diff_a < diff_b) {
			if(diff_a < diff_c) {
				// x is a
				axes[0] = 0; axes[1] = 1; axes[2] = 2;
			} else {
				// x is c
				axes[0] = 2; axes[1] = 0; axes[2] = 1;
			} // if-else
		} else {
			if(diff_b < diff_c) {
				// x is b
				axes[0] = 2; axes[1] = 0; axes[2] = 1;
			} else {
				// x is c
				axes[0] = 2; axes[1] = 0; axes[2] = 1;
			} // if-else
		} // if-else
	} // NumericFormFactor::find_axes_orientation()
	
	
	/**
	 * Function to read the shape definition input file in HDF5 format.
	 */
	template<typename float_t, typename complex_t>
	unsigned int NumericFormFactor<float_t, complex_t>::read_shapes_hdf5(const char* filename,
			std::vector<float_t> &shape_def, MPI::Intracomm& comm) {
		unsigned int num_triangles = 0;
		double* temp_shape_def = NULL;
	
		h5_shape_reader(filename, &temp_shape_def, &num_triangles/*, comm*/);
		for(int i = 0; i < num_triangles * 7; ++ i)	shape_def.push_back((float_t)temp_shape_def[i]);
	
		return num_triangles;
	} // NumericFormFactor::read_shapes_hdf5()
	
	
	/**
	 * Write a slice to output file
	 */
	/*template<typename complex_t>
	void write_slice_to_file(complex_t *ff, int nqx, int nqy, int nqz, char* filename, int axis, int slice) {
		if(ff == NULL) return;
	
		std::cout << "** Writing slice to file ...";
		switch(axis) {
			case 0:
				break;
			default:
				std::cout << "Given axis slice writing not implemented yet!" << std::endl;
		} // switch
	
		if(slice >= nqx || slice < 0) {
			std::cout << "Given slice does not exist!" << std::endl;
			return;
		} // if
	
		std::ofstream slice_file;
		//char* outfilename = "output_ff.dat";
		char outfilename[50];
		sprintf(outfilename, "ff_p%d.dat", MPI::COMM_WORLD.Get_size());
		std::cout << " " << outfilename << " ";
		slice_file.open(outfilename);
		if(!slice_file.is_open()) {
			std::cerr << "Error opening slice file for writing." << std::endl;
			return;
		} // if
	
		for(int y = 0; y < nqy; ++ y) {
			for(int z = 0; z < nqz; ++ z) {
				unsigned long int index = nqx * nqy * z + nqx * y + slice;
				slice_file << ff[index].x << "," << ff[index].y << "\t";
			} // for z
			slice_file << std::endl;
		} // for y
		slice_file.close();
		std::cout << " done." << std::endl;
	} // write_slice_to_file() */


} // namespace hig

#endif // _NUMERIC_FF_HPP_
