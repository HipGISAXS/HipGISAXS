/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_mic.cpp
 *  Created: Apr 02, 2013
 *  Modified: Fri 05 Apr 2013 10:41:57 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */
	
#include <iostream>
#include <cmath>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "woo/timer/woo_boostchronotimers.hpp"

#include "constants.hpp"
#include "parameters.hpp"
#include "parameters_mic.hpp"
#include "definitions_mic.hpp"
#include "mic_complex_numeric.hpp"

#include "ff_num_mic.hpp"
	
namespace hig {
	
	NumericFormFactorM::NumericFormFactorM() { }

	NumericFormFactorM::~NumericFormFactorM() { }

	// TODO: ...
	bool NumericFormFactorM::init() { return true; }

	/**
	 * The main host function called from outside, as part of the API for a single node.
	 */
	unsigned int NumericFormFactorM::compute_form_factor(int rank,
						float_vec_t &shape_def_vec, complex_t* &ff,
						float_t* qx, int nqx, float_t* qy, int nqy, complex_t* qz, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
						#ifdef FINDBLOCK
							, const int block_x, const int block_y, const int block_z, const int block_t
						#endif
						) {
		double kernel_time = 0.0, reduce_time = 0.0, total_kernel_time = 0.0, total_reduce_time = 0.0,
				temp_mem_time = 0.0, total_mem_time = 0.0;
		#ifdef _OPENMP
			if(rank == 0)
				std::cout << "++ Number of Host OpenMP threads: " << omp_get_max_threads() << std::endl;
		#endif

		int num_mic = 0;
		#ifdef __INTEL_OFFLOAD
			num_mic = _Offload_number_of_devices();
			std::cout << "++         Number of Target MICs: " << num_mic << std::endl;
			if(num_mic == 0) {
				std::cerr << "error: no Target MIC found!" << std::endl;
				return 0;
			} // if
		#else
			std::cerr << "warning: offloading to MIC not set! using host." << std::endl;
		#endif

		unsigned int num_triangles = shape_def_vec.size() / 7;
		if(num_triangles < 1) return 0;

		float_t* shape_def = &shape_def_vec[0];
	
		unsigned int total_qpoints = nqx * nqy * nqz;
	
		// allocate memory for the final FF 3D matrix
		ff = new (std::nothrow) complex_t[total_qpoints];	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(complex_t) << " b" << std::endl;
			return 0;
		} // if
		memset(ff, 0, total_qpoints * sizeof(complex_t));

		scomplex_t* qz_flat = new (std::nothrow) scomplex_t[nqz];
		if(qz_flat == NULL) {
			std::cerr << "Memory allocation failed for qz_flat. Size = "
						<< nqz * sizeof(scomplex_t) << " b" << std::endl;
			return 0;
		} // if
		for(int i = 0; i < nqz; ++ i) {
			//qz_flat[i].x = qz[i].real(); qz_flat[i].y = qz[i].imag();
			qz_flat[i] = make_sC(qz[i].real(), qz[i].imag());
		} // for

		unsigned int host_mem_usage = ((unsigned int) nqx + nqy) * sizeof(float_t) +	// qx, qy
										(unsigned int) nqz * sizeof(complex_t) +		// qz
										(unsigned int) nqz * sizeof(scomplex_t) +		// qz_flat
										(unsigned int) shape_def_vec.size() * sizeof(float_t) +	// shape_def
										total_qpoints * sizeof(complex_t);				// ff

		// allocate memory buffers on the target and transfer data
		#pragma offload_transfer target(mic:0) \
								in(qx: length(nqx) MIC_ALLOC) \
								in(qy: length(nqy) MIC_ALLOC) \
								in(qz_flat: length(nqz) MIC_ALLOC) \
								in(shape_def: length(7 * num_triangles) MIC_ALLOC)
	
		unsigned int target_mem_usage = ((unsigned int) nqx + nqy) * sizeof(float_t) +
										(unsigned int) nqz * sizeof(scomplex_t) +
										(unsigned int) num_triangles * 7 * sizeof(float_t);

		unsigned int matrix_size = (unsigned int) nqx * nqy * nqz * num_triangles;

		// get some information about the target and display ...
		// TODO ...
		
		// do hyperblocking to use less memory
		unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
		compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
							b_nqx, b_nqy, b_nqz, b_num_triangles
							#ifdef FINDBLOCK
								, block_x, block_y, block_z, block_t
							#endif
							);
		if(rank == 0) {
			std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
						<< " x " << b_nqz << " x " << b_num_triangles << std::endl;
		} // if
	
		unsigned int blocked_3d_matrix_size = (unsigned int) b_nqx * b_nqy * b_nqz;
		unsigned int blocked_matrix_size = (unsigned int) blocked_3d_matrix_size * b_num_triangles;

		// do target memory usage estimation ...
		// TODO ...
		
		size_t estimated_host_mem_need = host_mem_usage + blocked_matrix_size * sizeof(complex_t);
		if(rank == 0) {
			std::cout << "++    Estimated host memory need: " << (float) estimated_host_mem_need / 1024 / 1024
						<< " MB" << std::endl;
		} // if

		// allocate ff and fq buffers on host first
		scomplex_t *fq_buffer = new (std::nothrow) scomplex_t[blocked_matrix_size]();
		scomplex_t *ff_buffer = new (std::nothrow) scomplex_t[blocked_3d_matrix_size]();
		if(fq_buffer == NULL || ff_buffer == NULL) {
			std::cerr << "Memory allocation failed for f buffers. blocked_matrix_size = "
						<< blocked_matrix_size << std::endl
						<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
			delete[] ff;
			return 0;
		} // if

		host_mem_usage += (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(complex_t);

		if(rank == 0) {
			std::cout << "++              Host memory used: " << (float) host_mem_usage / 1024 / 1024
						<< " MB" << std::endl << std::flush;
		} // if

		// allocate ff and fq buffers on target
		#pragma offload_transfer target(mic:0) \
							nocopy(fq_buffer: length(blocked_matrix_size) MIC_ALLOC) \
							nocopy(ff_buffer: length(blocked_3d_matrix_size) MIC_ALLOC)

		// display actual target memory used ...
		// TODO ...
	
		// compute the number of sub-blocks, along each of the 4 dimensions
		// formulate loops over each dimension, to go over each sub block
		unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
		unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
		unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
		unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);
	
		unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
		unsigned int curr_b_num_triangles = b_num_triangles;
		unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;
	
		if(rank == 0) {
			std::cout << "++         Number of hyperblocks: " << num_blocks
						<< " [" << nb_x << " x " << nb_y << " x " << nb_z << " x " << nb_t << "]"
						<< std::endl;
		} // if
	
		unsigned int block_num = 0;

		if(rank == 0) std::cout << "-- Computing form factor on MIC ... " << std::flush;

		// compute for each hyperblock
		curr_b_nqx = b_nqx;
		for(unsigned int ib_x = 0; ib_x < nb_x; ++ ib_x) {
			if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
			curr_b_nqy = b_nqy;
			for(unsigned int ib_y = 0; ib_y < nb_y; ++ ib_y) {
				if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
				curr_b_nqz = b_nqz;
				for(unsigned int ib_z = 0; ib_z < nb_z; ++ ib_z) {
					if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
					curr_b_num_triangles = b_num_triangles;
					for(unsigned int ib_t = 0; ib_t < nb_t; ++ ib_t) {
						if(ib_t == nb_t - 1)
							curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

						#ifdef VERBOSITY_3
							if(rank == 0) {
								std::cout << "-- Processing hyperblock " << ++ block_num << " of "
										<< num_blocks << ": Kernel... " << std::flush;
							} // if
						#endif

						// call the main kernel offloaded to MIC
						form_factor_kernel(qx, qy, qz_flat, shape_def,
								curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
								b_nqx, b_nqy, b_nqz, b_num_triangles,
								ib_x, ib_y, ib_z, ib_t,
								fq_buffer);

						#ifdef VERBOSITY_3
							if(rank == 0) {
								std::cout << "done [" << temp_kernel_time << "s]. Reduction... "
										<< std::flush;
							} // if
						#endif

						// call the reduction kernel offloaded to MIC
						reduction_kernel(curr_b_nqx, curr_b_nqy, curr_b_nqz,
								curr_b_num_triangles, blocked_matrix_size,
								b_nqx, b_nqy, b_nqz, num_triangles,
								nqx, nqy, nqz,
								ib_x, ib_y, ib_z, ib_t,
								fq_buffer, ff_buffer);

						#ifdef VERBOSITY_3
							if(rank == 0) {
								std::cout << "done [" << temp_reduce_time << "s]." << std::endl;
							} // if
						#endif

						move_to_main_ff(ff_buffer, curr_b_nqx, curr_b_nqy, curr_b_nqz,
										b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
										ib_x, ib_y, ib_z, ff);
					} // for ib_t
				} // for ib_z
			} // for ib_y
		} // for ib_x

		// free f buffers on target
		#pragma offload_transfer target(mic:0) \
								nocopy(fq_buffer: length(0) MIC_FREE) \
								nocopy(ff_buffer: length(0) MIC_FREE)
		// and host
		delete[] ff_buffer;
		delete[] fq_buffer;

		#pragma offload_transfer target(mic:0) \
								nocopy(shape_def: length(0) MIC_FREE) \
								nocopy(qx: length(0) MIC_FREE) \
								nocopy(qy: length(0) MIC_FREE) \
								nocopy(qz_flat: length(0) MIC_FREE)
	
		delete[] qz_flat;

		if(rank == 0) std::cout << "done." << std::endl;

		pass_kernel_time = total_kernel_time;
		red_time = total_reduce_time;
		mem_time = total_mem_time;
	
		return num_triangles;
	} // NumericFormFactorM::compute_form_factor()
		

	/**
	 * the main Form Factor kernel function - for one hyperblock.
	 */
	void NumericFormFactorM::form_factor_kernel(float_t* qx, float_t* qy, scomplex_t* qz_flat,
					float_t* shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					scomplex_t* fq_buffer) {
		if(fq_buffer == NULL || qx == NULL || qy == NULL || qz_flat == NULL) return;
	
		#pragma offload target(mic:0) \
						in(shape_def: length(0) MIC_REUSE) \
						in(qx: length(0) MIC_REUSE) \
						in(qy: length(0) MIC_REUSE) \
						in(qz_flat: length(0) MIC_REUSE) \
						out(fq_buffer: length(0) MIC_REUSE)
		{
			omp_set_num_threads(MIC_OMP_NUM_THREADS_);
			#pragma omp parallel
			{
				#pragma omp for nowait //schedule(auto)
				for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
					unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
					float_t s = shape_def[shape_off];
					float_t nx = shape_def[shape_off + 1];
					float_t ny = shape_def[shape_off + 2];
					float_t nz = shape_def[shape_off + 3];
					float_t x = shape_def[shape_off + 4];
					float_t y = shape_def[shape_off + 5];
					float_t z = shape_def[shape_off + 6];
	
					unsigned int xy_size = curr_nqx * curr_nqy;
					unsigned int matrix_off = xy_size * curr_nqz * i_t;
					unsigned int start_z = b_nqz * ib_z;
					unsigned int start_y = b_nqy * ib_y;
					unsigned int start_x = b_nqx * ib_x;
	
					for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
							++ i_z, ++ global_i_z) {
						unsigned int off_start = matrix_off + xy_size * i_z;
						scomplex_t temp_z = qz_flat[global_i_z];
						scomplex_t qz2 = temp_z * temp_z;
						scomplex_t qzn = temp_z * nz;
						scomplex_t qzt = temp_z * z;
	
						for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
								++ i_y, ++ global_i_y) {
							unsigned int xy_off_start = curr_nqx * i_y;
							float_t temp_y = qy[global_i_y];
							float_t qy2 = temp_y * temp_y;
							float_t qyn = temp_y * ny;
							float_t qyt = temp_y * y;
	
							for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
									++ i_x, ++ global_i_x) {
								unsigned int off = off_start + xy_off_start + i_x;
								float_t temp_x = qx[global_i_x];
								scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
								scomplex_t qt = temp_x * x + qyt + qzt;
								scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

								fq_buffer[off] = compute_fq(s, qt, qn);
							} // for z
						} // for y
					} // for x
				} // for t
			} // pragma omp parallel
		} // pragma offload
	} // NumericFormFactorM::form_factor_kernel()
	
	
	/**
	 * Computational kernel function.
	 */

	// single precision
	__attribute__((target(mic:0)))
	float2_t NumericFormFactorM::compute_fq(float s, float2_t qt, float2_t qn) {
		float2_t v1 = qn * make_sC(cosf(qt.x), sinf(qt.x));
		float v2 = s * exp(qt.y);
		return v1 * v2;
	} // NumericFormFactorM::compute_fq()
	

	// double precision
	__attribute__((target(mic:0)))
	double2_t NumericFormFactorM::compute_fq(double s, double2_t qt, double2_t qn) {
		double2_t v1 = qn * make_sC(cos(qt.x), sin(qt.x));
		double v2 = s * exp(qt.y);
		return v1 * v2;
	} // NumericFormFactorM::compute_fq()
	
	
	/**
	 * Reduction kernel
	 */
	void NumericFormFactorM::reduction_kernel(
			unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
			unsigned int curr_b_num_triangles, unsigned int blocked_matrix_size,
			unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
			unsigned int num_triangles, unsigned int nqx, unsigned int nqy, unsigned int nqz,
			unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
			scomplex_t* fq_buffer, scomplex_t* ff_buffer) {
		if(fq_buffer == NULL || ff_buffer == NULL) return;
		
		#pragma offload target(mic:0) \
						in(fq_buffer: length(0) MIC_REUSE) \
						out(ff_buffer: length(curr_b_nqx * curr_b_nqy * curr_b_nqz) MIC_REUSE)
		{
			unsigned int curr_b_nqxyz = curr_b_nqx * curr_b_nqy * curr_b_nqz;
			unsigned int curr_b_nqxy = curr_b_nqx * curr_b_nqy;
			scomplex_t temp_complex = make_sC((float_t) 0.0, (float_t) -1.0);

			omp_set_num_threads(MIC_OMP_NUM_THREADS_);
	
			// reduction over all triangles
			for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
				#pragma omp parallel
				{
					#pragma omp for nowait //schedule(auto)
					for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
						for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
							scomplex_t total = make_sC((float_t) 0.0, (float_t) 0.0);
							for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
								unsigned int i_fq = curr_b_nqxyz * i_t + curr_b_nqxy * i_z +
													curr_b_nqx * i_y + i_x;
								total = total + fq_buffer[i_fq];
							} // for i_t
							unsigned int i_ff = curr_b_nqxy * i_z + curr_b_nqx * i_y + i_x;
							//ff_buffer[i_ff] = total * temp_complex;
							ff_buffer[i_ff] = make_sC(total.y, - total.x);
						} // for i_z
					} // for i_y
				} // pragma omp parallel
			} // for i_x
		} // pragma offload
	} // NumericFormFactorM::reduction_kernel()
	
	
	/**
	 * The main host function called from outside, as part of the API for a single node.
	 * Double buffered version.
	 */
	unsigned int NumericFormFactorM::compute_form_factor_db(int rank,
						float_vec_t &shape_def_vec, complex_t* &ff,
						float_t* qx, int nqx, float_t* qy, int nqy, complex_t* qz, int nqz,
						float_t& pass_kernel_time, float_t& red_time, float_t& mem_time
						#ifdef FINDBLOCK
							, const int block_x, const int block_y, const int block_z, const int block_t
						#endif
						) {
		double kernel_time = 0.0, reduce_time = 0.0, total_kernel_time = 0.0, total_reduce_time = 0.0,
				temp_mem_time = 0.0, total_mem_time = 0.0;
		woo::BoostChronoTimer kerneltimer;

		#ifdef _OPENMP
			if(rank == 0)
				std::cout << "++ Number of Host OpenMP threads: " << omp_get_max_threads() << std::endl;
		#endif

		int num_mic = 0;
		#ifdef __INTEL_OFFLOAD
			num_mic = _Offload_number_of_devices();
			std::cout << "++         Number of Target MICs: " << num_mic << std::endl;
			if(num_mic == 0) {
				std::cerr << "error: no Target MIC found!" << std::endl;
				return 0;
			} // if
		#else
			std::cerr << "warning: offloading to MIC not set! using host." << std::endl;
		#endif

		unsigned int num_triangles = shape_def_vec.size() / 7;
		if(num_triangles < 1) return 0;

		float_t* shape_def = &shape_def_vec[0];
	
		unsigned int total_qpoints = nqx * nqy * nqz;
	
		// allocate memory for the final FF 3D matrix
		ff = new (std::nothrow) complex_t[total_qpoints];	// allocate and initialize to 0
		if(ff == NULL) {
			std::cerr << "Memory allocation failed for ff. Size = "
						<< total_qpoints * sizeof(complex_t) << " b" << std::endl;
			return 0;
		} // if
		memset(ff, 0, total_qpoints * sizeof(complex_t));
		unsigned long int host_mem_usage = total_qpoints * sizeof(complex_t);

		//scomplex_t* qz_flat = new (std::nothrow) scomplex_t[nqz];
		scomplex_t* qz_flat = (scomplex_t*) _mm_malloc(nqz * sizeof(scomplex_t), 64);
		if(qz_flat == NULL) {
			std::cerr << "Memory allocation failed for qz_flat. Size = "
						<< nqz * sizeof(scomplex_t) << " b" << std::endl;
			return 0;
		} // if
		for(int i = 0; i < nqz; ++ i) {
			//qz_flat[i].x = qz[i].real(); qz_flat[i].y = qz[i].imag();
			qz_flat[i] = make_sC(qz[i].real(), qz[i].imag());
		} // for

		host_mem_usage += ((unsigned int) nqx + nqy) * sizeof(float_t) +	// qx, qy
							(unsigned int) nqz * sizeof(complex_t) +        // qz
							(unsigned int) nqz * sizeof(scomplex_t) +       // qz_flat
							(unsigned int) shape_def_vec.size() * sizeof(float_t) + // shape_def
							total_qpoints * sizeof(complex_t);              // ff

		// allocate memory buffers on the target and transfer data asynchronously
		#pragma offload_transfer target(mic:0) \
								in(qx: length(nqx) MIC_ALLOC) \
								in(qy: length(nqy) MIC_ALLOC) \
								in(qz_flat: length(nqz) MIC_ALLOC) \
								in(shape_def: length(7 * num_triangles) MIC_ALLOC) \
								signal(shape_def)
		
		unsigned long int target_mem_usage = ((unsigned int) nqx + nqy) * sizeof(float_t) +
											(unsigned int) nqz * sizeof(scomplex_t) +
											(unsigned int) num_triangles * 7 * sizeof(float_t);

		unsigned int matrix_size = (unsigned int) nqx * nqy * nqz * num_triangles;

		// do hyperblocking to use less memory and computation - mem transfer overlap
		unsigned int b_nqx = 0, b_nqy = 0, b_nqz = 0, b_num_triangles = 0;
		compute_hyperblock_size(nqx, nqy, nqz, num_triangles,
								b_nqx, b_nqy, b_nqz, b_num_triangles
								#ifdef FINDBLOCK
									, block_x, block_y, block_z, block_t
								#endif
								);
		// compute the number of sub-blocks, along each of the 4 dimensions
		// formulate loops over each dimension, to go over each sub block
		unsigned int nb_x = (unsigned int) ceil((float) nqx / b_nqx);
		unsigned int nb_y = (unsigned int) ceil((float) nqy / b_nqy);
		unsigned int nb_z = (unsigned int) ceil((float) nqz / b_nqz);
		unsigned int nb_t = (unsigned int) ceil((float) num_triangles / b_num_triangles);
	
		unsigned int curr_b_nqx = b_nqx, curr_b_nqy = b_nqy, curr_b_nqz = b_nqz;
		unsigned int curr_b_num_triangles = b_num_triangles;
		unsigned int num_blocks = nb_x * nb_y * nb_z * nb_t;
		// for the progress indicator
		unsigned int progress_delta = 10; // display progress at 10% intervals
		unsigned int progress_step = (unsigned int) ceil((float) num_blocks / progress_delta);
	
		if(rank == 0) {
			std::cout << "++               Hyperblock size: " << b_nqx << " x " << b_nqy
						<< " x " << b_nqz << " x " << b_num_triangles << std::endl;
			std::cout << "++         Number of hyperblocks: " << num_blocks
						<< " [" << nb_x << " x " << nb_y << " x " << nb_z << " x " << nb_t << "]"
						<< std::endl;
		} // if
	
		unsigned int blocked_3d_matrix_size = (unsigned int) b_nqx * b_nqy * b_nqz;
		unsigned long int blocked_matrix_size = (unsigned int) blocked_3d_matrix_size * b_num_triangles;

		// allocate ff and fq buffers on host first (double buffering)
		scomplex_t *fq_buffer0, *fq_buffer1, *ff_buffer0, *ff_buffer1;
		//fq_buffer0 = new (std::nothrow) scomplex_t[blocked_matrix_size]();
		//fq_buffer1 = new (std::nothrow) scomplex_t[blocked_matrix_size]();
		//ff_buffer0 = new (std::nothrow) scomplex_t[blocked_3d_matrix_size]();
		//ff_buffer1 = new (std::nothrow) scomplex_t[blocked_3d_matrix_size]();
		#ifdef MIC_PADDING
			// to make buffers multiple of 64B:
			unsigned int pad_fq = (unsigned int) ceil(blocked_matrix_size * sizeof(scomplex_t) / 64.0f) *
													64 - blocked_matrix_size * sizeof(scomplex_t);
			unsigned int pad_ff = (unsigned int) ceil(blocked_3d_matrix_size * sizeof(scomplex_t) / 64.0f) *
													64 - blocked_3d_matrix_size * sizeof(scomplex_t);
			unsigned int fq_bytes = blocked_3d_matrix_size * sizeof(scomplex_t) + pad_fq;
			unsigned int ff_bytes = blocked_matrix_size * sizeof(scomplex_t) + pad_ff;
			unsigned int fq_size = fq_bytes / sizeof(scomplex_t);
			unsigned int ff_size = ff_bytes / sizeof(scomplex_t);
			if(rank == 0) {
				std::cout << "++                Buffer padding: " << pad_ff << " B" << std::endl;
			} // if
			fq_buffer0 = (scomplex_t*) _mm_malloc(fq_bytes, 64);
			fq_buffer1 = (scomplex_t*) _mm_malloc(fq_bytes, 64);
			ff_buffer0 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
			ff_buffer1 = (scomplex_t*) _mm_malloc(ff_bytes, 64);
			host_mem_usage += 2 * (fq_bytes + ff_bytes);
			target_mem_usage += 2 * (fq_bytes + ff_bytes);
		#else
			fq_buffer0 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_matrix_size, 64);
			fq_buffer1 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_matrix_size, 64);
			ff_buffer0 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
			ff_buffer1 = (scomplex_t*) _mm_malloc(sizeof(scomplex_t) * blocked_3d_matrix_size, 2147483648);
			host_mem_usage += 2 * (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(scomplex_t);
			target_mem_usage += 2 * (blocked_matrix_size + blocked_3d_matrix_size) * sizeof(scomplex_t);
		#endif
		if(fq_buffer0 == NULL || fq_buffer1 == NULL || ff_buffer0 == NULL || ff_buffer1 == NULL) {
			std::cerr << "Memory allocation failed for f buffers. blocked_matrix_size = "
						<< blocked_matrix_size << std::endl
						<< "Host memory usage = " << (float) host_mem_usage / 1024 / 1024 << " MB"
						<< std::endl;
			return 0;
		} // if

		if(rank == 0) {
			std::cout << "++             Host memory usage: " << (float) host_mem_usage / 1024 / 1024
						<< " MB" << std::endl << std::flush;
			#ifdef __INTEL_OFFLOAD
				std::cout << "++           Target memory usage: " << (float) target_mem_usage / 1024 / 1024
							<< " MB" << std::endl << std::flush;
			#endif
		} // if

		// allocate ff and fq buffers on target
		#ifdef MIC_PADDING
			#pragma offload_transfer target(mic:0) \
							nocopy(fq_buffer0: length(fq_size) MIC_ALLOC) \
							nocopy(fq_buffer1: length(fq_size) MIC_ALLOC) \
							nocopy(ff_buffer0: length(ff_size) MIC_ALLOC) \
							nocopy(ff_buffer1: length(ff_size) MIC_ALLOC)
		#else
			#pragma offload_transfer target(mic:0) \
							nocopy(fq_buffer0: length(blocked_matrix_size) MIC_ALLOC) \
							nocopy(fq_buffer1: length(blocked_matrix_size) MIC_ALLOC) \
							nocopy(ff_buffer0: length(blocked_3d_matrix_size) MIC_ALLOC) \
							nocopy(ff_buffer1: length(blocked_3d_matrix_size) MIC_ALLOC)
		#endif

		if(rank == 0) {
			#ifdef __INTEL_OFFLOAD
				std::cout << "-- Computing form factor on MIC (DB) ... " << std::flush;
			#else
				std::cout << "-- Computing form factor on CPU (fallback) ... " << std::flush;
			#endif
		} // if

		int curr_buffer_i = 0;	// buffer index
		unsigned int prev_b_nqx, prev_b_nqy, prev_b_nqz, prev_b_num_triangles;
		int prev_ib_x, prev_ib_y, prev_ib_z, prev_ib_t;

		unsigned int hblock_counter = 0;

		// wait for input transfer to complete
		#pragma offload_wait target(mic:0) wait(shape_def)

		kerneltimer.start();

		// compute for each hyperblock
		curr_b_nqx = b_nqx;
		for(int ib_x = 0; ib_x < nb_x; ++ ib_x) {
			if(ib_x == nb_x - 1) curr_b_nqx = nqx - b_nqx * ib_x;
			curr_b_nqy = b_nqy;
			for(int ib_y = 0; ib_y < nb_y; ++ ib_y) {
				if(ib_y == nb_y - 1) curr_b_nqy = nqy - b_nqy * ib_y;
				curr_b_nqz = b_nqz;
				for(int ib_z = 0; ib_z < nb_z; ++ ib_z) {
					if(ib_z == nb_z - 1) curr_b_nqz = nqz - b_nqz * ib_z;
					curr_b_num_triangles = b_num_triangles;
					for(int ib_t = 0; ib_t < nb_t; ++ ib_t) {
						if(ib_t == nb_t - 1)
							curr_b_num_triangles = num_triangles - b_num_triangles * ib_t;

						// double buffering
						if(curr_buffer_i == 0) {
							// call the main kernel offloaded to MIC
							#ifndef FF_MIC_OPT
								#pragma offload target(mic:0) \
										in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
										in(b_nqx, b_nqy, b_nqz, num_triangles) \
										in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
										in(shape_def: length(0) MIC_REUSE) \
										in(qx: length(0) MIC_REUSE) \
										in(qy: length(0) MIC_REUSE) \
										in(qz_flat: length(0) MIC_REUSE) \
										out(fq_buffer0: length(0) MIC_REUSE)
								form_factor_kernel_db(qx, qy, qz_flat, shape_def,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										fq_buffer0);

								// call the reduction kernel offloaded to MIC
								#ifdef MIC_PADDING
									#pragma offload target(mic:0) signal(ff_buffer0) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(fq_buffer0: length(0) MIC_REUSE) \
											out(ff_buffer0: length(ff_size) MIC_REUSE)
								#else
									#pragma offload target(mic:0) signal(ff_buffer0) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(fq_buffer0: length(0) MIC_REUSE) \
											out(ff_buffer0: \
													length(curr_b_nqx * curr_b_nqy * curr_b_nqz) MIC_REUSE)
								#endif
								reduction_kernel_db(curr_b_nqx, curr_b_nqy, curr_b_nqz,
										curr_b_num_triangles, blocked_matrix_size,
										b_nqx, b_nqy, b_nqz, num_triangles,
										nqx, nqy, nqz,
										ib_x, ib_y, ib_z, ib_t,
										fq_buffer0, ff_buffer0);
							#else
								#ifdef MIC_PADDING
									#pragma offload target(mic:0) signal(ff_buffer0) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(shape_def: length(0) MIC_REUSE) \
											in(qx: length(0) MIC_REUSE) \
											in(qy: length(0) MIC_REUSE) \
											in(qz_flat: length(0) MIC_REUSE) \
											out(fq_buffer0: length(0) MIC_REUSE) \
											out(ff_buffer0: length(ff_size) MIC_REUSE)
								#else
									#pragma offload target(mic:0) signal(ff_buffer0) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(shape_def: length(0) MIC_REUSE) \
											in(qx: length(0) MIC_REUSE) \
											in(qy: length(0) MIC_REUSE) \
											in(qz_flat: length(0) MIC_REUSE) \
											out(fq_buffer0: length(0) MIC_REUSE) \
											out(ff_buffer0: \
													length(curr_b_nqx * curr_b_nqy * curr_b_nqz) MIC_REUSE)
								#endif
								form_factor_kernel_opt(qx, qy, qz_flat, shape_def,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										fq_buffer0, ff_buffer0);
							#endif

							if(ib_x + ib_y + ib_z + ib_t != 0) {
								// wait for transfer of 1 to finish before moving
								#pragma offload_wait target(mic:0) wait(ff_buffer1)

								// move computed ff block from buffer to final ff
								move_to_main_ff(ff_buffer1, prev_b_nqx, prev_b_nqy, prev_b_nqz,
												b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
												prev_ib_x, prev_ib_y, prev_ib_z, ff);
							} // if
						} else {
							// call the main kernel offloaded to MIC
							#ifndef FF_MIC_OPT
								#pragma offload target(mic:0) \
										in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
										in(b_nqx, b_nqy, b_nqz, num_triangles) \
										in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
										in(shape_def: length(0) MIC_REUSE) \
										in(qx: length(0) MIC_REUSE) \
										in(qy: length(0) MIC_REUSE) \
										in(qz_flat: length(0) MIC_REUSE) \
										out(fq_buffer1: length(0) MIC_REUSE)
								form_factor_kernel_db(qx, qy, qz_flat, shape_def,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										fq_buffer1);

								// call the reduction kernel offloaded to MIC
								#ifdef MIC_PADDING
									#pragma offload target(mic:0) signal(ff_buffer1) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(fq_buffer1: length(0) MIC_REUSE) \
											out(ff_buffer1: length(ff_size) MIC_REUSE)
								#else
									#pragma offload target(mic:0) signal(ff_buffer1) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(fq_buffer1: length(0) MIC_REUSE) \
											out(ff_buffer1: \
													length(curr_b_nqx * curr_b_nqy * curr_b_nqz) MIC_REUSE)
								#endif
								reduction_kernel_db(curr_b_nqx, curr_b_nqy, curr_b_nqz,
										curr_b_num_triangles, blocked_matrix_size,
										b_nqx, b_nqy, b_nqz, num_triangles,
										nqx, nqy, nqz,
										ib_x, ib_y, ib_z, ib_t,
										fq_buffer1, ff_buffer1);
							#else
								#ifdef MIC_PADDING
									#pragma offload target(mic:0) signal(ff_buffer1) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(shape_def: length(0) MIC_REUSE) \
											in(qx: length(0) MIC_REUSE) \
											in(qy: length(0) MIC_REUSE) \
											in(qz_flat: length(0) MIC_REUSE) \
											out(fq_buffer1: length(0) MIC_REUSE) \
											out(ff_buffer1: length(ff_size) MIC_REUSE)
								#else
									#pragma offload target(mic:0) signal(ff_buffer1) \
											in(curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles) \
											in(b_nqx, b_nqy, b_nqz, num_triangles) \
											in(nqx, nqy, nqz, ib_x, ib_y, ib_z, ib_t) \
											in(shape_def: length(0) MIC_REUSE) \
											in(qx: length(0) MIC_REUSE) \
											in(qy: length(0) MIC_REUSE) \
											in(qz_flat: length(0) MIC_REUSE) \
											out(fq_buffer1: length(0) MIC_REUSE) \
											out(ff_buffer1: \
													length(curr_b_nqx * curr_b_nqy * curr_b_nqz) MIC_REUSE)
								#endif
								form_factor_kernel_opt(qx, qy, qz_flat, shape_def,
										curr_b_nqx, curr_b_nqy, curr_b_nqz, curr_b_num_triangles,
										b_nqx, b_nqy, b_nqz, b_num_triangles,
										ib_x, ib_y, ib_z, ib_t,
										fq_buffer1, ff_buffer1);
							#endif

							if(ib_x + ib_y + ib_z + ib_t != 0) {
								// wait for transfer of 0 to finish before moving
								#pragma offload_wait target(mic:0) wait(ff_buffer0)

								// move computed ff block from buffer to final ff
								move_to_main_ff(ff_buffer0, prev_b_nqx, prev_b_nqy, prev_b_nqz,
												b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
												prev_ib_x, prev_ib_y, prev_ib_z, ff);
							} // if
						} // if-else

						// flip
						curr_buffer_i = 1 - curr_buffer_i;
						prev_b_nqx = curr_b_nqx; prev_b_nqy = curr_b_nqy; prev_b_nqz = curr_b_nqz;
						prev_b_num_triangles = curr_b_num_triangles;
						prev_ib_x = ib_x; prev_ib_y = ib_y; prev_ib_z = ib_z; prev_ib_t = ib_t;

						++ hblock_counter;

						if(rank == 0) {
							float progress = ((float) hblock_counter / num_blocks) * 100;
							double intpart, fracpart;
							fracpart = modf(progress, &intpart);
							if(((unsigned int)intpart) % progress_delta == 0 && fracpart < 100.0 / num_blocks)
								std::cout << intpart << "% ";
						} // if
					} // for ib_t
				} // for ib_z
			} // for ib_y
		} // for ib_x

		// move the last part
		if(curr_buffer_i == 0) {
			// wait for transfer of 1 to finish before moving
			#pragma offload_wait target(mic:0) wait(ff_buffer1)

			// move computed ff block from buffer to final ff
			move_to_main_ff(ff_buffer1, prev_b_nqx, prev_b_nqy, prev_b_nqz,
							b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
							prev_ib_x, prev_ib_y, prev_ib_z, ff);
		} else {
			// wait for transfer of 0 to finish before moving
			#pragma offload_wait target(mic:0) wait(ff_buffer0)

			// move computed ff block from buffer to final ff
			move_to_main_ff(ff_buffer0, prev_b_nqx, prev_b_nqy, prev_b_nqz,
							b_nqx, b_nqy, b_nqz, nqx, nqy, nqz,
							prev_ib_x, prev_ib_y, prev_ib_z, ff);
		} // if-else

		kerneltimer.stop();

		// free f buffers on target
		#pragma offload_transfer target(mic:0) \
				nocopy(fq_buffer0: length(0) MIC_FREE) \
				nocopy(fq_buffer1: length(0) MIC_FREE) \
				nocopy(ff_buffer0: length(0) MIC_FREE) \
				nocopy(ff_buffer1: length(0) MIC_FREE)
		// and host
		//delete[] ff_buffer0;
		//delete[] ff_buffer1;
		//delete[] fq_buffer0;
		//delete[] fq_buffer1;
		_mm_free(ff_buffer0);
		_mm_free(ff_buffer1);
		_mm_free(fq_buffer0);
		_mm_free(fq_buffer1);

		#pragma offload_transfer target(mic:0) \
				nocopy(shape_def: length(0) MIC_FREE) \
				nocopy(qx: length(0) MIC_FREE) \
				nocopy(qy: length(0) MIC_FREE) \
				nocopy(qz_flat: length(0) MIC_FREE)
	
		//delete[] qz_flat;
		_mm_free(qz_flat);

		if(rank == 0) std::cout << "done." << std::endl;

		pass_kernel_time = kerneltimer.elapsed_msec();
		red_time = total_reduce_time;
		mem_time = total_mem_time;
	
		return num_triangles;
	} // NumericFormFactorM::compute_form_factor_db()
		

	/***
	 * the following kernels are used in double buffering case
	 * TODO ... clean this
	 */

	__attribute__((target(mic:0)))
	void NumericFormFactorM::form_factor_kernel_db(float_t* qx, float_t* qy, scomplex_t* qz_flat,
					float_t* shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					scomplex_t* fq_buffer) {

		#ifdef __INTEL_OFFLOAD
			omp_set_num_threads(MIC_OMP_NUM_THREADS_);
		#endif
		#pragma omp parallel
		{
			#pragma omp for nowait //schedule(auto)
			for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
				unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
				float_t s = shape_def[shape_off];
				float_t nx = shape_def[shape_off + 1];
				float_t ny = shape_def[shape_off + 2];
				float_t nz = shape_def[shape_off + 3];
				float_t x = shape_def[shape_off + 4];
				float_t y = shape_def[shape_off + 5];
				float_t z = shape_def[shape_off + 6];

				unsigned int xy_size = curr_nqx * curr_nqy;
				unsigned int matrix_off = xy_size * curr_nqz * i_t;
				unsigned int start_z = b_nqz * ib_z;
				unsigned int start_y = b_nqy * ib_y;
				unsigned int start_x = b_nqx * ib_x;

				for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
						++ i_z, ++ global_i_z) {
					unsigned int off_start = matrix_off + xy_size * i_z;
					scomplex_t temp_z = qz_flat[global_i_z];
					scomplex_t qz2 = temp_z * temp_z;
					scomplex_t qzn = temp_z * nz;
					scomplex_t qzt = temp_z * z;

					for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
							++ i_y, ++ global_i_y) {
						unsigned int xy_off_start = curr_nqx * i_y;
						float_t temp_y = qy[global_i_y];
						float_t qy2 = temp_y * temp_y;
						float_t qyn = temp_y * ny;
						float_t qyt = temp_y * y;

						for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
								++ i_x, ++ global_i_x) {
							unsigned int off = off_start + xy_off_start + i_x;
							float_t temp_x = qx[global_i_x];
							scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
							scomplex_t qt = temp_x * x + qyt + qzt;
							scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

							fq_buffer[off] = compute_fq(s, qt, qn);
						} // for z
					} // for y
				} // for x
			} // for t
		} // pragma omp parallel
	} // NumericFormFactorM::form_factor_kernel_db()
	
	
	__attribute__((target(mic:0)))
	void NumericFormFactorM::reduction_kernel_db(
				unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
				unsigned int curr_b_num_triangles, unsigned int blocked_matrix_size,
				unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
				unsigned int num_triangles, unsigned int nqx, unsigned int nqy, unsigned int nqz,
				unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
				scomplex_t* fq_buffer, scomplex_t* ff_buffer) {

		unsigned int curr_b_nqxyz = curr_b_nqx * curr_b_nqy * curr_b_nqz;
		unsigned int curr_b_nqxy = curr_b_nqx * curr_b_nqy;
		scomplex_t temp_complex = make_sC((float_t) 0.0, (float_t) -1.0);

		#ifdef __INTEL_OFFLOAD
			omp_set_num_threads(MIC_OMP_NUM_THREADS_);
		#endif

		// reduction over all triangles
		for(unsigned int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
			#pragma omp parallel
			{
				#pragma omp for nowait //schedule(auto)
				for(unsigned int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
					for(unsigned int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
						scomplex_t total = make_sC((float_t) 0.0, (float_t) 0.0);
						for(unsigned int i_t = 0; i_t < curr_b_num_triangles; ++ i_t) {
							unsigned int i_fq = curr_b_nqxyz * i_t + curr_b_nqxy * i_z +
												curr_b_nqx * i_y + i_x;
							total = total + fq_buffer[i_fq];
						} // for i_t
						unsigned int i_ff = curr_b_nqxy * i_z + curr_b_nqx * i_y + i_x;
						//ff_buffer[i_ff] = total * temp_complex;
						ff_buffer[i_ff] = make_sC(total.y, - total.x);
					} // for i_z
				} // for i_y
			} // pragma omp parallel
		} // for i_x
	} // NumericFormFactorM::reduction_kernel_db()
	
	
	/***
	 * the following kernels are with optimizations (trials)
	 * TODO ... clean this
	 */

	__attribute__((target(mic:0)))
	void NumericFormFactorM::form_factor_kernel_opt(float_t* qx, float_t* qy, scomplex_t* qz_flat,
					float_t* shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					scomplex_t* fq_buffer, scomplex_t* ff_buffer) {

		#ifdef __INTEL_OFFLOAD
			omp_set_num_threads(MIC_OMP_NUM_THREADS_);
		#endif

		unsigned int curr_nqxy = curr_nqx * curr_nqy;
		#pragma omp parallel
		{
			#pragma omp for nowait //schedule(auto)
			for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
				unsigned int shape_off = (ib_t * b_num_triangles + i_t) * 7;
				float_t s = shape_def[shape_off];
				float_t nx = shape_def[shape_off + 1];
				float_t ny = shape_def[shape_off + 2];
				float_t nz = shape_def[shape_off + 3];
				float_t x = shape_def[shape_off + 4];
				float_t y = shape_def[shape_off + 5];
				float_t z = shape_def[shape_off + 6];

				unsigned int xy_size = curr_nqxy;
				unsigned int matrix_off = xy_size * curr_nqz * i_t;
				unsigned int start_z = b_nqz * ib_z;
				unsigned int start_y = b_nqy * ib_y;
				unsigned int start_x = b_nqx * ib_x;

				for(unsigned int i_z = 0, global_i_z = start_z; i_z < curr_nqz;
						++ i_z, ++ global_i_z) {
					unsigned int off_start = matrix_off + xy_size * i_z;
					scomplex_t temp_z = qz_flat[global_i_z];
					scomplex_t qz2 = temp_z * temp_z;
					scomplex_t qzn = temp_z * nz;
					scomplex_t qzt = temp_z * z;

					for(unsigned int i_y = 0, global_i_y = start_y; i_y < curr_nqy;
							++ i_y, ++ global_i_y) {
						unsigned int xy_off_start = curr_nqx * i_y;
						float_t temp_y = qy[global_i_y];
						float_t qy2 = temp_y * temp_y;
						float_t qyn = temp_y * ny;
						float_t qyt = temp_y * y;

						for(unsigned int i_x = 0, global_i_x = start_x; i_x < curr_nqx;
								++ i_x, ++ global_i_x) {
							unsigned int off = off_start + xy_off_start + i_x;
							float_t temp_x = qx[global_i_x];
							scomplex_t q2 = temp_x * temp_x + qy2 + qz2;
							scomplex_t qt = temp_x * x + qyt + qzt;
							scomplex_t qn = (temp_x * nx + qyn + qzn) / q2;

							fq_buffer[off] = compute_fq(s, qt, qn);
							/*unsigned int i_ff = curr_nqxy * i_z + curr_nqx * i_y + i_x;
							scomplex_t res = compute_fq(s, qt, qn);
							res = make_sC(res.y, - res.x);
							#pragma omp atomic
							ff_buffer[i_ff].x += res.x;
							#pragma omp atomic
							ff_buffer[i_ff].y += res.y;*/
						} // for z
					} // for y
				} // for x
			} // for t
		} // pragma omp parallel
		
		// doing reduction here itself
		// a separate reduction helps avoid use of atomics otherwise (which are SLOW!)
		unsigned int curr_nqxyz = curr_nqx * curr_nqy * curr_nqz;
		//scomplex_t temp_complex = make_sC((float_t) 0.0, (float_t) -1.0);
		#pragma omp parallel
		{
			#pragma omp for nowait collapse(3) //schedule(auto)
			for(unsigned int i_x = 0; i_x < curr_nqx; ++ i_x) {
				for(unsigned int i_y = 0; i_y < curr_nqy; ++ i_y) {
					for(unsigned int i_z = 0; i_z < curr_nqz; ++ i_z) {
						scomplex_t total = make_sC((float_t) 0.0, (float_t) 0.0);
						for(unsigned int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
							unsigned int i_fq = curr_nqxyz * i_t + curr_nqxy * i_z +
												curr_nqx * i_y + i_x;
							total = total + fq_buffer[i_fq];
						} // for i_t
						unsigned int i_ff = curr_nqxy * i_z + curr_nqx * i_y + i_x;
						//ff_buffer[i_ff] = total * temp_complex;
						ff_buffer[i_ff] = make_sC(total.y, - total.x);
					} // for i_z
				} // for i_y
			} // for i_x
		} // pragma omp parallel
	} // NumericFormFactorM::form_factor_kernel_opt()
	
	
	/**
	 * Function to move a computed block of ff to its right place in ff.
	 */
	void NumericFormFactorM::move_to_main_ff(scomplex_t* ff_buffer,
				unsigned int curr_b_nqx, unsigned int curr_b_nqy, unsigned int curr_b_nqz,
				unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
				unsigned int nqx, unsigned int nqy, unsigned int nqz,
				unsigned int ib_x, unsigned int ib_y, unsigned int ib_z,
				complex_t* ff) {
		unsigned int temp1 = nqx * nqy;
		unsigned int temp2 = curr_b_nqx * curr_b_nqy;
		unsigned int base_i = nqx * nqy * ib_z * b_nqz + nqx * ib_y * b_nqy + ib_x * b_nqx;
		// make these copy operations contiguous ... (very low priority)
		for(int i_z = 0; i_z < curr_b_nqz; ++ i_z) {
			unsigned int start_i = base_i + temp1 * i_z;
			for(int i_y = 0; i_y < curr_b_nqy; ++ i_y) {
				unsigned int super_i = start_i + nqx * i_y;
				for(int i_x = 0; i_x < curr_b_nqx; ++ i_x) {
					unsigned int final_i = super_i + i_x;
					unsigned int block_i = temp2 * i_z + curr_b_nqx * i_y + i_x;
					ff[final_i] += complex_t(ff_buffer[block_i].x, ff_buffer[block_i].y);
				} // for i_x
			} // for i_y
		} // for i_z
	} // move_to_main_ff()
	
	
	/**
	 * Function to compute the decomposition block size
	 * TODO: Improve it later ...
	 */
	void NumericFormFactorM::compute_hyperblock_size(int nqx, int nqy, int nqz, int num_triangles,
			unsigned int& b_nqx, unsigned int& b_nqy, unsigned int& b_nqz, unsigned int& b_num_triangles
			#ifdef FINDBLOCK
				, const int block_x, const int block_y, const int block_z, const int block_t
			#endif
			) {
		b_nqx = (unsigned int) nqx; b_nqy = (unsigned int) nqy; b_nqz = (unsigned int) nqz;
		b_num_triangles = (unsigned int) num_triangles;
	
		#ifdef FINDBLOCK
			b_nqx = (b_nqx > block_x) ? block_x : b_nqx;
			b_nqy = (b_nqy > block_y) ? block_y : b_nqy;
			b_nqz = (b_nqz > block_z) ? block_z : b_nqz;
			b_num_triangles = (b_num_triangles > block_t) ? block_t : b_num_triangles;
		#else
			b_nqx = (b_nqx > MIC_BLOCK_X_) ? MIC_BLOCK_X_ : b_nqx;
			b_nqy = (b_nqy > MIC_BLOCK_Y_) ? MIC_BLOCK_Y_ : b_nqy;
			b_nqz = (b_nqz > MIC_BLOCK_Z_) ? MIC_BLOCK_Z_ : b_nqz;
			b_num_triangles = (b_num_triangles > MIC_BLOCK_T_) ? MIC_BLOCK_T_ : b_num_triangles;
		#endif
	} // NumericFormFactorM::compute_hyperblock_size()


} // namespace hig
