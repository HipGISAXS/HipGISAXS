/***
  *  Project:
  *
  *  File: ff_num_gpu_kernels.cuh
  *  Created: Apr 11, 2013
  *  Modified: Thu 11 Apr 2013 02:34:55 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */


	/***
	 * special case when b_nqx == 1
	 */
	__global__ void form_factor_kernel_fused_nqx1(
						const float_t* qx, const float_t* qy, const cucomplex_t* qz,
						const float_t* shape_def, const short int* axes,
						const unsigned int curr_nqx, const unsigned int curr_nqy,
						const unsigned int curr_nqz, const unsigned int curr_num_triangles,
						const unsigned int b_nqx, const unsigned int b_nqy,
						const unsigned int b_nqz, const unsigned int b_num_triangles,
						const unsigned int ib_x, const unsigned int ib_y,
						const unsigned int ib_z, const unsigned int ib_t,
						cucomplex_t* ff) {
		unsigned int i_y = blockDim.x * blockIdx.x + threadIdx.x;
		unsigned int i_z = blockDim.y * blockIdx.y + threadIdx.y;
		unsigned int i_thread = blockDim.x * threadIdx.y + threadIdx.x;
		unsigned int num_threads = blockDim.x * blockDim.y;

		// sizes are:	shared_shape_def = T_PROP_SIZE_ * blockDim.x
		// 				shared_qx = curr_nqx	// the whole qx
		// 				shared_qy = blockDim.y
		// 				shared_qz = blockDim.z
		// reversed to fix alignment:
		float_t *shared_shape_def = (float_t*) dynamic_shared;
		cucomplex_t *shared_qz = (cucomplex_t*) &shared_shape_def[T_PROP_SIZE_ * curr_num_triangles];
		float_t *shared_qy = (float_t*) &shared_qz[blockDim.y];
		float_t *shared_qx = (float_t*) &shared_qy[blockDim.x];

		unsigned int i_shared, base_offset, num_loads;

		// load triangles
		unsigned int shape_def_size = T_PROP_SIZE_ * curr_num_triangles;
		num_loads = __float2int_ru(__fdividef(__int2float_ru(shape_def_size), num_threads));
		base_offset = T_PROP_SIZE_ * b_num_triangles * ib_t;
		for(int l = 0; l < num_loads; ++ l) {
			i_shared = num_threads * l + i_thread;
			if(i_shared < shape_def_size) shared_shape_def[i_shared] = shape_def[base_offset + i_shared];
		} // for

		// load qz
		unsigned int i_qz = b_nqz * ib_z + i_z;
		if(threadIdx.x == 0 && i_z < curr_nqz)
			shared_qz[threadIdx.y] = qz[i_qz];	// M: spread about access ...

		// load qy
		unsigned int i_qy = b_nqy * ib_y + i_y;
		if(threadIdx.y == 0 && i_y < curr_nqy)
			shared_qy[threadIdx.x] = qy[i_qy];	// M: spread about access ...

		// load qx
		if(i_thread == 0) shared_qx[0] = qx[ib_x];

		__syncthreads();	// sync to make sure all data is loaded and available

		cucomplex_t ff_tot = make_cuC((float_t) 0.0, (float_t) 0.0);
		if(i_y < curr_nqy && i_z < curr_nqz) {
			float_t temp_y = shared_qy[threadIdx.x];
			float_t qy2 = temp_y * temp_y;
			cucomplex_t temp_z = shared_qz[threadIdx.y];
			float_t temp_x = shared_qx[0];

			for(int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
				unsigned int shape_off = T_PROP_SIZE_ * i_t;
				float_t s = shared_shape_def[shape_off];
				float_t nx = shared_shape_def[shape_off + 1];
				float_t ny = shared_shape_def[shape_off + 2];
				float_t nz = shared_shape_def[shape_off + 3];
				float_t x = shared_shape_def[shape_off + 4];
				float_t y = shared_shape_def[shape_off + 5];
				float_t z = shared_shape_def[shape_off + 6];

				cucomplex_t qz2, qzn, qzt;
				cucomplex_t qn_d, qt_d;

				float_t qyn = temp_y * ny;
				float_t qyt = temp_y * y;
				compute_z(temp_z, nz, z, qz2, qzn, qzt);

				compute_x(temp_x, qy2, qz2, nx, qyn, qzn, x, qyt, qzt, qn_d, qt_d);
				cucomplex_t fq_temp = compute_fq(s, qt_d, qn_d);
				ff_tot = ff_tot + fq_temp;
			} // for
			ff[curr_nqy * i_z + i_y] = ff_tot;
		} // if
	} // form_factor_kernel_fused()


