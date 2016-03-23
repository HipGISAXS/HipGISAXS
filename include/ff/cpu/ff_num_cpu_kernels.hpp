/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */


//  #ifndef __SSE3__	// when not using sse3 optimizations

	/**
	 * the main Form Factor kernel with fused reduction function - for one hyperblock.
	 */
	void NumericFormFactorC::form_factor_kernel_fused_unroll4(real_t* qx, real_t* qy, complex_t* qz,
					real_vec_t& shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
					unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					real_t* rot,
					complex_t* ff) {
		if(ff == NULL || qx == NULL || qy == NULL || qz == NULL) return;

		unsigned long int xy_size = (unsigned long int) curr_nqx * curr_nqy;
		unsigned int start_z = b_nqz * ib_z;
		unsigned int start_y = b_nqy * ib_y;
		unsigned int start_x = b_nqx * ib_x;
		unsigned int start_t = b_num_triangles * ib_t * CPU_T_PROP_SIZE_;

		int unroll_rem = curr_num_triangles % 4;
		int unroll_off = curr_num_triangles - unroll_rem;

		#pragma omp parallel
		{
			#pragma omp for collapse(3)
			for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
				for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
					for(int i_x = 0; i_x < curr_nqx; ++ i_x) {
						unsigned long int super_i = (unsigned long int) nqx * nqy * (ib_z * b_nqz + i_z) +
													nqx * (ib_y * b_nqy + i_y) + ib_x * b_nqx + i_x;
						complex_t temp_qz = qz[start_z + i_z];
						real_t temp_qy = qy[start_y + i_y];
						real_t temp_qx = qx[start_x + i_x];

						// compute rotation stuff ... FIXME double check transposition etc ...
						complex_t temp_x = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
						complex_t temp_y = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
						complex_t temp_z = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;

						complex_t qz2 = temp_z * temp_z;
						complex_t qy2 = temp_y * temp_y;
						complex_t q2 = temp_x * temp_x + qy2 + qz2;

						for(int i_t = 0; i_t < curr_num_triangles; i_t += 4) {
							unsigned int shape_off1 = start_t + i_t * CPU_T_PROP_SIZE_;
							real_t s1 = shape_def[shape_off1];
							real_t nx1 = shape_def[shape_off1 + 1];
							real_t ny1 = shape_def[shape_off1 + 2];
							real_t nz1 = shape_def[shape_off1 + 3];
							real_t x1 = shape_def[shape_off1 + 4];
							real_t y1 = shape_def[shape_off1 + 5];
							real_t z1 = shape_def[shape_off1 + 6];
							unsigned int shape_off2 = start_t + (i_t + 1) * CPU_T_PROP_SIZE_;
							real_t s2 = shape_def[shape_off2];
							real_t nx2 = shape_def[shape_off2 + 1];
							real_t ny2 = shape_def[shape_off2 + 2];
							real_t nz2 = shape_def[shape_off2 + 3];
							real_t x2 = shape_def[shape_off2 + 4];
							real_t y2 = shape_def[shape_off2 + 5];
							real_t z2 = shape_def[shape_off2 + 6];
							unsigned int shape_off3 = start_t + (i_t + 2) * CPU_T_PROP_SIZE_;
							real_t s3 = shape_def[shape_off3];
							real_t nx3 = shape_def[shape_off3 + 1];
							real_t ny3 = shape_def[shape_off3 + 2];
							real_t nz3 = shape_def[shape_off3 + 3];
							real_t x3 = shape_def[shape_off3 + 4];
							real_t y3 = shape_def[shape_off3 + 5];
							real_t z3 = shape_def[shape_off3 + 6];
							unsigned int shape_off4 = start_t + (i_t + 3) * CPU_T_PROP_SIZE_;
							real_t s4 = shape_def[shape_off4];
							real_t nx4 = shape_def[shape_off4 + 1];
							real_t ny4 = shape_def[shape_off4 + 2];
							real_t nz4 = shape_def[shape_off4 + 3];
							real_t x4 = shape_def[shape_off4 + 4];
							real_t y4 = shape_def[shape_off4 + 5];
							real_t z4 = shape_def[shape_off4 + 6];
	
							complex_t qzn1 = temp_z * nz1;
							complex_t qzt1 = temp_z * z1;
							complex_t qzn2 = temp_z * nz2;
							complex_t qzt2 = temp_z * z2;
							complex_t qzn3 = temp_z * nz3;
							complex_t qzt3 = temp_z * z3;
							complex_t qzn4 = temp_z * nz4;
							complex_t qzt4 = temp_z * z4;
	
							complex_t qyn1 = temp_y * ny1;
							complex_t qyt1 = temp_y * y1;
							complex_t qyn2 = temp_y * ny2;
							complex_t qyt2 = temp_y * y2;
							complex_t qyn3 = temp_y * ny3;
							complex_t qyt3 = temp_y * y3;
							complex_t qyn4 = temp_y * ny4;
							complex_t qyt4 = temp_y * y4;
	
							complex_t qt1 = temp_x * x1 + qyt1 + qzt1;
							complex_t qn1 = (temp_x * nx1 + qyn1 + qzn1) / q2;
							complex_t qt2 = temp_x * x2 + qyt2 + qzt2;
							complex_t qn2 = (temp_x * nx2 + qyn2 + qzn2) / q2;
							complex_t qt3 = temp_x * x3 + qyt3 + qzt3;
							complex_t qn3 = (temp_x * nx3 + qyn3 + qzn3) / q2;
							complex_t qt4 = temp_x * x4 + qyt4 + qzt4;
							complex_t qn4 = (temp_x * nx4 + qyn4 + qzn4) / q2;

							complex_t fq1 = compute_fq(s1, qt1, qn1);
							complex_t fq2 = compute_fq(s2, qt2, qn2);
							complex_t fq3 = compute_fq(s3, qt3, qn3);
							complex_t fq4 = compute_fq(s4, qt4, qn4);

							ff[super_i] += fq1 + fq2 + fq3 + fq4;
						} // for t

						unsigned int shape_off;
						real_t s, nx, ny, nz, x, y, z;
						complex_t qzn, qzt;
						complex_t qyn, qyt;
						complex_t qn, qt;
	
						switch(unroll_rem) {
							case 3:
								shape_off = start_t + (unroll_off + 2) * CPU_T_PROP_SIZE_;
								s = shape_def[shape_off];
								nx = shape_def[shape_off + 1];
								ny = shape_def[shape_off + 2];
								nz = shape_def[shape_off + 3];
								x = shape_def[shape_off + 4];
								y = shape_def[shape_off + 5];
								z = shape_def[shape_off + 6];
	
								qzn = temp_z * nz;
								qzt = temp_z * z;
								qyn = temp_y * ny;
								qyt = temp_y * y;
								qt = temp_x * x + qyt + qzt;
								qn = (temp_x * nx + qyn + qzn) / q2;

								ff[super_i] += compute_fq(s, qt, qn);

							case 2:
								shape_off = start_t + (unroll_off + 1) * CPU_T_PROP_SIZE_;
								s = shape_def[shape_off];
								nx = shape_def[shape_off + 1];
								ny = shape_def[shape_off + 2];
								nz = shape_def[shape_off + 3];
								x = shape_def[shape_off + 4];
								y = shape_def[shape_off + 5];
								z = shape_def[shape_off + 6];
	
								qzn = temp_z * nz;
								qzt = temp_z * z;
								qyn = temp_y * ny;
								qyt = temp_y * y;
								qt = temp_x * x + qyt + qzt;
								qn = (temp_x * nx + qyn + qzn) / q2;

								ff[super_i] += compute_fq(s, qt, qn);

							case 1:
								shape_off = start_t + (unroll_off) * CPU_T_PROP_SIZE_;
								s = shape_def[shape_off];
								nx = shape_def[shape_off + 1];
								ny = shape_def[shape_off + 2];
								nz = shape_def[shape_off + 3];
								x = shape_def[shape_off + 4];
								y = shape_def[shape_off + 5];
								z = shape_def[shape_off + 6];
	
								qzn = temp_z * nz;
								qzt = temp_z * z;
								qyn = temp_y * ny;
								qyt = temp_y * y;
								qt = temp_x * x + qyt + qzt;
								qn = (temp_x * nx + qyn + qzn) / q2;

								ff[super_i] += compute_fq(s, qt, qn);

							case 0:
								// nothing to do
								break;
						} // switch()
					} // for x
				} // for y
			} // for z
		} // pragma omp parallel
	} // NumericFormFactorC::form_factor_kernel_fused_unroll4()

//	#endif
	
	/**
	 * special case of nqx == 1 of Form Factor kernel with fused reduction
	 */
	void NumericFormFactorC::form_factor_kernel_fused_nqx1(
					const real_t* __restrict__ qx, const real_t* __restrict__ qy,
					const complex_t* __restrict__ qz,
//					#ifndef __SSE3__
						real_vec_t& shape_def,
//					#else
//						real_t* __restrict__ shape_def,
//					#endif
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
					unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					real_t* rot,
					complex_t* ff) {
		if(ff == NULL || qx == NULL || qy == NULL || qz == NULL) return;

		if(nqx != 1) {
			std::cout << "error: wrong form factor function called!" << std::endl;
			return;
		} // if
	
		unsigned long int xy_size = (unsigned long int) curr_nqy;
		unsigned int start_z = b_nqz * ib_z;
		unsigned int start_y = b_nqy * ib_y;
		unsigned int start_x = 0;
		unsigned int start_t = b_num_triangles * ib_t * CPU_T_PROP_SIZE_;

		//#ifdef PROFILE_PAPI
			//int papi_events[3] = { PAPI_FP_OPS, PAPI_SP_OPS, PAPI_DP_OPS };
			//int papi_events[3] = { PAPI_FML_INS, PAPI_FAD_INS, PAPI_FDV_INS };
			//long long papi_counter_vals[3];
		//#endif

//		#ifndef __SSE3__	// fallback

			#pragma omp parallel
			{
				#pragma omp for collapse(2)
				for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
					for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
						//#ifdef PROFILE_PAPI
							//if(ib_y + ib_z + ib_t + i_z + i_y == 0) PAPI_start_counters(papi_events, 3);
						//#endif

						unsigned long int super_i = (unsigned long int) nqy * (ib_z * b_nqz + i_z) +
													(ib_y * b_nqy + i_y);
						complex_t temp_z = qz[start_z + i_z];
						real_t temp_y = qy[start_y + i_y];
						real_t temp_x = qx[0];

						// compute rotation stuff ... FIXME double check transposition etc ...
						complex_t mqx = rot[0] * temp_x + rot[1] * temp_y + rot[2] * temp_z;
						complex_t mqy = rot[3] * temp_x + rot[4] * temp_y + rot[5] * temp_z;
						complex_t mqz = rot[6] * temp_x + rot[7] * temp_y + rot[8] * temp_z;

						//complex_t qz2 = temp_z * temp_z;
						//real_t qy2 = temp_y * temp_y;
						//real_t qx2 = temp_x * temp_x;
						complex_t qz2 = mqz * mqz;
						complex_t qy2 = mqy * mqy;
						complex_t qx2 = mqx * mqx;
						complex_t q2 = qx2 + qy2 + qz2;
						complex_t q2_inv = ((real_t) 1.0) / q2;

						complex_t total(0.0, 0.0);

						for(int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
							unsigned int shape_off = start_t + i_t * CPU_T_PROP_SIZE_;
							real_t s = shape_def[shape_off];
							real_t nx = shape_def[shape_off + 1];
							real_t ny = shape_def[shape_off + 2];
							real_t nz = shape_def[shape_off + 3];
							real_t x = shape_def[shape_off + 4];
							real_t y = shape_def[shape_off + 5];
							real_t z = shape_def[shape_off + 6];

							//complex_t qzn = temp_z * nz;
							//complex_t qzt = temp_z * z;
							//real_t qyn = temp_y * ny;
							//real_t qyt = temp_y * y;
							//real_t qxn = temp_x * nx;
							//real_t qxt = temp_x * x;
							complex_t qzn = mqz * nz;
							complex_t qzt = mqz * z;
							complex_t qyn = mqy * ny;
							complex_t qyt = mqy * y;
							complex_t qxn = mqx * nx;
							complex_t qxt = mqx * x;
							complex_t qt = qxt + qyt + qzt;
							complex_t qn = (qxn + qyn + qzn) * q2_inv;
							complex_t fq = compute_fq(s, qt, qn);

							total += fq;
						} // for t

						ff[super_i] += total;

						//#ifdef PROFILE_PAPI
							//if(ib_y + ib_z + ib_t + i_z + i_y == 0) {
							//	PAPI_stop_counters(papi_counter_vals, 3);
							//	std::cout << "==== FP_OPS: " << papi_counter_vals[0] << std::endl;
							//	std::cout << "==== SP_OPS: " << papi_counter_vals[1] << std::endl;
							//	std::cout << "==== DP_OPS: " << papi_counter_vals[2] << std::endl;
							//} // if
						//#endif
					} // for y
				} // for z
			} // pragma omp parallel

/*		#elif defined INTEL_SB_AVX	// avx specific optimizations for intel sandy bridge (edison)

			// FIXME: assuming only float for now (4 bytes per float)
			unsigned int shape_padding = (8 - (num_triangles & 7)) & 7;
			unsigned int vec_size = 8; // 32 / sizeof(real_t);
			unsigned int padded_num_triangles = num_triangles + shape_padding;
			start_t = b_num_triangles * ib_t;

			// shape padding guarantees that padded_num_triangles is a multiple of 8
			// FIXME: assuming that b_num_triangles is a multiple of 8 as well ...

			#pragma omp parallel
			{
				#pragma omp for collapse(2)
				for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
					for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
						//#ifdef PROFILE_PAPI
						//	if(ib_y + ib_z + ib_t + i_z + i_y == 0) PAPI_start_counters(papi_events, 3);
						//#endif

						// ////////////////////////////////////////////////////
						// intrinsic wrapper naming (only for floating-point)
						// avx_xxx_abc  => a = r|c, b = p|s, c = s|d
						// avx_xxx_abcd => a = r|c, b = r|c, c = p|s, d = s|d
						// r = real,             c = complex,
						// p = packed (vector),  s = scalar,
						// s = single-precision, d = double-precision
						// ////////////////////////////////////////////////////

						// TODO this kernel has not been updated for the rotation stuff ...

						unsigned long int super_i = (unsigned long int) nqy * (ib_z * b_nqz + i_z) +
													(ib_y * b_nqy + i_y);

						// rotation stuff ... TODO: optimize later
						real_t temp_qx = qx[start_x], temp_qy = qy[start_y + i_y];
						complex_t temp_qz = qz[start_z + i_z];
						complex_t mqx = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
						complex_t mqy = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
						complex_t mqz = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;


						avx_m256c_t temp_z = avx_set1_cps(mqz);
						avx_m256c_t temp_y = avx_set1_cps(mqy);
						avx_m256c_t temp_x = avx_set1_cps(mqx);

						avx_m256c_t qx2 = avx_mul_ccps(temp_x, temp_x);
						avx_m256c_t qy2 = avx_mul_ccps(temp_y, temp_y);
						avx_m256c_t qz2 = avx_mul_ccps(temp_z, temp_z);
						avx_m256c_t qxy2 = avx_add_ccps(qx2, qy2);
						avx_m256c_t q2 = avx_add_ccps(qxy2, qz2);
						avx_m256c_t q2_inv = avx_rcp_cps(q2);

						// unrolling twice, using function call fusion!

						avx_m256c_t total1 = avx_setzero_cps();
						avx_m256c_t total2 = avx_setzero_cps();

						for(int i_t = 0; i_t < curr_num_triangles; i_t += 2 * vec_size) {
							unsigned int shape_off = start_t + i_t;
							avx_m256_t s1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t s2 = avx_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							avx_m256_t nx1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t nx2 = avx_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							avx_m256_t ny1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t ny2 = avx_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							avx_m256_t nz1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t nz2 = avx_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							avx_m256_t x1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t x2 = avx_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							avx_m256_t y1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t y2 = avx_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							avx_m256_t z1 = avx_load_rps(& shape_def[shape_off]);
							avx_m256_t z2 = avx_load_rps(& shape_def[shape_off + vec_size]);

							avx_m256c_t qxn1 = avx_mul_crps(temp_x, nx1);
							avx_m256c_t qxn2 = avx_mul_crps(temp_x, nx2);
							avx_m256c_t qxt1 = avx_mul_crps(temp_x, x1);
							avx_m256c_t qxt2 = avx_mul_crps(temp_x, x2);
							avx_m256c_t qyn1 = avx_mul_crps(temp_y, ny1);
							avx_m256c_t qyn2 = avx_mul_crps(temp_y, ny2);
							avx_m256c_t qyt1 = avx_mul_crps(temp_y, y1);
							avx_m256c_t qyt2 = avx_mul_crps(temp_y, y2);
							avx_m256c_t qzn1 = avx_mul_crps(temp_z, nz1);
							avx_m256c_t qzn2 = avx_mul_crps(temp_z, nz2);
							avx_m256c_t qzt1 = avx_mul_crps(temp_z, z1);
							avx_m256c_t qzt2 = avx_mul_crps(temp_z, z2);
							avx_m256c_t qxyt1 = avx_add_ccps(qxt1, qyt1);
							avx_m256c_t qxyt2 = avx_add_ccps(qxt2, qyt2);
							avx_m256c_t qxyn1 = avx_add_ccps(qxn1, qyn1);
							avx_m256c_t qxyn2 = avx_add_ccps(qxn2, qyn2);
							avx_m256c_t qt1 = avx_add_ccps(qxyt1, qzt1);
							avx_m256c_t qt2 = avx_add_ccps(qxyt2, qzt2);
							avx_m256c_t temp_qn1 = avx_add_ccps(qxyn1, qzn1);
							avx_m256c_t temp_qn2 = avx_add_ccps(qxyn2, qzn2);
							avx_m256c_t qn1 = avx_mul_ccps(temp_qn1, q2_inv);
							avx_m256c_t qn2 = avx_mul_ccps(temp_qn2, q2_inv);
							avx_m256c_t temp1, temp2;
							avx_sincos_rps_dual(qt1.xvec, qt2.xvec, &temp1.yvec, &temp2.yvec,
												&temp1.xvec, &temp2.xvec);
							//avx_sincos_rps(qt1.xvec, &temp1.yvec, &temp1.xvec);
							//avx_sincos_rps(qt2.xvec, &temp2.yvec, &temp2.xvec);
							avx_m256_t temp_v21, temp_v22;
							avx_exp_rps_dual(qt1.yvec, qt2.yvec, &temp_v21, &temp_v22);
							//temp_v21 = avx_exp_rps(qt1.yvec);
							//temp_v22 = avx_exp_rps(qt2.yvec);
							avx_m256_t v21 = avx_mul_rrps(s1, temp_v21);
							avx_m256_t v22 = avx_mul_rrps(s2, temp_v22);
							avx_m256c_t v11 = avx_mul_ccps(qn1, temp1);
							avx_m256c_t v12 = avx_mul_ccps(qn2, temp2);
							avx_m256c_t fq1 = avx_mul_crps(v11, v21);
							avx_m256c_t fq2 = avx_mul_crps(v12, v22);

							total1 = avx_add_ccps(total1, fq1);
							total2 = avx_add_ccps(total2, fq2);
						} // for t

						total1 = avx_hadd_ccps(total1, total1);
						total2 = avx_hadd_ccps(total2, total2);
						total1 = avx_hadd_ccps(total1, total1);
						total2 = avx_hadd_ccps(total2, total2);

						total1 = avx_add_ccps(total1, total2);
						avx_addstore_css(&(ff[super_i]), total1);

						//#ifdef PROFILE_PAPI
						//	if(ib_y + ib_z + ib_t + i_z + i_y == 0) {
						//		PAPI_stop_counters(papi_counter_vals, 3);
						//		std::cout << "==== FP_OPS: " << papi_counter_vals[0] << std::endl;
						//		std::cout << "==== SP_OPS: " << papi_counter_vals[1] << std::endl;
						//		std::cout << "==== DP_OPS: " << papi_counter_vals[2] << std::endl;
						//	} // if
						//#endif
					} // for y
				} // for z
			} // pragma omp parallel

		#else		// sse3 specific optimizations (also for amd magny cours (hopper))

			// FIXME: assuming only float for now (4 bytes per float) do for double also
			unsigned int shape_padding = (4 - (num_triangles & 3)) & 3;
			unsigned int vec_size = 4; //16 / sizeof(real_t);	// FIXME: for now assuming SP (float) only ...
			unsigned int padded_num_triangles = num_triangles + shape_padding;
			start_t = b_num_triangles * ib_t;

			// shape padding guarantees that padded_num_triangles is a multiple of 4
			// FIXME: assuming that b_num_triangles is a multiple of 4 as well ...

			#pragma omp parallel
			{
				#pragma omp for collapse(2)
				for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
					for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
						//#ifdef PROFILE_PAPI
						//	if(ib_y + ib_z + ib_t + i_z + i_y == 0) PAPI_start_counters(papi_events, 3);
						//#endif

						// ////////////////////////////////////////////////////
						// intrinsic wrapper naming (only for floating-point)
						// _mm_xxx_abc  => a = r|c, b = p|s, c = s|d
						// _mm_xxx_abcd => a = r|c, b = r|c, c = p|s, d = s|d
						// r = real,             c = complex,
						// p = packed (vector),  s = scalar,
						// s = single-precision, d = double-precision
						// ////////////////////////////////////////////////////

						unsigned long int super_i = (unsigned long int) nqy * (ib_z * b_nqz + i_z) +
													(ib_y * b_nqy + i_y);

						// rotation stuff ... TODO: optimize later
						real_t temp_qx = qx[start_x], temp_qy = qy[start_y + i_y];
						complex_t temp_qz = qz[start_z + i_z];
						complex_t mqx = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
						complex_t mqy = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
						complex_t mqz = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;

						//sse_m128c_t temp_z = sse_set1_cps(qz[start_z + i_z]);
						//sse_m128_t temp_y = sse_set1_rps(qy[start_y + i_y]);
						//sse_m128_t temp_x = sse_set1_rps(qx[start_x]);

						sse_m128c_t temp_z = sse_set1_cps(mqz);
						sse_m128c_t temp_y = sse_set1_cps(mqy);
						sse_m128c_t temp_x = sse_set1_cps(mqx);

						//sse_m128_t qx2 = sse_mul_rrps(temp_x, temp_x);
						//sse_m128_t qy2 = sse_mul_rrps(temp_y, temp_y);
						//sse_m128c_t qz2 = sse_mul_ccps(temp_z, temp_z);
						//sse_m128_t qxy2 = sse_add_rrps(qx2, qy2);
						//sse_m128c_t q2 = sse_add_rcps(qxy2, qz2);
						//sse_m128c_t q2_inv = sse_rcp_cps(q2);

						sse_m128c_t qx2 = sse_mul_ccps(temp_x, temp_x);
						sse_m128c_t qy2 = sse_mul_ccps(temp_y, temp_y);
						sse_m128c_t qz2 = sse_mul_ccps(temp_z, temp_z);
						sse_m128c_t qxy2 = sse_add_ccps(qx2, qy2);
						sse_m128c_t q2 = sse_add_ccps(qxy2, qz2);
						sse_m128c_t q2_inv = sse_rcp_cps(q2);

						// unrolling twice, using function call fusion!

						sse_m128c_t total1 = sse_setzero_cps();
						sse_m128c_t total2 = sse_setzero_cps();

						for(int i_t = 0; i_t < curr_num_triangles; i_t += 2 * vec_size) {
							unsigned int shape_off = start_t + i_t;
							sse_m128_t s1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t s2 = sse_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							sse_m128_t nx1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t nx2 = sse_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							sse_m128_t ny1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t ny2 = sse_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							sse_m128_t nz1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t nz2 = sse_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							sse_m128_t x1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t x2 = sse_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							sse_m128_t y1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t y2 = sse_load_rps(& shape_def[shape_off + vec_size]);
							shape_off += padded_num_triangles;
							sse_m128_t z1 = sse_load_rps(& shape_def[shape_off]);
							sse_m128_t z2 = sse_load_rps(& shape_def[shape_off + vec_size]);

							//sse_m128_t qxn1 = sse_mul_rrps(temp_x, nx1);
							//sse_m128_t qxn2 = sse_mul_rrps(temp_x, nx2);
							//sse_m128_t qxt1 = sse_mul_rrps(temp_x, x1);
							//sse_m128_t qxt2 = sse_mul_rrps(temp_x, x2);
							//sse_m128_t qyn1 = sse_mul_rrps(temp_y, ny1);
							//sse_m128_t qyn2 = sse_mul_rrps(temp_y, ny2);
							//sse_m128_t qyt1 = sse_mul_rrps(temp_y, y1);
							//sse_m128_t qyt2 = sse_mul_rrps(temp_y, y2);
							sse_m128c_t qxn1 = sse_mul_crps(temp_x, nx1);
							sse_m128c_t qxn2 = sse_mul_crps(temp_x, nx2);
							sse_m128c_t qxt1 = sse_mul_crps(temp_x, x1);
							sse_m128c_t qxt2 = sse_mul_crps(temp_x, x2);
							sse_m128c_t qyn1 = sse_mul_crps(temp_y, ny1);
							sse_m128c_t qyn2 = sse_mul_crps(temp_y, ny2);
							sse_m128c_t qyt1 = sse_mul_crps(temp_y, y1);
							sse_m128c_t qyt2 = sse_mul_crps(temp_y, y2);

							sse_m128c_t qzn1 = sse_mul_crps(temp_z, nz1);
							sse_m128c_t qzn2 = sse_mul_crps(temp_z, nz2);
							sse_m128c_t qzt1 = sse_mul_crps(temp_z, z1);
							sse_m128c_t qzt2 = sse_mul_crps(temp_z, z2);

							//sse_m128_t qxyt1 = sse_add_rrps(qxt1, qyt1);
							//sse_m128_t qxyt2 = sse_add_rrps(qxt2, qyt2);
							//sse_m128_t qxyn1 = sse_add_rrps(qxn1, qyn1);
							//sse_m128_t qxyn2 = sse_add_rrps(qxn2, qyn2);
							sse_m128c_t qxyt1 = sse_add_ccps(qxt1, qyt1);
							sse_m128c_t qxyt2 = sse_add_ccps(qxt2, qyt2);
							sse_m128c_t qxyn1 = sse_add_ccps(qxn1, qyn1);
							sse_m128c_t qxyn2 = sse_add_ccps(qxn2, qyn2);

							//sse_m128c_t qt1 = sse_add_rcps(qxyt1, qzt1);
							//sse_m128c_t qt2 = sse_add_rcps(qxyt2, qzt2);
							//sse_m128c_t temp_qn1 = sse_add_rcps(qxyn1, qzn1);
							//sse_m128c_t temp_qn2 = sse_add_rcps(qxyn2, qzn2);
							sse_m128c_t qt1 = sse_add_ccps(qxyt1, qzt1);
							sse_m128c_t qt2 = sse_add_ccps(qxyt2, qzt2);
							sse_m128c_t temp_qn1 = sse_add_ccps(qxyn1, qzn1);
							sse_m128c_t temp_qn2 = sse_add_ccps(qxyn2, qzn2);

							sse_m128c_t qn1 = sse_mul_ccps(temp_qn1, q2_inv);
							sse_m128c_t qn2 = sse_mul_ccps(temp_qn2, q2_inv);
							sse_m128c_t temp1, temp2;
							sse_sincos_rps_dual(qt1.xvec, qt2.xvec, &temp1.yvec, &temp2.yvec,
												&temp1.xvec, &temp2.xvec);
							//sse_sincos_rps(qt1.xvec, &temp1.yvec, &temp1.xvec);
							//sse_sincos_rps(qt2.xvec, &temp2.yvec, &temp2.xvec);
							sse_m128_t temp_v21, temp_v22;
							sse_exp_rps_dual(qt1.yvec, qt2.yvec, &temp_v21, &temp_v22);
							//temp_v21 = sse_exp_rps(qt1.yvec);
							//temp_v22 = sse_exp_rps(qt2.yvec);
							sse_m128_t v21 = sse_mul_rrps(s1, temp_v21);
							sse_m128_t v22 = sse_mul_rrps(s2, temp_v22);
							sse_m128c_t v11 = sse_mul_ccps(qn1, temp1);
							sse_m128c_t v12 = sse_mul_ccps(qn2, temp2);
							sse_m128c_t fq1 = sse_mul_crps(v11, v21);
							sse_m128c_t fq2 = sse_mul_crps(v12, v22);

							total1 = sse_add_ccps(total1, fq1);
							total2 = sse_add_ccps(total2, fq2);
						} // for t

						total1 = sse_hadd_ccps(total1, total1);
						total2 = sse_hadd_ccps(total2, total2);
						total1 = sse_hadd_ccps(total1, total1);
						total2 = sse_hadd_ccps(total2, total2);

						total1 = sse_add_ccps(total1, total2);
						sse_addstore_css(&(ff[super_i]), total1);

						//#ifdef PROFILE_PAPI
						//	if(ib_y + ib_z + ib_t + i_z + i_y == 0) {
						//		PAPI_stop_counters(papi_counter_vals, 3);
						//		std::cout << "==== FP_OPS: " << papi_counter_vals[0] << std::endl;
						//		std::cout << "==== SP_OPS: " << papi_counter_vals[1] << std::endl;
						//		std::cout << "==== DP_OPS: " << papi_counter_vals[2] << std::endl;
						//	} // if
						//#endif
					} // for y
				} // for z
			} // pragma omp parallel

		#endif // SSE3 AVX etc */
	} // NumericFormFactorC::form_factor_kernel_fused_nqx1()


/*	#ifdef INTEL_SB_AVX

		inline avx_m256c_t NumericFormFactorC::avx_compute_fq(avx_m256_t s,
												avx_m256c_t qt, avx_m256c_t qn) {
			avx_m256c_t temp;
			//temp.xvec = avx_cos_rps(qt.xvec);
			//temp.yvec = avx_sin_rps(qt.xvec);
			avx_sincos_rps(qt.xvec, &temp.yvec, &temp.xvec);
			avx_m256_t v2 = avx_mul_rrps(s, avx_exp_rps(qt.yvec));
			avx_m256c_t v1 = avx_mul_ccps(qn, temp);
			return avx_mul_crps(v1, v2);
		} // NumericFormFactorC::avx_compute_fq()

	#elif defined __SSE3__

		inline sse_m128c_t NumericFormFactorC::sse_compute_fq(sse_m128_t s,
												sse_m128c_t qt, sse_m128c_t qn) {
			sse_m128c_t temp;
			//temp.xvec = sse_cos_rps(qt.xvec);
			//temp.yvec = sse_sin_rps(qt.xvec);
			sse_sincos_rps(qt.xvec, &temp.yvec, &temp.xvec);
			sse_m128_t v2 = sse_mul_rrps(s, sse_exp_rps(qt.yvec));
			sse_m128c_t v1 = sse_mul_ccps(qn, temp);
			return sse_mul_crps(v1, v2);
		} // NumericFormFactorC::sse_compute_fq()

	#endif
*/


//	#ifndef __SSE3__	// not using sse3 optimizations

	/**
	 * special case of nqx == 1 of Form Factor kernel with fused reduction, and loop unrolling
	 */
	void NumericFormFactorC::form_factor_kernel_fused_nqx1_unroll4(
					real_t* qx, real_t* qy, complex_t* qz,
					real_vec_t& shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz,
					unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					real_t* rot,
					complex_t* ff) {
		if(ff == NULL || qx == NULL || qy == NULL || qz == NULL) return;

		if(nqx != 1) {
			std::cout << "error: wrong form factor function called!" << std::endl;
			return;
		} // if
	
		unsigned long int xy_size = (unsigned long int) curr_nqy;
		unsigned int start_z = b_nqz * ib_z;
		unsigned int start_y = b_nqy * ib_y;
		unsigned int start_x = 0;
		unsigned int start_t = b_num_triangles * ib_t * CPU_T_PROP_SIZE_;

		//int unroll_rem = curr_num_triangles % 4;
		int unroll_rem = curr_num_triangles & 3;
		int unroll_off = curr_num_triangles - unroll_rem;

		#pragma omp parallel
		{
			#pragma omp for collapse(2)
			for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
				for(int i_y = 0; i_y < curr_nqy; ++ i_y) {
					unsigned long int super_i = (unsigned long int) nqy * (ib_z * b_nqz + i_z) +
													(ib_y * b_nqy + i_y);
					complex_t temp_qz = qz[start_z + i_z];
					real_t temp_qy = qy[start_y + i_y];
					real_t temp_qx = qx[0];

					// compute rotation stuff ... FIXME double check transposition etc ...
					complex_t temp_x = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
					complex_t temp_y = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
					complex_t temp_z = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;

					complex_t qz2 = temp_z * temp_z;
					complex_t qy2 = temp_y * temp_y;
					complex_t q2 = temp_x * temp_x + qy2 + qz2;
					complex_t q2_inv = ((real_t) 1.0) / q2;

					complex_t total(0.0, 0.0);

					for(int i_t = 0; i_t < curr_num_triangles; i_t += 4) {
						unsigned int shape_off1 = start_t + i_t * CPU_T_PROP_SIZE_;
						real_t s1 = shape_def[shape_off1];
						real_t nx1 = shape_def[shape_off1 + 1];
						real_t ny1 = shape_def[shape_off1 + 2];
						real_t nz1 = shape_def[shape_off1 + 3];
						real_t x1 = shape_def[shape_off1 + 4];
						real_t y1 = shape_def[shape_off1 + 5];
						real_t z1 = shape_def[shape_off1 + 6];
						unsigned int shape_off2 = start_t + (i_t + 1) * CPU_T_PROP_SIZE_;
						real_t s2 = shape_def[shape_off2];
						real_t nx2 = shape_def[shape_off2 + 1];
						real_t ny2 = shape_def[shape_off2 + 2];
						real_t nz2 = shape_def[shape_off2 + 3];
						real_t x2 = shape_def[shape_off2 + 4];
						real_t y2 = shape_def[shape_off2 + 5];
						real_t z2 = shape_def[shape_off2 + 6];
						unsigned int shape_off3 = start_t + (i_t + 2) * CPU_T_PROP_SIZE_;
						real_t s3 = shape_def[shape_off3];
						real_t nx3 = shape_def[shape_off3 + 1];
						real_t ny3 = shape_def[shape_off3 + 2];
						real_t nz3 = shape_def[shape_off3 + 3];
						real_t x3 = shape_def[shape_off3 + 4];
						real_t y3 = shape_def[shape_off3 + 5];
						real_t z3 = shape_def[shape_off3 + 6];
						unsigned int shape_off4 = start_t + (i_t + 3) * CPU_T_PROP_SIZE_;
						real_t s4 = shape_def[shape_off4];
						real_t nx4 = shape_def[shape_off4 + 1];
						real_t ny4 = shape_def[shape_off4 + 2];
						real_t nz4 = shape_def[shape_off4 + 3];
						real_t x4 = shape_def[shape_off4 + 4];
						real_t y4 = shape_def[shape_off4 + 5];
						real_t z4 = shape_def[shape_off4 + 6];

						complex_t qzn1 = temp_z * nz1;
						complex_t qzt1 = temp_z * z1;
						complex_t qzn2 = temp_z * nz2;
						complex_t qzt2 = temp_z * z2;
						complex_t qzn3 = temp_z * nz3;
						complex_t qzt3 = temp_z * z3;
						complex_t qzn4 = temp_z * nz4;
						complex_t qzt4 = temp_z * z4;

						complex_t qyn1 = temp_y * ny1;
						complex_t qyt1 = temp_y * y1;
						complex_t qyn2 = temp_y * ny2;
						complex_t qyt2 = temp_y * y2;
						complex_t qyn3 = temp_y * ny3;
						complex_t qyt3 = temp_y * y3;
						complex_t qyn4 = temp_y * ny4;
						complex_t qyt4 = temp_y * y4;

						complex_t qt1 = temp_x * x1 + qyt1 + qzt1;
						complex_t qn1 = (temp_x * nx1 + qyn1 + qzn1) * q2_inv;
						complex_t qt2 = temp_x * x2 + qyt2 + qzt2;
						complex_t qn2 = (temp_x * nx2 + qyn2 + qzn2) * q2_inv;
						complex_t qt3 = temp_x * x3 + qyt3 + qzt3;
						complex_t qn3 = (temp_x * nx3 + qyn3 + qzn3) * q2_inv;
						complex_t qt4 = temp_x * x4 + qyt4 + qzt4;
						complex_t qn4 = (temp_x * nx4 + qyn4 + qzn4) * q2_inv;

						complex_t fq1 = compute_fq(s1, qt1, qn1);
						complex_t fq2 = compute_fq(s2, qt2, qn2);
						complex_t fq3 = compute_fq(s3, qt3, qn3);
						complex_t fq4 = compute_fq(s4, qt4, qn4);

						//ff[super_i] += fq1 + fq2 + fq3 + fq4;
						total += fq1 + fq2 + fq3 + fq4;
					} // for t

					unsigned int shape_off;
					real_t s, nx, ny, nz, x, y, z;
					complex_t qzn, qzt;
					complex_t qyn, qyt;
					complex_t qn, qt;

					switch(unroll_rem) {
						case 3:
							shape_off = start_t + (unroll_off + 2) * CPU_T_PROP_SIZE_;
							s = shape_def[shape_off];
							nx = shape_def[shape_off + 1];
							ny = shape_def[shape_off + 2];
							nz = shape_def[shape_off + 3];
							x = shape_def[shape_off + 4];
							y = shape_def[shape_off + 5];
							z = shape_def[shape_off + 6];

							qzn = temp_z * nz;
							qzt = temp_z * z;
							qyn = temp_y * ny;
							qyt = temp_y * y;
							qt = temp_x * x + qyt + qzt;
							qn = (temp_x * nx + qyn + qzn) * q2_inv;

							total += compute_fq(s, qt, qn);

						case 2:
							shape_off = start_t + (unroll_off + 1) * CPU_T_PROP_SIZE_;
							s = shape_def[shape_off];
							nx = shape_def[shape_off + 1];
							ny = shape_def[shape_off + 2];
							nz = shape_def[shape_off + 3];
							x = shape_def[shape_off + 4];
							y = shape_def[shape_off + 5];
							z = shape_def[shape_off + 6];

							qzn = temp_z * nz;
							qzt = temp_z * z;
							qyn = temp_y * ny;
							qyt = temp_y * y;
							qt = temp_x * x + qyt + qzt;
							qn = (temp_x * nx + qyn + qzn) * q2_inv;

							total += compute_fq(s, qt, qn);

						case 1:
							shape_off = start_t + (unroll_off) * CPU_T_PROP_SIZE_;
							s = shape_def[shape_off];
							nx = shape_def[shape_off + 1];
							ny = shape_def[shape_off + 2];
							nz = shape_def[shape_off + 3];
							x = shape_def[shape_off + 4];
							y = shape_def[shape_off + 5];
							z = shape_def[shape_off + 6];

							qzn = temp_z * nz;
							qzt = temp_z * z;
							qyn = temp_y * ny;
							qyt = temp_y * y;
							qt = temp_x * x + qyt + qzt;
							qn = (temp_x * nx + qyn + qzn) * q2_inv;

							total += compute_fq(s, qt, qn);

						case 0:
							// nothing to do
							break;
					} // switch()

					ff[super_i] += total;
				} // for y
			} // for z
		} // pragma omp parallel
	} // NumericFormFactorC::form_factor_kernel_fused_nqx1_unroll4()

//	#endif
