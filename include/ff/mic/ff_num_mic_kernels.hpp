/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_mic_kernels.hpp
 *  Created: Apr 22, 2013
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

	#ifdef FF_NUM_MIC_SWAP
	
	// special cases of the kernel when nqx == 1
	
	// for compiler vectorization
	__attribute__((target(mic:0)))
	void NumericFormFactorM::form_factor_kernel_loopswap_nqx1(real_t* qx, real_t* qy, scomplex_t* qz_flat,
					real_t* shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					real_t* rot,
					scomplex_t* ff_buffer) {

		const unsigned int start_z = b_nqz * ib_z;
		const unsigned int start_y = b_nqy * ib_y;
		const unsigned int start_t = b_num_triangles * ib_t; // (due to data re-org)

		// FIXME: for now assuming float only (4 bytes per float => 16 floats)
		const unsigned int vec_size = 16;
		const unsigned int shape_padding = (vec_size - (num_triangles & 15)) & 15;
		const unsigned int padded_num_triangles = num_triangles + shape_padding;

		#pragma omp parallel
		{
			#pragma omp for collapse(2) nowait
			for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
				for(int i_y = 0; i_y < curr_nqy; ++ i_y) {

					scomplex_t temp_qz = qz_flat[start_z + i_z];
					real_t temp_qy = qy[start_y + i_y];
					real_t temp_qx = qx[0];

					scomplex_t temp_x = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
					scomplex_t temp_y = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
					scomplex_t temp_z = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;

					scomplex_t qz2 = temp_z * temp_z;
					scomplex_t qy2 = temp_y * temp_y;
					scomplex_t qx2 = temp_x * temp_x;

					scomplex_t q2 = qx2 + qy2 + qz2;
					scomplex_t q2_inv = (real_t) 1.0 / q2;

					scomplex_t total = make_sC((real_t) 0.0, (real_t) 0.0);

					// TODO: do blocking for cache ... ?
					for(int i_t = 0; i_t < curr_num_triangles; ++ i_t) {
						unsigned int shape_off = start_t + i_t;
						real_t s = shape_def[shape_off];
						shape_off += padded_num_triangles;
						real_t nx = shape_def[shape_off];
						shape_off += padded_num_triangles;
						real_t ny = shape_def[shape_off];
						shape_off += padded_num_triangles;
						real_t nz = shape_def[shape_off];
						shape_off += padded_num_triangles;
						real_t x = shape_def[shape_off];
						shape_off += padded_num_triangles;
						real_t y = shape_def[shape_off];
						shape_off += padded_num_triangles;
						real_t z = shape_def[shape_off];

						scomplex_t qzn = temp_z * nz;
						scomplex_t qzt = temp_z * z;
						scomplex_t qyn = temp_y * ny;
						scomplex_t qyt = temp_y * y;
						scomplex_t qxn = temp_x * nx;
						scomplex_t qxt = temp_x * x;
						scomplex_t qt = qxt + qyt + qzt;
						scomplex_t qn = (qxn + qyn + qzn) * q2_inv;

						total = total + compute_fq_nqx1(s, qt, qn);
					} // for t

					unsigned int i_ff = curr_nqy * i_z + i_y;
					ff_buffer[i_ff] = make_sC(total.y, - total.x);
				} // for y
			} // for z
		} // pragma omp parallel
	} // NumericFormFactorM::form_factor_kernel_loopswap_nqx1()
	
	// single precision, not vectorized
	__attribute__((target(mic:0)))
	inline float2_t NumericFormFactorM::compute_fq_nqx1(float s, float2_t qt, float2_t qn) {
		float2_t v1 = qn * make_sC(cosf(qt.x), sinf(qt.x));
		float v2 = s * exp(qt.y);
		return v1 * v2;
	} // NumericFormFactorM::compute_fq()
	

	#ifdef __MIC__

	// with manual vectorization:
	__attribute__((target(mic:0)))
	void NumericFormFactorM::form_factor_kernel_loopswap_vec_nqx1(
					real_t* qx, real_t* qy, scomplex_t* qz_flat,
					real_t* shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					real_t* rot,
					scomplex_t* ff_buffer) {

		const unsigned int start_z = b_nqz * ib_z;
		const unsigned int start_y = b_nqy * ib_y;
		const unsigned int start_t = b_num_triangles * ib_t; // (due to data re-org)

		// FIXME: for now assuming float only (4 bytes per float => 16 floats)
		const unsigned int vec_size = 16;
		const unsigned int shape_padding = (vec_size - (num_triangles & 15)) & 15;
		const unsigned int padded_num_triangles = num_triangles + shape_padding;

		// shape padding guarantees that padded_num_triangles is a multiple of 16
		// FIXME: assuming that b_num_triangles is a multiple of 16 as well ...

/*		#pragma omp parallel
		{
			#pragma omp for collapse(2) nowait
			for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
				for(int i_y = 0; i_y < curr_nqy; ++ i_y) {

					// //////////////////////////////////////////////////////////
					// mic intrinsic wrappers naming convention (floating only)
					// mic_xxx_abc	=> a = r|c, b = p|s, c = s|d
					// mic_xxx_abcd	=> a = r|c, b = r|c, c = p|s, d = s|d
					// r = real, 				c = complex,
					// p = packed (vector),		s = scalar,
					// s = single-precision,	d = double-precision
					// //////////////////////////////////////////////////////////

					mic_m512c_t temp_z = mic_set1_cps(qz_flat[start_z + i_z]);
					mic_m512_t temp_y = mic_set1_rps(qy[start_y + i_y]);
					mic_m512_t temp_x = mic_set1_rps(qx[0]);

					mic_m512c_t qz2 = mic_mul_ccps(temp_z, temp_z);
					mic_m512_t qy2 = mic_mul_rrps(temp_y, temp_y);
					mic_m512_t qx2 = mic_mul_rrps(temp_x, temp_x);

					mic_m512c_t q2 = mic_add_rcps(mic_add_rrps(qx2, qy2), qz2);
					mic_m512c_t q2_inv = mic_rcp_cps(q2);

					mic_m512c_t total = mic_setzero_cps();

					// TODO: do blocking for cache ... ?
					// TODO: do prefetching ... ?
					for(int i_t = 0; i_t < curr_num_triangles; i_t += vec_size) {
						// load 16 floats at a time:
						unsigned int shape_off = start_t + i_t;
						mic_m512_t s = mic_load_rps(& shape_def[shape_off]);
						shape_off += padded_num_triangles;
						mic_m512_t nx = mic_load_rps(& shape_def[shape_off]);
						shape_off += padded_num_triangles;
						mic_m512_t ny = mic_load_rps(& shape_def[shape_off]);
						shape_off += padded_num_triangles;
						mic_m512_t nz = mic_load_rps(& shape_def[shape_off]);
						shape_off += padded_num_triangles;
						mic_m512_t x = mic_load_rps(& shape_def[shape_off]);
						shape_off += padded_num_triangles;
						mic_m512_t y = mic_load_rps(& shape_def[shape_off]);
						shape_off += padded_num_triangles;
						mic_m512_t z = mic_load_rps(& shape_def[shape_off]);

						mic_m512c_t qzn = mic_mul_crps(temp_z, nz);
						mic_m512c_t qzt = mic_mul_crps(temp_z, z);
						mic_m512_t qyn = mic_mul_rrps(temp_y, ny);
						mic_m512_t qyt = mic_mul_rrps(temp_y, y);
						mic_m512_t qxn = mic_mul_rrps(temp_x, nx);
						mic_m512_t qxt = mic_mul_rrps(temp_x, x);
						mic_m512c_t qt = mic_add_rcps(mic_add_rrps(qxt, qyt), qzt);
						mic_m512c_t temp_qn = mic_add_rcps(mic_add_rrps(qxn, qyn), qzn);
						mic_m512c_t qn = mic_mul_ccps(temp_qn, q2_inv);
						mic_m512c_t fq = compute_fq_vec_nqx1(s, qt, qn);

						total = mic_add_ccps(total, fq);
					} // for t

					unsigned int i_ff = curr_nqy * i_z + i_y;	// location in current buffer

					scomplex_t ff = mic_reduce_add_cps(total);
					ff_buffer[i_ff] = make_sC(ff.y, - ff.x);
				} // for y
			} // for z
		} // pragma omp parallel
*/
		#pragma omp parallel
		{
			#pragma omp for collapse(2) nowait
			for(int i_z = 0; i_z < curr_nqz; ++ i_z) {
				for(int i_y = 0; i_y < curr_nqy; ++ i_y) {

					// //////////////////////////////////////////////////////////
					// mic intrinsic wrappers naming convention (floating only)
					// mic_xxx_abc	=> a = r|c, b = p|s, c = s|d
					// mic_xxx_abcd	=> a = r|c, b = r|c, c = p|s, d = s|d
					// r = real, 				c = complex,
					// p = packed (vector),		s = scalar,
					// s = single-precision,	d = double-precision
					// //////////////////////////////////////////////////////////

					// TODO: optimize rot computation ... 
					scomplex_t temp_qz = qz_flat[start_z + i_z];
					real_t temp_qy = qy[start_y + i_y];
					real_t temp_qx = qx[0];
					scomplex_t temp_sx = rot[0] * temp_qx + rot[1] * temp_qy + rot[2] * temp_qz;
					scomplex_t temp_sy = rot[3] * temp_qx + rot[4] * temp_qy + rot[5] * temp_qz;
					scomplex_t temp_sz = rot[6] * temp_qx + rot[7] * temp_qy + rot[8] * temp_qz;

					mic_m512c_t temp_z = mic_set1_cps(temp_sz);
					mic_m512c_t temp_y = mic_set1_cps(temp_sy);
					mic_m512c_t temp_x = mic_set1_cps(temp_sx);

					mic_m512c_t qz2 = mic_mul_ccps(temp_z, temp_z);
					mic_m512c_t qy2 = mic_mul_ccps(temp_y, temp_y);
					mic_m512c_t qx2 = mic_mul_ccps(temp_x, temp_x);

					mic_m512c_t q2 = mic_add_ccps(mic_add_ccps(qx2, qy2), qz2);
					mic_m512c_t q2_inv = mic_rcp_cps(q2);

					mic_m512c_t total1 = mic_setzero_cps();
					mic_m512c_t total2 = mic_setzero_cps();

					// TODO: do blocking for cache ... ?
					// TODO: do prefetching ... ?
					for(int i_t = 0; i_t < curr_num_triangles; i_t += vec_size * 2) {
						// load 16 floats at a time:
						unsigned int shape_off = start_t + i_t;
						mic_m512_t s1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t s2 = mic_load_rps(& shape_def[shape_off + vec_size]);
						shape_off += padded_num_triangles;
						mic_m512_t nx1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t nx2 = mic_load_rps(& shape_def[shape_off + vec_size]);
						shape_off += padded_num_triangles;
						mic_m512_t ny1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t ny2 = mic_load_rps(& shape_def[shape_off + vec_size]);
						shape_off += padded_num_triangles;
						mic_m512_t nz1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t nz2 = mic_load_rps(& shape_def[shape_off + vec_size]);
						shape_off += padded_num_triangles;
						mic_m512_t x1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t x2 = mic_load_rps(& shape_def[shape_off + vec_size]);
						shape_off += padded_num_triangles;
						mic_m512_t y1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t y2 = mic_load_rps(& shape_def[shape_off + vec_size]);
						shape_off += padded_num_triangles;
						mic_m512_t z1 = mic_load_rps(& shape_def[shape_off]);
						mic_m512_t z2 = mic_load_rps(& shape_def[shape_off + vec_size]);

						mic_m512c_t qzn1 = mic_mul_crps(temp_z, nz1);
						mic_m512c_t qzn2 = mic_mul_crps(temp_z, nz2);
						mic_m512c_t qzt1 = mic_mul_crps(temp_z, z1);
						mic_m512c_t qzt2 = mic_mul_crps(temp_z, z2);
						mic_m512c_t qyn1 = mic_mul_crps(temp_y, ny1);
						mic_m512c_t qyn2 = mic_mul_crps(temp_y, ny2);
						mic_m512c_t qyt1 = mic_mul_crps(temp_y, y1);
						mic_m512c_t qyt2 = mic_mul_crps(temp_y, y2);
						mic_m512c_t qxn1 = mic_mul_crps(temp_x, nx1);
						mic_m512c_t qxn2 = mic_mul_crps(temp_x, nx2);
						mic_m512c_t qxt1 = mic_mul_crps(temp_x, x1);
						mic_m512c_t qxt2 = mic_mul_crps(temp_x, x2);
						mic_m512c_t qt1 = mic_add_ccps(mic_add_ccps(qxt1, qyt1), qzt1);
						mic_m512c_t qt2 = mic_add_ccps(mic_add_ccps(qxt2, qyt2), qzt2);
						mic_m512c_t temp_qn1 = mic_add_ccps(mic_add_ccps(qxn1, qyn1), qzn1);
						mic_m512c_t temp_qn2 = mic_add_ccps(mic_add_ccps(qxn2, qyn2), qzn2);
						mic_m512c_t qn1 = mic_mul_ccps(temp_qn1, q2_inv);
						mic_m512c_t qn2 = mic_mul_ccps(temp_qn2, q2_inv);

						mic_m512c_t temp1;
						mic_m512c_t temp2;
						mic_sincos_rps(qt1.xvec, & temp1.yvec, & temp1.xvec);
						mic_sincos_rps(qt2.xvec, & temp2.yvec, & temp2.xvec);
						mic_m512c_t v11 = mic_mul_ccps(qn1, temp1);
						mic_m512c_t v12 = mic_mul_ccps(qn2, temp2);
						mic_m512_t v21 = mic_mul_rrps(s1, mic_exp_rps(qt1.yvec));
						mic_m512_t v22 = mic_mul_rrps(s2, mic_exp_rps(qt2.yvec));
						mic_m512c_t fq1 =  mic_mul_crps(v11, v21);
						mic_m512c_t fq2 =  mic_mul_crps(v12, v22);

						total1 = mic_add_ccps(total1, fq1);
						total2 = mic_add_ccps(total2, fq2);
					} // for t

					unsigned int i_ff = curr_nqy * i_z + i_y;	// location in current buffer
					scomplex_t ff = mic_reduce_add_cps(total1);
					ff = ff + mic_reduce_add_cps(total2);
					ff_buffer[i_ff] = make_sC(ff.y, - ff.x);
				} // for y
			} // for z
		} // pragma omp parallel
	} // NumericFormFactorM::form_factor_kernel_loopswap_vec_nqx1()
	

	// single precision, vectorized
	__attribute__((target(mic:0)))
	inline mic_m512c_t NumericFormFactorM::compute_fq_vec_nqx1(mic_m512_t s, mic_m512c_t qt, mic_m512c_t qn) {
		mic_m512c_t temp;
		//temp.xvec = mic_cos_rps(qt.xvec);
		//temp.yvec = mic_sin_rps(qt.xvec);
		mic_sincos_rps(qt.xvec, & temp.yvec, & temp.xvec);
		mic_m512c_t v1 = mic_mul_ccps(qn, temp);
		mic_m512_t v2 = mic_mul_rrps(s, mic_exp_rps(qt.yvec));
		return mic_mul_crps(v1, v2);
	} // NumericFormFactorM::compute_fq()

	#endif // __MIC__
	#endif // FF_NUM_MIC_SWAP
