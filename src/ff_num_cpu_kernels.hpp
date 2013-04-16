/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: ff_num_cpu.hpp
 *  Created: Nov 05, 2011
 *  Modified: Wed 03 Apr 2013 11:18:54 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */
	

	/**
	 * the main Form Factor kernel with fused reduction function - for one hyperblock.
	 */
	void NumericFormFactorC::form_factor_kernel_fused_unroll4(float_t* qx, float_t* qy, complex_t* qz,
					float_vec_t& shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
					complex_t* ff) {
		if(ff == NULL || qx == NULL || qy == NULL || qz == NULL) return;
	
		unsigned long int xy_size = (unsigned long int) curr_nqx * curr_nqy;
		unsigned int start_z = b_nqz * ib_z;
		unsigned int start_y = b_nqy * ib_y;
		unsigned int start_x = b_nqx * ib_x;
		unsigned int start_t = b_num_triangles * ib_t * 7;

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
						complex_t temp_z = qz[start_z + i_z];
						float_t temp_y = qy[start_y + i_y];
						float_t temp_x = qx[start_x + i_x];
						complex_t qz2 = temp_z * temp_z;
						float_t qy2 = temp_y * temp_y;
						complex_t q2 = temp_x * temp_x + qy2 + qz2;

						for(int i_t = 0; i_t < curr_num_triangles; i_t += 4) {
							unsigned int shape_off1 = start_t + i_t * 7;
							float_t s1 = shape_def[shape_off1];
							float_t nx1 = shape_def[shape_off1 + 1];
							float_t ny1 = shape_def[shape_off1 + 2];
							float_t nz1 = shape_def[shape_off1 + 3];
							float_t x1 = shape_def[shape_off1 + 4];
							float_t y1 = shape_def[shape_off1 + 5];
							float_t z1 = shape_def[shape_off1 + 6];
							unsigned int shape_off2 = start_t + (i_t + 1) * 7;
							float_t s2 = shape_def[shape_off2];
							float_t nx2 = shape_def[shape_off2 + 1];
							float_t ny2 = shape_def[shape_off2 + 2];
							float_t nz2 = shape_def[shape_off2 + 3];
							float_t x2 = shape_def[shape_off2 + 4];
							float_t y2 = shape_def[shape_off2 + 5];
							float_t z2 = shape_def[shape_off2 + 6];
							unsigned int shape_off3 = start_t + (i_t + 2) * 7;
							float_t s3 = shape_def[shape_off3];
							float_t nx3 = shape_def[shape_off3 + 1];
							float_t ny3 = shape_def[shape_off3 + 2];
							float_t nz3 = shape_def[shape_off3 + 3];
							float_t x3 = shape_def[shape_off3 + 4];
							float_t y3 = shape_def[shape_off3 + 5];
							float_t z3 = shape_def[shape_off3 + 6];
							unsigned int shape_off4 = start_t + (i_t + 3) * 7;
							float_t s4 = shape_def[shape_off4];
							float_t nx4 = shape_def[shape_off4 + 1];
							float_t ny4 = shape_def[shape_off4 + 2];
							float_t nz4 = shape_def[shape_off4 + 3];
							float_t x4 = shape_def[shape_off4 + 4];
							float_t y4 = shape_def[shape_off4 + 5];
							float_t z4 = shape_def[shape_off4 + 6];
	
							complex_t qzn1 = temp_z * nz1;
							complex_t qzt1 = temp_z * z1;
							complex_t qzn2 = temp_z * nz2;
							complex_t qzt2 = temp_z * z2;
							complex_t qzn3 = temp_z * nz3;
							complex_t qzt3 = temp_z * z3;
							complex_t qzn4 = temp_z * nz4;
							complex_t qzt4 = temp_z * z4;
	
							float_t qyn1 = temp_y * ny1;
							float_t qyt1 = temp_y * y1;
							float_t qyn2 = temp_y * ny2;
							float_t qyt2 = temp_y * y2;
							float_t qyn3 = temp_y * ny3;
							float_t qyt3 = temp_y * y3;
							float_t qyn4 = temp_y * ny4;
							float_t qyt4 = temp_y * y4;
	
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
						float_t s, nx, ny, nz, x, y, z;
						complex_t qzn, qzt;
						complex_t qyn, qyt;
						complex_t qn, qt;
	
						switch(unroll_rem) {
							case 3:
								shape_off = start_t + (unroll_off + 2) * 7;
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
								shape_off = start_t + (unroll_off + 1) * 7;
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
								shape_off = start_t + (unroll_off) * 7;
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
	} // NumericFormFactorC::form_factor_kernel_fused_unroll()
	

	/**
	 * special case of nqx == 1 of Form Factor kernel with fused reduction
	 */
	void NumericFormFactorC::form_factor_kernel_fused_nqx1_unroll4(float_t* qx, float_t* qy, complex_t* qz,
					float_vec_t& shape_def,
					unsigned int curr_nqx, unsigned int curr_nqy, unsigned int curr_nqz,
					unsigned int curr_num_triangles,
					unsigned int b_nqx, unsigned int b_nqy, unsigned int b_nqz, unsigned int b_num_triangles,
					unsigned int nqx, unsigned int nqy, unsigned int nqz, unsigned int num_triangles,
					unsigned int ib_x, unsigned int ib_y, unsigned int ib_z, unsigned int ib_t,
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
		unsigned int start_t = b_num_triangles * ib_t * 7;

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
					complex_t temp_z = qz[start_z + i_z];
					float_t temp_y = qy[start_y + i_y];
					float_t temp_x = qx[0];
					complex_t qz2 = temp_z * temp_z;
					float_t qy2 = temp_y * temp_y;
					complex_t q2 = temp_x * temp_x + qy2 + qz2;
					complex_t q2_inv = ((float_t) 1.0) / q2;

					complex_t total(0.0, 0.0);

					for(int i_t = 0; i_t < curr_num_triangles; i_t += 4) {
						unsigned int shape_off1 = start_t + i_t * 7;
						float_t s1 = shape_def[shape_off1];
						float_t nx1 = shape_def[shape_off1 + 1];
						float_t ny1 = shape_def[shape_off1 + 2];
						float_t nz1 = shape_def[shape_off1 + 3];
						float_t x1 = shape_def[shape_off1 + 4];
						float_t y1 = shape_def[shape_off1 + 5];
						float_t z1 = shape_def[shape_off1 + 6];
						unsigned int shape_off2 = start_t + (i_t + 1) * 7;
						float_t s2 = shape_def[shape_off2];
						float_t nx2 = shape_def[shape_off2 + 1];
						float_t ny2 = shape_def[shape_off2 + 2];
						float_t nz2 = shape_def[shape_off2 + 3];
						float_t x2 = shape_def[shape_off2 + 4];
						float_t y2 = shape_def[shape_off2 + 5];
						float_t z2 = shape_def[shape_off2 + 6];
						unsigned int shape_off3 = start_t + (i_t + 2) * 7;
						float_t s3 = shape_def[shape_off3];
						float_t nx3 = shape_def[shape_off3 + 1];
						float_t ny3 = shape_def[shape_off3 + 2];
						float_t nz3 = shape_def[shape_off3 + 3];
						float_t x3 = shape_def[shape_off3 + 4];
						float_t y3 = shape_def[shape_off3 + 5];
						float_t z3 = shape_def[shape_off3 + 6];
						unsigned int shape_off4 = start_t + (i_t + 3) * 7;
						float_t s4 = shape_def[shape_off4];
						float_t nx4 = shape_def[shape_off4 + 1];
						float_t ny4 = shape_def[shape_off4 + 2];
						float_t nz4 = shape_def[shape_off4 + 3];
						float_t x4 = shape_def[shape_off4 + 4];
						float_t y4 = shape_def[shape_off4 + 5];
						float_t z4 = shape_def[shape_off4 + 6];

						complex_t qzn1 = temp_z * nz1;
						complex_t qzt1 = temp_z * z1;
						complex_t qzn2 = temp_z * nz2;
						complex_t qzt2 = temp_z * z2;
						complex_t qzn3 = temp_z * nz3;
						complex_t qzt3 = temp_z * z3;
						complex_t qzn4 = temp_z * nz4;
						complex_t qzt4 = temp_z * z4;

						float_t qyn1 = temp_y * ny1;
						float_t qyt1 = temp_y * y1;
						float_t qyn2 = temp_y * ny2;
						float_t qyt2 = temp_y * y2;
						float_t qyn3 = temp_y * ny3;
						float_t qyt3 = temp_y * y3;
						float_t qyn4 = temp_y * ny4;
						float_t qyt4 = temp_y * y4;

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
					float_t s, nx, ny, nz, x, y, z;
					complex_t qzn, qzt;
					complex_t qyn, qyt;
					complex_t qn, qt;

					switch(unroll_rem) {
						case 3:
							shape_off = start_t + (unroll_off + 2) * 7;
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
							shape_off = start_t + (unroll_off + 1) * 7;
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
							shape_off = start_t + (unroll_off) * 7;
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
