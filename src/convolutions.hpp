/***
  *  $Id: convolutions.hpp 47 2012-08-23 21:05:16Z asarje $
  *
  *  Project:
  *
  *  File: convolutions.hpp
  *  Created: Jul 03, 2012
  *  Modified: Thu 23 Aug 2012 01:57:49 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _CONVOLUTIONS_HPP_
#define _CONVOLUTIONS_HPP_

#include <cmath>
#include <omp.h>

#include "utilities.hpp"

// use some other library like mkl or gsl or scsl or boost
// ...

namespace hig {

	// this class is stateless
	class Convolutions {
		public:
			enum shape_param_t {
				conv_null = 0,	/* undefined, or cleared */
				conv_full,		/* return full matrix */
				conv_same,		/* return matrix with size of input */
				conv_valid,		/* return matrix containing only values computed
								   without padded zeros */
				conv_error		/* some error occured */
			};

		private:
			Convolutions() { }
			Convolutions(const Convolutions&);
			Convolutions& operator=(const Convolutions&);

		public:
			static Convolutions& instance() {
				static Convolutions conv;
				return conv;
			} // instance()

			bool convolution_2d(shape_param_t type,
								unsigned int a_xsize, unsigned int a_ysize, const float_t *a,
								unsigned int b_xsize, unsigned int b_ysize, const float_t *b,
								unsigned int& c_xsize, unsigned int& c_ysize, float_t* &c) {

				switch(type) {
					case conv_valid:
						// brute force:
						// compute the size of output matrix
						// accordingly compute the new sizes of the input matrices - ????
						// compute amount of padding needed	- ????
						// copy input matrices to new 0-padded matrices - ????
						// execute loops to compute c
						// return
						return compute_conv_2d_valid(a_xsize, a_ysize, a, b_xsize, b_ysize, b,
														c_xsize, c_ysize, c);
						break;

					case conv_full:
						return compute_conv_2d_full(a_xsize, a_ysize, a, b_xsize, b_ysize, b,
														c_xsize, c_ysize, c);
						break;

					case conv_same:
						std::cerr << "uh-oh: convolution type 'same' is not "
									<< "yet implemented!" << std::endl
									<< "use convolution type 'valid' or 'full' instead" << std::endl;
						return false;

					default:
					case conv_null:
					case conv_error:
						std::cerr << "error: convolution type is not defined properly" << std::endl;
						return false;
				} // switch

				return true;
			} // convolution_2d()


/*			bool compute_conv_2d_valid1(unsigned int a_xsize, unsigned int a_ysize, const double *a,
									unsigned int b_xsize, unsigned int b_ysize, const double *b,
									unsigned int& c_xsize, unsigned int& c_ysize, double* &c) {

				unsigned int cmax_xsize = max(a_xsize + b_xsize - 1, a_xsize, b_xsize);
				unsigned int cmax_ysize = max(a_ysize + b_ysize - 1, a_ysize, b_ysize);
				c_xsize = max((int)a_xsize - max(0, (int)b_xsize - 1), 0);
				c_ysize = max((int)a_ysize - max(0, (int)b_ysize - 1), 0);

				if(c_xsize <= 0 || c_ysize <= 0) {
					c = NULL;
					return true;
				} // if

				c = new (std::nothrow) double[c_xsize * c_ysize];
				bool skip = false, enter = false;

				int cx = 0, cy = 0;
				for(int l2 = 0; l2 < cmax_ysize; ++ l2) {
					for(int l1 = 0; l1 < cmax_xsize; ++ l1) {
						int j1_min = max(0, l1 - (int)b_xsize);
						int j1_max = min(cx, (int)a_xsize - 1);
						int j2_min = max(0, l2 - (int)b_ysize);
						int j2_max = min(cy, (int)a_ysize - 1);
						double sum = 0.0;
						int ax = 0, ay = 0, bx = 0, by = 0;
						skip = false; enter = false;
						for(int j2 = j2_min; j2 < j2_max; ++ j2) {
							ay = j2;
							by = cy - ay;
							// skip sum for zero paddings
							if(ay < 0 || by < 0 || ay >= a_ysize || by >= b_ysize) {
								skip = true;
								continue;
							} // if
							for(int j1 = j1_min; j1 < j1_max; ++ j1) {
								ax = j1;
								bx = cx - ax;
								// skip sum for zero paddings
								if(ax < 0 || bx < 0 || ax >= a_xsize || bx >= b_xsize) {
									skip = true;
									continue;
								} // if
								//enter = true;
								std::cout << "hiya!" << std::endl;
								sum += a[a_xsize * ay + ax] * b[b_xsize * by + bx];
							} // for j1
							if(skip) break;
						} // for j2
						if((!skip) && enter) {
							c[c_xsize * cy + cx] = sum;
							if(cx < c_xsize - 1) {
								++ cx;
							} else {
								++ cy; cx = 0;
							} // if-else
							if(cy >= c_ysize) {
								std::cerr << "error: i am overflowing the c convolution matrix "
										<< "with cy = " << cy << ", cx = " << cx
										<< ", (" << cmax_ysize << ", " << cmax_xsize << ")" << std::endl;
								break; 			///////////////////////////////////// ..............
								return false;
							} // if
						} // if
					} // for l1
				} // for l2
				return true;
			} // compute_conv_valid()


			bool compute_conv_2d_valid2(unsigned int a_xsize, unsigned int a_ysize, const double *a,
										unsigned int b_xsize, unsigned int b_ysize, const double *b,
										unsigned int& c_xsize, unsigned int& c_ysize, double* &c) {
				unsigned int cmax_xsize = max(a_xsize + b_xsize - 1, a_xsize, b_xsize);
				unsigned int cmax_ysize = max(a_ysize + b_ysize - 1, a_ysize, b_ysize);
				c_xsize = max((int)a_xsize - max(0, (int)b_xsize - 1), 0);
				c_ysize = max((int)a_ysize - max(0, (int)b_ysize - 1), 0);

				if(c_xsize <= 0 || c_ysize <= 0) {
					c = NULL;
					return true;
				} // if

				c = new (std::nothrow) double[c_xsize * c_ysize];
				bool skip = false, enter = false;

				int cx = 0, cy = 0;	// indices for c (since it is 'valid')
				for(int n2 = 0; n2 < cmax_ysize; ++ n2) {
					for(int n1 = 0; n1 < cmax_xsize; ++ n1) {
						double sum = 0.0;
						skip = false;
						for(int aby = 0; aby < max(a_ysize, b_ysize); ++ aby) {
							int ay = aby;
							int by = cy - ay + 1;
							if(ay < 0 || ay >= a_ysize || by < 0 || by >= b_ysize) {
								skip = true;
//								break;
								continue;
							} // if
							for(int abx = 0; abx < max(a_xsize, b_xsize); ++ abx) {
								int ax = abx;
								int bx = cx - ax + 1;
								if(ax < 0 || ax >= a_xsize || bx < 0 || bx >= b_xsize) {
									skip = true;
//									break;
									continue;
								} // if
								sum += a[a_xsize * ay + ax] * b[b_xsize * by + bx];
							} // for x
							if(skip) break;
						} // for y
//						if(!skip) {
							int c_index = c_xsize * cy + cx;
							c[c_index] = sum;
							if(cx < c_xsize - 1) ++ cx;
							else { ++ cy; cx = 0; }
//						} // if
					} // for cx
				} // for cy

				return true;
			} // compute_conv_2d_valid()
*/

			bool compute_conv_2d_valid(unsigned int a_xsize, unsigned int a_ysize, const float_t *a,
										unsigned int b_xsize, unsigned int b_ysize, const float_t *b,
										unsigned int& c_xsize, unsigned int& c_ysize, float_t* &c) {
				unsigned int cmax_xsize = max(a_xsize + b_xsize - 1, a_xsize, b_xsize);
				unsigned int cmax_ysize = max(a_ysize + b_ysize - 1, a_ysize, b_ysize);
				c_xsize = max((int)a_xsize - max(0, (int)b_xsize - 1), 0);	// a_xsize - b_xsize + 1
				c_ysize = max((int)a_ysize - max(0, (int)b_ysize - 1), 0);	// a_ysize - b_ysize + 1

				if(c_xsize <= 0 || c_ysize <= 0) {
					c = NULL;
					return true;
				} // if

				// assuming a_xsize > c_xsize and a_ysize > c_ysize:
				int c_xmin = ((int)cmax_xsize - (int)c_xsize) >> 1;
				int c_xmax = c_xmin + (int)c_xsize;
				int c_ymin = ((int)cmax_ysize - (int)c_ysize) >> 1;
				int c_ymax = c_ymin + (int)c_ysize;

				//std::cout << "c_xsize: " << c_xsize << ", c_ysize: " << c_ysize << std::endl;
				//std::cout << "c_xmin: " << c_xmin << ", c_ymin: " << c_ymin << std::endl;
				//std::cout << "c_xmax: " << c_xmax << ", c_ymax: " << c_ymax << std::endl;

				c = new (std::nothrow) float_t[c_xsize * c_ysize];
				int cx = 0, cy = 0;
				//std::cout << "  ";
				int count = 0;
				for(int i = c_ymin; i < c_ymax; ++ i) {
					for(int j = c_xmin; j < c_xmax; ++ j) {
						float_t sum = 0.0;
						//for(int m = 0; m <= std::min(i, (int)a_ysize); ++ m) {
						//	for(int n = 0; n <= std::min(j, (int)a_xsize); ++ n) {
						for(int m = 0; m <= a_ysize; ++ m) {
							for(int n = 0; n <= a_xsize; ++ n) {
								int ax = n;
								int ay = m;
								int bx = j - n;
								int by = i - m;
								if(bx < b_xsize && by < b_ysize && bx >= 0 && by >= 0)
									sum += a[a_xsize * ay + ax] * b[b_xsize * by + bx];
							} // for n
						} // for m
						if(cy * c_xsize + cx >= c_xsize * c_ysize) {
							std::cerr << "error: overflowing convolution buffer: ("
									<< cy * c_xsize + cx << " / " << c_xsize * c_ysize << ") "
									<< "cx: " << cx << ", cy: " << cy
									<< ", c_xsize: " << c_xsize << ", c_ysize: " << c_ysize
									<< std::endl;
							//return false;
						} // if
						c[cy * c_xsize + cx] = sum;
						++ cx;
						if(cx >= c_xsize) { cx = 0; ++ cy; }
						if(cy > c_ysize) std::cout << "error: overflowing convolutions" << std::endl;

						++ count;
						if(count % 1000 == 0) std::cout << "." << std::flush;
					} // for j
				} // for i
				std::cout << std::endl;
				return true;
			} // compute_conv_2d_valid()


			bool compute_conv_2d_full(unsigned int a_xsize, unsigned int a_ysize, const float_t *a,
										unsigned int b_xsize, unsigned int b_ysize, const float_t *b,
										unsigned int& c_xsize, unsigned int& c_ysize, float_t* &c) {
				unsigned int cmax_xsize = max(a_xsize + b_xsize - 1, a_xsize, b_xsize);
				unsigned int cmax_ysize = max(a_ysize + b_ysize - 1, a_ysize, b_ysize);
				c_xsize = cmax_xsize;
				c_ysize = cmax_ysize;
//				c_xsize = max((int)a_xsize - max(0, (int)b_xsize - 1), 0);	// a_xsize - b_xsize + 1
//				c_ysize = max((int)a_ysize - max(0, (int)b_ysize - 1), 0);	// a_ysize - b_ysize + 1

				if(c_xsize <= 0 || c_ysize <= 0) {
					c = NULL;
					return true;
				} // if

				// assuming a_xsize > c_xsize and a_ysize > c_ysize:
//				int c_xmin = ((int)cmax_xsize - (int)c_xsize) >> 1;
//				int c_xmax = c_xmin + (int)c_xsize;
//				int c_ymin = ((int)cmax_ysize - (int)c_ysize) >> 1;
//				int c_ymax = c_ymin + (int)c_ysize;

				c = new (std::nothrow) float_t[c_xsize * c_ysize];
//				int cx = 0, cy = 0;
				for(int i = 0; i < c_ysize; ++ i) {
					for(int j = 0; j < c_xsize; ++ j) {
						float_t sum = 0.0;
						for(int m = 0; m < a_ysize; ++ m) {
							for(int n = 0; n < a_xsize; ++ n) {
								int ax = n;
								int ay = m;
								int bx = j - n;
								int by = i - m;
								if(bx >= 0 && by >= 0 && bx < b_xsize && by < b_ysize)
									sum += a[a_xsize * ay + ax] * b[b_xsize * by + bx];
							} // for n
						} // for m
						c[i * c_xsize + j] = sum;
//						++ cx;
//						if(cx >= c_xsize) { cx = 0; ++ cy; }
					} // for j
				} // for i
				return true;
			} // compute_conv_2d_valid()

	}; // class Convolutions
} // namespace

#endif /* _CONVOLUTIONS_HPP_ */
