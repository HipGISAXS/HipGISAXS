/**
 * $Id: reduction.hpp 38 2012-08-09 23:01:20Z asarje $
 */

#ifndef _MISC_H_
#define _MISC_H_

#ifndef GPUR
void reduction_kernel(unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned long int,
						unsigned int, unsigned int, unsigned int,
						unsigned int,
						unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuComplex*, cuComplex*);

void reduction_kernel(unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned long int,
						unsigned int, unsigned int, unsigned int,
						unsigned int,
						unsigned int, unsigned int, unsigned int,
						unsigned int, unsigned int, unsigned int, unsigned int,
						cuDoubleComplex*, cuDoubleComplex*);
#else
__global__ void reduction_kernel(cuComplex*,
									unsigned int, unsigned int, unsigned int,
									unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									cuComplex*);

__global__ void reduction_kernel(cuDoubleComplex*,
									unsigned int, unsigned int, unsigned int,
									unsigned int,
									unsigned int, unsigned int, unsigned int, unsigned int,
									cuDoubleComplex*);
#endif // GPUR

#endif // _MISC_H_
