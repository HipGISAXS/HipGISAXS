/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: reduction.hpp
 *  Created: Aug 25, 2012
 *  Modified: Thu 26 Sep 2013 08:18:58 AM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Elaine Chan <erchan@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _MISC_H_
#define _MISC_H_

namespace hig {

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

} // namespace hig

#endif // _MISC_H_
