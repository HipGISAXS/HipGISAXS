/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: reduction.hpp
 *  Created: Aug 25, 2012
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

#ifndef __HIG_REDUCTION_HPP__
#define __HIG_REDUCTION_HPP__

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

#endif // __HIG_REDUCTION_HPP__
