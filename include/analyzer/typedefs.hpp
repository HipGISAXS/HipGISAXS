/**
 *  Project: AnalyzeHipGISAXS (High-Performance GISAXS Data Analysis)
 *
 *  File: typdefs.hpp
 *  Created: Dec 26, 2013
 *  Modified: Mon 27 Jan 2014 08:22:24 AM PST
 *  Description: Defines the commom data types used.
 *
 *  Author: Slim Chourou <stchourou@lbl.gov>
 *  Developers: Slim Chourou <stchourou@lbl.gov>
 *              Abhinav Sarje <asarje@lbl.gov>
 *              Alexander Hexemer <ahexemer@lbl.gov>
 *              Xiaoye Li <xsli@lbl.gov>
 *
 *  Licensing: The AnalyzeHipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#ifndef _TYPEDEFS_HPP_
#define _TYPEFEDS_HPP_

#include <vector>
#include <complex>
#include <string>

namespace hig{

  /* strings */
  typedef std::string string_t;
  typedef std::vector<std::string> string_vec_t;
  typedef std::vector<std::string>::const_iterator   string_vec_cit;

  /* numbers  */

  /* These types are taken form HipGISAXS. Will remove after merge */
#ifdef DOUBLEP                                          // double precision
  typedef double                                          float_t;
#ifdef USE_GPU
  typedef cuDoubleComplex                 cucomplex_t;
#endif
#if defined USE_MIC
  typedef double2_t                               scomplex_t;
#endif
#else                                                           // single precision
  typedef float                                           float_t;
#ifdef USE_GPU
  typedef cuFloatComplex                  cucomplex_t;
#endif
#if defined USE_MIC
  typedef float2_t                                scomplex_t;
#endif
#endif

  typedef std::complex<float_t>                   complex_t;
  typedef std::vector<float_t>                    float_vec_t;
  typedef std::vector<complex_t>                  complex_vec_t;
  typedef std::vector<float_vec_t>                float_mat_t;
  /***************************************/

} /* namespace hig */

#endif /* TYPEDEFS_HPP_ */
