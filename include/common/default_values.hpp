/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: default_values.hpp
 *  Created: Jul 12, 2013
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

#ifndef __DEFAULT_VALUES_HPP__
#define __DEFAULT_VALUES_HPP__

#include <common/typedefs.hpp>

namespace hig {

  /**
   * DEFAULT VALUES FOR OPTIONAL PARAMETERS
   */

  /**
   * shape and shape parameter parameter defaults
   */
  const real_t    DEFAULT_P1_           = 0.0;      // p1/mean
  const real_t    DEFAULT_P2_           = 0.0;      // p2/sd
  const int       DEFAULT_NVALUES_      = 1;        // num values
  const real_t    DEFAULT_ORIGIN_VEC_X_ = 0.0;
  const real_t    DEFAULT_ORIGIN_VEC_Y_ = 0.0;
  const real_t    DEFAULT_ORIGIN_VEC_Z_ = 0.0;
  const real_t    DEFAULT_Z_TILT_       = 0.0;
  const real_t    DEFAULT_XY_ROTATION_  = 0.0;

  // TODO: ...

  /**
   * lattice parameter defaults
   */
  // ...

  /**
   * grain parameter defaults
   */
  // ...

  /**
   * ensemble parameter defaults
   */
  // ...

  /**
   * structure parameter defaults
   */
  // ...

  /**
   * layer parameter defaults
   */
  // ...

  /**
   * scattering parameter defaults
   */
  // ...

  /**
   * compute parameter defaults
   */
  // ...

  /**
   * detector parameters
   */
  const std::string DEFAULT_ORIGIN_("bottom-left");     // not yet implemented
  const int         DEFAULT_TOTAL_PIXELS_Y_ = 5000;     // horizontal axis
  const int         DEFAULT_TOTAL_PIXELS_Z_ = 5000;     // vertical axis
  const real_t      DEFAULT_PIXEL_SIZE_     = 0.172;    // the detector pixel size
  const int         DEFAULT_SDD_            = 4000;     // not used currently
  const real_t      DEFAULT_DIRECT_BEAM_Y_  = 0.0;      // direct beam coordinates on detector
  const real_t      DEFAULT_DIRECT_BEAM_Z_  = 0.0;      // direct beam coordinates on detector


} // namespace hig

#endif // __DEFAULT_VALUES_HPP__
