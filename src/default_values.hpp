/***
  *  Project:
  *
  *  File: default_values.hpp
  *  Created: Jul 12, 2013
  *  Modified: Mon 15 Jul 2013 03:11:25 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __DEFAULT_VALUES_HPP__
#define __DEFAULT_VALUES_HPP__

#include "typedefs.hpp"

namespace hig {

	/**
	 * DEFAULT VALUES FOR OPTIONAL PARAMETERS
	 */

	/**
	 * shape and shape parameter parameter defaults
	 */
	const float_t		DEFAULT_P1_					= 0.0;			// p1/mean
	const float_t		DEFAULT_P2_					= 0.0;			// p2/sd
	const int			DEFAULT_NVALUES_			= 1;			// num values
	const float_t		DEFAULT_ORIGIN_VEC_X_		= 0.0;
	const float_t		DEFAULT_ORIGIN_VEC_Y_		= 0.0;
	const float_t		DEFAULT_ORIGIN_VEC_Z_		= 0.0;
	const float_t		DEFAULT_Z_TILT_				= 0.0;
	const float_t		DEFAULT_XY_ROTATION_		= 0.0;

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
	const std::string	DEFAULT_ORIGIN_("bottom-left");				// not yet implemented
	const int			DEFAULT_TOTAL_PIXELS_Y_ 	= 1000;			// horizontal axis
	const int			DEFAULT_TOTAL_PIXELS_Z_ 	= 1000;			// vertical axis
	const float_t		DEFAULT_PIXEL_SIZE_ 		= 0.172;		// the detector pixel size
	const int			DEFAULT_SDD_ 				= 4000;			// not used currently
	const float_t		DEFAULT_DIRECT_BEAM_Y_ 		= 0.0;			// direct beam coordinates on detector
	const float_t		DEFAULT_DIRECT_BEAM_Z_ 		= 0.0;			// direct beam coordinates on detector


} // namespace hig

#endif // __DEFAULT_VALUES_HPP__
