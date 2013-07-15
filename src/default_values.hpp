/***
  *  Project:
  *
  *  File: default_values.hpp
  *  Created: Jul 12, 2013
  *  Modified: Fri 12 Jul 2013 02:32:19 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef __DEFAULT_VALUES_HPP__
#define __DEFAULT_VALUES_HPP__

#include "typedefs.hpp"

namespace hig {

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
