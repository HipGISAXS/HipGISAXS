/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: tokens.hpp
 *  Created: Jun 05, 2012
 *  Modified: Tue 16 Jul 2013 11:52:26 AM PDT
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

#ifndef _TOKENS_HPP_
#define _TOKENS_HPP_

namespace hig {

	/**
	 * fundamental datatypes used
	 */

	enum ValueType {
		null_value = 0,					/* an empty/non-existant/null value */
		int_value,						/* integral number */
		real_value,						/* real number */
		string_value,					/* quoted strings */
		vector2_value,					/* a pair of two real numbers */
		vector3_value					/* a tuple of three real numbers */
	}; // enum ValueType


	/**
	 * set of valid tokens for an input component
	 */

	enum TokenType {
		/* sanity */
		null_token = 0,					/* when there is nothing yet */
		error_token,					/* there is something that is not supposed to be there */

		/* fundamental raw tokens: single character (they make bigger stuff) */
		character_token,				/* a single alphabet character */
		digit_token,					/* a single digit [0-9\-] */
		negative_token,					/* the negative sign '-' (for numbers only) */
		white_space_token,				/* a white space character: ' ' '\n' '\t' '\r' '\v' etc. */
		object_begin_token,				/* '{' */
		object_end_token,				/* '}' */
		array_begin_token,				/* '[' */
		array_end_token,				/* ']' */
		string_begin_end_token,			/* '"' */
		assignment_token,				/* '=' */
		separator_token,				/* ',' */
		comment_token,					/* sign '#' for a comment start */

		/**
		 * hig tokens: valid keywords and values
		 */

		hipgisaxs_token,				/* the main keyword representing the input object */

		/* value datatypes */
		number_token,					/* a number: integeral or real */
		string_token,					/* a quoted string */

		/* miscellaneous, used in multiple places */
		key_token,						/* key value for various: shape, layer, structure, etc. */
		min_token,						/* min value for various: shape param, .. */
		max_token,						/* max value for various: shape param, .. */
		step_token,						/* step from min to max for various */
		rot_token,						/* rotation value for various */
		type_token,						/* various types */
		stat_token,						/* statistic type: gaussian, uniform, random, etc. */

		/* refractive index components */
		refindex_token,					/* refractive index of various */
		refindex_delta_token,			/* delta value for a refractive index */
		refindex_beta_token,			/* beta value for a refractive index */

		/* shape name and parameters*/
		shape_token,					/* a shape */
		shape_name_token,				/* name of a shape: box, truncpyr, .., or [filename] */
		shape_originvec_token,			/* reference origin point for a shape */
		shape_ztilt_token,				/* angle of shape tilt off the z-axis */
		shape_xyrot_token,				/* rotation of a shape in xy-plane */
		shape_param_token,				/* parameter of a shape */
		shape_param_p1_token,			/* mean for a shape's statistic */
		shape_param_p2_token,			/* deviation for a shape's statistic */
		shape_param_nvalues_token,		/* ... */

		/* layer parameters */
		layer_token,
		layer_order_token,
		layer_thickness_token,

		/* structure components */
		struct_token,
		struct_grain_token,
		struct_grain_skey_token,
		struct_grain_lkey_token,
		struct_grain_lattice_token,
		struct_grain_lattice_a_token,
		struct_grain_lattice_b_token,
		struct_grain_lattice_c_token,
		struct_grain_lattice_hkl_token,
		struct_grain_lattice_abangle_token,
		struct_grain_lattice_caratio_token,
		struct_grain_transvec_token,
		struct_grain_scaling_token,
		struct_grain_repetition_token,
		struct_ensemble_token,
		struct_ensemble_spacing_token,
		struct_ensemble_maxgrains_token,
		struct_ensemble_distribution_token,
		struct_ensemble_orient_token,
		struct_ensemble_orient_stat_token,
		struct_ensemble_orient_rot1_token,
		struct_ensemble_orient_rot2_token,
		struct_ensemble_orient_rot3_token,
		struct_ensemble_orient_rot_axis_token,
		struct_ensemble_orient_rot_angles_token,

		/* compute parameters */
		compute_token,
		compute_path_token,
		compute_runname_token,
		compute_method_token,
		compute_resolution_token,
		compute_outregion_token,
		compute_outregion_minpoint_token,
		compute_outregion_maxpoint_token,
		compute_nslices_token,
		compute_structcorr_token,			/* defined grain/ensemble correlations */

		/* experiment instrumentation - scatter and detector */
		instrument_token,
		instrument_scatter_token,
		instrument_scatter_expt_token,
		instrument_scatter_alphai_token,
		instrument_scatter_inplanerot_token,
		instrument_scatter_tilt_token,
		instrument_scatter_photon_token,
		instrument_scatter_photon_value_token,
		instrument_scatter_photon_unit_token,
		instrument_scatter_polarize_token,
		instrument_scatter_coherence_token,
		instrument_scatter_spotarea_token,
		instrument_scatter_smearing_token,
		instrument_detector_token,
		instrument_detector_origin_token,
		instrument_detector_totpix_token,
		instrument_detector_sdd_token,
		instrument_detector_pixsize_token,
		instrument_detector_dirbeam_token
	}; // enum TokenType


	/**
	 * token class storing details of a token
	 */

	class Token {
		public:
			Token() : type_(null_token), svalue_(""), dvalue_(0.0) { }
			~Token() { }

			TokenType type_;		/* token type */
			std::string svalue_;	/* token's actual string value - if non-numeric */
			float_t dvalue_;		/* token's actual numeric value */
	}; // class Token

} // namespace hig

#endif /* _TOKENS_HPP_ */
