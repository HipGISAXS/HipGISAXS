/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: enums.hpp
 *  Created: Jun 11, 2012
 *  Modified: Tue 16 Jul 2013 11:49:46 AM PDT
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

#ifndef _ENUMS_HPP_
#define _ENUMS_HPP_

namespace hig {

	// these are values of string_token:

	enum ShapeName {
		shape_null,
		shape_error,
		shape_box,				/* box */
		shape_cylinder,			/* cylinder */
		shape_sphere,			/* sphere */
		shape_truncpyr,			/* truncated pyramid */
		shape_trunccone,		/* truncated cone */
		shape_prism3,			/* 3-fold prism */
		shape_prism6,			/* 6-fold prism */
		shape_prism3x,			/* triangular grating in x direction */
		shape_sawtooth,			/* sawtooth (prism along x) (same as above) */
		shape_random_cylinders,
		shape_horizontal_cylinder,
		// the following have not yet been defined or implemented ... :
		shape_sawtooth_down,
		shape_sawtooth_up,
		shape_pyramid,
		// adding my own
		shape_cone,
		shape_prism5,
		// custom (numerical)
		shape_custom			/* a custom shape */
	}; // enum ShapeName

	enum ShapeParamType {
		param_null,
		param_error,
		param_radius,			/* radius */
		param_xsize,			/* length in x dimension */
		param_ysize,			/* length in y dimension */
		param_height,			/* length in z dimension */
		param_edge,				/* length of an edge (for symmetric objects) */
		param_baseangle			/* base angle */
	}; // enum ShapeParamType

	enum LatticeType {
		lattice_null,
		lattice_error,
		lattice_bcc,			/* bcc */
		lattice_cubic,			/* default cubic */
		lattice_fcc,			/* fcc */
		lattice_fco,			/* fco */
		lattice_hcp,			/* hcp */
		lattice_hex				/* hex */
	}; // enum LatticeType

	enum StatisticType {
		stat_null,
		stat_error,
		stat_gaussian,			/* gaussian distribution */
		stat_none,				/* default */
		stat_random,			/* random distribution */
		stat_uniform			/* uniform distribution */
	}; // enum StatisticType

	enum OutputRegionType {
		region_null,
		region_error,
		region_angles,			/* angles region */
		region_pixels,			/* pixels region */
		region_qspace			/* default qspace region */
	}; // enum OutputRegionType

	enum ShapeFileType {
		shape_file_null,
		shape_file_error,
		shape_file_shape,		/* shape file format */
		shape_file_hdf5,		/* shape file in HDF5 format */
		shape_file_object		/* raw object file (e.g. from maya) */
	}; // enum ShapeFileType

	enum StructCorrelationType {
		structcorr_error,		/* error type */
		structcorr_null,		/* default, no correlation */
		structcorr_nGnE,		/* another name for default */
		structcorr_nGE,			/* non-corr grain, corr ensemble */
		structcorr_GnE,			/* corr grain, non-corr ensemble */
		structcorr_GE,			/* both correlated */
	}; // enum StructCorrelationType

} // namespace hig

#endif /* _ENUMS_HPP_ */
