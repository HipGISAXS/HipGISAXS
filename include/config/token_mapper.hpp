/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: token_mapper.hpp
 *  Created: Jun 05, 2012
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

#ifndef _TOKEN_MAPPER_HPP_
#define _TOKEN_MAPPER_HPP_

#include <unordered_map>
#include <string>

#include <config/tokens.hpp>
#include <common/enums.hpp>

namespace hig {

	class TokenMapper {
		private:
			std::unordered_map <std::string, TokenType> KeyWords_;
			std::unordered_map <std::string, ShapeName> ShapeKeyWords_;
			std::unordered_map <std::string, ShapeParamType> ShapeParamKeyWords_;
			std::unordered_map <std::string, StatisticType> StatisticKeyWords_;
			std::unordered_map <std::string, LatticeType> LatticeKeyWords_;
			std::unordered_map <std::string, OutputRegionType> OutputRegionKeyWords_;
			std::unordered_map <std::string, StructCorrelationType> StructCorrelationKeyWords_;

		public:
			static TokenMapper& instance() {
				static TokenMapper token_mapper;
				return token_mapper;
			} // instance()

			
			TokenType get_keyword_token(const std::string& str) {
				if(KeyWords_.count(str) > 0) return KeyWords_[str];
				else return error_token;
			} // get_keyword_token()

			
			ShapeName get_shapename_token(const std::string& str) {
				if(ShapeKeyWords_.count(str) > 0) return ShapeKeyWords_[str];
				else if(has_extension(str, std::string(".shp")) ||
						has_extension(str, std::string(".hd5")) ||
						has_extension(str, std::string(".dat")) ||
						has_extension(str, std::string(".obj")))
					return shape_custom;
				else
					return shape_error;
			} // get_shapename_token()

			
			ShapeParamType get_shapeparam_token(const std::string& str) {
				if(ShapeParamKeyWords_.count(str) > 0) return ShapeParamKeyWords_[str];
				else return param_error;
			} // get_shapeparam_token()

			
			StatisticType get_stattype_token(const std::string& str) {
				if(StatisticKeyWords_.count(str) > 0) return StatisticKeyWords_[str];
				else return stat_error;
			} // get_stattype_token()


			LatticeType get_lattice_type(const std::string& str) {
				if(LatticeKeyWords_.count(str) > 0) return LatticeKeyWords_[str];
				else return lattice_error;
			} // get_lattice_type_token()


			OutputRegionType get_output_region_type(const std::string& str) {
				if(OutputRegionKeyWords_.count(str) > 0) return OutputRegionKeyWords_[str];
				else return region_error;
			} // get_output_region_type()


			StructCorrelationType get_structcorr_type(const std::string& str) {
				if(StructCorrelationKeyWords_.count(str) > 0) return StructCorrelationKeyWords_[str];
				else return structcorr_error;
			} // get_structcorr_type()

			
			bool key_exists(const std::string& str) {
				if(KeyWords_.count(str) > 0 ||
						ShapeKeyWords_.count(str) > 0 ||
						ShapeParamKeyWords_.count(str) > 0 ||
						StatisticKeyWords_.count(str) > 0 ||
						LatticeKeyWords_.count(str) > 0 ||
						OutputRegionKeyWords_.count(str) > 0)
					return true;
				return false;
			} // keyword_exists()


			bool keyword_token_exists(TokenType token) {
				// a sequential search for value 'token'
				std::unordered_map <std::string, TokenType>::iterator i = KeyWords_.begin();
				while(i != KeyWords_.end()) { if((*i).second == token) return true; i ++; }
				return false;
			} // token_exists()


		private:

			/* constructor */
			TokenMapper() {

				/* language keywords */

				KeyWords_[std::string("a")]					= struct_grain_lattice_a_token;
				KeyWords_[std::string("abangle")]			= struct_grain_lattice_abangle_token;
				KeyWords_[std::string("alphai")]			= instrument_scatter_alphai_token;
				KeyWords_[std::string("angles")]			= struct_ensemble_orient_rot_angles_token;
				KeyWords_[std::string("axis")]				= struct_ensemble_orient_rot_axis_token;
				KeyWords_[std::string("b")]					= struct_grain_lattice_b_token;
				KeyWords_[std::string("beta")]				= refindex_beta_token;
				KeyWords_[std::string("c")]					= struct_grain_lattice_c_token;
				KeyWords_[std::string("caratio")]			= struct_grain_lattice_caratio_token;
				KeyWords_[std::string("coherence")]			= instrument_scatter_coherence_token;
				KeyWords_[std::string("computation")]		= compute_token;
				KeyWords_[std::string("delta")]				= refindex_delta_token;
				KeyWords_[std::string("detector")]			= instrument_detector_token;
				KeyWords_[std::string("directbeam")]		= instrument_detector_dirbeam_token;
				KeyWords_[std::string("distribution")]		= struct_ensemble_distribution_token;
				KeyWords_[std::string("ensemble")]			= struct_ensemble_token;
				KeyWords_[std::string("expt")]				= instrument_scatter_expt_token;
				KeyWords_[std::string("fitparam")]			= fit_param_token;
				KeyWords_[std::string("fitting")]			= fit_token;
				KeyWords_[std::string("grain")]				= struct_grain_token;
				KeyWords_[std::string("hipGisaxsInput")]	= hipgisaxs_token;
				KeyWords_[std::string("hkl")]				= struct_grain_lattice_hkl_token;
				KeyWords_[std::string("init")]				= fit_param_init_token;
				KeyWords_[std::string("inplanerot")]		= instrument_scatter_inplanerot_token;
				KeyWords_[std::string("instrumentation")]	= instrument_token;
				KeyWords_[std::string("iratio")]			= struct_iratio_token;
				KeyWords_[std::string("key")]				= key_token;
				KeyWords_[std::string("lattice")]			= struct_grain_lattice_token;
				KeyWords_[std::string("layer")]				= layer_token;
				KeyWords_[std::string("layer:key")]			= struct_grain_lkey_token;
				KeyWords_[std::string("max")]				= max_token;
				KeyWords_[std::string("maxgrains")]			= struct_ensemble_maxgrains_token;
				KeyWords_[std::string("maxpoint")]			= compute_outregion_maxpoint_token;
				KeyWords_[std::string("method")]			= compute_method_token;
				KeyWords_[std::string("min")]				= min_token;
				KeyWords_[std::string("minpoint")]			= compute_outregion_minpoint_token;
				KeyWords_[std::string("name")]				= shape_name_token;
				KeyWords_[std::string("nslices")]			= compute_nslices_token;
				KeyWords_[std::string("nvalues")]			= shape_param_nvalues_token;
				KeyWords_[std::string("order")]				= layer_order_token;
				KeyWords_[std::string("orientations")]		= struct_ensemble_orient_token;
				KeyWords_[std::string("origin")]			= instrument_detector_origin_token;
				KeyWords_[std::string("originvec")]			= shape_originvec_token;
				KeyWords_[std::string("outputregion")]		= compute_outregion_token;
				KeyWords_[std::string("p1")]				= shape_param_p1_token;		// mean
				KeyWords_[std::string("p2")]				= shape_param_p2_token;		// std dev
				KeyWords_[std::string("param")]				= shape_param_token;
				KeyWords_[std::string("pathprefix")]		= compute_path_token;
				KeyWords_[std::string("photon")]			= instrument_scatter_photon_token;
				KeyWords_[std::string("pixelsize")]			= instrument_detector_pixsize_token;
				KeyWords_[std::string("polarization")]		= instrument_scatter_polarize_token;
				KeyWords_[std::string("range")]				= fit_param_range_token;
				KeyWords_[std::string("repetition")]		= struct_grain_repetition_token;
				KeyWords_[std::string("refindex")]			= refindex_token;
				KeyWords_[std::string("resolution")]		= compute_resolution_token;
				KeyWords_[std::string("rot")]				= rot_token;
				KeyWords_[std::string("rot1")]				= struct_ensemble_orient_rot1_token;
				KeyWords_[std::string("rot2")]				= struct_ensemble_orient_rot2_token;
				KeyWords_[std::string("rot3")]				= struct_ensemble_orient_rot3_token;
				KeyWords_[std::string("runname")]			= compute_runname_token;
				KeyWords_[std::string("scaling")]			= struct_grain_scaling_token;
				KeyWords_[std::string("scattering")]		= instrument_scatter_token;
				KeyWords_[std::string("sdd")]				= instrument_detector_sdd_token;
				KeyWords_[std::string("shape")]				= shape_token;
				KeyWords_[std::string("shape:key")]			= struct_grain_skey_token;
				KeyWords_[std::string("smearing")]			= instrument_scatter_smearing_token;
				KeyWords_[std::string("spacing")]			= struct_ensemble_spacing_token;
				KeyWords_[std::string("spotarea")]			= instrument_scatter_spotarea_token;
				KeyWords_[std::string("stat")]				= stat_token;
				KeyWords_[std::string("step")]				= step_token;
				KeyWords_[std::string("structcorrelation")]	= compute_structcorr_token;
				KeyWords_[std::string("structure")]			= struct_token;
				KeyWords_[std::string("thickness")]			= layer_thickness_token;
				KeyWords_[std::string("tilt")]				= instrument_scatter_tilt_token;
				KeyWords_[std::string("totalpixels")]		= instrument_detector_totpix_token;
				KeyWords_[std::string("transvec")]			= struct_grain_transvec_token;
				KeyWords_[std::string("type")]				= type_token;
				KeyWords_[std::string("unit")]				= instrument_scatter_photon_unit_token;
				KeyWords_[std::string("value")]				= instrument_scatter_photon_value_token;
				KeyWords_[std::string("variable")]			= fit_param_variable_token;
				KeyWords_[std::string("xyrotation")]		= shape_xyrot_token;
				KeyWords_[std::string("ztilt")]				= shape_ztilt_token;
			
				/* shape name keywords */

				ShapeKeyWords_[std::string("box")]			= shape_box;
				ShapeKeyWords_[std::string("cylinder")]		= shape_cylinder;
				ShapeKeyWords_[std::string("hcylinder")]	= shape_horizontal_cylinder;
				ShapeKeyWords_[std::string("randcylinders")]= shape_random_cylinders;
				ShapeKeyWords_[std::string("sphere")]		= shape_sphere;
				ShapeKeyWords_[std::string("truncpyr")]		= shape_truncpyr;
				ShapeKeyWords_[std::string("trunccone")]	= shape_trunccone;
				ShapeKeyWords_[std::string("prism3")]		= shape_prism3;
				ShapeKeyWords_[std::string("prism6")]		= shape_prism6;
				ShapeKeyWords_[std::string("prism3x")]		= shape_prism3x;
				ShapeKeyWords_[std::string("sawtooth")]		= shape_sawtooth;
				ShapeKeyWords_[std::string("custom")]		= shape_custom;
			
				/* shape parameter type keywords */

				ShapeParamKeyWords_[std::string("radius")]		= param_radius;
				ShapeParamKeyWords_[std::string("xsize")]		= param_xsize;
				ShapeParamKeyWords_[std::string("ysize")]		= param_ysize;
				ShapeParamKeyWords_[std::string("height")]		= param_height;
				ShapeParamKeyWords_[std::string("edge")]		= param_edge;
				ShapeParamKeyWords_[std::string("baseangle")]	= param_baseangle;
			
				/* statistic type keywords */

				StatisticKeyWords_[std::string("none")]		= stat_none;
				StatisticKeyWords_[std::string("uniform")]	= stat_uniform;
				StatisticKeyWords_[std::string("random")]	= stat_random;
				StatisticKeyWords_[std::string("gaussian")]	= stat_gaussian;

				/* lattice type keywords */

				LatticeKeyWords_[std::string("bcc")]	= lattice_bcc;
				LatticeKeyWords_[std::string("cubic")]	= lattice_cubic;
				LatticeKeyWords_[std::string("fcc")]	= lattice_fcc;
				LatticeKeyWords_[std::string("fco")]	= lattice_fco;
				LatticeKeyWords_[std::string("hcp")]	= lattice_hcp;
				LatticeKeyWords_[std::string("hex")]	= lattice_hex;

				/* output region type keywords */

				OutputRegionKeyWords_[std::string("angles")]	= region_angles;
				OutputRegionKeyWords_[std::string("pixels")]	= region_pixels;
				OutputRegionKeyWords_[std::string("qspace")]	= region_qspace;

				/* structure grain/ensemble correlation keywords */

				StructCorrelationKeyWords_[std::string("nGnE")]	= structcorr_nGnE;
				StructCorrelationKeyWords_[std::string("nGE")]	= structcorr_nGE;
				StructCorrelationKeyWords_[std::string("GnE")]	= structcorr_GnE;
				StructCorrelationKeyWords_[std::string("GE")]	= structcorr_GE;

		} // TokenMapper()

		// singleton
		TokenMapper(const TokenMapper&);
		TokenMapper& operator=(const TokenMapper&);

		bool has_extension(const std::string& s, const std::string& e) const {
			unsigned int p = s.rfind(e);
			if(p != s.length() - e.length())	// s is not a suffix of
				return false;
			return true;
		} // has_extension()

	}; // class TokenMapper

} // namespace hig

#endif /* _TOKEN_MAPPER_HPP_ */
