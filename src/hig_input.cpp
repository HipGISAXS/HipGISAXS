/***
  *  $Id: hig_input.cpp 47 2012-08-23 21:05:16Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: hig_input.cpp
  *  Created: Jun 11, 2012
  *  Modified: Mon 08 Apr 2013 04:05:50 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <cfloat>

#include "hig_input.hpp"
#include "utilities.hpp"
#include "parameters.hpp"
#include "hig_file_reader.hpp"
//#include "object2hdf5.h"


namespace hig {


	HiGInput::HiGInput() {
		InputReader::instance();
		TokenMapper::instance();
		HiGFileReader::instance();
	} // HiGInput::HiGInput()


	void HiGInput::init() {
		shapes_.clear();
		layers_.clear();
		structures_.clear();
		scattering_.init();
		detector_.init();
		compute_.init();
		struct_in_layer_ = false;

		shape_def_.clear();
	} // init();


	bool HiGInput::construct_input_config(char* filename) {
		//if(!input_reader_.init(filename)) {
		if(!InputReader::instance().read_input(filename)) {
			std::cerr << "fatal error: some error happened in opening or reading "
						<< "input config file. aborting"
						<< std::endl;
			return false;
		} // if

		curr_keyword_ = null_token; past_keyword_ = null_token;
		curr_token_ = InputReader::instance().get_next_token();
		past_token_.type_ = null_token;
		while(curr_token_.type_ != null_token) {
			if(curr_token_.type_ == error_token) {
				std::cerr << "aborting due to fatal error" << std::endl;
				return false;
			} // if
			if(!process_curr_token()) {
				std::cerr << "aborting due to fatal error" << std::endl;
				return false;
			} // if
			past_token_ = curr_token_;
			curr_token_ = InputReader::instance().get_next_token();
		} // while

		return true;
	} // HiGInput::read_input_config()


	bool HiGInput::process_curr_token() {
		TokenType parent = null_token;
		TokenType gparent = null_token;

		// process the token, do some basic syntax checking (improve with AST later) ...
		switch(curr_token_.type_) {

			case error_token:
				std::cerr << "aborting due to error" << std::endl;
				return false;

			case null_token:
				std::cerr << "error: something went wrong - should have already stopped!"
							<< std::endl;
				return false;

			case white_space_token:
				std::cerr << "error: something went wrong - "
							<< "seeing whitespace when not supposed to!" << std::endl;
				return false;

			case comment_token:	// do nothing
				//std::cerr << "error: something went wrong - "
				//			<< "seeing comments when not supposed to!" << std::endl;
				return true;

			case object_begin_token:	// should be preceeded by '='
				if(past_token_.type_ != assignment_token && past_token_.type_ != comment_token) {
					std::cerr << "fatal error: unexpected object begin token '{'"
								<< std::endl;
					return false;
				} // if
				keyword_stack_.push(curr_keyword_);
				break;

			case object_end_token:		// preceeded by number or string or '}'
				if(past_token_.type_ != number_token &&
						past_token_.type_ != string_token &&
						past_token_.type_ != object_end_token &&
						past_token_.type_ != array_end_token &&
						past_token_.type_ != object_begin_token &&
						past_token_.type_ != comment_token) {
					std::cerr << "fatal error: unexpected object close token '}'" << std::endl;
					return false;
				} // if
				if(keyword_stack_.size() < 1) {
					std::cerr << "fatal error: unexpected object close token '}'. "
								<< "no matching object open token found" << std::endl;
					return false;
				} // if

				parent = get_curr_parent();
				switch(parent) {
					case shape_token:
						shapes_[curr_shape_.key()] = curr_shape_;	// insert the shape
						curr_shape_.clear();
						break;

					case shape_param_token:
						curr_shape_param_.set();
						//curr_shape_param_.print();
						curr_shape_.insert_param(curr_shape_param_.type_name(), curr_shape_param_);
						curr_shape_param_.clear();
						break;

					case refindex_token:
						// find out which ref index is this for
						gparent = get_curr_grandparent();
						switch(gparent) {
							case layer_token:	// nothing to do :-/ ??
								break;

							case struct_grain_token:	// nothing to do :-/ ??
								break;

							default:
								std::cerr << "error: wrong place for a refindex" << std::endl;
								return false;
						} // switch
						break;

					case layer_token:			/* insert the current layer, and its key map */
						layers_[curr_layer_.order()] = curr_layer_;
						layer_key_map_[curr_layer_.key()] = curr_layer_.order();
						curr_layer_.clear();
						break;

					case struct_token:			/* insert the current structure */
						structures_[curr_structure_.key()] = curr_structure_;
						curr_structure_.clear();
						break;

					case struct_grain_lattice_token:	// nothing to do :-/
					case struct_grain_token:	// nothing to do :-/
					case struct_ensemble_orient_stat_token:	// nothing to do :-/
					case struct_ensemble_orient_rot1_token:	// nothing to do :-/
					case struct_ensemble_orient_rot2_token:	// nothing to do :-/
					case struct_ensemble_orient_rot3_token:	// nothing to do :-/
					case struct_ensemble_orient_token:	// nothing to do :-/
					case struct_ensemble_token:	// nothing to do :-/
					case instrument_scatter_alphai_token:	// nothing to do :-/
					case instrument_scatter_inplanerot_token:	// nothing to do :-/
					case instrument_scatter_tilt_token:	// nothing to do :-/
					case instrument_scatter_photon_token:	// nothing to do :-/
					case instrument_scatter_token:	// nothing to do :-/
					case instrument_detector_token:	// nothing to do :-/
					case instrument_token:	// nothing to do :-/
					case compute_outregion_token:	// nothing to do :-/
					case compute_token:	// nothing to do :-/
					case hipgisaxs_token:	// nothing to do :-/
						break;

					default:
						std::cerr << "error: something is wrong with one of your objects"
									<< std::endl;
						return false;
				} // switch
				if(keyword_stack_.size() < 1) {
					std::cerr << "something is really wrong. keyword_stack_ is empty when "
					   			<< "object end was found" << std::endl;
					return false;
				} // if
				past_keyword_ = curr_keyword_;
				curr_keyword_ = keyword_stack_.top();
				keyword_stack_.pop();
				break;

			case array_begin_token:	// should be preceeded by '='
				if(past_token_.type_ != assignment_token) {
					std::cerr << "fatal error: unexpected array begin token '['"
								<< std::endl;
					return false;
				} // if
				keyword_stack_.push(curr_keyword_);
				break;

			case array_end_token:	// preceeded by number_token or array_begin_token
				if(past_token_.type_ != number_token &&
						past_token_.type_ != array_begin_token &&
						past_token_.type_ != comment_token) {
					std::cerr << "fatal error: unexpected array close token ']'"
								<< std::endl;
					return false;
				} // if
				if(keyword_stack_.size() < 1) {
					std::cerr << "fatal error: unexpected array close token ']', "
								<< "no matching array open token found" << std::endl;
					return false;
				} // if

				parent = keyword_stack_.top();
				switch(parent) {
					case shape_originvec_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in originvec" << std::endl;
							return false;
						} // if
						curr_shape_.originvec(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						break;

					case struct_grain_lattice_a_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in lattice vector a"
										<< std::endl;
							return false;
						} // if
						curr_structure_.lattice_vec_a(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						curr_structure_.lattice_abc_set(true);
						break;

					case struct_grain_lattice_b_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in lattice vector b"
										<< std::endl;
							return false;
						} // if
						curr_structure_.lattice_vec_b(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						curr_structure_.lattice_abc_set(true);
						break;

					case struct_grain_lattice_c_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in lattice vector c"
										<< std::endl;
							return false;
						} // if
						curr_structure_.lattice_vec_c(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						curr_structure_.lattice_abc_set(true);
						break;

					case struct_grain_transvec_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in grain transvec"
										<< std::endl;
							return false;
						} // if
						curr_structure_.grain_transvec(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						break;

					case struct_grain_repetition_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in grain repetition"
										<< std::endl;
							return false;
						} // if
						curr_structure_.grain_repetition(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						break;

					case struct_ensemble_spacing_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in ensemble spacing"
										<< std::endl;
							return false;
						} // if
						curr_structure_.ensemble_spacing(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						break;

					case struct_ensemble_maxgrains_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: less than 3 values in ensemble maxgrains"
										<< std::endl;
							return false;
						} // if
						curr_structure_.ensemble_maxgrains(
								curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						break;

					case struct_ensemble_orient_rot_angles_token:
						if(curr_vector_.size() != 2) {
							std::cerr << "error: values in orientation rotation angles should be 2"
										<< std::endl;
							return false;
						} // if
						// find out which rot is this for
						gparent = get_curr_grandparent();
						switch(gparent) {
							case struct_ensemble_orient_rot1_token:
								curr_structure_.grain_orientation_rot1_angles(
										curr_vector_[0], curr_vector_[1]);
								break;
							case struct_ensemble_orient_rot2_token:
								curr_structure_.grain_orientation_rot2_angles(
										curr_vector_[0], curr_vector_[1]);
								break;
							case struct_ensemble_orient_rot3_token:
								curr_structure_.grain_orientation_rot3_angles(
										curr_vector_[0], curr_vector_[1]);
								break;
							default:
								std::cerr << "error: something wrong in the rot angles" << std::endl;
								return false;
						} // switch
						break;

					case instrument_scatter_smearing_token:
						if(curr_vector_.size() != 3) {
							std::cerr << "error: scattering smearing vector size should be 3"
										<< std::endl;
							return false;
						} // if
						scattering_.smearing(curr_vector_[0], curr_vector_[1], curr_vector_[2]);
						break;

					case instrument_detector_totpix_token:
						if(curr_vector_.size() != 2) {
							std::cerr << "error: totalpixels vector size should be 2"
										<< std::endl;
							return false;
						} // if
						detector_.total_pixels(curr_vector_[0], curr_vector_[1]);
						break;

					case instrument_detector_dirbeam_token:
						if(curr_vector_.size() != 2) {
							std::cerr << "error: detector direct beam vector size should be 2"
										<< std::endl;
							return false;
						} // if
						detector_.direct_beam(curr_vector_[0], curr_vector_[1]);
						break;

					case compute_resolution_token:
						if(curr_vector_.size() != 2) {
							std::cerr << "error: resolution vector size should be 2"
										<< std::endl;
							return false;
						} // if
						compute_.resolution(curr_vector_[0], curr_vector_[1]);
						break;

					case compute_outregion_minpoint_token:
						if(curr_vector_.size() != 2) {
							std::cerr << "error: output region min point vector size should be 2"
										<< std::endl;
							return false;
						} // if
						compute_.output_region_minpoint(curr_vector_[0], curr_vector_[1]);
						break;

					case compute_outregion_maxpoint_token:
						if(curr_vector_.size() != 2) {
							std::cerr << "error: output region max point vector size should be 2"
										<< std::endl;
							return false;
						} // if
						compute_.output_region_maxpoint(curr_vector_[0], curr_vector_[1]);
						break;

					default:
						std::cerr << "error: found array value in place of non-array type" << std::endl;
						return false;
				} // switch
				curr_vector_.clear();
				keyword_stack_.pop();
				past_keyword_ = curr_keyword_;
				curr_keyword_ = keyword_stack_.top();
				break;

			case assignment_token:	// should be preceeded by a 'keyword',
									// followed by '{' (object) or '[' (array)
									// or string or number
				if(!preceeded_by_keyword()) {
					std::cerr << "error: misplaced assignment token '='" << std::endl;
					return false;
				} // if
				break;

			case number_token:		// preceeded by '=' or '[' or number_token
				if(past_token_.type_ != assignment_token &&
						past_token_.type_ != array_begin_token &&
						past_token_.type_ != number_token &&
						past_token_.type_ != comment_token &&
						past_token_.type_ != white_space_token) {
					std::cerr << "error: unexpected number '"
								<< curr_token_.dvalue_ << "'" << std::endl;
					return false;
				} // if
				if(!process_number(curr_token_.dvalue_)) {
					std::cerr << "error: could not process number '"
								<< curr_token_.dvalue_ << "'" << std::endl;
					return false;
				} // if
				break;

			case string_token:		// preceeded by '='
				if(past_token_.type_ != assignment_token &&
						past_token_.type_ != comment_token) {
					std::cerr << "error: stray string found '"
								<< curr_token_.svalue_ << "'" << std::endl;
					return false;
				} // if
				if(!process_string(curr_token_.svalue_)) {
					std::cerr << "error: could not process string "
								<< curr_token_.svalue_ << std::endl;
					return false;
				} // if
				break;

			case separator_token:	// should be preceeded by
									// array_end or string or number or object_end
				if(past_token_.type_ != array_end_token &&
						past_token_.type_ != object_end_token &&
						past_token_.type_ != string_token &&
						past_token_.type_ != number_token &&
						past_token_.type_ != comment_token) {
					std::cerr << "error: stray seperator token ',' found" << std::endl;
					return false;
				} // if
				break;

			default:				// this is for keyword tokens
									// read_oo_input makes sure there are no illegal tokens
									// this is always preceeded by ',' or '{'
				if(curr_token_.type_ != hipgisaxs_token &&
						past_token_.type_ != object_begin_token &&
						past_token_.type_ != separator_token &&
						past_token_.type_ != comment_token) {
					std::cerr << "error: keyword '" << curr_token_.svalue_
								<< "' not placed properly" << std::endl;
					return false;
				} // if
				past_keyword_ = curr_keyword_;
				curr_keyword_ = curr_token_.type_;
				if(!process_curr_keyword()) {
					std::cerr << "error: could not process current keyword '" << curr_token_.svalue_
								<< "'" << std::endl;
					return false;
				} // if
				break;
		} // switch

		return true;
	} // HiGInput::process_curr_token()


	bool HiGInput::process_curr_keyword() {
		// do some syntax error checkings
		switch(curr_keyword_) {

			case hipgisaxs_token:	// this will always be the first token
				if(past_token_.type_ != null_token || keyword_stack_.size() != 0) {
					std::cerr << "fatal error: 'hipGisaxsInput' token is not at the beginning!"
								<< std::endl;
					return false;
				} // if

				init();		// initialize everything
				break;

			case key_token:
			case min_token:
			case max_token:
			case step_token:
			case rot_token:
			case type_token:
			case stat_token:
				break;

			case refindex_token:
			case refindex_delta_token:
			case refindex_beta_token:
				break;

			case shape_token:
				curr_shape_.init();
				break;

			case shape_name_token:
			case shape_originvec_token:
			case shape_ztilt_token:
			case shape_xyrot_token:
				break;

			case shape_param_token:
				curr_shape_param_.init();
				break;

			case shape_param_p1_token:
			case shape_param_p2_token:
			case shape_param_nvalues_token:
				break;

			case layer_token:
				curr_layer_.init();
				break;

			case layer_order_token:
			case layer_thickness_token:
				break;

			case struct_token:
				curr_structure_.init();
				break;

			case struct_grain_token:
			case struct_grain_skey_token:
			case struct_grain_lkey_token:
			case struct_grain_lattice_token:
			case struct_grain_lattice_a_token:
			case struct_grain_lattice_b_token:
			case struct_grain_lattice_c_token:
			case struct_grain_lattice_hkl_token:
			case struct_grain_lattice_abangle_token:
			case struct_grain_lattice_caratio_token:
			case struct_grain_transvec_token:
			case struct_grain_scaling_token:
			case struct_grain_repetition_token:
				break;

			case struct_ensemble_token:
			case struct_ensemble_spacing_token:
			case struct_ensemble_maxgrains_token:
			case struct_ensemble_distribution_token:
			case struct_ensemble_orient_token:
			case struct_ensemble_orient_stat_token:
			case struct_ensemble_orient_rot1_token:
			case struct_ensemble_orient_rot2_token:
			case struct_ensemble_orient_rot3_token:
			case struct_ensemble_orient_rot_axis_token:
			case struct_ensemble_orient_rot_angles_token:
				break;

			case compute_token:
			case compute_path_token:
			case compute_runname_token:
			case compute_method_token:
			case compute_resolution_token:
			case compute_nslices_token:
			case compute_outregion_token:
			case compute_outregion_maxpoint_token:
			case compute_outregion_minpoint_token:
				break;

			case instrument_token:
			case instrument_scatter_token:
			case instrument_scatter_expt_token:
			case instrument_scatter_alphai_token:
			case instrument_scatter_inplanerot_token:
			case instrument_scatter_tilt_token:
			case instrument_scatter_photon_token:
			case instrument_scatter_photon_value_token:
			case instrument_scatter_photon_unit_token:
			case instrument_scatter_polarize_token:
			case instrument_scatter_coherence_token:
			case instrument_scatter_spotarea_token:
			case instrument_scatter_smearing_token:
				break;

			case instrument_detector_token:
			case instrument_detector_origin_token:
			case instrument_detector_totpix_token:
			case instrument_detector_sdd_token:
			case instrument_detector_pixsize_token:
			case instrument_detector_dirbeam_token:
				break;

			default:
				std::cerr << "error: non keyword token in keyword's position"
							<< std::endl;
				return false;
		} // switch()

		return true;
	} // HiGInput::process_curr_keyword()


	inline TokenType HiGInput::get_curr_parent() {
		if(keyword_stack_.size() < 1) return null_token;
		return keyword_stack_.top();
	} // HiGInput::get_curr_parent()


	inline TokenType HiGInput::get_curr_grandparent() {
		if(keyword_stack_.size() < 1) return null_token;
		TokenType temp = keyword_stack_.top();
		keyword_stack_.pop();
		if(keyword_stack_.size() < 1) { keyword_stack_.push(temp); return null_token; }
		TokenType gparent = keyword_stack_.top();
		keyword_stack_.push(temp);
		return gparent;
	} // HiGInput::get_curr_grandparent()


	bool HiGInput::process_number(const float_t& num) {
		TokenType parent = null_token;
		TokenType gparent = null_token;

		switch(curr_keyword_) {
			 
			case min_token:
				// find out which min is this for
				// shape param, scattering alphai, inplanerot, tilt
				parent = get_curr_parent();
				switch(parent) {
					case shape_param_token:
						curr_shape_param_.min(num);
						break;

					case instrument_scatter_alphai_token:
						scattering_.alphai_min(num);
						break;

					case instrument_scatter_inplanerot_token:
						scattering_.inplane_rot_min(num);
						break;

					case instrument_scatter_tilt_token:
						scattering_.tilt_min(num);
						break;

					default:
						std::cerr << "'min' token appears in wrong place" << std::endl;
						return false;
				} // switch
				break;

			case max_token:
				// find out which max is this for
				// shape param, scattering alphai, inplanerot, tilt
				parent = get_curr_parent();
				switch(parent) {
					case shape_param_token:
						curr_shape_param_.max(num);
						break;

					case instrument_scatter_alphai_token:
						scattering_.alphai_max(num);
						break;

					case instrument_scatter_inplanerot_token:
						scattering_.inplane_rot_max(num);
						break;

					case instrument_scatter_tilt_token:
						scattering_.tilt_max(num);
						break;

					default:
						std::cerr << "'max' token appears in wrong place" << std::endl;
						return false;
				} // switch
				break;

			case step_token:
				// find out which step is this for
				// scattering alphai, inplanerot, tilt
				parent = get_curr_parent();
				switch(parent) {
					case instrument_scatter_alphai_token:
						scattering_.alphai_step(num);
					case instrument_scatter_inplanerot_token:
						scattering_.inplane_rot_step(num);
						break;

					case instrument_scatter_tilt_token:
						scattering_.tilt_step(num);
						break;

					default:
						std::cerr << "'step' token appears in a wrong place" << std::endl;
						return false;
				} // switch
				break;

			//case rot_token:
			//	break;

			case refindex_delta_token:
				// find out which ref index is this for
				// layer, grain
				parent = get_curr_parent();
				gparent = get_curr_grandparent();
				switch(gparent) {
					case layer_token:
						curr_layer_.refindex_delta(num);
						break;

					case struct_grain_token:
						curr_structure_.grain_refindex_delta(num);
						break;

					default:
						std::cerr << "'refindex' token appears in a wrong place" << std::endl;
						return false;
				} // switch
				break;

			case refindex_beta_token:
				// find out which ref index is this for
				// layer, grain
				parent = get_curr_parent();
				gparent = get_curr_grandparent();
				switch(gparent) {
					case layer_token:
						curr_layer_.refindex_beta(num);
						break;

					case struct_grain_token:
						curr_structure_.grain_refindex_beta(num);
						break;

					default:
						std::cerr << "'refindex' token appears in a wrong place" << std::endl;
						return false;
				} // switch
				break;


			case shape_originvec_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in origin vector" << std::endl;
					return false;
				} // if
				break;

			case shape_ztilt_token:
				curr_shape_.ztilt(num);
				break;

			case shape_xyrot_token:
				curr_shape_.xyrotation(num);
				break;


			case shape_param_p1_token:
				curr_shape_param_.p1(num);
				break;

			case shape_param_p2_token:
				curr_shape_param_.p2(num);
				break;

			case shape_param_nvalues_token:
				curr_shape_param_.nvalues(num);
				break;


			case layer_order_token:
				curr_layer_.order(num);
				break;

			case layer_thickness_token:
				curr_layer_.thickness(num);
				break;

			case struct_grain_lattice_a_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in lattice vector a" << std::endl;
					return false;
				} // if
				break;

			case struct_grain_lattice_b_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in lattice vector b" << std::endl;
					return false;
				} // if
				break;

			case struct_grain_lattice_c_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in lattice vector c" << std::endl;
					return false;
				} // if
				break;

			case struct_grain_lattice_abangle_token:
				curr_structure_.lattice_abangle(num);
				break;
				
			case struct_grain_lattice_caratio_token:
				curr_structure_.lattice_caratio(num);
				break;
				
			case struct_grain_transvec_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in trans vector" << std::endl;
					return false;
				} // if
				break;

			case struct_grain_scaling_token:
				curr_structure_.grain_scaling(num);
				break;

			case struct_grain_repetition_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in repetition vector" << std::endl;
					return false;
				} // if
				break;

			case struct_ensemble_spacing_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in spacing vector" << std::endl;
					return false;
				} // if
				break;

			case struct_ensemble_maxgrains_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in maxgrains vector" << std::endl;
					return false;
				} // if
				break;

			case struct_ensemble_orient_rot_angles_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 2) {
					std::cerr << "error: more than 2 values in angles vector" << std::endl;
					return false;
				} // if
				break;

			case compute_outregion_maxpoint_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 2) {
					std::cerr << "error: more than 2 values in maxpoint" << std::endl;
					return false;
				} // if
				break;

			case compute_outregion_minpoint_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 2) {
					std::cerr << "error: more than 2 values in minpoint" << std::endl;
					return false;
				} // if
				break;

			case compute_resolution_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 2) {
					std::cerr << "error: more than 2 values in resolution vector" << std::endl;
					return false;
				} // if
				break;

			case compute_nslices_token:
				compute_.nslices(num);
				break;


			case instrument_scatter_photon_value_token:
				scattering_.photon_value(num);
				break;

			case instrument_scatter_coherence_token:
				scattering_.coherence(num);
				break;

			case instrument_scatter_spotarea_token:
				scattering_.spot_area(num);
				break;

			case instrument_scatter_smearing_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 3) {
					std::cerr << "error: more than 3 values in scatter smearing" << std::endl;
					return false;
				} // if
				break;

			case instrument_detector_totpix_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 2) {
					std::cerr << "error: more than 2 values in totalpixels vector" << std::endl;
					return false;
				} // if
				break;

			case instrument_detector_sdd_token:
				detector_.sd_distance(num);
				break;

			case instrument_detector_pixsize_token:
				detector_.pixel_size(num);
				break;

			case instrument_detector_dirbeam_token:
				curr_vector_.push_back(num);
				if(curr_vector_.size() > 2) {
					std::cerr << "error: more than 2 values in directbeam vector" << std::endl;
					return false;
				} // if
				break;

			default:
				std::cerr << "fatal error: found a number '" << num
							<< "' where it should not be" << std::endl;
				return false;
		} // switch

		return true;
	} // HiGInput::process_number()


	bool HiGInput::process_string(const std::string& str) {
		TokenType parent = null_token;
		ShapeName shp = shape_null;

		switch(curr_keyword_) {

			case key_token:
				// find out which key is this for
				parent = get_curr_parent();
				switch(parent) {
					case shape_token:
						curr_shape_.key(str);
						break;

					case layer_token:
						curr_layer_.key(str);
						break;

					case struct_token:
						curr_structure_.key(str);
						break;

					default:
						std::cerr << "error: extraneous key" << std::endl;
						return false;
				} // switch
				break;

			case shape_name_token:
				shp = TokenMapper::instance().get_shapename_token(str);
				if(shp == shape_error) {
					std::cerr << "error: shape name '" << str
								<< "' is an unknown shape, and is not a shape file" << std::endl;
					return false;
				} // if
				curr_shape_.name_str(str);
				curr_shape_.name(shp);
				break;

			case type_token:
				// find out which type is this for
				parent = get_curr_parent();
				switch(parent) {
					case shape_param_token:
						curr_shape_param_.type(TokenMapper::instance().get_shapeparam_token(str));
						curr_shape_param_.type_name(str);
						break;

					case struct_grain_lattice_token:
						curr_structure_.lattice_type(
								TokenMapper::instance().get_lattice_type(str));
						break;

					case compute_outregion_token:
						compute_.output_region_type(
								TokenMapper::instance().get_output_region_type(str));
						break;

					default:
						std::cerr << "error: 'type' token in wrong place" << std::endl;
						return false;
				} // switch
				break;

			case stat_token:
				// find out which stat is this for
				parent = get_curr_parent();
				switch(parent) {
					case shape_param_token:
						curr_shape_param_.stat(TokenMapper::instance().get_stattype_token(str));
						break;

					case struct_ensemble_orient_token:
						curr_structure_.ensemble_orientation_stat(str);
						break;

					default:
						std::cerr << "error: 'stat' token in wrong place" << std::endl;
						return false;
				} // switch
				break;

			case struct_grain_skey_token:
				curr_structure_.grain_shape_key(str);
				break;

			case struct_grain_lkey_token:
				curr_structure_.grain_layer_key(str);
				struct_in_layer_ = true;
				break;

			case struct_ensemble_distribution_token:
				curr_structure_.ensemble_distribution(str);
				break;

			case struct_ensemble_orient_rot_axis_token:
				// find out which of the 3 rot is this for
				parent = get_curr_parent();
				switch(parent) {
					case struct_ensemble_orient_rot1_token:
						curr_structure_.grain_orientation_rot1_axis(str.c_str()[0]);
						break;

					case struct_ensemble_orient_rot2_token:
						curr_structure_.grain_orientation_rot2_axis(str.c_str()[0]);
						break;

					case struct_ensemble_orient_rot3_token:
						curr_structure_.grain_orientation_rot3_axis(str.c_str()[0]);
						break;

					default:
						std::cerr << "error: 'axis' token in wrong place" << std::endl;
						return false;
				} // switch
				break;

			case struct_grain_lattice_hkl_token:
				curr_structure_.lattice_hkl(str);
				break;

			case instrument_scatter_expt_token:
				scattering_.expt(str);
				break;

			case instrument_scatter_photon_unit_token:
				scattering_.photon_unit(str);
				break;

			case instrument_scatter_polarize_token:
				scattering_.polarization(str);
				break;

			case instrument_detector_origin_token:
				detector_.origin(str);
				break;

			case compute_path_token:
				compute_.pathprefix(str);
				break;

			case compute_runname_token:
				compute_.runname(str);
				break;

			case compute_method_token:
				compute_.method(str);
				break;

			default:
				std::cerr << "fatal error: found a string '"
							<< str << "' where it should not be" << std::endl;
				return false;
		} // switch

		return true;
	} // HiGInput::process_string()



	/**
	 * input accessor and modifier functions
	 */


	/* shapes */

	// technically some of these are not completely correct. correct them later ...	// this function can go into the shape class ...
	bool HiGInput::compute_shape_domain(Shape &shape, vector3_t& min_dim, vector3_t& max_dim, int mpi_rank) {
		min_dim[0] = min_dim[1] = min_dim[2] = 0.0;
		max_dim[0] = max_dim[1] = max_dim[2] = 0.0;
		shape_param_iterator_t param = shape.param_begin();
		std::string shape_filename;
		switch(shape.name()) {
			case shape_box:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							std::cerr << "warning: ignoring the radius value provided for a box shape"
										<< std::endl;
							break;
						case param_xsize:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							break;
						case param_ysize:
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_height:
							max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[2] = -max_dim[2];
							break;
						case param_edge:
							max_dim[0] = max_dim[1] = max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = min_dim[1] = min_dim[2] = -max_dim[0];
							break;
						case param_baseangle:
							// do nothing
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_cylinder:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for cylinder shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for cylinder shape"
										<< std::endl;
							break;
						case param_height:
							max_dim[2] = 2.0 * max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							std::cerr << "warning: ignoring the edge values given for cylinder shape"
										<< std::endl;
							break;
						case param_baseangle:
							// do nothing
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_horizontal_cylinder:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for hcylinder shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for hcylinder shape"
										<< std::endl;
							break;
						case param_height:
							max_dim[2] = 2.0 * max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							std::cerr << "warning: ignoring the edge values given for hcylinder shape"
										<< std::endl;
							break;
						case param_baseangle:
							// do nothing
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_random_cylinders:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for cylinder shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for cylinder shape"
										<< std::endl;
							break;
						case param_height:
							max_dim[2] = 2.0 * max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							std::cerr << "warning: ignoring the edge values given for cylinder shape"
										<< std::endl;
							break;
						case param_baseangle:
							// do nothing
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_sphere:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							max_dim[2] = 2.0 * max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for cylinder shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for cylinder shape"
										<< std::endl;
							break;
						case param_height:
							std::cerr << "warning: ignoring the height values given for cylinder shape"
										<< std::endl;
							break;
						case param_edge:
							std::cerr << "warning: ignoring the edge values given for cylinder shape"
										<< std::endl;
							break;
						case param_baseangle:
							// do nothing
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_truncpyr:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							std::cerr << "warning: ignoring the radius values given for truncpyr shape"
										<< std::endl;
							break;
						case param_xsize:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							break;
						case param_ysize:
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_height:
							max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_baseangle:
							// defines the angle at the base -> xy size at height z
							// do this later ...
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_trunccone:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for trunccone shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for trunccone shape"
										<< std::endl;
							break;
						case param_height:
							max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							std::cerr << "warning: ignoring the edge values given for trunccone shape"
										<< std::endl;
							break;
						case param_baseangle:
							// defines the angle at the base -> xy radius at height z
							// if == 90 is cylinder, if > 90, then max x and y will change ...
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_prism3:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							std::cerr << "warning: ignoring the radius values given for prism3 shape"
										<< std::endl;
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for prism3 shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for prism3 shape"
										<< std::endl;
							break;
						case param_height:
							max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							max_dim[0] = max_dim[1] = max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = min_dim[1] = min_dim[2] = -max_dim[0];
							break;
						case param_baseangle:
							std::cerr << "warning: ignoring the baseangle values given for prism3 shape"
										<< std::endl;
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_prism6:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							std::cerr << "warning: ignoring the radius values given for prism6 shape"
										<< std::endl;
							break;
						case param_xsize:
							std::cerr << "warning: ignoring the xsize values given for prism6 shape"
										<< std::endl;
							break;
						case param_ysize:
							std::cerr << "warning: ignoring the ysize values given for prism6 shape"
										<< std::endl;
							break;
						case param_height:
							max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							max_dim[0] = max_dim[1] = max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = min_dim[1] = min_dim[2] = -max_dim[0];
							break;
						case param_baseangle:
							std::cerr << "warning: ignoring the baseangle values given for prism3 shape"
										<< std::endl;
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_prism3x:
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
							std::cerr << "warning: ignoring the radius values given for prism3x shape"
										<< std::endl;
							break;
						case param_xsize:
							max_dim[0] = max((*param).second.max(), (*param).second.min());
							min_dim[0] = -max_dim[0];
							break;
						case param_ysize:
							max_dim[1] = max((*param).second.max(), (*param).second.min());
							min_dim[1] = -max_dim[1];
							break;
						case param_height:
							max_dim[2] = max((*param).second.max(), (*param).second.min());
							min_dim[2] = 0.0;
							break;
						case param_edge:
							std::cerr << "warning: ignoring the edge values given for prism3x shape"
										<< std::endl;
							break;
						case param_baseangle:
							std::cerr << "warning: ignoring the baseangle values given for prism3x shape"
										<< std::endl;
							break;
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_sawtooth:
				std::cerr << "uh-oh: this shape has not been implemented yet" << std::endl;
				return false;
				while(param != shape.param_end()) {
					switch((*param).second.type()) {
						case param_radius:
						case param_xsize:
						case param_ysize:
						case param_height:
						case param_edge:
						case param_baseangle:
						default:
							std::cerr << "error: invalid parameter found in a shape" << std::endl;
							return false;
					} // switch
					++ param;
				} // while
				return true;

			case shape_custom:
				shape_filename = shape.filename();
				read_shape_definition(shape_filename.c_str());
				compute_shapedef_minmax(min_dim, max_dim, mpi_rank);
				break;

			case shape_null:
				std::cerr << "error: null shape encountered" << std::endl;
				return false;

			case shape_error:
				std::cerr << "error: the shape is an error" << std::endl;
				return false;

			default:
				std::cerr << "error: unknown shape name stored" << std::endl;
				return false;
		} // switch

		return true;
	} // HiGInput::compute_shape_domain()


	bool HiGInput::compute_shapedef_minmax(vector3_t& min_dim, vector3_t& max_dim, int mpi_rank) {
		float_t min_a = shape_def_[4], max_a = shape_def_[4];
		float_t min_b = shape_def_[5], max_b = shape_def_[5];
		float_t min_c = shape_def_[6], max_c = shape_def_[6];

		for(int i = 0; i + 6 < shape_def_.size(); i += 7) {
			min_a = (min_a > shape_def_[i + 4]) ? shape_def_[i + 4] : min_a ;
			max_a = (max_a < shape_def_[i + 4]) ? shape_def_[i + 4] : max_a ;
			min_b = (min_b > shape_def_[i + 5]) ? shape_def_[i + 5] : min_b ;
			max_b = (max_b < shape_def_[i + 5]) ? shape_def_[i + 5] : max_b ;
			min_c = (min_c > shape_def_[i + 6]) ? shape_def_[i + 6] : min_c ;
			max_c = (max_c < shape_def_[i + 6]) ? shape_def_[i + 6] : max_c ;
		} // for

		//std::cout << "------ min = " << min_a << ", " << min_b << ", " << min_c << std::endl;
		//std::cout << "------ max = " << max_a << ", " << max_b << ", " << max_c << std::endl;

		float_t diff_a = max_a - min_a;
        float_t diff_b = max_b - min_b;
        float_t diff_c = max_c - min_c;

		//std::cout << "++ diff_a = " << diff_a << ", diff_b = " << diff_b
		//			<< ", diff_c = " << diff_c << std::endl;

		vector3_t axes;
        // axes[i] = j
        // i: x=0 y=1 z=2
        // j: 0=a 1=b 2=c

#ifndef AXIS_ROT		// no rotation of shape axes
        axes[0] = 0; axes[1] = 1; axes[2] = 2;
        min_dim[0] = min_a; min_dim[1] = min_b; min_dim[2] = min_c;
        max_dim[0] = max_a; max_dim[1] = max_b; max_dim[2] = max_c;
#else
        // the smallest one is x, other two are y and z
        if(diff_a < diff_b) {
            if(diff_a < diff_c) {
                // x is a
                axes[0] = 0; axes[1] = 1; axes[2] = 2;
                min_dim[0] = min_a; min_dim[1] = min_b; min_dim[2] = min_c;
                max_dim[0] = max_a; max_dim[1] = max_b; max_dim[2] = max_c;
            } else {
                // x is c
                axes[0] = 2; axes[1] = 0; axes[2] = 1;
                min_dim[0] = min_c; min_dim[1] = min_a; min_dim[2] = min_b;
                max_dim[0] = max_c; max_dim[1] = max_a; max_dim[2] = max_b;
            } // if-else
        } else {
            if(diff_b < diff_c) {
                // x is b
                axes[0] = 1; axes[1] = 0; axes[2] = 2;
                min_dim[0] = min_b; min_dim[1] = min_a; min_dim[2] = min_c;
                max_dim[0] = max_b; max_dim[1] = max_a; max_dim[2] = max_c;
            } else {
                // x is c
                axes[0] = 2; axes[1] = 0; axes[2] = 1;
                min_dim[0] = min_c; min_dim[1] = min_a; min_dim[2] = min_b;
                max_dim[0] = max_c; max_dim[1] = max_a; max_dim[2] = max_b;
            } // if-else
        } // if-else
#endif
		if(mpi_rank == 0) {
	        std::cout << "**              Shape size range: ("
						<< min_dim[0] << ", " << min_dim[1] << ", " << min_dim[2] << ") x ("
						<< max_dim[0] << ", " << max_dim[1] << ", "	<< max_dim[2] << ")" << std::endl;
	        /*std::cout << "++ Final shape min point: " << min_dim[0] << ", "
    	                << min_dim[1] << ", " << min_dim[2] << std::endl;
        	std::cout << "++ Final shape max point: " << max_dim[0] << ", "
            	        << max_dim[1] << ", " << max_dim[2] << std::endl;
	        std::cout << "++ Final shape dimensions: "
    	                << fabs(max_dim[0] - min_dim[0]) << " x "
        	            << fabs(max_dim[1] - min_dim[1]) << " x "
            	        << fabs(max_dim[2] - min_dim[2]) << std::endl;*/
		} // if

		return true;
	} // HiGInput::compute_shapedef_minmax()


	unsigned int HiGInput::read_shape_definition(const char* shape_file) {
		ShapeFileType file_type = shape_filetype(shape_file);

		if(file_type == shape_file_shape) {
			return read_shape_file_shape(shape_file);
		} else if(file_type == shape_file_hdf5) {
			return read_shape_file_hdf5(shape_file);
		} else if(file_type == shape_file_object) {
			return read_shape_file_object(shape_file);
		} else {
			std::cerr << "error: unknown shape file extension in '" << shape_file << "'" << std::endl;
			return false;
		} // if-else
		
		return true;
	} // HiGInput::read_shape_definition()


	ShapeFileType HiGInput::shape_filetype(const char* filename) {

		// ... implement this later ...
		return shape_file_hdf5;
		return shape_file_shape;
	} // HiGInput::shape_filetype()


	unsigned int HiGInput::read_shape_file_object(const char* filename) {

		// ... implement this later ...
		std::cerr << "uh-oh: given shape file type reader not yet implemented" << std::endl;
		return false;

		return true;
	} // HiGInput::read_shape_file_shape()


	unsigned int HiGInput::read_shape_file_hdf5(const char* filename) {
		unsigned int num_triangles = 0;
		double* temp_shape_def;

		// ASSUMING THERE ARE 7 ENTRIES FOR EACH TRIANGLE ... IMPROVE/GENERALIZE ...
		HiGFileReader::instance().hdf5_shape_reader(filename, temp_shape_def, num_triangles);
		shape_def_.clear();
#ifndef KERNEL2
		shape_def_.reserve(7 * num_triangles);
		unsigned int max_size = shape_def_.max_size();
		if(7 * num_triangles > max_size) {
			std::cerr << "error: number of triangles more than what can be handled currently ["
						<< max_size / 7 << "]" << std::endl;
			return 0;
		} // if
		for(unsigned int i = 0; i < 7 * num_triangles; ++ i) {
			shape_def_.push_back((float_t)temp_shape_def[i]);
		} // for
#else	// KERNEL2
		shape_def_.reserve(T_PROP_SIZE_ * num_triangles);
		unsigned int max_size = shape_def_.max_size();
		if(T_PROP_SIZE_ * num_triangles > max_size) {
			std::cerr << "error: number of triangles more than what can be handled currently ["
						<< max_size / T_PROP_SIZE_ << "]" << std::endl;
			return 0;
		} // if
		for(unsigned int i = 0, count = 0; i < T_PROP_SIZE_ * num_triangles; ++ i) {
			if((i + 1) % T_PROP_SIZE_ == 0) shape_def_.push_back((float_t)0.0); // for padding
			else shape_def_.push_back((float_t)temp_shape_def[count ++]);
		} // for
#endif // KERNEL2
		return num_triangles;
	} // HiGInput::read_shape_file_shape()


	unsigned int HiGInput::read_shape_file_shape(const char* filename) {
		//std::cerr << "uh-oh: given shape file type reader not yet implemented" << std::endl;
		//return false;

		unsigned int num_triangles = 0;
		HiGFileReader::instance().shape_shape_reader(filename, shape_def_, num_triangles);

		return true;
	} // HiGInput::read_shape_file_shape()


	/* grains */

	bool HiGInput::construct_lattice_vectors() {
		if(structures_.size() < 1) return false;
		for(structure_iterator_t i = structures_.begin(); i != structures_.end(); ++ i) {
			if(!(*i).second.construct_lattice_vectors()) return false;
		} // for

		return true;
	} // HiGInput::construct_lattice_vectors()


	/* layers */

	unsigned int HiGInput::num_layers() const {
		// -1 is substrate layer
		// 0 is vacuum
		// dont count these two
		//if(has_substrate_layer() && layers_.count(0) != 0) return layers_.size() - 2;
		//if(has_substrate_layer() || layers_.count(0) != 0) return layers_.size() - 1;
		//return layers_.size();
		if(has_substrate_layer() && has_vacuum_layer())
			return layers_.size() - 2;
		if((!has_substrate_layer()) && has_vacuum_layer() ||
				has_substrate_layer() && (!has_vacuum_layer())) 
			return layers_.size() - 1;
		return layers_.size();
	} // HiGInput::num_layers()


	bool HiGInput::is_single_layer() const {
		//if(has_substrate_layer() && has_vacuum_layer() && layers_.size() == 3 ||
		//		(!has_substrate_layer()) && has_vacuum_layer() && layers_.size() == 2 ||
		//		has_substrate_layer() && (!has_vacuum_layer()) && layers_.size() == 2 ||
		//		(!has_substrate_layer()) && (!has_vacuum_layer()) && layers_.size == 1)
		//	return true;
		//else
		//	return false;
		return (num_layers() == 1);
	} // HiGInput::is_single_layer()


	bool HiGInput::has_vacuum_layer() const {
		return (layers_.count(0) != 0);
	} // HiGInput::has_vacuum_layer()


	bool HiGInput::has_substrate_layer() const {	// substrate is the one with order -1
													// it extends to infinity
		//layer_iterator_t i = layers_.begin();
		/*auto i = layers_.begin();
		while(i != layers_.end() && (*i).second.order() != -1) i ++;
		if(i == layers_.end()) return false;
		return true;*/
		return (layers_.count(-1) != 0);
	} // HiGInput::substrate_layer()


	Layer& HiGInput::substrate_layer() {	// substrate is the one with order -1
											// it extends to infinity
		//layer_iterator_t i = layers_.begin();
		/*auto i = layers_.begin();
		while(i != layers_.end() && (*i).second.order() != -1) ++ i;
		if(i == layers_.end()) return (*i).second;	// what to send here ... ? NULL doesnt work ...
		return (*i).second;*/
		return layers_[-1];
	} // HiGInput::substrate_layer()


	RefractiveIndex HiGInput::substrate_refindex() {
		if(has_substrate_layer())
			return layers_[-1].refindex();
		else
			return RefractiveIndex(0, 0);
	} // HiGInput::substrate_refindex()


	Layer& HiGInput::single_layer() {	// if there is exactly 1 layer
										// excluding substrate and vacuum
		layer_iterator_t i = layers_.begin();
		if(has_substrate_layer() && has_vacuum_layer() && layers_.size() == 3 ||
				(!has_vacuum_layer()) && has_substrate_layer() && layers_.size() == 2 ||
				(!has_substrate_layer()) && has_vacuum_layer() && layers_.size() == 2 ||
				(!has_substrate_layer()) && (!has_vacuum_layer()) && layers_.size() == 1) {
			while(i != layers_.end() && (*i).first != 0) ++ i;
			if((*i).first == 0) ++i;
			if(i != layers_.end()) return (*i).second;
			else return layers_[0];
		//} else if(has_substrate_layer() && layers_.size() == 2) {
		//	if((*i).first == -1) { ++ i; return (*i).second; }
		//	else return (*i).second;
		} else {
			std::cerr << "error: single_layer requested on multiple layers" << std::endl;
			return (*i).second;
		} // if-else
	} // HiGInput::single_layer()


	int HiGInput::min_layer_order() {
		layer_iterator_t iter = layers_.begin();
		while((*iter).second.order() < 0) ++ iter;			// since layers are sorted
		return (*iter).second.order();
	} // HiGInput::min_layer_order()


/*	struct CompareLayers {
		  bool operator() (const std::pair<std::string, Layer>& a, const std::pair<std::string, Layer>& b) { return (a.order() < b.order()); }
	} comp_layer;
*/

	bool HiGInput::construct_layer_profile() {
		// if there is substrate, its the first one, followed by others
		// insert layer 0, with refindex 0, 0
		
		if(has_vacuum_layer()) return true;

		Layer vacuum;
		vacuum.key(std::string("vacuum"));
		vacuum.refindex_delta(0.0);
		vacuum.refindex_beta(0.0);
		vacuum.thickness(0.0);
		vacuum.order(0);
		vacuum.z_val(0.0);
		layers_[vacuum.order()] = vacuum;
		// compute the cumulative z value
		float_t curr_z = 0.0;
		for(layer_iterator_t i = layers_.begin(); i != layers_.end(); i ++) {
			if((*i).second.order() == -1) { (*i).second.z_val(0.0); continue; }
			if((*i).second.order() == 0) continue;
			curr_z = curr_z - (*i).second.thickness();
			(*i).second.z_val(curr_z);
		} // for

		return true;
	} // HiGInput::construct_layer_profile()


	/* structures */

	unsigned int HiGInput::num_structures() const {
		return structures_.size();
	} // HiGInput::num_structures()


	float_t HiGInput::layers_z_min() {		// check ... (*i).first is layer order, not position ...
		float_t min_val = FLT_MAX;
		for(layer_iterator_t i = layers_.begin(); i != layers_.end(); i ++) {
			if(min_val > (*i).first && (*i).first >= 0) min_val = (*i).first;
		} // for
		return min_val;
	} // HiGInput::layers_z_min()


	bool HiGInput::compute_domain_size(vector3_t& min_vec, vector3_t& max_vec,
										float_t& z_min_0, float_t& z_max_0, int mpi_rank) {
		float_t ma = FLT_MAX;
		float_t mi = -FLT_MAX;

		vector3_t max_l(mi, mi, layers_z_min());
		vector3_t min_l(ma, ma, ma);
		z_max_0 = mi;
		z_min_0 = ma;

		// iterate over structures
		for(structure_iterator_t s = structures_.begin(); s != structures_.end(); s ++) {

			Shape curr_shape = shapes_[(*s).second.grain_shape_key()];
			vector3_t shape_min(0.0, 0.0, 0.0), shape_max(0.0, 0.0, 0.0);
			compute_shape_domain(curr_shape, shape_min, shape_max, mpi_rank);

			/*std::cout << "++ Shape min point: " << shape_min[0] << ", " << shape_min[1]
						<< ", " << shape_min[2] << std::endl;
			std::cout << "++ Shape max point: " << shape_max[0] << ", " << shape_max[1]
						<< ", " << shape_max[2] << std::endl;
			std::cout << "++ Shape dimensions: " << shape_max[0] - shape_min[0] << " x "
						<< shape_max[1] - shape_min[1] << " x "
						<< shape_max[2] - shape_min[2] << std::endl;*/

			/* determine the structure's position in the sample configuration */
			float_t zc_l = layer_origin_z((*s).second);

			vector3_t n = (*s).second.grain_repetition();
			-- n[0]; -- n[1]; -- n[2];

			// get the lattice vectors a, b, c t, and translate
			Lattice *curr_lattice = (Lattice*)(*s).second.lattice();
			vector3_t a = curr_lattice->a();
			vector3_t b = curr_lattice->b();
			vector3_t c = curr_lattice->c();
			vector3_t t = curr_lattice->t();
			vector3_t transvec = (*s).second.grain_transvec();
			t[0] += transvec[0]; t[1] += transvec[1]; t[2] += transvec[2];

			vector3_t a_max(0.0, 0.0, 0.0), b_max(0.0, 0.0, 0.0), c_max(0.0, 0.0, 0.0);
			vector3_t a_min(0.0, 0.0, 0.0), b_min(0.0, 0.0, 0.0), c_min(0.0, 0.0, 0.0);

			// a along x
			if(a[0] > 0) {
				a_max[0] = n[0] * a[0] + transvec[0] + shape_max[0];
				a_min[0] = transvec[0] + shape_min[0];
			} else {
				a_max[0] = transvec[0] + shape_max[0];
				a_min[0] = n[0] * a[0] + transvec[0] + shape_min[0];
			} // if-else
			// a along y
			if(a[1] > 0) {
				a_max[1] = n[0] * a[1] + transvec[1] + shape_max[1];
				a_min[1] = transvec[1] + shape_min[1];
			} else {
				a_max[1] = transvec[1] + shape_max[1];
				a_min[1] = n[0] * a[1] + transvec[1] + shape_min[1];
			} // if-else
			// a along z
			if(a[2] > 0) {
				a_max[2] = n[0] * a[2] + zc_l + shape_max[2];
				a_min[2] = zc_l + shape_min[2];
			} else {
				a_max[2] = zc_l + shape_max[2];
				a_min[2] = n[0] * a[2] + zc_l + shape_min[2];
			} // if-else
			
			// b along x
			if(b[0] > 0) {
				b_max[0] = n[1] * b[0] + transvec[0] + shape_max[0];
				b_min[0] = transvec[0] + shape_min[0];
			} else {
				b_max[0] = transvec[0] + shape_max[0];
				b_min[0] = n[1] * b[0] + transvec[0] + shape_min[0];
			} // if-else
			// b along y
			if(b[1] > 0) {
				b_max[1] = n[1] * b[1] + transvec[1] + shape_max[1];
				b_min[1] = transvec[1] + shape_min[1];
			} else {
				b_max[1] = transvec[1] + shape_max[1];
				b_min[1] = n[1] * b[1] + transvec[1] + shape_min[1];
			} // if-else
			// b along z
			if(b[2] > 0) {
				b_max[2] = n[1] * b[2] + zc_l + shape_max[2];
				b_min[2] = zc_l + shape_min[2];
			} else {
				b_max[2] = zc_l + shape_max[2];
				b_min[2] = n[1] * b[2] + zc_l + shape_min[2];
			} // if-else
			
			// c along x
			if(c[0] > 0) {
				c_max[0] = n[2] * c[0] + transvec[0] + shape_max[0];
				c_min[0] = transvec[0] + shape_min[0];
			} else {
				c_max[0] = transvec[0] + shape_max[0];
				c_min[0] = n[2] * c[0] + transvec[0] + shape_min[0];
			} // if-else
			// c along y
			if(c[1] > 0) {
				c_max[1] = n[2] * c[1] + transvec[1] + shape_max[1];
				c_min[1] = transvec[1] + shape_min[1];
			} else {
				c_max[1] = transvec[1] + shape_max[1];
				c_min[1] = n[2] * c[1] + transvec[1] + shape_min[1];
			} // if-else
			// c along z
			if(c[2] > 0) {
				c_max[2] = n[2] * c[2] + zc_l + shape_max[2];
				c_min[2] = zc_l + shape_min[2];
			} else {
				c_max[2] = zc_l + shape_max[2];
				c_min[2] = n[2] * c[2] + zc_l + shape_min[2];
			} // if-else

			vector3_t d_min, d_max;
			d_min[0] = t[0] + min(a_min[0], b_min[0], c_min[0]);
			d_max[0] = t[0] + max(a_max[0], b_max[0], c_max[0]);
			d_min[1] = t[1] + min(a_min[1], b_min[1], c_min[1]);
			d_max[1] = t[1] + max(a_max[1], b_max[1], c_max[1]);
			d_min[2] = t[2] + min(a_min[2], b_min[2], c_min[2]);
			d_max[2] = t[2] + max(a_max[2], b_max[2], c_max[2]);

			// case with structure elements in vacuum
			Layer curr_layer = layers_[layer_key_map_[(*s).second.grain_layer_key()]];
			if(curr_layer.order() == 0) {
				z_min_0 = (d_min[2] < z_min_0) ? d_min[2] : z_min_0;
				z_max_0 = (d_max[2] > max_l[2]) ? d_max[2] : z_max_0;
			} // if

			// compute values of min_l and max_l
			max_l[0] = max(max_l[0], d_max[0]);
			max_l[1] = max(max_l[1], d_max[1]);
			max_l[2] = max(max_l[2], d_max[2]);
			min_l[0] = min(min_l[0], d_min[0]);
			min_l[1] = min(min_l[1], d_min[1]);
			min_l[2] = min(min_l[2], d_min[2]);
		} // for

		max_vec[0] = max_l[0];
		max_vec[1] = max_l[1];
		max_vec[2] = max_l[2];
		min_vec[0] = min_l[0];
		min_vec[1] = min_l[1];
		min_vec[2] = min_l[2];

		z_min_0 = min(min_l[2], z_min_0);
		z_max_0 = max(max_l[2], z_max_0);

		return true;
	} // HiGInput::construct_domain_size()


	/** print functions for testing only
	 */


	void HiGInput::print_all() {
		std::cout << "HipGISAXS Inputs: " << std::endl;
		print_shapes();
		print_layers();
		print_structures();
		print_scattering_params();
		print_detector_params();
		print_compute_params();
	} // HiGInput::print_all()


	void HiGInput::print_shapes() {
		std::cout << "Shapes:" << std::endl;
		for(shape_iterator_t i = shapes_.begin(); i != shapes_.end(); i ++) {
			(*i).second.print();
		} // for
	} // HiGInput::print_shapes()


	void HiGInput::print_layers() {
		std::cout << "Layers:" << std::endl;
		for(layer_iterator_t i = layers_.begin(); i != layers_.end(); i ++) {
			(*i).second.print();
		} // for
	} // HiGInput::print_layers()


	void HiGInput::print_structures() {
		std::cout << "Structures:" << std::endl;
		for(structure_iterator_t i = structures_.begin(); i != structures_.end(); i ++) {
			(*i).second.print();
		} // for
	} // HiGInput::print_structures()


	void HiGInput::print_scattering_params() {
		scattering_.print();
	} // HiGInput::print_scattering_params()


	void HiGInput::print_detector_params() {
		detector_.print();
	} // HiGInput::print_detector_params()


	void HiGInput::print_compute_params() {
		compute_.print();
	} // HiGInput::print_compute_params()

} // namespace hig

