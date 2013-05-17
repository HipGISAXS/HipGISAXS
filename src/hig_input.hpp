/***
  *  $Id: hig_input.hpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: hig_input.hpp
  *  Created: Jun 11, 2012
  *  Modified: Fri 17 May 2013 11:50:18 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _HIG_INPUT_HPP_
#define _HIG_INPUT_HPP_

#include <stack>
#include <vector>
#include <unordered_map>
#include <string>

#include "globals.hpp"
#include "tokens.hpp"
#include "token_mapper.hpp"
#include "shape.hpp"
#include "layer.hpp"
#include "structure.hpp"
#include "inst_scattering.hpp"
#include "inst_detector.hpp"
#include "compute_params.hpp"
#include "read_oo_input.hpp"

namespace hig {

	// later on create a syntax tree out of the input reading
	// for that create a class with generic 'object type' and parent, children pointers
	// ...

	class HiGInput {

		private:
			/*containers */

			shape_list_t shapes_;
			layer_list_t layers_;
			layer_key_t layer_key_map_;
			structure_list_t structures_;
			ScatteringParams scattering_;
			DetectorParams detector_;
			ComputeParams compute_;
			bool struct_in_layer_;

			std::vector<float_t> shape_def_;	/* shape definition from a file */
			// there may be multiple shape files ... do this later ...

			/* helpers */

			Token curr_token_;
			Token past_token_;
			TokenType curr_keyword_;
			TokenType past_keyword_;

			std::stack <TokenType> keyword_stack_;	// for keyword tokens
													// keyword tokens get pushed on '{' and '['
													// and popped on '}' and ']'
			Shape curr_shape_;
			ShapeParam curr_shape_param_;
			Layer curr_layer_;
			Structure curr_structure_;
			std::vector <float_t> curr_vector_;		// to store values in a vector while parsing it


			/**
			 * methods
			 */

			/* singleton */

			HiGInput();
			HiGInput(const HiGInput&);
			HiGInput& operator=(const HiGInput&);

			void init();

			/* setters */

			TokenType get_curr_parent();
			TokenType get_curr_grandparent();

			bool process_curr_keyword();
			bool process_curr_token();

			bool process_number(const float_t&);
			bool process_string(const std::string&);

			unsigned int read_shape_definition(const char* shape_file);
			//unsigned int read_shape_definition(std::string shape_file) {
			//	return read_shape_definition(shape_file.c_str());
			//} // read_shape_definition()


			/* getters */

			ShapeFileType shape_filetype(const char*);

			/* misc */

			inline bool preceeded_by_keyword() {
				return TokenMapper::instance().keyword_token_exists(past_token_.type_);
			} // HiGInput::preceeded_by_keyword()


			/* computers */

			bool compute_shape_domain(Shape&, vector3_t&, vector3_t&, int);
			bool compute_shapedef_minmax(vector3_t&, vector3_t&, int);

			/* iterators */

			template <typename type_t>
			class HiGIterators {
				// TODO ...
			}; // class HiGIterators

			/* testers */

			void print_shapes();
			void print_layers();
			void print_structures();
			void print_scattering_params();
			void print_detector_params();
			void print_compute_params();

		public:
			// TODO: ...
			//typedef HiGIterators <Shape> shape_iterator_t;
			//typedef std::unordered_map <std::string, Structure>::iterator structure_iterator_t;
			//typedef structure_list_t::iterator structure_iterator_t;

			static HiGInput& instance() {
				static HiGInput hig_input;
				return hig_input;
			} // instance()

			bool construct_input_config(char* filename);
			bool construct_lattice_vectors();
			bool construct_layer_profile();

			bool compute_domain_size(vector3_t&, vector3_t&, float_t&, float_t&, int);

			const std::string& path() const { return compute_.pathprefix(); }
			const std::string& runname() const { return compute_.runname(); }

			void photon_energy(float_t& value, std::string& unit) const {
				value = scattering_.photon_energy().value_;
				unit = scattering_.photon_energy().unit_;
			} // photon_energy()

			unsigned int num_layers() const;
			bool is_single_layer() const;
			int min_layer_order();
			bool has_vacuum_layer() const;
			bool has_substrate_layer() const;
			Layer& substrate_layer();		// the one with order -1
			RefractiveIndex substrate_refindex();
			Layer& single_layer();			// if there is exactly 1 layer
											// excluding substrate
			float_t layers_z_min();
			unsigned int num_structures() const;

			float_t scattering_spot_area() const { return scattering_.spot_area_; }
			float_t scattering_min_alpha_i() const { return scattering_.alpha_i_.min_; }
			void scattering_alphai(float_t& min, float_t& max, float_t& step) {
				min = scattering_.alpha_i_.min_;
				max = scattering_.alpha_i_.max_;
				step = scattering_.alpha_i_.step_; }
			void scattering_inplanerot(float_t& min, float_t& max, float_t& step) {
				min = scattering_.inplane_rot_.min_;
				max = scattering_.inplane_rot_.max_;
				step = scattering_.inplane_rot_.step_; }
			void scattering_tilt(float_t& min, float_t& max, float_t& step) {
				min = scattering_.tilt_.min_;
				max = scattering_.tilt_.max_;
				step = scattering_.tilt_.step_; }
			std::string experiment() const { return scattering_.expt_; }
			vector2_t detector_total_pixels() const { return detector_.total_pixels_; }
			vector2_t detector_direct_beam() const { return detector_.direct_beam_; }
			float_t detector_pixel_size() const { return detector_.pixel_size_; }
			float_t detector_sd_distance() const { return detector_.sd_distance_; }
			vector2_t param_output_minpoint() { return compute_.output_region_.minpoint_; }
			vector2_t param_output_maxpoint() { return compute_.output_region_.maxpoint_; }
			OutputRegionType param_output_type() { return compute_.output_region_.type_; }
			vector2_t param_resolution() const { return compute_.resolution_; }
			const std::string& param_pathprefix() const { return compute_.pathprefix_; }
			unsigned int param_nslices() const { return compute_.nslices_; }
			StructCorrelationType param_structcorrelation() const { return compute_.correlation_; }

			Shape* shape(Structure& s) { return &(shapes_[s.grain_shape_key()]); }
			const Lattice* lattice(Structure& s) { return s.lattice(); }
			ShapeName shape_name(Structure& s) { return shapes_[s.grain_shape_key()].name(); }
			float_t  shape_tau(Structure& s) { return shapes_[s.grain_shape_key()].ztilt(); }
			float_t shape_eta(Structure& s) { return shapes_[s.grain_shape_key()].xyrotation(); }
			vector3_t shape_originvec(Structure& s) { return shapes_[s.grain_shape_key()].originvec(); }
			std::string shape_filename(Structure& s) { return shapes_[s.grain_shape_key()].filename(); }
			shape_param_list_t& shape_params(Structure& s) {
				return shapes_[s.grain_shape_key()].param_list(); }
			bool struct_in_layer() { return struct_in_layer_; }

			unsigned int read_shape_file_shape(const char*);
			unsigned int read_shape_file_hdf5(const char*);
			unsigned int read_shape_file_object(const char*);

			int structure_layer_order(Structure& s) { return layer_key_map_[s.grain_layer_key()]; }
			float_t layer_z_val(int order) { return layers_[order].z_val(); }
			float_t layer_z_val_min() {
				layer_iterator_t begin = layers_.begin();
				return (*begin).second.z_val();
			} // layer_z_val()
			float_t layer_origin_z(Structure& s) {
				int order = structure_layer_order(s);
				vector3_t transvec = s.grain_transvec();
				if(order >= 0) {
					float_t layer_z_val = layers_[order].z_val();
					return layer_z_val + transvec[2];
				} else {
					float_t layer_z_val = (*(layers_.begin())).second.z_val();
					return layer_z_val - transvec[2];
				} // if-else
			} // layer_origin_z()

			// implement better iterators for structures, shapes and layers ...
			structure_iterator_t structure_begin() { return structures_.begin(); }
			structure_iterator_t structure_end() { return structures_.end(); }

			/* printing for testing */
			void print_all();

	}; // class HiGInput

} // namespace hig

#endif /* _HIG_INPUT_HPP_ */
