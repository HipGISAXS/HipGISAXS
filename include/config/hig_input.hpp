/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hig_input.hpp
 *  Created: Jun 11, 2012
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

#ifndef __HIG_INPUT_HPP__
#define __HIG_INPUT_HPP__

#include <stack>
#include <vector>
#include <unordered_map>
#include <string>

#include <common/globals.hpp>
#include <common/constants.hpp>
#include <config/tokens.hpp>
#include <config/token_mapper.hpp>
#include <config/input.hpp>
#include <model/shape.hpp>
#include <model/layer.hpp>
#include <model/unitcell.hpp>
#include <model/structure.hpp>
#include <model/inst_scattering.hpp>
#include <model/inst_detector.hpp>
#include <model/compute_params.hpp>
#include <file/read_oo_input.hpp>

#include <config/temp_helpers.hpp>
#include <model/fitting_params.hpp>

namespace hig {

  // TODO: later on create a syntax tree out of the input reading
  // for that create a class with generic 'object type' and parent, children pointers
  // ...

  class HiGInput : public Input {

    private:
      /*containers */

      std::vector<real_t> shape_def_;  /* shape definition from a file */
      // there may be multiple shape files ... do this later ...

      /* helpers */
      bool struct_in_layer_;

      Token curr_token_;
      Token past_token_;
      TokenType curr_keyword_;
      TokenType past_keyword_;

      std::stack <TokenType> keyword_stack_;  // for keyword tokens
                          // keyword tokens get pushed on '{' and '['
                          // and popped on '}' and ']'
      Shape curr_shape_;
      ShapeParam curr_shape_param_;
      Unitcell curr_unitcell_;
      Unitcell::element_list_t curr_element_list_;
      string_t curr_element_shape_key_;
      Layer curr_layer_;
      Structure curr_structure_;
      std::vector <real_t> curr_vector_;    // to store values in a vector while parsing it


      analysis_algo_list_t analysis_algos_;            // list of algorithms
      std::map <std::string, std::string> param_key_map_;      // maps keys to param strings
      std::map <std::string, ParamSpace> param_space_key_map_;  // maps keys to param space

      // TODO: ...
      std::vector <FitReferenceData> reference_data_;
      bool reference_data_set_;

      /* helpers */
      std::map <std::string, FitParam> param_data_key_map_;  // temporary, to be merged above ...
      FitParam curr_fit_param_;
      AnalysisAlgorithmData curr_fit_algo_;
      AnalysisAlgorithmParamData curr_fit_algo_param_;
      FitReferenceData curr_ref_data_;

      /**
       * methods
       */


      /* setters */

      TokenType get_curr_parent();
      TokenType get_curr_grandparent();

      bool process_curr_keyword();
      bool process_curr_token();

      bool process_number(const real_t&);
      bool process_string(const std::string&);

      unsigned int read_shape_definition(const char* shape_file);
      //unsigned int read_shape_definition(std::string shape_file) {
      //  return read_shape_definition(shape_file.c_str());
      //} // read_shape_definition()


      /* getters */

      ShapeFileType shape_filetype(const char*);

      /* misc */

      inline bool preceeded_by_keyword() {
        return TokenMapper::instance().keyword_token_exists(past_token_.type_);
      } // HiGInput::preceeded_by_keyword()


      /* computers */

      bool compute_shape_domain(Shape&, vector3_t&, vector3_t&);
      bool compute_shapedef_minmax(vector3_t&, vector3_t&);

      /* iterators */

      template <typename type_t>
      class HiGIterators {
        // TODO ...
      }; // class HiGIterators

      /* testers */

      void print_shapes();
      void print_unitcells();
      void print_layers();
      void print_structures();
      void print_scattering_params();
      void print_detector_params();
      void print_compute_params();

      void print_fit_params();
      void print_ref_data();
      void print_fit_algos();

    public:

      HiGInput();
      //HiGInput(const HiGInput&);
      HiGInput& operator=(const HiGInput&);
      ~HiGInput() { }

      void init();
      // TODO: ...
      //typedef HiGIterators <Shape> shape_iterator_t;
      //typedef std::unordered_map <std::string, Structure>::iterator structure_iterator_t;
      //typedef structure_list_t::iterator structure_iterator_t;

      bool construct_input_config(const char* filename);
      bool construct_lattice_vectors();
      bool construct_layer_profile();

      bool compute_domain_size(vector3_t&, vector3_t&, real_t&, real_t&);

      const std::string& path() const { return compute_.pathprefix(); }
      const std::string& runname() const { return compute_.runname(); }
      bool saveff() const { return compute_.saveff(); }
      bool savesf() const { return compute_.savesf(); }

      void photon_energy(real_t& value, std::string& unit) const {
        value = scattering_.photon_energy().value_;
        unit = scattering_.photon_energy().unit_;
      } // photon_energy()

      int num_of_layers() const  { return layers_.size(); }
      unsigned int num_layers() const;
      bool is_single_layer() const;
      int min_layer_order();
      bool has_vacuum_layer() const;
      bool has_substrate_layer() const;
      Layer& substrate_layer();    // the one with order -1
      RefractiveIndex substrate_refindex();
      Layer& single_layer();      // if there is exactly 1 layer
                      // excluding substrate
      real_t layers_z_min();
      unsigned int num_structures() const;

      real_t scattering_spot_area() const { return scattering_.spot_area_; }
      real_t scattering_min_alpha_i() const { return scattering_.alpha_i_.min_; }
      void scattering_alphai(real_t& min, real_t& max, real_t& step) {
        min = scattering_.alpha_i_.min_;
        max = scattering_.alpha_i_.max_;
        step = scattering_.alpha_i_.step_; }
      void scattering_inplanerot(real_t& min, real_t& max, real_t& step) {
        min = scattering_.inplane_rot_.min_;
        max = scattering_.inplane_rot_.max_;
        step = scattering_.inplane_rot_.step_; }
      void scattering_tilt(real_t& min, real_t& max, real_t& step) {
        min = scattering_.tilt_.min_;
        max = scattering_.tilt_.max_;
        step = scattering_.tilt_.step_; }
      std::string experiment() const { return scattering_.expt_; }
      real_t scattering_smearing() const { return scattering_.smearing_; }
      vector2_t detector_total_pixels() const { return detector_.total_pixels_; }
      vector2_t detector_direct_beam() const { return detector_.direct_beam_; }
      real_t detector_pixel_size() const { return detector_.pixel_size_; }
      real_t detector_sd_distance() const { return detector_.sd_distance_; }
      vector2_t param_output_minpoint() { return compute_.output_region_.minpoint_; }
      vector2_t param_output_maxpoint() { return compute_.output_region_.maxpoint_; }
      OutputRegionType param_output_type() { return compute_.output_region_.type_; }
      std::vector<int> param_resolution() const { return compute_.resolution_; }
      const std::string& param_pathprefix() const { return compute_.pathprefix_; }
      unsigned int param_nslices() const { return compute_.nslices_; }
      StructCorrelationType param_structcorrelation() const { return compute_.correlation_; }
      std::string palette() const { return compute_.palette_; }

      const Lattice* lattice(Structure& s) { return s.lattice(); }
//#ifdef __INTEL_COMPILER
//      const Unitcell* unitcell(Structure& s) {
//		  if(unitcells_.count(s.grain_unitcell_key()) > 0) return &(unitcells_[s.grain_unitcell_key()]);
//		  else return nullptr;
//	  } // unitcell()
//#else
      const Unitcell* unitcell(Structure& s) { return &(unitcells_.at(s.grain_unitcell_key())); }
//#endif
      //Shape* shape(Structure& s) { return &(shapes_[s.grain_shape_key()]); }
      //ShapeName shape_name(Structure& s) { return shapes_[s.grain_shape_key()].name(); }
      //vector3_t shape_originvec(Structure& s) { return shapes_[s.grain_shape_key()].originvec(); }
      //std::string shape_filename(Structure& s) { return shapes_[s.grain_shape_key()].filename(); }
      //shape_param_list_t& shape_params(Structure& s) {
      //  return shapes_[s.grain_shape_key()].param_list(); }
      ShapeName shape_name(const string_t& key) { return shapes_[key].name(); }
      real_t shape_zrot(const string_t& key) { return shapes_[key].zrot(); }
      real_t shape_yrot(const string_t& key) { return shapes_[key].yrot(); }
      real_t shape_xrot(const string_t& key) { return shapes_[key].xrot(); }
      vector3_t shape_originvec(const string_t& key) { return shapes_[key].originvec(); }
      std::string shape_filename(const string_t& key) { return shapes_[key].filename(); }
      shape_param_list_t& shape_params(const string_t& key) {
        return shapes_[key].param_list(); }
      bool struct_in_layer() { return struct_in_layer_; }

      unsigned int read_shape_file_data(const char*);
      unsigned int read_shape_file_object(const char*);
      #ifdef USE_PARALLEL_HDF5
      unsigned int read_shape_file_hdf5(const char*);
      #endif

      layer_citerator_t layers_begin() const { return layers_.begin(); }
      layer_citerator_t layers_end() const { return layers_.end(); }

      int structure_layer_order(Structure& s) { return layer_key_map_[s.grain_layer_key()]; }
      real_t layer_z_val(int order) { return layers_[order].z_val(); }
      real_t layer_z_val_min() {
        layer_iterator_t begin = layers_.begin();
        return (*begin).second.z_val();
      } // layer_z_val()
      real_t layer_origin_z(Structure& s) {
        int order = structure_layer_order(s);
        vector3_t transvec = s.grain_transvec();
        if(order >= 0) {
          real_t layer_z_val = layers_[order].z_val();
          return layer_z_val + transvec[2];
        } else {
          real_t layer_z_val = (*(layers_.begin())).second.z_val();
          return layer_z_val - transvec[2];
        } // if-else
      } // layer_origin_z()

      // implement better iterators for structures, shapes and layers ...
      structure_iterator_t structure_begin() { return structures_.begin(); }
      structure_iterator_t structure_end() { return structures_.end(); }

      /* fitting related */
      bool update_params(const map_t&);
      // return list of parameter keys
      std::vector <std::string> fit_param_keys() const {
        std::vector <std::string> key_list;
        for(std::map <std::string, ParamSpace>::const_iterator i = param_space_key_map_.begin();
            i != param_space_key_map_.end(); ++ i)
          key_list.push_back((*i).first);
        return key_list;
      } // get_fit_param_keys()
      // return list of min-max for all parameters
      std::vector <std::pair <real_t, real_t> > fit_param_limits() const {
        std::vector <std::pair <real_t, real_t> > plimits;
        for(std::map <std::string, ParamSpace>::const_iterator i = param_space_key_map_.begin();
            i != param_space_key_map_.end(); ++ i)
          plimits.push_back(std::pair<real_t, real_t>((*i).second.min_, (*i).second.max_));
        return plimits;
      } // get_fit_param_limits()
      // return list of step values for all parameters
      real_vec_t fit_param_step_values() const {
        real_vec_t steps;
        for(std::map <std::string, ParamSpace>::const_iterator i = param_space_key_map_.begin();
            i != param_space_key_map_.end(); ++ i)
          steps.push_back((*i).second.step_);
        return steps;
      } // fit_param_step_values()
      // return mean value of given parameter
      real_t param_space_mean(const std::string& key) {
        return (param_space_key_map_[key].max_ - param_space_key_map_[key].min_) / 2.0 + param_space_key_map_[key].min_;
      } // param_space_mean()
      std::string reference_data_path(int i) const { return reference_data_[i].image_path(); }
      std::string reference_data_mask(int i) const { return reference_data_[i].image_mask(); }
      OutputRegionType reference_region_type(int i) const {
        return reference_data_[i].get_region_type();
      } // reference_region_type()
      real_t reference_region_min_x(int i) const { return reference_data_[i].region_min_x(); }
      real_t reference_region_min_y(int i) const { return reference_data_[i].region_min_y(); }
      real_t reference_region_max_x(int i) const { return reference_data_[i].region_max_x(); }
      real_t reference_region_max_y(int i) const { return reference_data_[i].region_max_y(); }
      int num_analysis_data() const { return 1; }    // temp
      int num_fit_params() const { return param_key_map_.size(); }
      std::vector <real_t> fit_param_init_values() const {
        std::vector<real_t> init_vec;
        std::cout << "Initial Vector: ";
        for(std::map<std::string, FitParam>::const_iterator i = param_data_key_map_.begin();
            i != param_data_key_map_.end(); ++ i) {
          init_vec.push_back((*i).second.init_);
          std::cout << (*i).second.init_ << " ";
        } // for
        std::cout << std::endl;
        return init_vec;
      } // fit_param_init_vector()
      int num_analysis_algos() const { return analysis_algos_.size(); }
      FittingAlgorithmName analysis_algo(int i) const { return analysis_algos_[i].name(); }
      bool analysis_algo_param(int i, const std::string pstr, real_t& val) const {
        return analysis_algos_[i].param(pstr, val);
      } // analysis_algo_param()
      real_t analysis_tolerance(int i) const { return analysis_algos_[i].tolerance(); }
      real_t analysis_regularization(int i) const { return analysis_algos_[i].regularization(); }

      /* printing for testing */
      void print_all();

  }; // class HiGInput

} // namespace hig

#endif /* __HIG_INPUT_HPP__ */
