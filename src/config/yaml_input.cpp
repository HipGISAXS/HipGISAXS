/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: yaml_input.cpp
 *  Created:  Friday the 13th (11-13-2015)
 *
 *  Author: Dinesh Kumar <dkumar@lbl.gov>
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

#include <cfloat>
#include <config/yaml_input.hpp>
#include <utils/utilities.hpp>
#include <utils/string_utils.hpp>
#include <common/parameters.hpp>


namespace hig {
  YAMLInput::YAMLInput() {
    shapes_.clear();
    unitcells_.clear();
    layers_.clear();
    structures_.clear();
    scattering_.init();
    detector_.init();
    compute_.init();
  }

  void YAMLInput::init() {
    shapes_.clear();
    layers_.clear();
    structures_.clear();
    scattering_.init();
    detector_.init();
    compute_.init();
  } // init()

  bool YAMLInput::read_input(const char* filename){
    YAML::Node f = YAML::LoadFile(filename);
    config_ = f["HipGisaxsInput"];
    if (config_ == NULL){
      std::cerr << "fatal error: some error happened in opening or reading "
            << "input config file. aborting"
            << std::endl;
      return false;
    } // if
  }

  bool YAMLInput::extract_shapes(){
    YAML::Node shapes = config_["shapes"];
    num_shapes_ = shapes.size();
    for (int i = 0; i < num_shapes_; i++){
      curr_shape_.clear();
      YAML::Node curr_shape = shapes[i];
   
      // get key and shape
      std::string key = curr_shape["key"].as<std::string>();
      std::string name = curr_shape["name"].as<std::string>();
      RefractiveIndex refindex = curr_shape["refindex"].as<RefractiveIndex>();
      YAML::Node params = curr_shape["params"];
      if (params.Type() != YAML::NodeType::Sequence){
        std::cerr << "Error: Ill formed YAML Node at shape params." << std::endl;
        return false;
      }
      for (int j = 0; j < params.size(); j ++){
        curr_shape_param_.clear();
        curr_shape_param_ = params[i].as<ShapeParam>();
      }
      curr_shape_param_.set();
      curr_shape_.insert_param(curr_shape_param_.type_name(), curr_shape_param_);
    }
    return true;
  }

  bool YAMLInput::extract_unitcells(){
    if (config_["unitcells"]){
      YAML::Node unitcells = config_["unitcells"];
      for (int i = 0; i < unitcells.size(); i++){
        YAML::Node unitcell = unitcells[i];
        curr_unitcell_.clear();
        std::string key = unitcell["key"].as<std::string>(); 
        if (!unitcell["elements"]){
          std::cerr << "error: ill-formed unitcell." << std::endl;
          return false;
        }
        YAML::Node elements = unitcell["elements"];
        for (int j = 0; j < elements.size(); j++){
          YAML::Node element = elements[j];
          std::string shape_key = element["shape:key"].as<std::string>();
          curr_element_list_[shape_key] = element["locations"].as<std::vector<vector3_t>>(); 
        } // for j
        curr_unitcell_.element_list(curr_element_list_);
        unitcells_[key] = curr_unitcell_;
      } // for i
      return true;
    }
    return false; // TODO create a defualt unitcell
  }

  bool YAMLInput::extract_layers(){
    bool found_substr = false;
    YAML::Node layers = config_["layers"];
    num_layers_ = layers.size();
    for (int i = 0; i < num_layers_; i++){
      curr_layer_.clear();
      YAML::Node curr_layer = layers[i];
      std::string key = curr_layer["key"].as<std::string>();
      if (key.compare("substr")) found_substr = true; 

      // order : must begin at 1, unless it is substrate
      int order = curr_layer["order"].as<int>();
      if ((!key.compare("substr")) || (order < 1)) {
        std::cerr << "error: order must begin at 1, unless it is substrate." << std::endl;
        return false;
      }
      curr_layer_.order(order);

      // film thickness 
      if(curr_layer["thickness"]){
        curr_layer_.thickness(curr_layer["thickness"].as<real_t>());
      } else {
        if (!key.compare("substr")){
          std::cerr << "error: layer thickness is missing for \"" << key <<"\"" << std::endl;
          return false;
        }
      }

      // refrective index
      if (curr_layer["refindex"]){
        curr_layer_.refindex(curr_layer["refindex"].as<RefractiveIndex>());
      } else {
        std::cerr << "error: refrective index is missing for layer \"" << key <<"\"" << std::endl;
        return false;
      }
      layers_[order] = curr_layer_;
    } //for i (layers)
    if (!found_substr){
      std::cerr << "Warning: substrate not defined. Will atumatically add Si-substrate" << std::endl;
    }
  }

  bool YAMLInput::extract_instrumentation(){
    YAML::Node node = config_["instrumentation"];
    // scattering parameters
    YAML::Node scattering = node["scattering"];
    scattering_.alphai_min(scattering["min"].as<real_t>());
    if (scattering["max"]){
      scattering_.alphai_max(scattering["max"].as<real_t>()); 
      scattering_.alphai_step(scattering["step"].as<real_t>());
    }
    if (!scattering["expt"]) scattering_.experiment(std::string("gisaxs"));
    else scattering_.experiment(scattering["expt"].as<std::string>());
  }

  bool YAMLInput::extract_compute_params(){
    // extract values 
    YAML::Node node = config_["computation"];
    compute_.pathprefix(node["pathprefix"].as<std::string>());
    compute_.runname(node["runname"].as<std::string>());
    compute_.method(node["method"].as<std::string>());
    YAML::Node output = node["outputregion"];
    compute_.output_region_type(TokenMapper::instance().
      get_output_region_type(output["type"].as<std::string>()));
    vectorf_t v = output["minpoint"].as<vectorf_t>();
    compute_.output_region_minpoint(v[0], v[1]);
    v = output["maxpoint"].as<vectorf_t>();
    compute_.output_region_maxpoint(v[0], v[1]);
    vectori_t res = node["resolution"].as<vectori_t>();
    compute_.resolution(res[0], res[1]);
  }
    
  /* extract structures */
  bool YAMLInput::extract_structures(){
    YAML::Node structs = config_["structures"];
    num_structs_ = structs.size();
    for (int i = 0; i < structs.size(); i++){
      curr_structure_.clear();
      YAML::Node curr_struct = structs[i];

      /* structure two members as sub-structures */
      std::string key = curr_struct["key"].as<std::string>();
      curr_structure_.key(key);
      curr_structure_.iratio(curr_struct["iratio"].as<real_t>());
      YAML::Node grain = curr_struct["grain"];
      if (grain["unitcell_key"])
        curr_structure_.grain_unitcell_key(grain["unitcell_key"].as<std::string>());
      else if (grain["shape_key"])
        curr_structure_.grain_shape_key(grain["shape_key"].as<std::string>());
      else {
        std::cerr << "error: grain must contain either a unitcell or a shape key" << std::endl;
        return false;
      }

      // get lattice type or lattice vectors
      switch (grain["lattice"].Type()){
        case YAML::NodeType::Scalar:
          curr_structure_.lattice_type(
            TokenMapper::instance().get_lattice_type(
              grain["lattice"].as<std::string>()));
          if (grain["hkl"])
              curr_structure_.lattice_hkl(grain["hkl"].as<std::string>());
          break;
        case YAML::NodeType::Map:
          curr_structure_.lattice_abc_set(true);
          curr_structure_.lattice_vec_a(grain["lattice"]["a"].as<vector3_t>()); 
          curr_structure_.lattice_vec_b(grain["lattice"]["b"].as<vector3_t>()); 
          curr_structure_.lattice_vec_c(grain["lattice"]["c"].as<vector3_t>()); 
          break;
        default:
          curr_structure_.lattice_type(
            TokenMapper::instance().get_lattice_type("cubic"));
          break;
      } // switch

      // Grain size (repetition)
      real_t t1;
      vector3_t rep; // recast as int in the setter
      std::vector<StatisticType> stat; stat.resize(3);
      switch (grain["repetition"].Type()){
        case YAML::NodeType::Scalar:
          t1 = grain["repetition"].as<real_t>();
          curr_structure_.grain_repetition(vector3_t(t1, t1, t1));
          break;
        case YAML::NodeType::Sequence:
          curr_structure_.grain_repetition(grain["repetition"].as<vector3_t>());
          break;
        case YAML::NodeType::Map:

          // min values
          rep[0] = grain["repetition"]["a"]["min"].as<real_t>();
          rep[1] = grain["repetition"]["b"]["min"].as<real_t>();
          rep[2] = grain["repetition"]["c"]["min"].as<real_t>();
          curr_structure_.grain_repetition_min(rep);

          // max values
          rep[0] = grain["repetition"]["a"]["max"].as<real_t>();
          rep[1] = grain["repetition"]["b"]["max"].as<real_t>();
          rep[2] = grain["repetition"]["c"]["max"].as<real_t>();
          curr_structure_.grain_repetition_max(rep);

          // stat values
          stat[0] = TokenMapper::instance().get_stattype_token(
            grain["repetition"]["a"]["stat"].as<std::string>());
          stat[1] = TokenMapper::instance().get_stattype_token(
            grain["repetition"]["b"]["stat"].as<std::string>());
          stat[2] = TokenMapper::instance().get_stattype_token(
            grain["repetition"]["c"]["stat"].as<std::string>());
          curr_structure_.grain_repetition_stat(stat);
          break;
        default:
          std::cerr << "error: ill-formed repetition input" << std::endl;
          return false;
      } 

      // Grain Scaling
      vector3_t mean, stdev;
      vectori_t nvalues;
      switch (grain["scaling"].Type()){
        case YAML::NodeType::Scalar:
          t1 = grain["scaling"].as<real_t>();
          curr_structure_.grain_scaling_mean(vector3_t(t1, t1, t1));
          break; 
        case YAML::NodeType::Sequence:
          curr_structure_.grain_scaling_mean(grain["scaling"].as<vectorf_t>());
          break;
        case YAML::NodeType::Map:
          // all lattice vectors must have mean 
          mean[0] = grain["scaling"]["a"]["mean"].as<real_t>();
          mean[1] = grain["scaling"]["b"]["mean"].as<real_t>();
          mean[2] = grain["scaling"]["c"]["mean"].as<real_t>();
          curr_structure_.grain_scaling_mean(mean);

          // stdandard deviation is optional, with 0 default
          if (grain["scaling"]["a"]["stddev"])
            stdev[0] = grain["scaling"]["a"]["stddev"].as<real_t>();
          else
            stdev[0] = 0.;
          if (grain["scaling"]["b"]["stddev"])
            stdev[1] = grain["scaling"]["b"]["stddev"].as<real_t>();
          else 
            stdev[1] = 0.;
          if (grain["scaling"]["c"]["stddev"])
            stdev[2] = grain["scaling"]["c"]["stddev"].as<real_t>();
          else 
            stdev[2] = 0.;
          curr_structure_.grain_scaling_stddev(stdev);


          // stat is optional with "gaussian" as defalut
          if (grain["scaling"]["a"]["stat"])
            stat[0] = TokenMapper::instance().get_stattype_token(
              grain["scaling"]["a"]["stat"].as<std::string>());
          else
            stat[0] = stat_gaussian;
          if (grain["scaling"]["b"]["stat"])
            stat[1] = TokenMapper::instance().get_stattype_token(
              grain["scaling"]["b"]["stat"].as<std::string>());
          else
            stat[1] = stat_gaussian;
          if (grain["scaling"]["c"]["stat"])
            stat[1] = TokenMapper::instance().get_stattype_token(
              grain["scaling"]["c"]["stat"].as<std::string>());
          else
            stat[1] = stat_gaussian;

          // nvalues is optional with 1 as default
          if (grain["scaling"]["a"]["nvalues"])
            nvalues[0] = grain["scaling"]["a"]["nvalues"].as<int>();
          else
            nvalues[0] = 1;
          if (grain["scaling"]["b"]["nvalues"])
            nvalues[1] = grain["scaling"]["b"]["nvalues"].as<int>();
          else
            nvalues[1] = 1;
          if (grain["scaling"]["c"]["nvalues"])
            nvalues[2] = grain["scaling"]["c"]["nvalues"].as<int>();
          else
            nvalues[2] = 1;
          curr_structure_.grain_scaling_nsamples(nvalues);
          break;  
        default:
          std::cerr << "error: lattice input in structure " << key << " is ill formed" << std::endl;
          return false;
      } // switch  for lattice vectors

      /* ensemble */ 
      if (curr_struct["ensemble"]){
        YAML::Node ensemble = curr_struct["ensemble"];
        // max grains
        if (ensemble["maxgrains"])
          curr_structure_.ensemble_maxgrains(ensemble["maxgrains"].as<vector3_t>());

        if (ensemble["orientations"]){
          YAML::Node orientations = ensemble["orientations"];
          // get the distribution type
          curr_structure_.ensemble_distribution(orientations["stat"].as<std::string>());          
          if (orientations["rot1"]){
            curr_structure_.grain_orientation_rot1_axis(orientations["axis"].as<char>());
            curr_structure_.grain_orientation_rot1_angles(orientations["angles"].as<vector2_t>());
            if (orientations["rot1"]["mean"])
              curr_structure_.grain_orientation_rot1_anglemean(orientations["rot1"]["mean"].as<real_t>());
            if (orientations["rot1"]["stddev"])
              curr_structure_.grain_orientation_rot1_anglesd(orientations["rot1"]["stddev"].as<real_t>());
          } // rot1
          if (orientations["rot2"]){
            curr_structure_.grain_orientation_rot2_axis(orientations["axis"].as<char>());
            curr_structure_.grain_orientation_rot2_angles(orientations["angles"].as<vector2_t>());
            if (orientations["rot2"]["mean"])
              curr_structure_.grain_orientation_rot2_anglemean(orientations["rot2"]["mean"].as<real_t>());
            if (orientations["rot2"]["stddev"])
              curr_structure_.grain_orientation_rot2_anglesd(orientations["rot2"]["stddev"].as<real_t>());
          } // rot2
          if (orientations["rot3"]){
            curr_structure_.grain_orientation_rot3_axis(orientations["axis"].as<char>());
            curr_structure_.grain_orientation_rot3_angles(orientations["angles"].as<vector2_t>());
            if (orientations["rot3"]["mean"])
              curr_structure_.grain_orientation_rot3_anglemean(orientations["rot3"]["mean"].as<real_t>());
            if (orientations["rot3"]["stddev"])
              curr_structure_.grain_orientation_rot3_anglesd(orientations["rot3"]["stddev"].as<real_t>());
          } // rot3
        } //orientation
      } // if curr_struct["ensemble"]
      structures_[key] = curr_structure_;
    } // for i (structs)
  }

  bool YAMLInput::construct_input_config(const char * input_yaml){

    // open the file
    if (!read_input(input_yaml)) return false;

    // shapes
    if (!extract_shapes()) return false;

    // layers
    if (!extract_layers()) return false;

     // unitcell
    if (!extract_unitcells()) return false;

     // structrues
     if (!extract_structures()) return false;

     // instrumentation (energy etc) 
     if (!extract_instrumentation()) return false;

     // compute parameters
     if (!extract_compute_params()) return false;


     return true;
  }

} // namespace hig

