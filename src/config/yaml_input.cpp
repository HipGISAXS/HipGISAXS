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
    config_ = f["hipGisaxsInput"];
    if (config_ == NULL){
      std::cerr << "fatal error: some error happened in opening or reading "
            << "input config file. aborting"
            << std::endl;
      return false;
    } // if
    return true;
  }

  bool YAMLInput::read_shape_param(const YAML::Node &node, ShapeParam &param){
    /* unpack required values */
    // parameter type
    if (node["type"]) {
      std::string type = node["type"].as<std::string>();
      param.type_name(type);
      param.type(hig::TokenMapper::instance().get_shapeparam_token(type));
    } else {
      return false;
    }

    // min value 
    if (node["min"]) param.min(node["min"].as<real_t>());
    else return false;

    /* optional parameters */
    bool is_dist = false;
    bool is_normal = false;

    // max
    if (node["max"]) {
      param.max(node["max"].as<real_t>());
      is_dist = true;
    }

    if (is_dist) {
      if (node["stat"]) {
        std::string stat = node["stat"].as<std::string>();
        param.stat(TokenMapper::instance().get_stattype_token(stat));
        if ((!stat.compare("normal")) && (!stat.compare("gaussian"))) is_normal = true;
      } else {
        std::cerr << "error: stat is not defined for \"" << param.type_name() << "\"" << std::endl;
        return false;
      }

      if (is_normal) {
        if (node["p1"]) param.p1(node["p1"].as<real_t>());
        else {
          std::cerr << "error: mean value not defined for \"" << param.type_name() << "\"" << std::endl;
          return false;
        }
        if (node["p2"]) param.p2(node["p2"].as<real_t>());
        else {
          std::cerr << "error: std deviation value not defined for \"" << param.type_name() << "\"" << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  bool YAMLInput::extract_shapes(){
    YAML::Node shapes = config_["shapes"];
    num_shapes_ = shapes.size();
    for (int i = 0; i < num_shapes_; i++){
      curr_shape_.clear();
      YAML::Node curr_shape = shapes[i];
   
      // get key and shape
      std::string key = curr_shape["key"].as<std::string>();
      curr_shape_.key(key);
      std::string name = curr_shape["name"].as<std::string>();
      curr_shape_.name_str(name);
      curr_shape_.name(TokenMapper::instance().get_shapename_token(name));
      RefractiveIndex refindex = curr_shape["refindex"].as<RefractiveIndex>();
      curr_shape_.refindex(refindex);
      YAML::Node params = curr_shape["params"];
      if (params.Type() != YAML::NodeType::Sequence){
        std::cerr << "Error: Ill formed YAML Node at shape params." << std::endl;
        return false;
      }
      for (int j = 0; j < params.size(); j++) {
        curr_shape_param_.clear();
        if (!read_shape_param(params[j], curr_shape_param_)) return false;
        curr_shape_param_.set();
        curr_shape_.insert_param(curr_shape_param_.type_name(), curr_shape_param_);
      }
      shapes_[key] = curr_shape_;
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
        curr_unitcell_.key(key);
        if (!unitcell["elements"]){
          std::cerr << "error: ill-formed unitcell." << std::endl;
          return false;
        }
        std::string shape_key;
        YAML::Node elements = unitcell["elements"];
        for (int j = 0; j < elements.size(); j++){
          YAML::Node element = elements[j];
          if (element["shape_key"]) shape_key = element["shape_key"].as<std::string>();
          else {
            std::cerr << "error: missing shape_key in unitcell locations." << std::endl;
            return false;
          }
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
      if (!key.compare("substr")) found_substr = true; 

      int order;
      // order : must begin at 1, unless it is substrate
      if (curr_layer["order"]) {
        order = curr_layer["order"].as<int>();
        curr_layer_.order(order);
        if ((key.compare("substr")) && (order < 1)) {
          std::cerr << "error: order must begin at 1, unless it is substrate." << std::endl;
          return false;
        }
      } else {
        std::cerr << "error: order not defined for layer: " << key << std::endl;
        return false;
      }
 
      // film thickness 
      if(curr_layer["thickness"]){
        curr_layer_.thickness(curr_layer["thickness"].as<real_t>());
      } else {
        if (key.compare("substr")){
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
      std::cerr << "error: substrate not defined." << std::endl;
      return false;
    }
    return true;
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
    if (node["path"]) {
      compute_.pathprefix(node["path"].as<std::string>());
    } else {
      compute_.pathprefix(std::string("."));
    }

    if (node["runname"]) {
      compute_.runname(node["runname"].as<std::string>());
    } else {
      compute_.runname(std::string("run_"));
    }
    
    // compute_.method(node["method"].as<std::string>());
    if (node["output"]) {
      YAML::Node output = node["output"];
      compute_.output_region_type(TokenMapper::instance().get_output_region_type(output["type"].as<std::string>()));
      vector2_t v = output["minpoint"].as<vector2_t>();
      compute_.output_region_minpoint(v[0], v[1]);
      v = output["maxpoint"].as<vector2_t>();
      compute_.output_region_maxpoint(v[0], v[1]);
    } else {
      std::cerr << "error: output is ill formed" << std::endl;
      return false;
    }

    if (node["resolution"]) {
      vectori_t res = node["resolution"].as<vectori_t>();
      compute_.resolution(res[0], res[1]);
    } else {
      std::cerr << "error: problem reading resolution." << std::endl;
      return false;
    }

    // incoming angle
    if (node["alphai"]) {
      scattering_.alphai_min(node["alphai"]["min"].as<real_t>());
    } else {
      std::cerr << "error: unable to parse incoming angle" << std::endl;
      return false;
    }

    if (node["alphai"]["max"])
      scattering_.alphai_max(node["alphai"]["max"].as<real_t>());

    if (node["alphai"]["step"])
      scattering_.alphai_step(node["alphai"]["step"].as<real_t>());


    if (node["photon"]){
      scattering_.photon_value(node["photon"]["energy"].as<real_t>());
      scattering_.photon_unit(node["photon"]["unit"].as<std::string>());
    } else {
      std::cerr << "error: unable to read input energy." << std::endl;
      return false;
    }
    return true;
  }
    
  /* extract structures */
  bool YAMLInput::extract_structures(){
    YAML::Node structs = config_["structures"];
    num_structs_ = structs.size();
    for (int i = 0; i < structs.size(); i++){
      curr_structure_.init();
      YAML::Node curr_struct = structs[i];

      /* structure two members as sub-structures */
      std::string key;
      if (curr_struct["key"]){
        key = curr_struct["key"].as<std::string>();
        curr_structure_.key(key);
      } else {
        std::cerr << "error: key is missing from the structcure." << std::endl;
        return false;
      }

      if (curr_struct["iratio"]) 
        curr_structure_.iratio(curr_struct["iratio"].as<real_t>());

      if (curr_struct["transvec"])
        curr_structure_.grain_transvec(curr_struct.as<vector3_t>());

      if (curr_struct["grain"]) {
        YAML::Node grain = curr_struct["grain"];
        if (grain["unitcell_key"])
          curr_structure_.grain_unitcell_key(grain["unitcell_key"].as<std::string>());
        else {
          std::cerr << "error: grain must contain a unitcell" << std::endl;
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
        if (grain["repetition"]) {
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
        } else {
          curr_structure_.grain_repetition(vector3_t(1, 1, 1));
        } 

        // Grain Scaling
        vector3_t mean, stdev;
        vectori_t nvalues;
        if (grain["scaling"]) {
          switch (grain["scaling"].Type()){
            case YAML::NodeType::Scalar:
              t1 = grain["scaling"].as<real_t>();
              curr_structure_.grain_scaling_mean(vector3_t(t1, t1, t1));
              break; 
            case YAML::NodeType::Sequence:
              curr_structure_.grain_scaling_mean(grain["scaling"].as<vector3_t>());
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
              std::cerr << "error: scaling input in structure " << key << " is ill formed" << std::endl;
              return false;
          } // switch  for lattice vectors
          curr_structure_.construct_lattice_vectors();
        } 
      } else {
        std::cerr << "error: grain not defined in structcure: " << key << "." << std::endl;
        return false;
      }

      /* ensemble optional */ 
      if (curr_struct["ensemble"]){
        YAML::Node ensemble = curr_struct["ensemble"];
        // max grains
        if (ensemble["maxgrains"])
          curr_structure_.ensemble_maxgrains(ensemble["maxgrains"].as<vector3_t>());

        if (ensemble["orientations"]){
          YAML::Node orientations = ensemble["orientations"];
          // get the distribution type
          curr_structure_.ensemble_orientation_stat(orientations["stat"].as<std::string>());          
          if (orientations["rot1"]){
            curr_structure_.grain_orientation_rot1_axis(orientations["rot1"]["axis"].as<char>());
            curr_structure_.grain_orientation_rot1_angles(orientations["rot1"]["angles"].as<vector2_t>());
            if (orientations["rot1"]["mean"])
              curr_structure_.grain_orientation_rot1_anglemean(orientations["rot1"]["mean"].as<real_t>());
            if (orientations["rot1"]["stddev"])
              curr_structure_.grain_orientation_rot1_anglesd(orientations["rot1"]["stddev"].as<real_t>());
          } // rot1
          if (orientations["rot2"]){
            curr_structure_.grain_orientation_rot2_axis(orientations["rot2"]["axis"].as<char>());
            curr_structure_.grain_orientation_rot2_angles(orientations["rot2"]["angles"].as<vector2_t>());
            if (orientations["rot2"]["mean"])
              curr_structure_.grain_orientation_rot2_anglemean(orientations["rot2"]["mean"].as<real_t>());
            if (orientations["rot2"]["stddev"])
              curr_structure_.grain_orientation_rot2_anglesd(orientations["rot2"]["stddev"].as<real_t>());
          } // rot2
          if (orientations["rot3"]){
            curr_structure_.grain_orientation_rot3_axis(orientations["rot3"]["axis"].as<char>());
            curr_structure_.grain_orientation_rot3_angles(orientations["rot3"]["angles"].as<vector2_t>());
            if (orientations["rot3"]["mean"])
              curr_structure_.grain_orientation_rot3_anglemean(orientations["rot3"]["mean"].as<real_t>());
            if (orientations["rot3"]["stddev"])
              curr_structure_.grain_orientation_rot3_anglesd(orientations["rot3"]["stddev"].as<real_t>());
          } // rot3
        } //orientation
      } // if curr_struct["ensemble"]
      structures_[key] = curr_structure_;
      curr_structure_.clear();
    } // for i (structs)
    return true;
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
     // if (!extract_instrumentation()) return false;

     // compute parameters
     if (!extract_compute_params()) return false;


     return true;
  }

} // namespace hig

