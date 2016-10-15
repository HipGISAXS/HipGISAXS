/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hig_input.hpp
 *  Created: Jun 11, 2012
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

#ifndef YAML_INPUT__HPP
#define YAML_INPUT__HPP

#include <vector>
#include <unordered_map>
#include <string>

#include <common/globals.hpp>
#include <common/constants.hpp>
#include <config/tokens.hpp>
#include <config/input.hpp>
#include <config/token_mapper.hpp>
#include <model/shape.hpp>
#include <model/layer.hpp>
#include <model/unitcell.hpp>
#include <model/structure.hpp>
#include <model/inst_scattering.hpp>
#include <model/inst_detector.hpp>
#include <model/compute_params.hpp>

#include <config/temp_helpers.hpp>
#include <yaml-cpp/yaml.h>

namespace YAML {

  /* Refrective index */
  template<> struct convert <hig::RefractiveIndex> {
    static bool decode(const Node &node, hig::RefractiveIndex &refindex){
      // delta
      if (node["delta"])
        refindex.delta(node["delta"].as<hig::real_t>());
      else 
        return false;
 
      // beta
      if (node["beta"])
        refindex.beta(node["beta"].as<hig::real_t>());
      else
        return false;
      return true;
    }
  };

  /* vector2_t */
  template<> struct convert <hig::vector2_t> {
    static bool decode(const Node &node, hig::vector2_t & vec) {
      if (node.IsSequence() && (node.size() == 2)){
        vec[0] = node[0].as<hig::real_t>();
        vec[1] = node[1].as<hig::real_t>();
        return true;
      }
      return false;
    }
  };

  /* vector3_t */
  template<> struct convert <hig::vector3_t> {
    static bool decode(const Node &node, hig::vector3_t & vec) {
      if (node.IsSequence() && node.size() == 3){
        vec[0] = node[0].as<hig::real_t>();
        vec[1] = node[1].as<hig::real_t>();
        vec[2] = node[2].as<hig::real_t>();
        return true;
      }
      return false;
    }
  };

  /* integer vector */
  template<> struct convert <std::vector<int> > {
    static bool decode(const Node &node, std::vector<int> & vec){
      if (node.IsSequence()) {
        for (int i = 0; i < node.size(); i++)
          vec.push_back(node[i].as<int>());
      } else {
        std::cerr << "error: trying to convert illegal type to int-vector." << std::endl;
        return false;
      }
      return true;
    }
  };

  /* fundamental of grian orientations */
  template<> struct convert <hig::Rotation> {
    static bool decode(const Node &node, hig::Rotation &rot) {
      if (!node.IsMap()) {
        std::cerr << "error: ill-formed grain orientations." << std::endl;
        return false;
      }

      hig::vector2_t angles;
      // parse rotation
      if (node["axis"]) rot.axis(node["axis"].as<char>());
      else {
        std::cerr << "error: axis is not defined in grain orinentation." << std::endl;
        return false;
      }
      if (node["min"])  angles[0] = node["min"].as<hig::real_t>();
      if (node["max"])  angles[1] = node["max"].as<hig::real_t>();
      else (angles[1] = angles[0]);
      // set angles
      rot.angles(angles);
      if (node["mean"]) rot.angle_mean(node["mean"].as<hig::real_t>());
      if (node["std"]) rot.angle_sd(node["std"].as<hig::real_t>());
      if (node["stat"]) rot.stat(node["stat"].as<std::string>());
      return true;
    }
  };
  /* shape parameter */
  template<> struct convert <hig::ShapeParam> {
    static bool decode(const Node &node, hig::ShapeParam & param) {
      if (node["type"]) {
        std::string type = node["type"].as<std::string>();
        param.type_name(type);
        param.type(hig::TokenMapper::instance().get_shapeparam_token(type));
      } else {
        return false;
      }
      
      // required value
      if (node["min"]) 
        param.min(node["min"].as<hig::real_t>());
      else
        return false;

      /* optional values */
      bool is_dist = false;
      bool is_normal = false;
      // max
      if (node["max"]) {
        param.max(node["max"].as<hig::real_t>());
        is_dist = true;
      }
     
      if (is_dist) {
        // stat
        if (node["stat"]){
          std::string stat = node["stat"].as<std::string>();
          param.stat(hig::TokenMapper::instance().get_stattype_token(stat));
          if (stat.compare("normal") || stat.compare("gaussian")) is_normal = true;
        } else {
          std::cerr << "error: stat not defined for \"" << param.type_name() << "\"" << std::endl;
          return false;
        }

        // number of samples
        if (node["nvalues"])
          param.nvalues(node["nvalues"].as<int>());
        else {
          std::cerr << "error: nvalues not defined for \"" << param.type_name() << "\"" << std::endl;
          return false;
        }

        if (is_normal){
          if (node["p1"])
            param.p1(node["p1"].as<hig::real_t>());
          else {
            std::cerr << "error: mean is not defined for the gaussian distribution" << std::endl;
            return false;
          }
          if (node["p2"]) 
            param.p2(node["p2"].as<hig::real_t>());
          else {
            std::cerr << "error: standard deviation not defined for gaussian distribution" << std::endl;
            return false;
          }
        }
      }
      return true;
    }
  };
}

namespace hig {

  class YAMLInput : public Input {
    private:
      typedef std::vector<int> vectori_t;

      /* YAML */
      YAML::Node config_;
 
      /* helpers */
      int num_shapes_;
      int num_layers_;
      int num_structs_;
      Shape curr_shape_;
      ShapeParam curr_shape_param_;
      Unitcell curr_unitcell_;
      Unitcell::element_list_t curr_element_list_;
      string_t curr_element_shape_key_;
      Layer curr_layer_;
      Structure curr_structure_;



      /**
       * methods
       */
      YAMLInput(const YAMLInput&);
      YAMLInput& operator=(const YAMLInput&);
    
      void init();
      bool read_input(const char *);
      bool read_shape_param(const YAML::Node &, ShapeParam &);
      bool decode_shapes();
      bool decode_layers();
      bool decode_unitcells();
      bool decode_instrumentation();
      bool decode_compute_params();
      bool decode_structures();
      
    public:
    
      YAMLInput();
      bool construct_input_config(const char *);
      bool update_params(map_t & params) { return false; }
      const std::string& path() const { return compute_.pathprefix(); }
      const std::string& runname() const { return compute_.runname(); }

      Shape & shape(std::string key) { return shapes_[key]; }
      Unitcell & unitcell(std::string key) { return unitcells_[key]; }
      const shape_list_t & shapes() const { return shapes_; }
      const layer_list_t & layers() const { return layers_;}
      const unitcell_list_t & unitcells() const { return unitcells_; }
      const structure_list_t & structures() const { return structures_; }
      const ScatteringParams & scattering() const { return scattering_; }
      const DetectorParams & detector() const { return detector_; }
      const ComputeParams & compute() const { return compute_; }
      const FittingParams & fitting() const { return fitting_; }

  }; // class YAMLInput

} // namespace hig

#endif /* YAML_INPUT__HPP */
