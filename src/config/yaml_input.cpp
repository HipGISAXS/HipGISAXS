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

  void YAMLInput::init() {
    shapes_.clear();
    layers_.clear();
    structures_.clear();
    scattering_.init();
    detector_.init();
    compute_.init();
    struct_in_layer_ = false;

    shape_def_.clear();

    analysis_algos_.clear();
    param_key_map_.clear();
    param_space_key_map_.clear();
    param_data_key_map_.clear();

    curr_fit_param_.clear();
    curr_fit_algo_.clear();
    curr_fit_algo_param_.clear();

    // temp ...
    reference_data_.push_back(*(new FitReferenceData()));
    reference_data_set_ = false;
  } // init();


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

  void YAMLInput::get_shapes(){
    YAML::Node shapes = config_["shapes"];
    num_shapes_ = shapes.size();
    for (int i = 0; i < num_shapes; i++){
      curr_shape_.clear();
      YAML::Node curr_shape = shapes[i];
   
      // get key and shape
      std::string key = curr_shape["key"].as<std::string>();
      std::string name = curr_shape["name"].as<std::string>();
      YAML::Node params = curr_shape["params"];
      if (params.Type() != YAML::NodeType::Sequence){
        std::cerr << "Error: Ill formed YAML Node at shape params." << std::endl;
        return false;
      }
      for (int j = 0; j < params.size(); j ++){
        curr_shape_param_.clear();
        std::string type = params["type"].as<std::sting>();
        curr_shape_param_.type(TokenMapper::instance().get_shapeparam_token(type));
        curr_shape_param_.type_name(type);
  
        // get min value        
        real_t minval = params["min"].as<real_t>();
        curr_shape_param_.min(minval);

        // get optional values if present
        if (params["max"]){
          real_t maxval = params["max"].as<real_t>();
          curr_shape_param_.max(num);
        } else { curr_shape_param_.max(maxval); }

        // get parameters if the shape is a distribution
        std::string stat = nullptr;
        if (params["stat"]){
          stat = params["stat"].as<std::string>(); 
          curr_shape_param_.stat(TokenMapper::instance().get_stattype_token(stat); 
        }
        if (params["mean"])
          curr_shape_param_.p1(params["mean"].as<real_t>());
        if (params["stddev"])
          curr_shape_param_.p2(params["stddev"].as<real_t>());
        if (params["nvalues"])
          curr_shape_param_.nvalues(params["nvalues"].as<int>());
      }
      curr_shape_param_.set();
      curr_shape_.insert_param(curr_shape_param_.type_name(), curr_shape_param_);
    }
  }

  void YAMLInput::get_layers(){
    bool found_substr = false;
    YAML::Node layers config_["layers"];
    num_layers_ = layers.size();
    for (int i = 0; i < num_layers_; i++){
      curr_layer_.clear();
      YAML::Node curr_layer = layers[i];
      std::string key = curr_layer["key"].as<std::string>();
      if key.compare("substr") found_substr = true; 
      curr_layer_.order(curr_layer["order"].as<int>());
      if(curr_layer["thickness"]){
        curr_layer_.thickness(curr_layer["thickness"].as<real_t>());
      }
      curr_layer_.refindex_delta(curr_layer_["thickness"]["delta"].as<real_t>());
      curr_layer_.refindex_beta(curr_layer_["thickness"]["beta"].as<real_t>());
    } //for i (layers)
  }

  void YAMLInput::get_instrumentation(){
    YAML::Node node = config_["instrumentation"];
    // scattering parameters
    YAML::Node scattering = node["scattering"];
    scattering_.alphai_min(scattering["min"].as<real_t>());
    if (scattering["max"]){
      scattering_.alphai_max(scattering["max"].as<real_t>()); 
      scattering_.alphai_step(scattering["step"].as<real_t>());
    }
  }

  void YAMLInput::compute_params(){
    // extract values 
    YAML::Node node = config_["computation"];
    comoute_.pathprefix(node["pathprefix"].as<std::string>());
    compute_.runname(node["runname"].as<std::string>());
    compute_.method(node["method"].as<std::string>());
    YAML::node output = node["outputregion"];
    compute_.output_region_type(output["type"].as<std::srting>());
    compute_.output_region_minpoint(output["minpoint"].as<vector2_t>());
    compute_.output_region_maxpoint(output["maxpoint"].as<vector2_t>());
    compute_.resolution(node["resolution"].as<vector2_t>());
  }
    
  /* extract structures */
  void::YAMLInput::get_structures(){
    YAML::Node structs = config_["structures"];
    num_strcts = structs.size();
    for (int i = 0; i < structs.size(); i++){
      curr_structure_.clear();
      YAML::Node curr_struct = structs[i];

      /* structure two members as sub-structures */
      std::string key = curr_struct["key"].as<std::string>());
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
      switch (grain["lattice"]){
        case YAML::NodeType::Scalar:
          curr_structure_.grain_lattice(
            TokenMapper::instance().get_lattice_type(
              grain["lattice"].as<std::string>()));
          if (grain["hkl"])
              curr_structure_.lattice_hkl(grain["hkl"].as<std::string>());
          break;
        case YAML::NodeType::Map:
          curr_structure_.lattice_abc_set(true);
          curr_structure_.lattic_vec_a(grain["lattice"]["a"].as<vector3_t>()); 
          curr_structure_.lattic_vec_b(grain["lattice"]["b"].as<vector3_t>()); 
          curr_structure_.lattic_vec_c(grain["lattice"]["c"].as<vector3_t>()); 
          break;
        default:
          curr_structure_.grain_lattice(
            TokenMapper::instance().get_lattice_type("cubic"));
          break;
      } // switch

      // Grain size (repetition)
      real_t t1;
      vector3_t rep; // recast as int in the setter
      switch (grain["repetition"]){
        case YAML::NodeType::Scalar:
          t1 = grain["repetition"].as<real_t>();
          curr_structure_.grain_repetition(vector3_t(t1, t1, t1));
          break;
        case YAML::NodeType::Sequence:
          curr_structure_.grain_repetition(grain["repetition"].as<vector3_t>());
          break;
        case YAML::NodeType::Map:
          rep[0] = grain["repetition"]["a"]["mean"].as<real_t>(); 
          rep[1] = grain["repetition"]["b"]["mean"].as<real_t>(); 
          rep[2] = grain["repetition"]["c"]["mean"].as<real_t>(); 
          curr_structure_.grain_repetition(rep);

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
          break;
        default:
          std::cerr << "error: ill-formed repetition input" << std::endl;
          return false;
      } 

      // Grain Scaling
      vector3_t mean, stdev, nvalues;
      switch (grain["scaling"]){
        case YAML::NodeType::Scalar:
          t1 = grain["scaling"].as<real_t>();
          curr_structure_.grain_scaling_mean(vector3_t(t1, t1, t1);
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
          if (stdev[1] = grain["scaling"]["b"]["stddev"])
            stdev[1] = grain["scaling"]["b"]["stddev"].as<real_t>();
          else 
            stdev[1] 0.;
          if (stdev[2] = grain["scaling"]["c"]["stddev"])
            stdev[2] = grain["scaling"]["c"]["stddev"].as<real_t>();
          else 
            stdev[2] = 0.;
          curr_structure_.grain_scaling_stddev(stdev);

          // nvalues is optional with 1 as default
          if (grain["scaling"]["a"]["nvalues"])
            nvalues[0] = grain["scaling"]["a"]["nvalues"].as<real_t>();
          else
            nvalues[0] = 1;
          if (grain["scaling"]["b"]["nvalues"])
            nvalues[1] = grain["scaling"]["b"]["nvalues"].as<real_t>();
          else
            nvalues[1] = 1;
          if (grain["scaling"]["c"]["nvalues"])
            nvalues[2] = grain["scaling"]["c"]["nvalues"].as<real_t>();
          else
            nvalues[2] = 1;
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
            curr_structure_.grain_orientation_rot1_axis(orientations["axis"].as<std::string()>());
            curr_structure_.grain_orientation_rot1_angles(orientations["angles"].as<vector2_t>());
            if (orientations["rot1"]["mean"])
              curr_structure_.grain_orientation_rot1_anglesmean(orientations["rot1"]["mean"].as<real_t>());
            if (orientations["rot1"]["stddev"])
              curr_structure_.grain_orientation_rot1_anglessd(orientations["rot1"]["stddev"].as<real_t>());
          } // rot1
          if (orientations["rot2"]){
            curr_structure_.grain_orientation_rot2_axis(orientations["axis"].as<std::string()>());
            curr_structure_.grain_orientation_rot2_angles(orientations["angles"].as<vector2_t>());
            if (orientations["rot2"]["mean"])
              curr_structure_.grain_orientation_rot2_anglesmean(orientations["rot2"]["mean"].as<real_t>());
            if (orientations["rot2"]["stddev"])
              curr_structure_.grain_orientation_rot2_anglessd(orientations["rot2"]["stddev"].as<real_t>());
          } // rot2
          if (orientations["rot3"]){
            curr_structure_.grain_orientation_rot3_axis(orientations["axis"].as<std::string()>());
            curr_structure_.grain_orientation_rot3_angles(orientations["angles"].as<vector2_t>());
            if (orientations["rot3"]["mean"])
              curr_structure_.grain_orientation_rot3_anglesmean(orientations["rot3"]["mean"].as<real_t>());
            if (orientations["rot3"]["stddev"])
              curr_structure_.grain_orientation_rot3_anglessd(orientations["rot3"]["stddev"].as<real_t>());
          } // rot3
        } //orientation
      } // if curr_struct["ensemble"]
    } // for i (structs)
  }

  /**
   * input accessor and modifier functions
   */


  /* shapes */
  unsigned int YAMLInput::read_shape_definition(const char* shape_file) {
    ShapeFileType file_type = shape_filetype(shape_file);

    if(file_type == shape_file_data) {
      return read_shape_file_data(shape_file);
    } else if(file_type == shape_file_hdf5) {
      #ifdef USE_PARALLEL_HDF5
        return read_shape_file_hdf5(shape_file);
      #else
        std::cerr << "error: use of parallel hdf5 format has not been enabled in your installation. "
                  << "Please reinstal with the support enabled." << std::endl;
        return false;
      #endif
    } else if(file_type == shape_file_object) {
      return read_shape_file_object(shape_file);
    } else {
      std::cerr << "error: unknown shape file extension in '" << shape_file << "'" << std::endl;
      return false;
    } // if-else
    
    return true;
  } // YAMLInput::read_shape_definition()


  ShapeFileType YAMLInput::shape_filetype(const char* filename) {

    std::istringstream file(filename);
    std::string s;
    while(std::getline(file, s, '.'));// std::cout << s << std::endl;
    if(s.compare("") == 0) return shape_file_null;
    if(s.compare("dat") == 0) return shape_file_data;
    if(s.compare("Dat") == 0) return shape_file_data;
    if(s.compare("DAT") == 0) return shape_file_data;
    if(s.compare("hd5") == 0) return shape_file_hdf5;
    if(s.compare("Hd5") == 0) return shape_file_hdf5;
    if(s.compare("HD5") == 0) return shape_file_hdf5;
    if(s.compare("hdf5") == 0) return shape_file_hdf5;
    if(s.compare("Hdf5") == 0) return shape_file_hdf5;
    if(s.compare("HDF5") == 0) return shape_file_hdf5;
    if(s.compare("obj") == 0) return shape_file_object;
    if(s.compare("Obj") == 0) return shape_file_object;
    if(s.compare("OBJ") == 0) return shape_file_object;
    return shape_file_error;
  } // YAMLInput::shape_filetype()


  unsigned int YAMLInput::read_shape_file_object(const char* filename) {
    //std::cerr << "uh-oh: given shape file type reader not yet implemented" << std::endl;
    //return false;
    unsigned int num_triangles = 0;
    double* temp_shape_def;

    // ASSUMING THERE ARE 7 ENTRIES FOR EACH TRIANGLE ... IMPROVE/GENERALIZE ...
    HiGFileReader::instance().object_shape_reader(filename, temp_shape_def, num_triangles);
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
      shape_def_.push_back((real_t)temp_shape_def[i]);
    } // for
#else  // KERNEL2
    shape_def_.reserve(T_PROP_SIZE_ * num_triangles);
    unsigned int max_size = shape_def_.max_size();
    if(T_PROP_SIZE_ * num_triangles > max_size) {
      std::cerr << "error: number of triangles more than what can be handled currently ["
            << max_size / T_PROP_SIZE_ << "]" << std::endl;
      return 0;
    } // if
    for(unsigned int i = 0, count = 0; i < T_PROP_SIZE_ * num_triangles; ++ i) {
      if((i + 1) % T_PROP_SIZE_ == 0) shape_def_.push_back((real_t)0.0); // for padding
      else shape_def_.push_back((real_t)temp_shape_def[count ++]);
    } // for
#endif // KERNEL2
    return num_triangles;
  } // YAMLInput::read_shape_file_object()


  #ifdef USE_PARALLEL_HDF5
  unsigned int YAMLInput::read_shape_file_hdf5(const char* filename) {
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
      shape_def_.push_back((real_t)temp_shape_def[i]);
    } // for
    #else  // KERNEL2
    shape_def_.reserve(T_PROP_SIZE_ * num_triangles);
    unsigned int max_size = shape_def_.max_size();
    if(T_PROP_SIZE_ * num_triangles > max_size) {
      std::cerr << "error: number of triangles more than what can be handled currently ["
            << max_size / T_PROP_SIZE_ << "]" << std::endl;
      return 0;
    } // if
    for(unsigned int i = 0, count = 0; i < T_PROP_SIZE_ * num_triangles; ++ i) {
      if((i + 1) % T_PROP_SIZE_ == 0) shape_def_.push_back((real_t)0.0); // for padding
      else shape_def_.push_back((real_t)temp_shape_def[count ++]);
    } // for
    #endif // KERNEL2
    return num_triangles;
  } // YAMLInput::read_shape_file_hdf5()
  #endif // USE_PARALLEL_HDF5


  unsigned int YAMLInput::read_shape_file_data(const char* filename) {
    //std::cerr << "uh-oh: given shape file type reader not yet implemented" << std::endl;
    //return false;

    unsigned int num_triangles = 0;
    HiGFileReader::instance().shape_shape_reader(filename, shape_def_, num_triangles);

    return true;
  } // YAMLInput::read_shape_file_data()


  /* grains */

  bool YAMLInput::construct_lattice_vectors() {
    if(structures_.size() < 1) return false;
    for(structure_iterator_t i = structures_.begin(); i != structures_.end(); ++ i) {
      if(!(*i).second.construct_lattice_vectors()) return false;
    } // for

    return true;
  } // YAMLInput::construct_lattice_vectors()


  /* layers */

  Layer& YAMLInput::substrate_layer() {  // substrate is the one with order -1
    return layers_[-1];
  } // YAMLInput::substrate_layer()


  bool YAMLInput::construct_layer_profile() {
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
    real_t curr_z = 0.0;
    for(layer_iterator_t i = layers_.begin(); i != layers_.end(); i ++) {
      if((*i).second.order() == -1) { (*i).second.z_val(0.0); continue; }
      if((*i).second.order() == 0) continue;
      curr_z = curr_z - (*i).second.thickness();
      (*i).second.z_val(curr_z);
    } // for

    return true;
  } // YAMLInput::construct_layer_profile()


  /* structures */

  unsigned int YAMLInput::num_structures() const {
    return structures_.size();
  } // YAMLInput::num_structures()


  real_t YAMLInput::layers_z_min() {    // check ... (*i).first is layer order, not position ...
    real_t min_val = FLT_MAX;
    for(layer_iterator_t i = layers_.begin(); i != layers_.end(); i ++) {
      if(min_val > (*i).first && (*i).first >= 0) min_val = (*i).first;
    } // for
    return min_val;
  } // YAMLInput::layers_z_min()


  bool YAMLInput::compute_domain_size(vector3_t& min_vec, vector3_t& max_vec,
                    real_t& z_min_0, real_t& z_max_0) {
    real_t ma = FLT_MAX;
    real_t mi = -FLT_MAX;

    vector3_t max_l(mi, mi, layers_z_min());
    vector3_t min_l(ma, ma, ma);
    z_max_0 = mi;
    z_min_0 = ma;

    // iterate over structures
    for(structure_iterator_t s = structures_.begin(); s != structures_.end(); s ++) {

      Unitcell *curr_unitcell = &unitcells_.at((*s).second.grain_unitcell_key());
      vector3_t element_min(REAL_ZERO_, REAL_ZERO_, REAL_ZERO_);
      vector3_t element_max(REAL_ZERO_, REAL_ZERO_, REAL_ZERO_);
      for(Unitcell::element_iterator_t e = curr_unitcell->element_begin();
          e != curr_unitcell->element_end(); ++ e) {

        Shape curr_shape = shapes_[(*e).first];
        vector3_t shape_min(0.0, 0.0, 0.0), shape_max(0.0, 0.0, 0.0);
        compute_shape_domain(curr_shape, shape_min, shape_max);

        element_min[0] = std::min(element_min[0], shape_min[0]);
        element_min[1] = std::min(element_min[1], shape_min[1]);
        element_min[2] = std::min(element_min[2], shape_min[2]);
        element_max[0] = std::max(element_max[0], shape_max[0]);
        element_max[1] = std::max(element_max[1], shape_max[1]);
        element_max[2] = std::max(element_max[2], shape_max[2]);
      } // for elements

      /* determine the structure's position in the sample configuration */
      real_t zc_l = layer_origin_z((*s).second);

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
        a_max[0] = n[0] * a[0] + transvec[0] + element_max[0];
        a_min[0] = transvec[0] + element_min[0];
      } else {
        a_max[0] = transvec[0] + element_max[0];
        a_min[0] = n[0] * a[0] + transvec[0] + element_min[0];
      } // if-else
      // a along y
      if(a[1] > 0) {
        a_max[1] = n[0] * a[1] + transvec[1] + element_max[1];
        a_min[1] = transvec[1] + element_min[1];
      } else {
        a_max[1] = transvec[1] + element_max[1];
        a_min[1] = n[0] * a[1] + transvec[1] + element_min[1];
      } // if-else
      // a along z
      if(a[2] > 0) {
        a_max[2] = n[0] * a[2] + zc_l + element_max[2];
        a_min[2] = zc_l + element_min[2];
      } else {
        a_max[2] = zc_l + element_max[2];
        a_min[2] = n[0] * a[2] + zc_l + element_min[2];
      } // if-else
      
      // b along x
      if(b[0] > 0) {
        b_max[0] = n[1] * b[0] + transvec[0] + element_max[0];
        b_min[0] = transvec[0] + element_min[0];
      } else {
        b_max[0] = transvec[0] + element_max[0];
        b_min[0] = n[1] * b[0] + transvec[0] + element_min[0];
      } // if-else
      // b along y
      if(b[1] > 0) {
        b_max[1] = n[1] * b[1] + transvec[1] + element_max[1];
        b_min[1] = transvec[1] + element_min[1];
      } else {
        b_max[1] = transvec[1] + element_max[1];
        b_min[1] = n[1] * b[1] + transvec[1] + element_min[1];
      } // if-else
      // b along z
      if(b[2] > 0) {
        b_max[2] = n[1] * b[2] + zc_l + element_max[2];
        b_min[2] = zc_l + element_min[2];
      } else {
        b_max[2] = zc_l + element_max[2];
        b_min[2] = n[1] * b[2] + zc_l + element_min[2];
      } // if-else
      
      // c along x
      if(c[0] > 0) {
        c_max[0] = n[2] * c[0] + transvec[0] + element_max[0];
        c_min[0] = transvec[0] + element_min[0];
      } else {
        c_max[0] = transvec[0] + element_max[0];
        c_min[0] = n[2] * c[0] + transvec[0] + element_min[0];
      } // if-else
      // c along y
      if(c[1] > 0) {
        c_max[1] = n[2] * c[1] + transvec[1] + element_max[1];
        c_min[1] = transvec[1] + element_min[1];
      } else {
        c_max[1] = transvec[1] + element_max[1];
        c_min[1] = n[2] * c[1] + transvec[1] + element_min[1];
      } // if-else
      // c along z
      if(c[2] > 0) {
        c_max[2] = n[2] * c[2] + zc_l + element_max[2];
        c_min[2] = zc_l + element_min[2];
      } else {
        c_max[2] = zc_l + element_max[2];
        c_min[2] = n[2] * c[2] + zc_l + element_min[2];
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
    } // for structures

    max_vec[0] = max_l[0];
    max_vec[1] = max_l[1];
    max_vec[2] = max_l[2];
    min_vec[0] = min_l[0];
    min_vec[1] = min_l[1];
    min_vec[2] = min_l[2];

    z_min_0 = min(min_l[2], z_min_0);
    z_max_0 = max(max_l[2], z_max_0);

    return true;
  } // YAMLInput::compute_domain_size()


  /**
   * fitting related functions
   */

  bool YAMLInput::update_params(const map_t& params) {
    //print_all();
    for(map_t::const_iterator p = params.begin(); p != params.end(); ++ p) {
      real_t new_val = (*p).second;
      // check if new_val is within the param space
      ParamSpace ps = param_space_key_map_.at((*p).first);  // if not exist, exception!!
      if(new_val < ps.min_ || new_val > ps.max_) {
        std::cerr << "warning: given parameter value out of range space. resetting to limit."
              << std::endl;
        new_val = std::max(std::min(new_val, ps.max_), ps.min_);
      } // if
      std::string param = param_key_map_.at((*p).first);  // if not exist, exception!!
      // get first component from the string
      std::string keyword, rem_param;
      if(!extract_first_keyword(param, keyword, rem_param)) return false;
      // get keyword name and key (if any)
      std::string keyword_name, keyword_key;
      if(!extract_keyword_name_and_key(keyword, keyword_name, keyword_key)) return false;
      std::string rem_param2;
      switch(TokenMapper::instance().get_keyword_token(keyword_name)) {
        case shape_token:
          #ifdef __INTEL_COMPILER
            if(shapes_.count(keyword_key) == 0 ||
                !shapes_[keyword_key].update_param(rem_param, new_val)) {
              std::cerr << "error: failed to update param '" << param << "'" << std::endl;
              return false;
            } // if
          #else
            if(!shapes_.at(keyword_key).update_param(rem_param, new_val)) {
              std::cerr << "error: failed to update param '" << param << "'" << std::endl;
              return false;
            } // if
          #endif // __INTEL_COMPILER
          break;

        case layer_token:
          #ifdef __INTEL_COMPILER
            if(layer_key_map_.count(keyword_key) == 0 ||
                layers_.count(layer_key_map_[keyword_key]) == 0 ||
                !(layers_[layer_key_map_[keyword_key]]).update_param(rem_param, new_val)) {
              std::cerr << "error: failed to update param '" << param << "'" << std::endl;
              return false;
            } // if
          #else
            if(!layers_.at(layer_key_map_.at(keyword_key)).update_param(rem_param, new_val)) {
              std::cerr << "error: failed to update param '" << param << "'" << std::endl;
              return false;
            } // if
          #endif
          break;

        case struct_token:
          #ifdef __INTEL_COMPILER
            if(structures_.count(keyword_key) == 0 ||
                !structures_[keyword_key].update_param(rem_param, new_val)) {
              std::cerr << "error: failed to update param '" << param << "'" << std::endl;
              return false;
            } // if
          #else
            if(!structures_.at(keyword_key).update_param(rem_param, new_val)) {
              std::cerr << "error: failed to update param '" << param << "'" << std::endl;
              return false;
            } // if
          #endif
          break;

        case instrument_token:
          extract_first_keyword(rem_param, keyword, rem_param2);
          switch(TokenMapper::instance().get_keyword_token(keyword)) {
            case instrument_scatter_token:
              if(!scattering_.update_param(rem_param2, new_val)) {
                std::cerr << "error: failed to update param '" << param << "'"
                      << std::endl;
                return false;
              } // if
              break;

            case instrument_detector_token:
              if(!detector_.update_param(rem_param2, new_val)) {
                std::cerr << "error: failed to update param '" << param << "'"
                      << std::endl;
                return false;
              } // if
              break;

            case error_token:
              std::cerr << "error: invalid keyword '" << keyword
                    << "' in parameter variable name '" << param << "'"
                    << std::endl;
              return false;

            default:
              std::cerr << "error: misplaced keyword '" << keyword
                    << "' in parameter variable name '" << param << "'"
                    << std::endl;
              return false;
          } // switch
          break;

        case compute_token:
          if(!compute_.update_param(rem_param, new_val)) {
            std::cerr << "error: failed to update param '" << param << "'" << std::endl;
            return false;
          } // if
          break;

        case error_token:
          std::cerr << "error: invalid keyword '" << keyword_name
                << "' in parameter variable name '" << param << "'"
                << std::endl;
          return false;

        default:
          std::cerr << "error: misplaced keyword '" << keyword_name
                << "' in parameter variable name '" << param << "'"
                << std::endl;
          return false;
      } // switch
    } // for
    return true;
  } // YAMLInput::update_params()


  /** print functions for testing only
   */


  void YAMLInput::print_all() {
    std::cout << "HipGISAXS Inputs: " << std::endl;
    print_shapes();
    print_unitcells();
    print_layers();
    print_structures();
    print_scattering_params();
    print_detector_params();
    print_compute_params();
    print_fit_params();
    print_ref_data();
    print_fit_algos();
  } // YAMLInput::print_all()


  void YAMLInput::print_shapes() {
    std::cout << "Shapes:" << std::endl;
    for(shape_iterator_t i = shapes_.begin(); i != shapes_.end(); i ++) {
      (*i).second.print();
    } // for
  } // YAMLInput::print_shapes()


  void YAMLInput::print_unitcells() {
    std::cout << "Unitcells:" << std::endl;
    for(unitcell_iterator_t u = unitcells_.begin(); u != unitcells_.end(); ++ u) {
      (*u).second.print();
    } // for
  } // YAMLInput::print_unitcells()


  void YAMLInput::print_layers() {
    std::cout << "Layers:" << std::endl;
    for(layer_iterator_t i = layers_.begin(); i != layers_.end(); i ++) {
      (*i).second.print();
    } // for
  } // YAMLInput::print_layers()


  void YAMLInput::print_structures() {
    std::cout << "Structures:" << std::endl;
    for(structure_iterator_t i = structures_.begin(); i != structures_.end(); i ++) {
      (*i).second.print();
    } // for
  } // YAMLInput::print_structures()


  void YAMLInput::print_scattering_params() {
    scattering_.print();
  } // YAMLInput::print_scattering_params()


  void YAMLInput::print_detector_params() {
    detector_.print();
  } // YAMLInput::print_detector_params()


  void YAMLInput::print_compute_params() {
    compute_.print();
  } // YAMLInput::print_compute_params()

  void YAMLInput::print_fit_params() {
    if(param_key_map_.empty()) return;
    std::cout << "Fit Parameters: " << std::endl;
    for(std::map <std::string, std::string>::const_iterator i = param_key_map_.begin();
        i != param_key_map_.end(); ++ i) {
      ParamSpace temp = param_space_key_map_.at((*i).first);
      FitParam temp2 = param_data_key_map_.at((*i).first);
      std::cout << "  " << (*i).first << ": [" << temp.min_ << " " << temp.max_ << "] "
        << temp2.key_ << " " << temp2.variable_ << " " << temp2.init_
        << " (" << (*i).second << ")" << std::endl;
    } // for
  } // YAMLInput::print_fit_params()

  void YAMLInput::print_ref_data() {
    if(!reference_data_set_) return;
    reference_data_[0].print();
  } // YAMLInput::print_ref_data()

  void YAMLInput::print_fit_algos() {
    if(analysis_algos_.empty()) return;
    std::cout << "Analysis Algorithms: " << std::endl;
    for(analysis_algo_list_t::const_iterator i = analysis_algos_.begin();
        i != analysis_algos_.end(); ++ i) {
      (*i).print();
    } // for
  } // YAMLInput::print_fit_algos()

} // namespace hig

