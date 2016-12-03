/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: structure.cpp
 *  Created: Jun 12, 2012
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

#include <iostream>
#include <set>

#include <model/structure.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>
#include <numerics/matrix.hpp>
#include <numerics/numeric_utils.hpp>
#include <common/constants.hpp>
#include <model/qgrid.hpp>


namespace hig {


  /** grain scaling methods
   */

  GrainScaling::GrainScaling() {
    init();
  } // GrainScaling::GrainScaling()

  void GrainScaling::init() {
    mean_[0] = mean_[1] = mean_[2] = 1;
    stddev_[0] = stddev_[1] = stddev_[2] = 0;
    for(int i = 0; i < 3; ++ i) {
      dist_.push_back(stat_none);
      nvals_.push_back(40);   // FIXME
    } // for
  } // GrainScaling::init()

  void GrainScaling::clear() {
      mean_[0] = mean_[1] = mean_[2] = 0;
      stddev_[0] = stddev_[1] = stddev_[2] = 0;
      dist_.clear();
      nvals_.clear();
  } // GrainScaling::clear()


  /** lattice functions
   */

  Lattice::Lattice() { ca_ = 1; gamma_ = 0; }
  Lattice::~Lattice() { }


  void Lattice::init() {
    a_[1] = a_[2] = b_[0] = b_[2] = c_[0] = c_[1] = 0.0;
    a_[0] = b_[1] = c_[2] = 1.;
    t_[0] = t_[1] = t_[2] = 0.0;
    abc_set_ = false;
    type_ = lattice_cubic;    // default values
    hkl_ = "100";
    abangle_ = 60;
    caratio_ = 1;
  } // Lattice::init()


  void Lattice::clear() {
    a_[0] = a_[1] = a_[2] = b_[0] = b_[1] = b_[2] = c_[0] = c_[1] = c_[2] = 0.0;
    t_[0] = t_[1] = t_[2] = 0.0;
    abc_set_ = false;
    type_ = lattice_null;
    hkl_.clear();
    abangle_ = 0.0;
    caratio_ = 0.0;
  } // Lattice::clear()


  bool Lattice::construct_vectors(vector3_t scaling) {
    //if(abc_set_) return true;  // a b c are already defined in the input

    real_t sqrt2 = sqrt(2.0);
    real_t sqrt3 = sqrt(3.0);

    if(!abc_set_) {
      switch(type_) {
        case lattice_cubic:          // CUBIC
          if (hkl_ == "100") {
            a_[0] = 1; a_[1] = 0; a_[2] = 0;
            b_[0] = 0; b_[1] = 1; b_[2] = 0;
            c_[0] = 0; c_[1] = 0; c_[2] = 1;
            t_[0] = 0; t_[1] = 0; t_[2] = 0;
          } else if( hkl_ == "110") {
            a_[0] = 1; a_[1] = 0; a_[2] = 0;
            b_[0] = 0.; b_[1] = 1.4142; b_[2] = 0;
            c_[0] = 0.; c_[1] = 0.7071; c_[2] = 0.7071;
            t_[0] = 0; t_[1] = 0; t_[2] = 0;
          } else {
            std::cerr << "error: invalid lattice type and hkl combination" << std::endl;
            return false;
          }
          break;

        case lattice_bcc:          // BCC
          if(hkl_ == "100") {
            a_[0] = 1.; a_[1] = 0.; a_[2] = 0.;
            b_[0] = 0.; b_[1] = 1.; b_[2] = 0.;
            c_[0] = 0.5; c_[1] =0.5; c_[2]= 0.5;
            t_[0] = 0.0; t_[1] = 0; t_[2] = 0.0;
          } else if(hkl_ == "110") {
            a_[0] = 1.0; a_[1] = 0.0; a_[2] = 0;
            b_[0] = 0.5; b_[1] = 0.7071; b_[2] = 0;
            c_[0] = 0.0; c_[1] = 0.7071; c_[2] = 0.7071;
            t_[0] = 0.0; t_[1] = 0; t_[2] = 0.0;
          } else {
            std::cerr << "error: invalid lattice type and hkl combination" << std::endl;
            return false;
          } // if-else
          break;

        case lattice_fcc:          // FCC
          if(hkl_ == "100") {
            a_[0] = 0.0; a_[1] = 0.5; a_[2] = 0.5;
            b_[0] = 0.5; b_[1] = 0.0; b_[2] = 0.5;
            c_[0] = 0.5; c_[1] = 0.5; c_[2] = 0.0;
            t_[0] = 0; t_[1] = 0; t_[2] = 0;
            //t_[0] = 0.5; t_[1] = 0; t_[2] = 0.5;
          } else if(hkl_ == "111") {
            a_[0] = 1.; a_[1] = 0.; a_[2] = 0;
            b_[0] = 0.5; b_[1] = 0.866; b_[2] = -0.;
            c_[0] = 1.; c_[1] = 0.5773; c_[2] = 0.866;
            //t_[0] = 1.0 / sqrt3; t_[1] = 0; t_[2] = 1 / sqrt3;
          } else {
            std::cerr << "error: invalid lattice type and khl combination" << std::endl;
            return false;
          } // if-else
          break;

        case lattice_fco:          // FCO
          a_[0] = 1; a_[1] = 0; a_[2] = 0;
          b_[0] = 0.5; b_[1] = tan(gamma_) / 2.0; b_[2] = 0;
          c_[0] = 0; c_[1] = 0; c_[2] = ca_;
          t_[0] = 0.5; t_[1] = 0; t_[2] = ca_ / 2;
          break;

        case lattice_hcp:          // HCP
          a_[0] = 1; a_[1] = 0; a_[2] = 0;
          b_[0] = 0.5; b_[1] = 0.866; b_[2] = 0;
          c_[0] = 0.5; c_[1] = 0.2887; c_[2] = 0.8165;
          t_[0] = 0; t_[1] = 0; t_[2] = 0;
          break;

        case lattice_hex:          // HEX
          a_[0] = 1; a_[1] = 0; a_[2] = 0;
          b_[0] = 0.5; b_[1] = sqrt3 / 2.0; b_[2] = 0;
          c_[0] = 0; c_[1] = 0; c_[2] = 1;
          t_[0] = 0; t_[1] = 0; t_[2] = 0;
          break;

        case lattice_null:          // means custom vectors are defined
          t_[0] = t_[1] = t_[2] = 0;
          break;

        case lattice_error:
          std::cerr << "error: lattice type already has some previous error" << std::endl;
          return false;

        default:
          std::cerr << "error: lattice type has unknown value" << std::endl;
          return false;
      } // switch
    } // if

    return true;
  } // Lattice::construct_vectors()

  void Lattice::bragg_angles(vector3_t repeats, vector3_t scaling, real_t k0, real_vec_t & angles){
    angles.clear();

    vector3_t a = a_ * scaling[0];
    vector3_t b = b_ * scaling[1];
    vector3_t c = c_ * scaling[1];

    // compute reciprocal lattice
    vector3_t bc = cross(b, c);
    real_t vol = std::abs(dot(a, bc));
    real_t t1 = 2. * PI_ / vol;
    vector3_t mra = cross(b, c) * t1;
    vector3_t mrb = cross(c, a) * t1;
    vector3_t mrc = cross(a, b) * t1;
    vector2_t qmin = QGrid::instance().qmin();
    vector2_t qmax = QGrid::instance().qmax();

    const real_t d_ang = PI_ / 180. * 0.5;
    std::set<real_t> angs;
    for (real_t h = 0; h < repeats[0]; h++){
      for (real_t k = 0; k < repeats[1]; k++){
        for (real_t l = 0; l < repeats[2]; l++){
          vector3_t G = mra * h + mrb * k + mrc * l;
          real_t Gp = sign(G[1]) * std::sqrt(G[0]*G[0] + G[1]*G[1]);
          if ((Gp < qmin[0]) || (Gp > qmax[0])) continue;
          if ((G[2] < qmin[1]) || (G[2] > qmax[1])) continue;
          if (G.norm() > 0 ){
            real_t G_proj = std::sqrt(G[0] * G[0] + G[1] * G[1]);
            if ( G_proj < TINY_ ) continue;
            real_t k_proj = std::sqrt(k0 * k0 - G[2] * G[2]);
            real_t t1 = std::atan2(G[1], G[0]);
            real_t t2 = std::acos((k0 * k0 + G_proj * G_proj - k_proj * k_proj)/(2. * k0 * G_proj));
            real_t a = t1 - t2;
            angs.insert(a);
          }
        }
      }
    }
    for (std::set<real_t>::iterator it = angs.begin(); it != angs.end(); it++){
      angles.push_back(*it);
    }
  }



  /** grainorientations functions
   */

  GrainOrientations::GrainOrientations() { }
  GrainOrientations::~GrainOrientations() { }


  void GrainOrientations::init() {
    stat_ = "single";
    rot1_.axis('x');
    rot1_.angles(0, 0);
    rot2_.axis('y');
    rot2_.angles(0, 0);
    rot3_.axis('z');
    rot3_.angles(0, 0);
    //rot2_.axis('n');
    //rot3_.axis('n');    // for null. all could be same as rot1 too ...
  } // GrainOrientations::init()


  void GrainOrientations::clear() {
    stat_.clear();
    rot1_.axis('n'); rot2_.axis('n'); rot3_.axis('n');
    rot1_.angles(0, 0); rot2_.angles(0, 0); rot3_.angles(0, 0);
  } // GrainOrientations::clear()


  bool GrainOrientations::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    std::string keyword_name, key;
    if(!extract_keyword_name_and_key(keyword, keyword_name, key)) return false;
    std::string keyword2, rem_str2;
    switch(TokenMapper::instance().get_keyword_token(keyword_name)) {
      case struct_ensemble_orient_rot1_token:
        if(key.compare(std::string(1, rot1_.axis())) != 0) {
          std::cerr << "error: invalid axis key given in param '" << str << "'" << std::endl;
          return false;
        } // if
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        if(keyword2.compare("min") == 0 || keyword2.compare("0") == 0) {
          rot1_.angles_min(new_val);
        } else if(keyword2.compare("max") == 0 || keyword2.compare("1") == 0) {
          rot1_.angles_max(new_val);
        } else if(keyword2.compare("anglelocation") == 0) {
          rot1_.angle_location(new_val);
        } else if(keyword2.compare("anglemean") == 0) {
          rot1_.angle_mean(new_val);
        } else if(keyword2.compare("anglescale") == 0) {
          rot1_.angle_scale(new_val);
        } else if(keyword2.compare("anglesd") == 0) {
          rot1_.angle_sd(new_val);
        } else {
          std::cerr << "error: invalid keyword '" << keyword2 << "' in param '"
                << str << "'" << std::endl;
          return false;
        } // if-else
        break;

      case struct_ensemble_orient_rot2_token:
        if(key.compare(std::string(1, rot2_.axis())) != 0) {
          std::cerr << "error: invalid axis key given in param '" << str << "'" << std::endl;
          return false;
        } // if
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        if(keyword2.compare("min") == 0 || keyword2.compare("0") == 0) {
          rot2_.angles_min(new_val);
        } else if(keyword2.compare("max") == 0 || keyword2.compare("1") == 0) {
          rot2_.angles_min(new_val);
        } else if(keyword2.compare("anglelocation") == 0) {
          rot2_.angle_location(new_val);
        } else if(keyword2.compare("anglemean") == 0) {
          rot2_.angle_mean(new_val);
        } else if(keyword2.compare("anglescale") == 0) {
          rot2_.angle_scale(new_val);
        } else if(keyword2.compare("anglesd") == 0) {
          rot2_.angle_sd(new_val);
        } else {
          std::cerr << "error: invalid keyword '" << keyword2 << "' in param '"
                << str << "'" << std::endl;
          return false;
        } // if-else
        break;

      case struct_ensemble_orient_rot3_token:
        if(key.compare(std::string(1, rot3_.axis())) != 0) {
          std::cerr << "error: invalid axis key given in param '" << str << "'" << std::endl;
          return false;
        } // if
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        if(keyword2.compare("min") == 0 || keyword2.compare("0") == 0) {
          rot3_.angles_min(new_val);
        } else if(keyword2.compare("max") == 0 || keyword2.compare("1") == 0) {
          rot3_.angles_max(new_val);
        } else if(keyword2.compare("anglelocation") == 0) {
          rot3_.angle_location(new_val);
        } else if(keyword2.compare("anglemean") == 0) {
          rot3_.angle_mean(new_val);
        } else if(keyword2.compare("anglescale") == 0) {
          rot3_.angle_scale(new_val);
        } else if(keyword2.compare("anglesd") == 0) {
          rot3_.angle_sd(new_val);
        } else {
          std::cerr << "error: invalid keyword '" << keyword2 << "' in param '"
                << str << "'" << std::endl;
          return false;
        } // if-else
        break;

      default:
        std::cerr << "error: unknown keyword encountered in param '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // GrainOrientations::update_param()



  /** grain functions
   */

  Grain::Grain() { }
  Grain::~Grain() { }


  void Grain::init() {
    layer_key_ = "air";
    in_layer_ = true; 
    layer_order_ = 0;
    lattice_.init();
    transvec_[0] = transvec_[1] = transvec_[2] = 0;
    repetition_[0] = repetition_[1] = repetition_[2] = 1;
    is_repetition_dist_ = false;
        scaling_.init();
  } // Grain::init()


  void Grain::clear() {
    unitcell_key_.clear();
    layer_key_.clear();
    in_layer_ = false;
        scaling_.clear();
    lattice_.clear();
    is_repetition_dist_ = false;
  } // Grain::clear()



  /** ensemble functions
   */

  Ensemble::Ensemble() { }
  Ensemble::~Ensemble() { }


  void Ensemble::init() {
    spacing_[0] = spacing_[1] = spacing_[2] = 0;
    maxgrains_[0] = maxgrains_[1] =  maxgrains_[2] = 1;
    distribution_ = "regular";
    orientations_.init();
  } // Ensemble::init()


  void Ensemble::clear() {
    spacing_[0] = spacing_[1] = spacing_[2] = 0;
    maxgrains_[0] = maxgrains_[1] =  maxgrains_[2] = 0;
    orientations_.clear();
  } // Ensemble::clear()



  /** structure functions
   */

  Structure::Structure() { iratio_ = 1.0; type_ = default_type; }
  Structure::~Structure() { }

  void Structure::grain_repetition_min(vector3_t v){
    grain_.xrepetition_min((unsigned) v[0]);
    grain_.yrepetition_min((unsigned) v[1]);
    grain_.zrepetition_min((unsigned) v[2]);
  }
  void Structure::grain_repetition_max(vector3_t v){
    grain_.xrepetition_max((unsigned) v[0]);
    grain_.yrepetition_max((unsigned) v[1]);
    grain_.zrepetition_max((unsigned) v[2]);
  }
  void Structure::grain_repetition_stat(std::vector<StatisticType> v){
    grain_.xrepetition_stat(v[0]);
    grain_.yrepetition_stat(v[1]);
    grain_.zrepetition_stat(v[2]);
  }

  void Structure::init() {
    key_ = "";
    grain_.init();
    ensemble_.init();
    iratio_ = 1.0;
    type_ = default_type;
  } // Structure::init()


  void Structure::clear() {
    key_.clear();
    grain_.clear();
    ensemble_.clear();
    iratio_ = 1.0;
  } // Structure::clear()


  // copy
  Structure::Structure(const Structure & rhs){
    key_ = rhs.key_;
    grain_ = rhs.grain_;
    ensemble_ = rhs.ensemble_;
    iratio_ = rhs.iratio_;
    type_ = rhs.type_;
    paracrystal_ = std::move(rhs.paracrystal_);
    percusyevick_ = std::move(rhs.percusyevick_);
  }

  //assignment
  Structure &
  Structure::operator=(const Structure & rhs){
    key_ = rhs.key_;
    grain_ = rhs.grain_;
    ensemble_ = rhs.ensemble_;
    iratio_ = rhs.iratio_;
    type_ = rhs.type_;
    paracrystal_ = std::move(rhs.paracrystal_);
    percusyevick_ = std::move(rhs.percusyevick_);
  }


  /** printing for testing
   */

  void Structure::print() {
    std::cout << " key_ = " << key_ << std::endl;
    std::cout << " grain_: " << std::endl
          << "  unitcell_key_ = " << grain_.unitcell_key_ << std::endl
          << "  layer_ley_ = " << grain_.layer_key_ << std::endl
          << "  scaling_a_ = " << grain_.scaling_.mean_[0] << std::endl
          << "  scaling_b_ = " << grain_.scaling_.mean_[1] << std::endl
          << "  scaling_c_ = " << grain_.scaling_.mean_[2] << std::endl
          << "  transvec_ = [" << grain_.transvec_[0] << ", "
          << grain_.transvec_[1] << ", " << grain_.transvec_[2]
          << "]" << std::endl
          << "  repetition_ = [" << grain_.repetition_[0] << ", "
          << grain_.repetition_[1] << ", " << grain_.repetition_[2]
          << "]" << std::endl
          << "  refindex_ = [" << grain_.refindex_.delta() << ", "
          << grain_.refindex_.beta() << "]" << std::endl
          << "  lattice_: " << std::endl
          << "   type_ = " << grain_.lattice_.type() << std::endl
          << "   abc_set_ = " << grain_.lattice_abc_set() << std::endl
          << "   a_ = [" << grain_.lattice_.a()[0] << ", " << grain_.lattice_.a()[1]
          << ", " << grain_.lattice_.a()[2] << "]" << std::endl
          << "   b_ = [" << grain_.lattice_.b()[0] << ", " << grain_.lattice_.b()[1]
          << ", " << grain_.lattice_.b()[2] << "]" << std::endl
          << "   c_ = [" << grain_.lattice_.c()[0] << ", " << grain_.lattice_.c()[1]
          << ", " << grain_.lattice_.c()[2] << "]" << std::endl
          << "   hkl_ = " << grain_.lattice_.hkl() << std::endl
          << "   abangle_ = " << grain_.lattice_.abangle() << std::endl
          << "   caratio_ = " << grain_.lattice_.caratio() << std::endl
          << "  ensemble_: " << std::endl
          << "   spacing_ = [" << ensemble_.spacing_[0] << ", " << ensemble_.spacing_[1]
          << ", " << ensemble_.spacing_[2] << "]" << std::endl
          << "   maxgrains_ = [" << ensemble_.maxgrains_[0] << ", " << ensemble_.maxgrains_[1]
          << ", " << ensemble_.maxgrains_[2] << "]" << std::endl
          << "   distribution_ = " << ensemble_.distribution_ << std::endl
          << "   orientations_: " << std::endl
          << "    stat_ = " << ensemble_.orientations_.stat() << std::endl
          << "    rot1_ = [" << ensemble_.orientations_.rot1().axis() << ", "
          << ensemble_.orientations_.rot1().angles()[0] << ", "
          << ensemble_.orientations_.rot1().angles()[1] << "]" << std::endl
          << "    rot2_ = [" << ensemble_.orientations_.rot2().axis() << ", "
          << ensemble_.orientations_.rot2().angles()[0] << ", "
          << ensemble_.orientations_.rot2().angles()[1] << "]" << std::endl
          << "    rot3_ = [" << ensemble_.orientations_.rot3().axis() << ", "
          << ensemble_.orientations_.rot3().angles()[0] << ", "
          << ensemble_.orientations_.rot3().angles()[1] << "]" << std::endl
          << std::endl;
  } // Structure::print()


  /** fitting updates
   */

  bool Structure::update_param(const std::string& str, real_t new_val) {
    std::string keyword, rem_str;
    if(!extract_first_keyword(str, keyword, rem_str)) return false;
    std::string keyword2, rem_str2;
    std::string keyword3, rem_str3;
    std::string keyword4, rem_str4;
    switch(TokenMapper::instance().get_keyword_token(keyword)) {
      case key_token:
        std::cerr << "warning: immutable param in '" << str << "'. ignoring." << std::endl;
        break;

      case struct_iratio_token:
        iratio_ = new_val;
        break;

      case struct_grain_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case refindex_token:
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            switch(TokenMapper::instance().get_keyword_token(keyword3)) {
              case refindex_delta_token:
                grain_refindex_delta(new_val);
                break;

              case refindex_beta_token:
                grain_refindex_beta(new_val);
                break;

              case error_token:
                std::cerr << "error: invalid keyword in '" << rem_str2
                      << "'" << std::endl;
                return false;

              default:
                std::cerr << "error: misplaced keyword in '" << rem_str2
                      << "'" << std::endl;
                return false;
            } // switch
            break;

          case struct_grain_lattice_token:
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            switch(TokenMapper::instance().get_keyword_token(keyword3)) {
              case struct_grain_lattice_a_token:
                if(!extract_first_keyword(rem_str3, keyword4, rem_str4)) return false;
                if(keyword4.compare("x") == 0 || keyword4.compare("0") == 0) {
                  //grain_.lattice_.a_[0] = new_val;
                  grain_.lattice_.ax(new_val);
                } else if(keyword4.compare("y") == 0 || keyword4.compare("1") == 0) {
                  //grain_.lattice_.a_[1] = new_val;
                  grain_.lattice_.ay(new_val);
                } else if(keyword4.compare("z") == 0 || keyword4.compare("2") == 0) {
                  //grain_.lattice_.a_[2] = new_val;
                  grain_.lattice_.az(new_val);
                } else {
                  std::cerr << "error: invalid keyword '" << keyword4 << "' in param '"
                        << rem_str2 << "'" << std::endl;
                  return false;
                } // if-else
                break;

              case struct_grain_lattice_b_token:
                if(!extract_first_keyword(rem_str3, keyword4, rem_str4)) return false;
                if(keyword4.compare("x") == 0 || keyword4.compare("0") == 0) {
                  //grain_.lattice_.b_[0] = new_val;
                  grain_.lattice_.bx(new_val);
                } else if(keyword4.compare("y") == 0 || keyword4.compare("1") == 0) {
                  //grain_.lattice_.b_[1] = new_val;
                  grain_.lattice_.by(new_val);
                } else if(keyword4.compare("z") == 0 || keyword4.compare("2") == 0) {
                  //grain_.lattice_.b_[2] = new_val;
                  grain_.lattice_.bz(new_val);
                } else {
                  std::cerr << "error: invalid keyword '" << keyword4 << "' in param '"
                        << rem_str2 << "'" << std::endl;
                  return false;
                } // if-else
                break;

              case struct_grain_lattice_c_token:
                if(!extract_first_keyword(rem_str3, keyword4, rem_str4)) return false;
                if(keyword4.compare("x") == 0 || keyword4.compare("0") == 0) {
                  //grain_.lattice_.c_[0] = new_val;
                  grain_.lattice_.cx(new_val);
                } else if(keyword4.compare("y") == 0 || keyword4.compare("1") == 0) {
                  //grain_.lattice_.c_[1] = new_val;
                  grain_.lattice_.cy(new_val);
                } else if(keyword4.compare("z") == 0 || keyword4.compare("2") == 0) {
                  //grain_.lattice_.c_[2] = new_val;
                  grain_.lattice_.cz(new_val);
                } else {
                  std::cerr << "error: invalid keyword '" << keyword4 << "' in param '"
                        << rem_str2 << "'" << std::endl;
                  return false;
                } // if-else
                break;

              case type_token:
              case struct_grain_lattice_hkl_token:
                std::cerr << "warning: immutable param in '" << rem_str2
                      << "'. ignoring." << std::endl;
                break;

              case struct_grain_lattice_abangle_token:
                lattice_abangle(new_val);
                break;

              case struct_grain_lattice_caratio_token:
                lattice_caratio(new_val);
                break;

              case error_token:
                std::cerr << "error: invalid keyword in '" << rem_str2
                      << "'" << std::endl;
                return false;

              default:
                std::cerr << "error: misplaced keyword in '" << rem_str2
                      << "'" << std::endl;
                return false;
            } // switch()
            break;

          case struct_grain_scaling_token:
            if(rem_str2 == "") {
              grain_.scaling_.mean_[0] = new_val;
              grain_.scaling_.mean_[1] = new_val;
              grain_.scaling_.mean_[2] = new_val;
              return true;
            } // if
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            switch(TokenMapper::instance().get_keyword_token(keyword3)) {
              case struct_grain_lattice_a_token:
                if(!extract_first_keyword(rem_str3, keyword4, rem_str4)) return false;
                switch(TokenMapper::instance().get_keyword_token(keyword4)) {
                  case mean_token:
                    grain_.scaling_.mean_[0] = new_val;
                    break;
                  case stddev_token:
                    grain_.scaling_.stddev_[0] = new_val;
                    break;
                  case nsamples_token:
                    grain_.scaling_.nvals_[0] = new_val;
                    break;
                  default:
                    std::cerr << "error: invalid keyword '" << keyword4 << "' in param '"
                              << rem_str3 << "'" << std::endl;
                  return false;
                } // if-else
                break;

              case struct_grain_lattice_b_token:
                if(!extract_first_keyword(rem_str3, keyword4, rem_str4)) return false;
                switch(TokenMapper::instance().get_keyword_token(keyword4)) {
                  case mean_token:
                    grain_.scaling_.mean_[1] = new_val;
                    break;
                  case stddev_token:
                    grain_.scaling_.stddev_[1] = new_val;
                    break;
                  case nsamples_token:
                    grain_.scaling_.nvals_[1] = new_val;
                    break;
                  default:
                    std::cerr << "error: invalid keyword '" << keyword4 << "' in param '"
                              << rem_str3 << "'" << std::endl;
                  return false;
                } // if-else
                break;

              case struct_grain_lattice_c_token:
                if(!extract_first_keyword(rem_str3, keyword4, rem_str4)) return false;
                switch(TokenMapper::instance().get_keyword_token(keyword4)) {
                  case mean_token:
                    grain_.scaling_.mean_[2] = new_val;
                    break;
                  case stddev_token:
                    grain_.scaling_.stddev_[2] = new_val;
                    break;
                  case nsamples_token:
                    grain_.scaling_.nvals_[2] = new_val;
                    break;
                  default:
                    std::cerr << "error: invalid keyword '" << keyword4 << "' in param '"
                              << rem_str3 << "'" << std::endl;
                  return false;
                } // if-else
                break;

              default:
                std::cerr << "error: invalid keyword '" << keyword3 << "' in param '"
                              << rem_str2 << "'" << std::endl;
                return false;
            } // switch
            break;

          case struct_grain_transvec_token:
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            if(keyword3.compare("x") == 0 || keyword3.compare("0") == 0) {
              grain_.transvec_[0] = new_val;
            } else if(keyword3.compare("y") == 0 || keyword3.compare("1") == 0) {
              grain_.transvec_[1] = new_val;
            } else if(keyword3.compare("z") == 0 || keyword3.compare("2") == 0) {
              grain_.transvec_[2] = new_val;
            } else {
              std::cerr << "error: invalid keyword '" << keyword3 << "' in param '"
                    << rem_str << "'" << std::endl;
              return false;
            } // if-else
            break;

          case struct_grain_repetition_token:
            //std::cerr << "warning: immutable param in '" << rem_str
            //      << "'. ignoring." << std::endl;
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            if(keyword3.compare("x") == 0 || keyword3.compare("0") == 0) {
              grain_.repetition_[0] = new_val;
            } else if(keyword3.compare("y") == 0 || keyword3.compare("1") == 0) {
              grain_.repetition_[1] = new_val;
            } else if(keyword3.compare("z") == 0 || keyword3.compare("2") == 0) {
              grain_.repetition_[2] = new_val;
            } else {
              std::cerr << "error: invalid keyword '" << keyword3 << "' in param '"
                    << rem_str << "'" << std::endl;
              return false;
            } // if-else
            break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << rem_str << "'" << std::endl;
            return false;
    
          default:
            std::cerr << "error: misplaced keyword in '" << rem_str << "'" << std::endl;
            return false;
        } // switch
        break;

      case struct_ensemble_token:
        if(!extract_first_keyword(rem_str, keyword2, rem_str2)) return false;
        switch(TokenMapper::instance().get_keyword_token(keyword2)) {
          case struct_ensemble_spacing_token:
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            if(keyword3.compare("x") == 0 || keyword3.compare("0") == 0) {
              ensemble_.spacing_[0] = new_val;
            } else if(keyword3.compare("y") == 0 || keyword3.compare("1") == 0) {
              ensemble_.spacing_[1] = new_val;
            } else if(keyword3.compare("z") == 0 || keyword3.compare("1") == 0) {
              ensemble_.spacing_[2] = new_val;
            } else {
              std::cerr << "error: invalid keyword '" << keyword3 << "' in param '"
                    << rem_str << "'" << std::endl;
              return false;
            } // if-else
            break;

          case struct_ensemble_maxgrains_token:
            if(!extract_first_keyword(rem_str2, keyword3, rem_str3)) return false;
            if(keyword3.compare("x") == 0 || keyword3.compare("0") == 0) {
              ensemble_.maxgrains_[0] = new_val;
            } else if(keyword3.compare("y") == 0 || keyword3.compare("1") == 0) {
              ensemble_.maxgrains_[1] = new_val;
            } else if(keyword3.compare("z") == 0 || keyword3.compare("1") == 0) {
              ensemble_.maxgrains_[2] = new_val;
            } else {
              std::cerr << "error: invalid keyword '" << keyword3 << "' in param '"
                    << rem_str << "'" << std::endl;
              return false;
            } // if-else
            break;

          case struct_ensemble_orient_token:
            if(!ensemble_.orientations_.update_param(rem_str2, new_val)) return false;
            break;

          case struct_ensemble_distribution_token:
            std::cerr << "warning: immutable param in '" << rem_str
                  << "'. ignoring." << std::endl;
            break;

          case error_token:
            std::cerr << "error: invalid keyword in '" << rem_str << "'" << std::endl;
            return false;
    
          default:
            std::cerr << "error: misplaced keyword in '" << rem_str << "'" << std::endl;
            return false;
        } // switch
        break;

      case error_token:
        std::cerr << "error: invalid keyword in '" << str << "'" << std::endl;
        return false;

      default:
        std::cerr << "error: misplaced keyword in '" << str << "'" << std::endl;
        return false;
    } // switch

    return true;
  } // Structure::update_param()


} // namespace hig
