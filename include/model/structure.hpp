/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: structure.hpp
 *  Created: Jun 09, 2012
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

#ifndef _STRUCTURE_HPP_
#define _STRUCTURE_HPP_

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include <common/globals.hpp>
#include <common/enums.hpp>
#include <model/common.hpp>

#ifdef DEBUG
#include <iostream>
#endif 

namespace hig {

  // make stuff private (with help of friend) ...

  class Paracrystal {
    private:
      int dims_;
      real_t dist_[2];
      real_t sdev_[2];
      real_t domsize_;
      Paracrystal(const Paracrystal &){}
      Paracrystal& operator= (const Paracrystal &){}
    public:
      Paracrystal() : dims_(1) {};

      // the puts
      void putDims(int d)          { dims_ = d;    }
      void putDistYMean(real_t d)  { dist_[0] = d; }
      void putDistXMean(real_t d)  { dist_[1] = d; }
      void putDistYStddev(real_t d){ sdev_[0] = d; }
      void putDistXStddev(real_t d){ sdev_[1] = d; }
      void putDomainSize(real_t d) { domsize_ = d; }
  
      // ... and the gets
      int    getDims()       { return dims_;    }
      real_t getDomSize()    { return domsize_; }
      real_t getDistYMean()  { return dist_[0]; }
      real_t getDistXMean()  { return dist_[1]; }
      real_t getDistYStddev(){ return sdev_[0]; }
      real_t getDistXStddev(){ return sdev_[1]; }
  };

  class PercusYevick {
    private:
      int dims_;
      real_t dia_;
      real_t volf_;
      PercusYevick(const PercusYevick &){}
      PercusYevick& operator= (const PercusYevick &){}

    public:
      PercusYevick() : dims_(2) {}

      // puts
      void putDims(int d)       { dims_ = d;}
      void putDiameter(real_t d){ dia_ = d; }
      void putVolfract(real_t f){ volf_ = f;}

      // gets
      int getDims() { return dims_; }
      real_t getDiameter(){ return dia_; }
      real_t getVolfract(){ return volf_;}
  };


  class Lattice {
    private:
      LatticeType type_;  /* overrides a_ b_ c_ vectors */
      vector3_t a_;
      vector3_t b_;
      vector3_t c_;
      vector3_t t_;    // to be computed
      bool abc_set_;
      std::string hkl_;
      real_t abangle_;
      real_t caratio_;
      real_t ca_;      // currently predefined
      real_t gamma_;    // currently predefined

      bool construct_vectors(vector3_t scaling);

    public:
      Lattice();
      ~Lattice();

      void init();
      void clear();

      LatticeType type() const { return type_; }
      bool abc_set() const { return abc_set_; }
      vector3_t a() { return a_; }
      vector3_t b() { return b_; }
      vector3_t c() { return c_; }
      vector3_t t() {return t_; }
      std::string hkl() const { return hkl_; }
      real_t abangle() const { return abangle_; }
      real_t caratio() const { return caratio_; }

      void type(LatticeType l) { type_ = l; }
      void hkl(std::string s) { hkl_ = s; }

      void a(vector3_t v) { a_ = v; }
      void b(vector3_t v) { b_ = v; }
      void c(vector3_t v) { c_ = v; }

      void a(real_t x, real_t y, real_t z) { a_[0] = x; a_[1] = y; a_[2] = z; }
      void b(real_t x, real_t y, real_t z) { b_[0] = x; b_[1] = y; b_[2] = z; }
      void c(real_t x, real_t y, real_t z) { c_[0] = x; c_[1] = y; c_[2] = z; }
      void abc_set(bool v) { abc_set_ = v; }
      void ax(real_t val) { a_[0] = val; }
      void ay(real_t val) { a_[1] = val; }
      void az(real_t val) { a_[2] = val; }
      void bx(real_t val) { b_[0] = val; }
      void by(real_t val) { b_[1] = val; }
      void bz(real_t val) { b_[2] = val; }
      void cx(real_t val) { c_[0] = val; }
      void cy(real_t val) { c_[1] = val; }
      void cz(real_t val) { c_[2] = val; }

      void abangle(real_t d) { abangle_ = d; }
      void caratio(real_t d) { caratio_ = d; }
      void bragg_angles(vector3_t, vector3_t, real_t, real_vec_t &);

      friend class Grain;

  }; // class Lattice


  class GrainScaling {
    private:
      vector3_t mean_;
      vector3_t stddev_;
      std::vector<StatisticType> dist_;
      std::vector<int> nvals_;
    public:
      GrainScaling();
      void init();
      void clear();
    friend class Grain;
    friend class Structure;
  }; // class GrainScaling

  class Rotation {
    private: 
      char axis_;    // x y or z
      std::string stat_;
      vector2_t angles_;
      real_t mean_;  // for gaussian
      real_t sd_;  // for gaussian
      bool mean_set_;

    public:
      Rotation() : axis_('n'), mean_(0), sd_(0), mean_set_(false) { }
      ~Rotation() { }
      void init();

      // getters 
      char axis() const { return axis_; }
      std::string stat() const { return stat_; }
      vector2_t angles() const { return angles_; }
      real_t angle_mean() const { return mean_; }
      real_t angle_sd() const { return sd_; }

      // setters 
      void stat(std::string stat) { stat_ = stat; }
      void axis(char c) { axis_ = c; }
      void angles(vector2_t v) { angles_ = v; }
      void angles(real_t a, real_t b) { angles_[0] = a; angles_[1] = b; }

      void angles_min(real_t val) { angles_[0] = val; if(!mean_set_) mean_ = val; }
      void angles_max(real_t val) { angles_[1] = val; }

      void angle_mean(real_t val) { mean_ = val; mean_set_ = true; }
      void angle_sd(real_t val) { sd_ = val; }
#ifdef DEBUG
      void print(){
        std::cout << "orientation = { axis: " << axis_ << ", stat: " << stat_ << ", angles: [ " 
                   << angles_[0] << " " << angles_[1] 
                   << " ], mean: " << mean_ << ", std: " << sd_ << " }" << std::endl; 
      }
#endif
    }; // class Rotation

  class GrainOrientations {

    private:
      std::string stat_;    // "single", "range", "random", "filename.ori" - change to enum?
      Rotation rot1_;      // rotation 1
      Rotation rot2_;      // rotation 2
      Rotation rot3_;      // rotation 3

    public:
      GrainOrientations();
      ~GrainOrientations();

      void init();
      void clear();

      std::string stat() const { return stat_; }
      Rotation rot1() const { return rot1_; }
      Rotation rot2() const { return rot2_; }
      Rotation rot3() const { return rot3_; }

      void stat(std::string s) { stat_ = s; }
      void rot1(const Rotation & rot) { rot1_ = rot; }
      void rot2(const Rotation & rot) { rot2_ = rot; }
      void rot3(const Rotation & rot) { rot3_ = rot; }

      void rot1_angles(vector2_t r) { rot1_.angles(r); }
      void rot2_angles(vector2_t r) { rot2_.angles(r); }
      void rot3_angles(vector2_t r) { rot3_.angles(r); }
      void rot1_angles(real_t a, real_t b) { rot1_.angles(a, b); }
      void rot2_angles(real_t a, real_t b) { rot2_.angles(a, b); }
      void rot3_angles(real_t a, real_t b) { rot3_.angles(a, b); }

      void rot1_axis(char c) { rot1_.axis(c); }
      void rot2_axis(char c) { rot2_.axis(c); }
      void rot3_axis(char c) { rot3_.axis(c); }

      void rot1_anglemean(real_t m) { rot1_.angle_mean(m); }
      void rot2_anglemean(real_t m) { rot2_.angle_mean(m); }
      void rot3_anglemean(real_t m) { rot3_.angle_mean(m); }
      void rot1_anglesd(real_t m) { rot1_.angle_sd(m); }
      void rot2_anglesd(real_t m) { rot2_.angle_sd(m); }
      void rot3_anglesd(real_t m) { rot3_.angle_sd(m); }

      bool update_param(const std::string&, real_t);

      friend class Ensemble;

  }; // class GrainOrientations


  class GrainRepetitions {
    class Repetition {
      private:
        StatisticType stat_;
        unsigned int min_;  // repetitions have to be integers
        unsigned int max_;
        real_t mean_;    // for gaussian
        real_t sd_;    // for gaussian

      public:
        Repetition(): stat_(stat_null), min_(0), max_(0), mean_(0), sd_(0) { }
        ~Repetition() { }

        void stat(StatisticType s) { stat_ = s; }
        void min(unsigned int v) { min_ = v; }
        void max(unsigned int v) { max_ = v; }
        void mean(real_t v) { mean_ = v; }
        void sd(real_t v) { sd_ = v; }

        friend class GrainRepetitions;

    }; // class Repetition

    private:
      Repetition xrepetition_;
      Repetition yrepetition_;
      Repetition zrepetition_;

    public:
      GrainRepetitions(): xrepetition_(), yrepetition_(), zrepetition_() { }
      ~GrainRepetitions() { }

      void xrepetition_stat(StatisticType s) { xrepetition_.stat(s); }
      void yrepetition_stat(StatisticType s) { yrepetition_.stat(s); }
      void zrepetition_stat(StatisticType s) { zrepetition_.stat(s); }
      void xrepetition_min(unsigned int v) { xrepetition_.min(v); }
      void yrepetition_min(unsigned int v) { yrepetition_.min(v); }
      void zrepetition_min(unsigned int v) { zrepetition_.min(v); }
      void xrepetition_max(unsigned int v) { xrepetition_.max(v); }
      void yrepetition_max(unsigned int v) { yrepetition_.max(v); }
      void zrepetition_max(unsigned int v) { zrepetition_.max(v); }
      void xrepetition_mean(real_t v) { xrepetition_.mean(v); }
      void yrepetition_mean(real_t v) { yrepetition_.mean(v); }
      void zrepetition_mean(real_t v) { zrepetition_.mean(v); }
      void xrepetition_sd(real_t v) { xrepetition_.sd(v); }
      void yrepetition_sd(real_t v) { yrepetition_.sd(v); }
      void zrepetition_sd(real_t v) { zrepetition_.sd(v); }

      /* getters */
      StatisticType xstat() const { return xrepetition_.stat_; }
      StatisticType ystat() const { return yrepetition_.stat_; }
      StatisticType zstat() const { return zrepetition_.stat_; }
      unsigned int xmin() const { return xrepetition_.min_; }
      unsigned int ymin() const { return yrepetition_.min_; }
      unsigned int zmin() const { return zrepetition_.min_; }
      unsigned int xmax() const { return xrepetition_.max_; }
      unsigned int ymax() const { return yrepetition_.max_; }
      unsigned int zmax() const { return zrepetition_.max_; }


      friend class Grain;
  }; // class GrainRepetitions

  class Grain {
    private:
      int layer_order_;
      std::string shape_key_;
      std::string unitcell_key_;
      std::string layer_key_;
      bool in_layer_;
      GrainScaling scaling_;  
      vector3_t transvec_;
      vector3_t repetition_;
      GrainRepetitions repetitiondist_;
      bool is_repetition_dist_;    // true if repetitiondist_ is defined
      RefractiveIndex refindex_;
      Lattice lattice_;

      bool construct_lattice_vectors() { return lattice_.construct_vectors(scaling_.mean_); }

    public:
      Grain();
      ~Grain();

      void init();
      void clear();

      bool lattice_abc_set() { return lattice_.abc_set(); }

      void shape_key(std::string s) { shape_key_ = s; }
      void unitcell_key(std::string s) { unitcell_key_ = s; }
      void layer_key(std::string s) { layer_key_ = s; in_layer_ = true; }
      void layer_order(int d) { layer_order_ = d; }

      void lattice_vec_a(vector3_t v) { lattice_.a(v); }
      void lattice_vec_b(vector3_t v) { lattice_.b(v); }
      void lattice_vec_c(vector3_t v) { lattice_.c(v); }
      void lattice_vec_a(real_t v, real_t w, real_t x) { lattice_.a(v, w, x); }
      void lattice_vec_b(real_t v, real_t w, real_t x) { lattice_.b(v, w, x); }
      void lattice_vec_c(real_t v, real_t w, real_t x) { lattice_.c(v, w, x); }
      void lattice_abc_set(bool v) { lattice_.abc_set(v); }

      void transvec(vector3_t v) { transvec_ = v; }
      void repetition(vector3_t v) { repetition_ = v; }
      void transvec(real_t v, real_t w, real_t x) {
        transvec_[0] = v, transvec_[1] = w, transvec_[2] = x; }
      void repetition(real_t v, real_t w, real_t x) {
        repetition_[0] = v, repetition_[1] = w, repetition_[2] = x; }

      void is_repetition_dist(bool b) { is_repetition_dist_ = b; }
      void xrepetition_stat(StatisticType s) { repetitiondist_.xrepetition_stat(s); }
      void yrepetition_stat(StatisticType s) { repetitiondist_.yrepetition_stat(s); }
      void zrepetition_stat(StatisticType s) { repetitiondist_.zrepetition_stat(s); }
      void xrepetition_min(unsigned int v) { repetitiondist_.xrepetition_min(v); }
      void yrepetition_min(unsigned int v) { repetitiondist_.yrepetition_min(v); }
      void zrepetition_min(unsigned int v) { repetitiondist_.zrepetition_min(v); }
      void xrepetition_max(unsigned int v) { repetitiondist_.xrepetition_max(v); }
      void yrepetition_max(unsigned int v) { repetitiondist_.yrepetition_max(v); }
      void zrepetition_max(unsigned int v) { repetitiondist_.zrepetition_max(v); }
      void xrepetition_mean(real_t v) { repetitiondist_.xrepetition_mean(v); }
      void yrepetition_mean(real_t v) { repetitiondist_.yrepetition_mean(v); }
      void zrepetition_mean(real_t v) { repetitiondist_.zrepetition_mean(v); }
      void xrepetition_sd(real_t v) { repetitiondist_.xrepetition_sd(v); }
      void yrepetition_sd(real_t v) { repetitiondist_.yrepetition_sd(v); }
      void zrepetition_sd(real_t v) { repetitiondist_.zrepetition_sd(v); }

      void refindex_delta(real_t d) { refindex_.delta(d); }
      void refindex_beta(real_t d) { refindex_.beta(d); }

      void lattice_abangle(real_t d) { lattice_.abangle(d); }
      void lattice_caratio(real_t d) { lattice_.caratio(d); }

      void scaling_dist(std::vector<StatisticType> s) { scaling_.dist_ = s; }
      void scaling_mean(vector3_t d) { scaling_.mean_ = d; }
      void scaling_stddev(vector3_t d) { scaling_.stddev_ = d; }
      void scaling_nsamples(std::vector<int> d) { scaling_.nvals_ = d; }

      void lattice_type(LatticeType l) { lattice_.type(l); }
      void lattice_hkl(std::string l) { lattice_.hkl(l); }


      friend class Structure;

  }; // class Grain

  class Ensemble {
    private:
      vector3_t spacing_;
      vector3_t maxgrains_;
      std::string distribution_;      // "regular", "random", "filename.spa" - change to enum?
      GrainOrientations orientations_;

    public:
      Ensemble();
      ~Ensemble();

      void init();
      void clear();

      void spacing(vector3_t v) { spacing_ = v; }
      void maxgrains(vector3_t v) { maxgrains_ = v; }
      void spacing(real_t a, real_t b, real_t c) {
        spacing_[0] = a; spacing_[1] = b; spacing_[2] = c; }
      void maxgrains(real_t a, real_t b, real_t c) {
        maxgrains_[0] = a; maxgrains_[1] = b; maxgrains_[2] = c; }

      void grain_orientation_rot1(const Rotation & rot) { orientations_.rot1(rot); }
      void grain_orientation_rot2(const Rotation & rot) { orientations_.rot2(rot); }
      void grain_orientation_rot3(const Rotation & rot) { orientations_.rot3(rot); }

      void grain_orientation_rot1_angles(vector2_t v) { orientations_.rot1_angles(v); }
      void grain_orientation_rot2_angles(vector2_t v) { orientations_.rot2_angles(v); }
      void grain_orientation_rot3_angles(vector2_t v) { orientations_.rot3_angles(v); }
      void grain_orientation_rot1_angles(real_t a, real_t b) { orientations_.rot1_angles(a, b); }
      void grain_orientation_rot2_angles(real_t a, real_t b) { orientations_.rot2_angles(a, b); }
      void grain_orientation_rot3_angles(real_t a, real_t b) { orientations_.rot3_angles(a, b); }

      void grain_orientation_rot1_axis(char c) { orientations_.rot1_axis(c); }
      void grain_orientation_rot2_axis(char c) { orientations_.rot2_axis(c); }
      void grain_orientation_rot3_axis(char c) { orientations_.rot3_axis(c); }

      void grain_orientation_rot1_mean(real_t c) { orientations_.rot1_anglemean(c); }
      void grain_orientation_rot2_mean(real_t c) { orientations_.rot2_anglemean(c); }
      void grain_orientation_rot3_mean(real_t c) { orientations_.rot3_anglemean(c); }
      void grain_orientation_rot1_sd(real_t c) { orientations_.rot1_anglesd(c); }
      void grain_orientation_rot2_sd(real_t c) { orientations_.rot2_anglesd(c); }
      void grain_orientation_rot3_sd(real_t c) { orientations_.rot3_anglesd(c); }

      void grain_orientation_stat(std::string s) { orientations_.stat(s); }

      void distribution(std::string s) { distribution_ = s; }

      friend class Structure;

  }; // class Ensemble

  class Structure {
    private:
      int dims_;
      std::string key_;
      Grain grain_;
      Ensemble ensemble_;
      real_t iratio_;
      StructureType type_;
      std::shared_ptr<Paracrystal> paracrystal_;
      std::shared_ptr<PercusYevick> percusyevick_;

    public:
      Structure();
      ~Structure();

      void init();
      void clear();

      /* setters */

      void key(std::string s) { key_ = s; }
      void iratio(real_t i) { iratio_ = i; }

      void lattice_vec_a(vector3_t v) { grain_.lattice_vec_a(v); }
      void lattice_vec_b(vector3_t v) { grain_.lattice_vec_b(v); }
      void lattice_vec_c(vector3_t v) { grain_.lattice_vec_c(v); }
      void lattice_vec_a(real_t a, real_t b, real_t c) { grain_.lattice_vec_a(a, b, c); }
      void lattice_vec_b(real_t a, real_t b, real_t c) { grain_.lattice_vec_b(a, b, c); }
      void lattice_vec_c(real_t a, real_t b, real_t c) { grain_.lattice_vec_c(a, b, c); }
      void lattice_abc_set(bool v) { grain_.lattice_abc_set(v); }

      void grain_shape_key(std::string s) { grain_.shape_key(s); }
      void grain_unitcell_key(std::string s) { grain_.unitcell_key(s); }
      void grain_layer_key(std::string s) { grain_.layer_key(s); }
      void grain_layer_order(int d){ grain_.layer_order(d); }

      void grain_transvec(vector3_t v) { grain_.transvec(v); }
      void grain_repetition(vector3_t v) { grain_.repetition(v); }
      void grain_transvec(real_t v, real_t w, real_t x) { grain_.transvec(v, w, x); }
      void grain_repetition(real_t v, real_t w, real_t x) { grain_.repetition(v, w, x); }

      void grain_repetition_min(vector3_t);
      void grain_repetition_max(vector3_t);
      void grain_repetition_stat(std::vector<StatisticType>);
      void grain_is_repetition_dist(bool b) { grain_.is_repetition_dist(b); }
      void grain_xrepetition_min(unsigned int v) { grain_.xrepetition_min(v); }
      void grain_yrepetition_min(unsigned int v) { grain_.yrepetition_min(v); }
      void grain_zrepetition_min(unsigned int v) { grain_.zrepetition_min(v); }
      void grain_xrepetition_max(unsigned int v) { grain_.xrepetition_max(v); }
      void grain_yrepetition_max(unsigned int v) { grain_.yrepetition_max(v); }
      void grain_zrepetition_max(unsigned int v) { grain_.zrepetition_max(v); }
      void grain_xrepetition_stat(StatisticType s) { grain_.xrepetition_stat(s); }
      void grain_yrepetition_stat(StatisticType s) { grain_.yrepetition_stat(s); }
      void grain_zrepetition_stat(StatisticType s) { grain_.zrepetition_stat(s); }


      void grain_orientation_rot1(const Rotation & rot) { ensemble_.grain_orientation_rot1(rot); }
      void grain_orientation_rot2(const Rotation & rot) { ensemble_.grain_orientation_rot2(rot); }
      void grain_orientation_rot3(const Rotation & rot) { ensemble_.grain_orientation_rot3(rot); }

      void grain_orientation_rot1_angles(vector2_t v) { ensemble_.grain_orientation_rot1_angles(v); }
      void grain_orientation_rot2_angles(vector2_t v) { ensemble_.grain_orientation_rot2_angles(v); }
      void grain_orientation_rot3_angles(vector2_t v) { ensemble_.grain_orientation_rot3_angles(v); }
      void grain_orientation_rot1_angles(real_t a, real_t b) {
        ensemble_.grain_orientation_rot1_angles(a, b); }
      void grain_orientation_rot2_angles(real_t a, real_t b) {
        ensemble_.grain_orientation_rot2_angles(a, b); }
      void grain_orientation_rot3_angles(real_t a, real_t b) {
        ensemble_.grain_orientation_rot3_angles(a, b); }

      void grain_orientation_rot1_axis(char c) { ensemble_.grain_orientation_rot1_axis(c); }
      void grain_orientation_rot2_axis(char c) { ensemble_.grain_orientation_rot2_axis(c); }
      void grain_orientation_rot3_axis(char c) { ensemble_.grain_orientation_rot3_axis(c); }

      void grain_orientation_rot1_anglemean(real_t a) { ensemble_.grain_orientation_rot1_mean(a); }
      void grain_orientation_rot2_anglemean(real_t a) { ensemble_.grain_orientation_rot2_mean(a); }
      void grain_orientation_rot3_anglemean(real_t a) { ensemble_.grain_orientation_rot3_mean(a); }
      void grain_orientation_rot1_anglesd(real_t a) { ensemble_.grain_orientation_rot1_sd(a); }
      void grain_orientation_rot2_anglesd(real_t a) { ensemble_.grain_orientation_rot2_sd(a); }
      void grain_orientation_rot3_anglesd(real_t a) { ensemble_.grain_orientation_rot3_sd(a); }

      void grain_refindex_delta(real_t d) { grain_.refindex_delta(d); }
      void grain_refindex_beta(real_t d) { grain_.refindex_beta(d); }

            // set grain scaling parameters
      void grain_scaling_a_stat (StatisticType s) { grain_.scaling_.dist_[0] = s; }
      void grain_scaling_b_stat (StatisticType s) { grain_.scaling_.dist_[1] = s; }
      void grain_scaling_c_stat (StatisticType s) { grain_.scaling_.dist_[2] = s; }
      void grain_scaling_a_mean (real_t num) { grain_.scaling_.mean_[0] = num; } 
      void grain_scaling_b_mean (real_t num) { grain_.scaling_.mean_[1] = num; } 
      void grain_scaling_c_mean (real_t num) { grain_.scaling_.mean_[2] = num; }
      void grain_scaling_a_stddev (real_t num) { grain_.scaling_.stddev_[0] = num; }
      void grain_scaling_b_stddev (real_t num) { grain_.scaling_.stddev_[1] = num; } 
      void grain_scaling_c_stddev (real_t num) { grain_.scaling_.stddev_[2] = num; }
      void grain_scaling_a_nsamples (int num) { grain_.scaling_.nvals_[0] = num; }
      void grain_scaling_b_nsamples (int num) { grain_.scaling_.nvals_[1] = num; }
      void grain_scaling_c_nsamples (int num) { grain_.scaling_.nvals_[2] = num; }
      void grain_scaling_mean(vector3_t v) { grain_.scaling_.mean_ = v; }
      void grain_scaling_stddev(vector3_t v) { grain_.scaling_.stddev_ = v; }
      void grain_scaling_stat(std::vector<StatisticType> v) { grain_.scaling_.dist_ = v; }
      void grain_scaling_nsamples(std::vector<int> v) { grain_.scaling_.nvals_ = v; }

      void ensemble_spacing(vector3_t v) { ensemble_.spacing(v); }
      void ensemble_maxgrains(vector3_t v) { ensemble_.maxgrains(v); }
      void ensemble_spacing(real_t v, real_t w, real_t x) { ensemble_.spacing(v, w, x); }
      void ensemble_maxgrains(real_t v, real_t w, real_t x) { ensemble_.maxgrains(v, w, x); }

      void ensemble_orientation_stat(std::string s) { ensemble_.grain_orientation_stat(s); }
      void ensemble_distribution(std::string s) { ensemble_.distribution(s); }

      void lattice_abangle(real_t d) { grain_.lattice_abangle(d); }
      void lattice_caratio(real_t d) { grain_.lattice_caratio(d); }

      void lattice_type(LatticeType l) { grain_.lattice_type(l); }
      void lattice_hkl(std::string l) { grain_.lattice_hkl(l); }

      // **** liquids / begin ****
      // paracrystal
      
      void putStructureType (StructureType d)    { type_ = d;                              }
      StructureType getStructureType() const     { return type_;                           }
      void paracrystal_init()                    { paracrystal_.reset(new Paracrystal());  
                                                   type_ = paracrystal_type;               }
      void paracrystal_putDimensions(int d)      { paracrystal_->putDims(d);               }
      void paracrystal_putDistYMean(real_t d)    { paracrystal_->putDistYMean(d);          }
      void paracrystal_putDistXMean(real_t d)    { paracrystal_->putDistXMean(d);          }
      void paracrystal_putDistYStddev(real_t d)  { paracrystal_->putDistYStddev(d);        }
      void paracrystal_putDistXStddev(real_t d)  { paracrystal_->putDistXStddev(d);        }
      void paracrystal_putDomainSize(real_t d)   { paracrystal_->putDomainSize(d);         }
      int  paracrystal_getDimensions() const     { return paracrystal_->getDims();         }
      real_t paracrystal_getDomsSize() const     { return paracrystal_->getDomSize();      }
      real_t paracrystal_getDistXMean() const    { return paracrystal_->getDistXMean();    }
      real_t paracrystal_getDistYMean() const    { return paracrystal_->getDistYMean();    }
      real_t paracrystal_getDistXStddev() const  { return paracrystal_->getDistXStddev();  }
      real_t paracrystal_getDistYStddev() const  { return paracrystal_->getDistYStddev();  }

      // percus-yevick
      void percusyevick_init()                { percusyevick_.reset(new PercusYevick());   
                                                type_ = percusyevick_type;                 }
      void percusyevick_putDimensions(int d)  { percusyevick_->putDims(d);                 }
      void percusyevick_putDiameter(real_t d) { percusyevick_->putDiameter(d);             }
      void percusyevick_putVolf(real_t d)     { percusyevick_->putVolfract(d);             }
      int  percusyevick_getDimensions() const { return percusyevick_->getDims();           }
      real_t percusyevick_getDiameter() const { return percusyevick_->getDiameter();       }
      real_t percusyevick_getVolfract() const { return percusyevick_->getVolfract();       }

      std::shared_ptr<Paracrystal> paracrystal() const   { return paracrystal_; }
      std::shared_ptr<PercusYevick> percusyevick() const { return percusyevick_;}
      // **** liquids / end ****

      /* computers */

      bool construct_lattice_vectors() { return grain_.construct_lattice_vectors(); }

      /* getters */

      std::string key() const { return key_; }
      real_t iratio() const { return iratio_; }

      const Lattice* lattice() const { return &(grain_.lattice_); }
      bool lattice_abc_set() { return grain_.lattice_abc_set(); }

      vector3_t grain_repetition() const { return grain_.repetition_; }
      bool grain_is_repetition_dist() const { return grain_.is_repetition_dist_; }
      const GrainRepetitions& grain_repetitiondist() const { return grain_.repetitiondist_; }
      std::string grain_orientation() const { return ensemble_.orientations_.stat(); }
      RefractiveIndex grain_refindex() const { return grain_.refindex_; }
      complex_t one_minus_n2() const { return grain_.refindex_.one_minus_n2(); }
      const std::string& grain_unitcell_key() const { return grain_.unitcell_key_; }
      const std::string& grain_layer_key() const { return grain_.layer_key_; }
      int layer_order() const { return grain_.layer_order_; }
      bool grain_in_layer() { return grain_.in_layer_; }
      vector3_t grain_transvec() const { return grain_.transvec_; }

      // scaling related params
      bool grain_scaling_is_dist () const {
        for (int i = 0; i < 3; i++)
          if ( grain_.scaling_.stddev_[i] > 0 )
            return true;
        return false;
      }

      vector3_t grain_scaling() const { return grain_.scaling_.mean_; }
      vector3_t grain_scaling_stddev() const { return grain_.scaling_.stddev_; }
      std::vector<StatisticType> grain_scaling_dist() const { return grain_.scaling_.dist_; }
      std::vector<int> grain_scaling_nvals() const { return grain_.scaling_.nvals_; }

      vector3_t ensemble_spacing() const { return ensemble_.spacing_; }
      vector3_t ensemble_maxgrains() const { return ensemble_.maxgrains_; }
      std::string ensemble_distribution() const { return ensemble_.distribution_; }

      vector2_t rotation_tau() { return ensemble_.orientations_.rot1().angles(); }
      vector2_t rotation_eta() { return ensemble_.orientations_.rot2().angles(); }
      vector2_t rotation_zeta() { return ensemble_.orientations_.rot3().angles(); }

      vector3_t rotation_rot1() const {
        char axis = ensemble_.orientations_.rot1().axis();
        if(axis == 'n') return vector3_t(0, 0, 0);
        vector2_t angs = ensemble_.orientations_.rot1().angles();
        return vector3_t((real_t) ((char)axis - 'x'), angs[0], angs[1]);
      } // rotation_rot1()
      vector3_t rotation_rot2() const {
        char axis = ensemble_.orientations_.rot2().axis();
        if(axis == 'n') return vector3_t(0, 0, 0);
        vector2_t angs = ensemble_.orientations_.rot2().angles();
        return vector3_t((real_t) ((char)axis - 'x'), angs[0], angs[1]);
      } // rotation_rot2()
      vector3_t rotation_rot3() const {
        char axis = ensemble_.orientations_.rot3().axis();
        if(axis == 'n') return vector3_t(0, 0, 0);
        vector2_t angs = ensemble_.orientations_.rot3().angles();
        return vector3_t((real_t) ((char)axis - 'x'), angs[0], angs[1]);
      } // rotation_rot3()

      real_t rotation_rot1_anglemean() const { return ensemble_.orientations_.rot1().angle_mean(); }
      real_t rotation_rot2_anglemean() const { return ensemble_.orientations_.rot2().angle_mean(); }
      real_t rotation_rot3_anglemean() const { return ensemble_.orientations_.rot3().angle_mean(); }
      real_t rotation_rot1_anglesd() const { return ensemble_.orientations_.rot1().angle_sd(); }
      real_t rotation_rot2_anglesd() const { return ensemble_.orientations_.rot2().angle_sd(); }
      real_t rotation_rot3_anglesd() const { return ensemble_.orientations_.rot3().angle_sd(); }

      /* modifiers (updates) */
      bool update_param(const std::string&, real_t);


      // copy constructor and assignment
      Structure(const Structure & );
      Structure & operator=(const Structure & );
      /* testers */

      void print();

  }; // class Structure

  typedef std::unordered_map <std::string, Structure> structure_list_t;
  typedef structure_list_t::iterator structure_iterator_t;
  typedef structure_list_t::const_iterator structure_citerator_t;

} // namespace hig

#endif /* _STRUCTURE_HPP_ */
