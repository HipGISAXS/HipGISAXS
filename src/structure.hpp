/***
  *  $Id: structure.hpp 33 2012-08-06 16:22:01Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: structure.hpp
  *  Created: Jun 09, 2012
  *  Modified: Mon 01 Oct 2012 11:18:40 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _STRUCTURE_HPP_
#define _STRUCTURE_HPP_

#include <string>

#include "globals.hpp"
#include "enums.hpp"
#include "common.hpp"

namespace hig {

	// make stuff private (with help of friend) ...

	class Lattice {
		private:
			LatticeType type_;	/* overrides a_ b_ c_ vectors */
			vector3_t a_;
			vector3_t b_;
			vector3_t c_;
			vector3_t t_;		// to be computed
			bool abc_set_;
			std::string hkl_;
			float_t abangle_;
			float_t caratio_;
			float_t ca_;			// currently predefined
			float_t gamma_;		// currently predefined

			bool construct_vectors(float_t scaling);

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
			float_t abangle() const { return abangle_; }
			float_t caratio() const { return caratio_; }

			void type(LatticeType l) { type_ = l; }
			void hkl(std::string s) { hkl_ = s; }

			void a(vector3_t v) { a_ = v; }
			void b(vector3_t v) { b_ = v; }
			void c(vector3_t v) { c_ = v; }

			void a(float_t x, float_t y, float_t z) { a_[0] = x; a_[1] = y; a_[2] = z; }
			void b(float_t x, float_t y, float_t z) { b_[0] = x; b_[1] = y; b_[2] = z; }
			void c(float_t x, float_t y, float_t z) { c_[0] = x; c_[1] = y; c_[2] = z; }
			void abc_set(bool v) { abc_set_ = v; }

			void abangle(float_t d) { abangle_ = d; }
			void caratio(float_t d) { caratio_ = d; }

			friend class Grain;

	}; // class Lattice

	class GrainOrientations {

		class Rotation {
			private: 
				char axis_;		// x y or z
				vector2_t angles_;

			public:
				Rotation() { axis_ = 'n'; }
				~Rotation() { }

				void init();

				char axis() { return axis_; }
				vector2_t angles() { return angles_; }

				void axis(char c) { axis_ = c; }
				void angles(vector2_t v) { angles_ = v; }
				void angles(float_t a, float_t b) { angles_[0] = a; angles_[1] = b; }

		}; // class Rotation

		private:
			std::string stat_;		// "single", "range", "random", "filename.ori" - change to enum?
			Rotation rot1_;			// tau
			Rotation rot2_;			// eta
			Rotation rot3_;			// zeta

		public:
			GrainOrientations();
			~GrainOrientations();

			void init();
			void clear();

			std::string stat() { return stat_; }
			Rotation rot1() { return rot1_; }
			Rotation rot2() { return rot2_; }
			Rotation rot3() { return rot3_; }

			void stat(std::string s) { stat_ = s; }

			void rot1_angles(vector2_t r) { rot1_.angles(r); }
			void rot2_angles(vector2_t r) { rot2_.angles(r); }
			void rot3_angles(vector2_t r) { rot3_.angles(r); }
			void rot1_angles(float_t a, float_t b) { rot1_.angles(a, b); }
			void rot2_angles(float_t a, float_t b) { rot2_.angles(a, b); }
			void rot3_angles(float_t a, float_t b) { rot3_.angles(a, b); }

			void rot1_axis(char c) { rot1_.axis(c); }
			void rot2_axis(char c) { rot2_.axis(c); }
			void rot3_axis(char c) { rot3_.axis(c); }

			friend class Ensemble;

	}; // class GrainOrientations

	class Grain {
		private:
			std::string shape_key_;
			std::string layer_key_;
			float_t scaling_;
			vector3_t transvec_;
			vector3_t repetition_;
			RefractiveIndex refindex_;
			Lattice lattice_;

			bool construct_lattice_vectors() { return lattice_.construct_vectors(scaling_); }

		public:
			Grain();
			~Grain();

			void init();
			void clear();

			bool lattice_abc_set() { return lattice_.abc_set(); }

			void shape_key(std::string s) { shape_key_ = s; }
			void layer_key(std::string s) { layer_key_ = s; }

			void lattice_vec_a(vector3_t v) { lattice_.a(v); }
			void lattice_vec_b(vector3_t v) { lattice_.b(v); }
			void lattice_vec_c(vector3_t v) { lattice_.c(v); }
			void lattice_vec_a(float_t v, float_t w, float_t x) { lattice_.a(v, w, x); }
			void lattice_vec_b(float_t v, float_t w, float_t x) { lattice_.b(v, w, x); }
			void lattice_vec_c(float_t v, float_t w, float_t x) { lattice_.c(v, w, x); }
			void lattice_abc_set(bool v) { lattice_.abc_set(v); }

			void transvec(vector3_t v) { transvec_ = v; }
			void repetition(vector3_t v) { repetition_ = v; }
			void transvec(float_t v, float_t w, float_t x) {
				transvec_[0] = v, transvec_[1] = w, transvec_[2] = x; }
			void repetition(float_t v, float_t w, float_t x) {
				repetition_[0] = v, repetition_[1] = w, repetition_[2] = x; }

			void refindex_delta(float_t d) { refindex_.delta(d); }
			void refindex_beta(float_t d) { refindex_.beta(d); }

			void lattice_abangle(float_t d) { lattice_.abangle(d); }
			void lattice_caratio(float_t d) { lattice_.caratio(d); }

			void scaling(float_t d) { scaling_ = d; }

			void lattice_type(LatticeType l) { lattice_.type(l); }
			void lattice_hkl(std::string l) { lattice_.hkl(l); }

			friend class Structure;

	}; // class Grain

	class Ensemble {
		private:
			vector3_t spacing_;
			vector3_t maxgrains_;
			std::string distribution_;			// "regular", "random", "filename.spa" - change to enum?
			GrainOrientations orientations_;

		public:
			Ensemble();
			~Ensemble();

			void init();
			void clear();

			void spacing(vector3_t v) { spacing_ = v; }
			void maxgrains(vector3_t v) { maxgrains_ = v; }
			void spacing(float_t a, float_t b, float_t c) {
				spacing_[0] = a; spacing_[1] = b; spacing_[2] = c; }
			void maxgrains(float_t a, float_t b, float_t c) {
				maxgrains_[0] = a; maxgrains_[1] = b; maxgrains_[2] = c; }

			void grain_orientation_rot1_angles(vector2_t v) { orientations_.rot1_angles(v); }
			void grain_orientation_rot2_angles(vector2_t v) { orientations_.rot2_angles(v); }
			void grain_orientation_rot3_angles(vector2_t v) { orientations_.rot3_angles(v); }
			void grain_orientation_rot1_angles(float_t a, float_t b) { orientations_.rot1_angles(a, b); }
			void grain_orientation_rot2_angles(float_t a, float_t b) { orientations_.rot2_angles(a, b); }
			void grain_orientation_rot3_angles(float_t a, float_t b) { orientations_.rot3_angles(a, b); }

			void grain_orientation_rot1_axis(char c) { orientations_.rot1_axis(c); }
			void grain_orientation_rot2_axis(char c) { orientations_.rot2_axis(c); }
			void grain_orientation_rot3_axis(char c) { orientations_.rot3_axis(c); }

			void grain_orientation_stat(std::string s) { orientations_.stat(s); }

			void distribution(std::string s) { distribution_ = s; }

			friend class Structure;

	}; // class Ensemble

	class Structure {
		private:
			std::string key_;
			Grain grain_;
			Ensemble ensemble_;

		public:
			Structure();
			~Structure();

			void init();
			void clear();

			/* setters */

			void key(std::string s) { key_ = s; }

			void lattice_vec_a(vector3_t v) { grain_.lattice_vec_a(v); }
			void lattice_vec_b(vector3_t v) { grain_.lattice_vec_b(v); }
			void lattice_vec_c(vector3_t v) { grain_.lattice_vec_c(v); }
			void lattice_vec_a(float_t a, float_t b, float_t c) { grain_.lattice_vec_a(a, b, c); }
			void lattice_vec_b(float_t a, float_t b, float_t c) { grain_.lattice_vec_b(a, b, c); }
			void lattice_vec_c(float_t a, float_t b, float_t c) { grain_.lattice_vec_c(a, b, c); }
			void lattice_abc_set(bool v) { grain_.lattice_abc_set(v); }

			void grain_shape_key(std::string s) { grain_.shape_key(s); }
			void grain_layer_key(std::string s) { grain_.layer_key(s); }

			void grain_transvec(vector3_t v) { grain_.transvec(v); }
			void grain_repetition(vector3_t v) { grain_.repetition(v); }
			void grain_transvec(float_t v, float_t w, float_t x) { grain_.transvec(v, w, x); }
			void grain_repetition(float_t v, float_t w, float_t x) { grain_.repetition(v, w, x); }

			void grain_orientation_rot1_angles(vector2_t v) { ensemble_.grain_orientation_rot1_angles(v); }
			void grain_orientation_rot2_angles(vector2_t v) { ensemble_.grain_orientation_rot2_angles(v); }
			void grain_orientation_rot3_angles(vector2_t v) { ensemble_.grain_orientation_rot3_angles(v); }
			void grain_orientation_rot1_angles(float_t a, float_t b) {
				ensemble_.grain_orientation_rot1_angles(a, b); }
			void grain_orientation_rot2_angles(float_t a, float_t b) {
				ensemble_.grain_orientation_rot2_angles(a, b); }
			void grain_orientation_rot3_angles(float_t a, float_t b) {
				ensemble_.grain_orientation_rot3_angles(a, b); }

			void grain_orientation_rot1_axis(char c) { ensemble_.grain_orientation_rot1_axis(c); }
			void grain_orientation_rot2_axis(char c) { ensemble_.grain_orientation_rot2_axis(c); }
			void grain_orientation_rot3_axis(char c) { ensemble_.grain_orientation_rot3_axis(c); }

			void grain_refindex_delta(float_t d) { grain_.refindex_delta(d); }
			void grain_refindex_beta(float_t d) { grain_.refindex_beta(d); }

			void grain_scaling(float_t d) { grain_.scaling(d); }

			void ensemble_spacing(vector3_t v) { ensemble_.spacing(v); }
			void ensemble_maxgrains(vector3_t v) { ensemble_.maxgrains(v); }
			void ensemble_spacing(float_t v, float_t w, float_t x) { ensemble_.spacing(v, w, x); }
			void ensemble_maxgrains(float_t v, float_t w, float_t x) { ensemble_.maxgrains(v, w, x); }

			void ensemble_orientation_stat(std::string s) { ensemble_.grain_orientation_stat(s); }
			void ensemble_distribution(std::string s) { ensemble_.distribution(s); }

			void lattice_abangle(float_t d) { grain_.lattice_abangle(d); }
			void lattice_caratio(float_t d) { grain_.lattice_caratio(d); }

			void lattice_type(LatticeType l) { grain_.lattice_type(l); }
			void lattice_hkl(std::string l) { grain_.lattice_hkl(l); }

			/* computers */

			bool construct_lattice_vectors() { return grain_.construct_lattice_vectors(); }

			/* getters */

			std::string key() const { return key_; }

			const Lattice* lattice() const { return &(grain_.lattice_); }
			bool lattice_abc_set() { return grain_.lattice_abc_set(); }

			vector3_t grain_repetition() const { return grain_.repetition_; }
			std::string grain_orientation() { return ensemble_.orientations_.stat(); }
			RefractiveIndex grain_refindex() { return grain_.refindex_; }
			const std::string& grain_shape_key() { return grain_.shape_key_; }
			const std::string& grain_layer_key() { return grain_.layer_key_; }
			vector3_t grain_transvec() { return grain_.transvec_; }

			vector3_t ensemble_spacing() const { return ensemble_.spacing_; }
			vector3_t ensemble_maxgrains() const { return ensemble_.maxgrains_; }
			std::string ensemble_distribution() const { return ensemble_.distribution_; }

			vector2_t rotation_tau() { return ensemble_.orientations_.rot1().angles(); }
			vector2_t rotation_eta() { return ensemble_.orientations_.rot2().angles(); }
			vector2_t rotation_zeta() { return ensemble_.orientations_.rot3().angles(); }

			/* testers */

			void print();

	}; // class Structure

	typedef std::unordered_map <std::string, Structure> structure_list_t;
	typedef structure_list_t::iterator structure_iterator_t;

} // namespace hig

#endif /* _STRUCTURE_HPP_ */
