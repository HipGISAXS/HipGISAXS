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
#include <unordered_map>

#include <common/globals.hpp>
#include <common/enums.hpp>
#include <model/common.hpp>

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
			void ax(float_t val) { a_[0] = val; }
			void ay(float_t val) { a_[1] = val; }
			void az(float_t val) { a_[2] = val; }
			void bx(float_t val) { b_[0] = val; }
			void by(float_t val) { b_[1] = val; }
			void bz(float_t val) { b_[2] = val; }
			void cx(float_t val) { c_[0] = val; }
			void cy(float_t val) { c_[1] = val; }
			void cz(float_t val) { c_[2] = val; }

			void abangle(float_t d) { abangle_ = d; }
			void caratio(float_t d) { caratio_ = d; }

			friend class Grain;

	}; // class Lattice

	class GrainOrientations {

		class Rotation {
			private: 
				char axis_;		// x y or z
				vector2_t angles_;
				float_t mean_;	// for gaussian
				float_t sd_;	// for gaussian
				bool mean_set_;

			public:
				Rotation() : axis_('n'), mean_(0), sd_(0), mean_set_(false) { }
				~Rotation() { }

				void init();

				char axis() { return axis_; }
				vector2_t angles() { return angles_; }

				void axis(char c) { axis_ = c; }
				void angles(vector2_t v) { angles_ = v; }
				void angles(float_t a, float_t b) { angles_[0] = a; angles_[1] = b; }

				void angles_min(float_t val) { angles_[0] = val; if(!mean_set_) mean_ = val; }
				void angles_max(float_t val) { angles_[1] = val; }

				void angle_mean(float_t val) { mean_ = val; mean_set_ = true; }
				void angle_sd(float_t val) { sd_ = val; }
				float_t angle_mean() { return mean_; }
				float_t angle_sd() { return sd_; }

		}; // class Rotation

		private:
			std::string stat_;		// "single", "range", "random", "filename.ori" - change to enum?
			Rotation rot1_;			// rotation 1
			Rotation rot2_;			// rotation 2
			Rotation rot3_;			// rotation 3

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

			void rot1_anglemean(float_t m) { rot1_.angle_mean(m); }
			void rot2_anglemean(float_t m) { rot2_.angle_mean(m); }
			void rot3_anglemean(float_t m) { rot3_.angle_mean(m); }
			void rot1_anglesd(float_t m) { rot1_.angle_sd(m); }
			void rot2_anglesd(float_t m) { rot2_.angle_sd(m); }
			void rot3_anglesd(float_t m) { rot3_.angle_sd(m); }

			bool update_param(const std::string&, float_t);

			friend class Ensemble;

	}; // class GrainOrientations


	class GrainRepetitions {
		class Repetition {
			private:
				StatisticType stat_;
				unsigned int min_;	// repetitions have to be integers
				unsigned int max_;
				float_t mean_;		// for gaussian
				float_t sd_;		// for gaussian

			public:
				Repetition(): stat_(stat_null), min_(0), max_(0), mean_(0), sd_(0) { }
				~Repetition() { }

				void stat(StatisticType s) { stat_ = s; }
				void min(unsigned int v) { min_ = v; }
				void max(unsigned int v) { max_ = v; }
				void mean(float_t v) { mean_ = v; }
				void sd(float_t v) { sd_ = v; }

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
			void xrepetition_mean(float_t v) { xrepetition_.mean(v); }
			void yrepetition_mean(float_t v) { yrepetition_.mean(v); }
			void zrepetition_mean(float_t v) { zrepetition_.mean(v); }
			void xrepetition_sd(float_t v) { xrepetition_.sd(v); }
			void yrepetition_sd(float_t v) { yrepetition_.sd(v); }
			void zrepetition_sd(float_t v) { zrepetition_.sd(v); }

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
			std::string shape_key_;
			std::string layer_key_;
			bool in_layer_;
			float_t scaling_;
			vector3_t transvec_;
			vector3_t repetition_;
			GrainRepetitions repetitiondist_;
			bool is_repetition_dist_;		// true if repetitiondist_ is defined
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
			void layer_key(std::string s) { layer_key_ = s; in_layer_ = true; }

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
			void xrepetition_mean(float_t v) { repetitiondist_.xrepetition_mean(v); }
			void yrepetition_mean(float_t v) { repetitiondist_.yrepetition_mean(v); }
			void zrepetition_mean(float_t v) { repetitiondist_.zrepetition_mean(v); }
			void xrepetition_sd(float_t v) { repetitiondist_.xrepetition_sd(v); }
			void yrepetition_sd(float_t v) { repetitiondist_.yrepetition_sd(v); }
			void zrepetition_sd(float_t v) { repetitiondist_.zrepetition_sd(v); }

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

			void grain_orientation_rot1_mean(float_t c) { orientations_.rot1_anglemean(c); }
			void grain_orientation_rot2_mean(float_t c) { orientations_.rot2_anglemean(c); }
			void grain_orientation_rot3_mean(float_t c) { orientations_.rot3_anglemean(c); }
			void grain_orientation_rot1_sd(float_t c) { orientations_.rot1_anglesd(c); }
			void grain_orientation_rot2_sd(float_t c) { orientations_.rot2_anglesd(c); }
			void grain_orientation_rot3_sd(float_t c) { orientations_.rot3_anglesd(c); }

			void grain_orientation_stat(std::string s) { orientations_.stat(s); }

			void distribution(std::string s) { distribution_ = s; }

			friend class Structure;

	}; // class Ensemble

	class Structure {
		private:
			std::string key_;
			Grain grain_;
			Ensemble ensemble_;
			float_t iratio_;

		public:
			Structure();
			~Structure();

			void init();
			void clear();

			/* setters */

			void key(std::string s) { key_ = s; }
			void iratio(float_t i) { iratio_ = i; }

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

			void grain_orientation_rot1_anglemean(float_t a) { ensemble_.grain_orientation_rot1_mean(a); }
			void grain_orientation_rot2_anglemean(float_t a) { ensemble_.grain_orientation_rot2_mean(a); }
			void grain_orientation_rot3_anglemean(float_t a) { ensemble_.grain_orientation_rot3_mean(a); }
			void grain_orientation_rot1_anglesd(float_t a) { ensemble_.grain_orientation_rot1_sd(a); }
			void grain_orientation_rot2_anglesd(float_t a) { ensemble_.grain_orientation_rot2_sd(a); }
			void grain_orientation_rot3_anglesd(float_t a) { ensemble_.grain_orientation_rot3_sd(a); }

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
			float_t iratio() const { return iratio_; }

			const Lattice* lattice() const { return &(grain_.lattice_); }
			bool lattice_abc_set() { return grain_.lattice_abc_set(); }

			vector3_t grain_repetition() const { return grain_.repetition_; }
			bool grain_is_repetition_dist() const { return grain_.is_repetition_dist_; }
			const GrainRepetitions& grain_repetitiondist() const { return grain_.repetitiondist_; }
			std::string grain_orientation() { return ensemble_.orientations_.stat(); }
			RefractiveIndex grain_refindex() { return grain_.refindex_; }
			const std::string& grain_shape_key() { return grain_.shape_key_; }
			const std::string& grain_layer_key() { return grain_.layer_key_; }
			bool grain_in_layer() { return grain_.in_layer_; }
			vector3_t grain_transvec() { return grain_.transvec_; }
			float_t grain_scaling() const { return grain_.scaling_; }

			vector3_t ensemble_spacing() const { return ensemble_.spacing_; }
			vector3_t ensemble_maxgrains() const { return ensemble_.maxgrains_; }
			std::string ensemble_distribution() const { return ensemble_.distribution_; }

			vector2_t rotation_tau() { return ensemble_.orientations_.rot1().angles(); }
			vector2_t rotation_eta() { return ensemble_.orientations_.rot2().angles(); }
			vector2_t rotation_zeta() { return ensemble_.orientations_.rot3().angles(); }

			vector3_t rotation_rot1() {
				char axis = ensemble_.orientations_.rot1().axis();
				if(axis == 'n') return vector3_t(0, 0, 0);
				vector2_t angs = ensemble_.orientations_.rot1().angles();
				return vector3_t((float_t) ((char)axis - 'x'), angs[0], angs[1]);
				//return vector3_t((float_t) axis, angs[0], angs[1]);
			} // rotation_rot1()
			vector3_t rotation_rot2() {
				char axis = ensemble_.orientations_.rot2().axis();
				if(axis == 'n') return vector3_t(0, 0, 0);
				vector2_t angs = ensemble_.orientations_.rot2().angles();
				return vector3_t((float_t) ((char)axis - 'x'), angs[0], angs[1]);
				//return vector3_t((float_t) axis, angs[0], angs[1]);
			} // rotation_rot2()
			vector3_t rotation_rot3() {
				char axis = ensemble_.orientations_.rot3().axis();
				if(axis == 'n') return vector3_t(0, 0, 0);
				vector2_t angs = ensemble_.orientations_.rot3().angles();
				return vector3_t((float_t) ((char)axis - 'x'), angs[0], angs[1]);
				//return vector3_t((float_t) axis, angs[0], angs[1]);
			} // rotation_rot3()

			float_t rotation_rot1_anglemean() { return ensemble_.orientations_.rot1().angle_mean(); }
			float_t rotation_rot2_anglemean() { return ensemble_.orientations_.rot2().angle_mean(); }
			float_t rotation_rot3_anglemean() { return ensemble_.orientations_.rot3().angle_mean(); }
			float_t rotation_rot1_anglesd() { return ensemble_.orientations_.rot1().angle_sd(); }
			float_t rotation_rot2_anglesd() { return ensemble_.orientations_.rot2().angle_sd(); }
			float_t rotation_rot3_anglesd() { return ensemble_.orientations_.rot3().angle_sd(); }

			/* modifiers (updates) */
			bool update_param(const std::string&, float_t);

			/* testers */

			void print();

	}; // class Structure

	typedef std::unordered_map <std::string, Structure> structure_list_t;
	typedef structure_list_t::iterator structure_iterator_t;

} // namespace hig

#endif /* _STRUCTURE_HPP_ */
