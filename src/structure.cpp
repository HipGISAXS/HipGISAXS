/***
  *  $Id: structure.cpp 42 2012-08-22 05:07:05Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: structure.cpp
  *  Created: Jun 12, 2012
  *  Modified: Fri 07 Dec 2012 01:34:09 PM PST
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>

#include "structure.hpp"

namespace hig {


	/** lattice functions
	 */


	Lattice::Lattice() { ca_ = 1; gamma_ = 0; }
	Lattice::~Lattice() { }


	void Lattice::init() {
		a_[0] = a_[1] = a_[2] = b_[0] = b_[1] = b_[2] = c_[0] = c_[1] = c_[2] = 0.0;
		t_[0] = t_[1] = t_[2] = 0.0;
		abc_set_ = false;
		type_ = lattice_cubic;		// default values
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


	bool Lattice::construct_vectors(float_t scaling) {
		if(abc_set_) return true;	// a b c are already defined in the input

		float_t sqrt2 = sqrt(2.0);
		float_t sqrt3 = sqrt(3.0);

		switch(type_) {
			case lattice_bcc:					// BCC
				if(hkl_ == "100") {
					a_[0] = 1; a_[1] = 0; a_[2] = 0;
					b_[0] = 0; b_[1] = 1; b_[2] = 0;
					c_[0] = 0; b_[1] = 0; c_[2] = 1;
					t_[0] = 0.5; t_[1] = 0.5; t_[2] = 0.5;
				} else if(hkl_ == "110") {
					a_[0] = 1; a_[1] = 0; a_[2] = 0;
					b_[0] = 0.5; b_[1] = 0.5 / sqrt2; b_[2] = 0;
					c_[0] = 0; b_[1] = 0; c_[2] = 1;
					t_[0] = 0.5; t_[1] = 0; t_[2] = 0.5;
				} else {
					std::cerr << "error: invalid lattice type and hkl combination" << std::endl;
					return false;
				} // if-else
				break;

			case lattice_cubic:					// CUBIC
				a_[0] = 1; a_[1] = 0; a_[2] = 0;
				b_[0] = 0; b_[1] = 1; b_[2] = 0;
				c_[0] = 0; c_[1] = 0; c_[2] = 1;
				t_[0] = 0; t_[1] = 0; t_[2] = 0;
				break;

			case lattice_fcc:					// FCC
				if(hkl_ == "100") {
					a_[0] = 1; a_[1] = 0; a_[2] = 0;
					b_[0] = 0.5; b_[1] = 0.5; b_[2] = 0;
					c_[0] = 0; b_[1] = 0; c_[2] = 1;
					t_[0] = 0.5; t_[1] = 0; t_[2] = 0.5;
				} else if(hkl_ == "111") {
					a_[0] = 1; a_[1] = 0; a_[2] = 0;
					b_[0] = 0.5; b_[1] = 1.0 / (2.0 * sqrt3); b_[2] = 0;
					c_[0] = 0; b_[1] = 0; c_[2] = 2.0 / sqrt3;
					t_[0] = sqrt3 / 3; t_[1] = 0; t_[2] = 1 / sqrt3;
				} else {
					std::cerr << "error: invalid lattice type and khl combination" << std::endl;
					return false;
				} // if-else
				break;

			case lattice_fco:					// FCO
				a_[0] = 1; a_[1] = 0; a_[2] = 0;
				b_[0] = 0.5; b_[1] = tan(gamma_) / 2.0; b_[2] = 0;
				c_[0] = 0; b_[1] = 0; c_[2] = ca_;
				t_[0] = 0.5; t_[1] = 0; t_[2] = ca_ / 2;
				break;

			case lattice_hcp:					// HCP
				a_[0] = 1; a_[1] = 0; a_[2] = 0;
				b_[0] = 0.5; b_[1] = 1 / (2.0 * sqrt3); b_[2] = 0;
				c_[0] = 0; b_[1] = 0; c_[2] = ca_ / sqrt3;
				t_[0] = sqrt3 / 3.0; t_[1] = 0; t_[2] = ca_ / (2 * sqrt3);
				break;

			case lattice_hex:					// HEX
				a_[0] = 1; a_[1] = 0; a_[2] = 0;
				b_[0] = 0.5; b_[1] = 1 / (2.0 * sqrt3); b_[2] = 0;		// check ...
				c_[0] = 0; c_[1] = 0; c_[2] = 1;
				t_[0] = 0; t_[1] = 0; t_[2] = 0;
				break;

			case lattice_null:					// means custom vectors are defined
				t_[0] = t_[1] = t_[2] = 0;
				break;

			case lattice_error:
				std::cerr << "error: lattice type already has some previous error" << std::endl;
				return false;

			default:
				std::cerr << "error: lattice type has unknown value" << std::endl;
				return false;
		} // switch

		// scale the vector values
		a_[0] *= scaling; a_[1] *= scaling; a_[2] *= scaling;
		b_[0] *= scaling; b_[1] *= scaling; b_[2] *= scaling;
		c_[0] *= scaling; c_[1] *= scaling; c_[2] *= scaling;
		t_[0] *= scaling; t_[1] *= scaling; t_[2] *= scaling;

		return true;
	} // Lattice::construct_vectors()



	/** grainorientations functions
	 */

	GrainOrientations::GrainOrientations() { }
	GrainOrientations::~GrainOrientations() { }


	void GrainOrientations::init() {
		stat_ = "single";
		rot1_.axis('x');
		rot1_.angles(0, 0);
		rot2_.axis('n');
		rot3_.axis('n');		// for null. all could be same as rot1 too ...
	} // GrainOrientations::init()


	void GrainOrientations::clear() {
		stat_.clear();
		rot1_.axis('n'); rot2_.axis('n'); rot3_.axis('n');
		rot1_.angles(0, 0); rot2_.angles(0, 0); rot3_.angles(0, 0);
	} // GrainOrientations::clear()



	/** grain functions
	 */

	Grain::Grain() { }
	Grain::~Grain() { }


	void Grain::init() {
		layer_key_ = "";
		in_layer_ = false;
		lattice_.init();
		transvec_[0] = transvec_[1] = transvec_[2] = 0;
		scaling_ = 1;
		repetition_[0] = repetition_[1] = repetition_[2] = 1;
	} // Grain::init()


	void Grain::clear() {
		shape_key_.clear();
		layer_key_.clear();
		in_layer_ = false;
		scaling_ = 0.0;
		lattice_.clear();
	} // Grain::clear()



	/** ensemble functions
	 */

	Ensemble::Ensemble() { }
	Ensemble::~Ensemble() { }


	void Ensemble::init() {
		spacing_[0] = spacing_[1] = spacing_[2] = 0;
		maxgrains_[0] = maxgrains_[1] =	maxgrains_[2] = 1;
		distribution_ = "regular";
		orientations_.init();
	} // Ensemble::init()


	void Ensemble::clear() {
		spacing_[0] = spacing_[1] = spacing_[2] = 0;
		maxgrains_[0] = maxgrains_[1] =	maxgrains_[2] = 0;
		orientations_.clear();
	} // Ensemble::clear()



	/** structure functions
	 */

	Structure::Structure() { }
	Structure::~Structure() { }


	void Structure::init() {
		key_ = "";
		grain_.init();
		ensemble_.init();
	} // Structure::init()


	void Structure::clear() {
		key_.clear();
		grain_.clear();
		ensemble_.clear();
	} // Structure::clear()


	/** printing for testing
	 */

	void Structure::print() {
		std::cout << " key_ = " << key_ << std::endl;
		std::cout << " grain_: " << std::endl
					<< "  shape_key_ = " << grain_.shape_key_ << std::endl
					<< "  layer_ley_ = " << grain_.layer_key_ << std::endl
					<< "  scaling_ = " << grain_.scaling_ << std::endl
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

} // namespace hig
