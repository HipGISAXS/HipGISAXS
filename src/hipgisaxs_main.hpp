/***
  *  $Id: hipgisaxs_main.hpp 47 2012-08-23 21:05:16Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: hipgisaxs_main.hpp
  *  Created: Jun 11, 2012
  *  Modified: Mon 17 Jun 2013 01:59:00 PM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _HIPGISAXS_MAIN_HPP_
#define _HIPGISAXS_MAIN_HPP_

#include <complex>
#include <mpi.h>

#include "typedefs.hpp"
#include "globals.hpp"
#include "hig_input.hpp"
#include "common.hpp"
#include "qgrid.hpp"
#include "ff.hpp"
#include "sf.hpp"
#include "image.hpp"

//#include "hipgisaxs_error.hpp"

namespace hig {

	class HipGISAXSError;

	/* This class provides high level functions for all the
	 * components underneath in the HipGISAXS package. */

	class HipGISAXS {
		private:
			//////////////////////////////////////////////////////////
			// temporary ... FIXME
			// ///////////////////////////////////////////////////////
			float_t* final_data_;
			float_t* reference_data_;
			//std::vector<float_t> reference_data_;
			/////////////////////////////////////////////////////////

			float_t freq_;
			float_t k0_;

			int num_layers_;
			int num_structures_;

			RefractiveIndex substrate_refindex_;		// can these be avoided? ...
			RefractiveIndex single_layer_refindex_;		// or can these be filled while reading input? ...
			float_t single_layer_thickness_;

			complex_t dns2_; 	// find a better name ...

			vector3_t vol_;			/* represents the illuminated volume */
			vector3_t cell_;		/* stores the domain size in each dimension */

			class SampleRotation {
				friend class HipGISAXS;

				private:
					vector3_t r1_, r2_, r3_;

					SampleRotation():
						r1_(0.0, 0.0, 0.0),
						r2_(0.0, 0.0, 0.0),
						r3_(0.0, 0.0, 0.0) {
					} //SampleRotation()
			} rotation_matrix_;

			unsigned int nqx_;			/* number of q-points along x */
			unsigned int nqy_;			/* number of q-points along y */
			unsigned int nqz_;			/* number of q-points along z */
			unsigned int nqz_extended_;	/* number of q-points along z in case of gisaxs */

			complex_t* fc_;				/* fresnel coefficients */
			FormFactor ff_;				/* form factor object */
			StructureFactor sf_;		/* structure factor object */

			bool init(MPI::Intracomm&);	/* global initialization for all runs */
			bool init_steepest_fit(MPI::Intracomm&, float_t);	/* init for steepest descent fitting ... */
			bool run_init(float_t, float_t, float_t, MPI::Intracomm&); 	/* initialization for a single run */
			bool run_gisaxs(float_t, float_t, float_t, float_t, float_t*&, MPI::Intracomm&, int c = 0);
										/* a single GISAXS run */

			/* wrapper over sf function */
			bool structure_factor(std::string, vector3_t&, Lattice*&, vector3_t&,
									vector3_t&, vector3_t&, vector3_t&,
									MPI::Intracomm&);

			/* wrapper over ff function */
			bool form_factor(ShapeName, std::string, shape_param_list_t&, vector3_t&,
							float_t, float_t, vector3_t&, vector3_t&, vector3_t&,
							MPI::Intracomm&);

			bool layer_qgrid_qz(float_t, complex_t);
			bool compute_propagation_coefficients(float_t, complex_t*&, complex_t*&,
									complex_t*&, complex_t*&, complex_t*&, complex_t*&,
									complex_t*&, complex_t*&, complex_t*&);
			bool compute_fresnel_coefficients_embedded(float_t);
			bool compute_fresnel_coefficients_top_buried(float_t);
			bool compute_rotation_matrix_x(float_t, vector3_t&, vector3_t&, vector3_t&);
			bool compute_rotation_matrix_y(float_t, vector3_t&, vector3_t&, vector3_t&);
			bool compute_rotation_matrix_z(float_t, vector3_t&, vector3_t&, vector3_t&);

			void save_gisaxs(float_t *final_data, std::string output);

			bool illuminated_volume(float_t, float_t, int, RefractiveIndex);
			bool spatial_distribution(structure_iterator_t, float_t, int, int&, int&, float_t*&);
			bool orientation_distribution(structure_iterator_t, float_t*, int, int, float_t*&);

			/* some functions just for testing and debugging */
			bool write_qgrid(char* filename);
			bool read_form_factor(const char* filename);
			void printfc(const char*, complex_t*, unsigned int);
			void printfr(const char*, float_t*, unsigned int);

		public:
			HipGISAXS();
			~HipGISAXS();

			bool construct_input(char* filename) {
				return HiGInput::instance().construct_input_config(filename);
			} // construct_input()

			/* loops over all configs and computes GISAXS for each */
			bool run_all_gisaxs(MPI::Intracomm&, int = 0, int = 0, int = 0);


			/***
			 * the following are to incorporate fitting
			 */

			////////////////////////////////////////////////
			// temporary stuff ... FIXME
			// /////////////////////////////////////////////

			// FIXME ... these are all just temporary fixes ...

			/* get map of parameters that are being fitted */
			//FitParamMap& get_fit_parameters_map() { return HiGInput::instance().get_fit_params_map(); }

			float_t operator[](unsigned int i) { return final_data_[i]; }

			bool update_fit_param(std::string key1, float_t val1, std::string key2, float_t val2) {
				return HiGInput::instance().update_fit_param(key1, val1, key2, val2);
			} // update_fit_param()

			bool update_fit_param(std::string key1, float_t val1, std::string key2, float_t val2,
									std::string key3, float_t val3, std::string key4, float_t val4) {
				return HiGInput::instance().update_fit_param(key1, val1, key2, val2, key3, val3, key4, val4);
			} // update_fit_param()

			/* update values of parameters being fitted */
			//bool update_fit_parameters(FitParamMap*);

			/* to compute gisaxs pattern, and error with base input */
			template <typename ErrorFunction>
			float_t compute_gisaxs_error(ErrorFunction error, MPI::Intracomm& comm) {
				if(!run_all_gisaxs(comm)) {
					std::cerr << "error: failed to compute gisaxs" << std::endl;
					return 0;
				} // if
				return error(&(reference_data_[0]), final_data_, nqx_, nqy_, nqz_);
			} // compute_gisaxs_error()

			bool read_reference_gisaxs_data(const char* filename) {
				//reference_data_.clear();
				//reference_data_.reserve(300000000);
				float_t *temp_data = new (std::nothrow) float_t[9000000];
				std::ifstream fin(filename);
				if(!fin.is_open()) {
					std::cerr << "error: failed to read reference data" << std::endl;
					return false; 
				} // if
				unsigned int i = 0;
				while(true) {
					float_t temp;
					fin >> temp;
					if(fin.eof() || !fin.good()) break;
					//reference_data_.push_back(temp);
					temp_data[i ++] = temp;
				} // while
				fin.close();
				//std::vector <float_t> hmmm(temp_data, temp_data + i);
				reference_data_ = temp_data;
				temp_data = NULL;
				return true;
			} // read_reference_gisaxs_data()


			/* temporary -- does the initial 1D fitting */
			bool fit_steepest_descent(float_t, float_t, float_t, float_t,
										float_t, float_t, float_t, unsigned int,
										MPI::Intracomm&, int = 0, int = 0, int = 0);
			float_t compute_cut_fit_error(float_t*, float_t*, float_t);

	}; // class HipGISAXS

} // namespace hig

#endif /* _HIPGISAXS_MAIN_HPP_ */
