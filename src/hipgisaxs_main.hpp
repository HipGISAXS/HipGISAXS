/***
  *  $Id: hipgisaxs_main.hpp 47 2012-08-23 21:05:16Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: hipgisaxs_main.hpp
  *  Created: Jun 11, 2012
  *  Modified: Fri 07 Dec 2012 07:46:51 PM PST
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

namespace hig {

	/* This class provides high level functions for all the
	 * components underneath in the HipGISAXS package. */

	class HipGISAXS {
		private:
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
			bool run_init(float_t, float_t, float_t, MPI::Intracomm&); 	/* initialization for a single run */
			bool run_gisaxs(float_t, float_t, float_t, float_t, float_t*&, MPI::Intracomm&, int c = 0);
										/* a single GISAXS run */

			/* wrapper over sf function */
			bool structure_factor(std::string, vector3_t&, Lattice*&, vector3_t&,
									vector3_t&, vector3_t&, vector3_t&,
									MPI::Intracomm&);

			/* wrapper over ff function */
			//template <typename float_t, typename complex_t>
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
			HipGISAXS(): freq_(0.0), k0_(0.0),
						num_layers_(0), num_structures_(0),
						nqx_(0), nqy_(0), nqz_(0), nqz_extended_(0), 
#ifdef KERNEL2
						ff_(2, 4, 4)
#else
						ff_(64)
#endif // KERNEL2
						{
				single_layer_refindex_.delta(0.0);
				single_layer_refindex_.beta(0.0);
				single_layer_thickness_ = 0.0;
				HiGInput::instance();
				QGrid::instance();
				//MPI::Init();
			} // HipGISAXS()


			~HipGISAXS() {
				//MPI::Finalize();
			} // ~HipGISAXS()


			bool construct_input(char* filename) {
				return HiGInput::instance().construct_input_config(filename);
			} // construct_input()

			/* loops over all configs and computes GISAXS for each */
			bool run_all_gisaxs(MPI::Intracomm&, int = 0, int = 0, int = 0);

			/* does the initial 1D fitting */
			bool fit_steepest_descent(unsigned int, float_t, float_t, float_t,
										float_t, float_t, float_t, unsigned int,
										MPI::Intracomm&, int = 0, int = 0, int = 0);
			float_t compute_cut_fit_error(float_t*, float_t*, float_t);

	}; // class HipGISAXS

} // namespace hig

#endif /* _HIPGISAXS_MAIN_HPP_ */
