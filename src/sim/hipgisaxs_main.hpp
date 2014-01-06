/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hipgisaxs_main.hpp
 *  Created: Jun 11, 2012
 *  Modified: Sat 28 Dec 2013 09:17:29 AM PST
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

#ifndef _HIPGISAXS_MAIN_HPP_
#define _HIPGISAXS_MAIN_HPP_

#include <complex>
#ifdef USE_MPI
#include <mpi.h>
#include "../woo/comm/multi_node_comm.hpp"
#endif

#include "../common/typedefs.hpp"
#include "../common/globals.hpp"
#include "../config/hig_input.hpp"
#include "../model/common.hpp"
#include "../model/qgrid.hpp"
#include "../ff/ff.hpp"
#include "../sf/sf.hpp"
#include "../image/image.hpp"

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
			}; // rotation_matrix_;

			unsigned int nqx_;			/* number of q-points along x */
			unsigned int nqy_;			/* number of q-points along y */
			unsigned int nqz_;			/* number of q-points along z */
			unsigned int nqz_extended_;	/* number of q-points along z in case of gisaxs */

//			complex_t* fc_;				/* fresnel coefficients */
//			FormFactor ff_;				/* form factor object */
//			StructureFactor sf_;		/* structure factor object */

			#ifdef USE_MPI
				woo::MultiNode multi_node_;	/* for multi node communication */
			#endif

			bool init();	/* global initialization for all runs */
			bool init_steepest_fit(float_t);	/* init for steepest descent fitting */
			bool run_init(float_t, float_t, float_t, SampleRotation&); 	/* init for a single run */
			bool run_gisaxs(float_t, float_t, float_t, float_t, float_t*&, const char*, int c = 0);
										/* a single GISAXS run */

			/* wrapper over sf function */
			bool structure_factor(StructureFactor&, std::string, vector3_t&, Lattice*&, vector3_t&,
									vector3_t&, vector3_t&, vector3_t&
									#ifdef USE_MPI
										, const char*
									#endif
									);

			/* wrapper over ff function */
			bool form_factor(FormFactor&, ShapeName, std::string, shape_param_list_t&, vector3_t&,
									float_t, float_t, vector3_t&, vector3_t&, vector3_t&
									#ifdef USE_MPI
										, const char*
									#endif
									);

			bool layer_qgrid_qz(float_t, complex_t);
			bool compute_propagation_coefficients(float_t, complex_t*&, complex_t*&,
									complex_t*&, complex_t*&, complex_t*&, complex_t*&,
									complex_t*&, complex_t*&, complex_t*&, complex_t*&);
			bool compute_fresnel_coefficients_embedded(float_t, complex_t*&);
			bool compute_fresnel_coefficients_top_buried(float_t, complex_t*&);

			bool compute_rotation_matrix_x(float_t, vector3_t&, vector3_t&, vector3_t&);
			bool compute_rotation_matrix_y(float_t, vector3_t&, vector3_t&, vector3_t&);
			bool compute_rotation_matrix_z(float_t, vector3_t&, vector3_t&, vector3_t&);

			void save_gisaxs(float_t *final_data, std::string output);

			bool illuminated_volume(float_t, float_t, int, RefractiveIndex);
			bool spatial_distribution(structure_iterator_t, float_t, int, int&, int&, float_t*&);
			bool orientation_distribution(structure_iterator_t, float_t*, int, int, float_t*&);

			/* some functions just for testing and debugging */
			bool write_qgrid(char* filename);
			bool read_form_factor(FormFactor&, const char* filename);
			void printfc(const char*, complex_t*, unsigned int);
			void printfr(const char*, float_t*, unsigned int);

		public:
			HipGISAXS(int, char**);
			~HipGISAXS();

			bool construct_input(char* filename) {
				return HiGInput::instance().construct_input_config(filename);
			} // construct_input()

			/* loops over all configs and computes GISAXS for each */
			bool run_all_gisaxs(int = 0, int = 0, int = 0);

			/* update parameters given a map:key->value */
			bool update_params(const map_t&);

			/* temporary: does the initial 1D fitting */
			bool fit_steepest_descent(float_t, float_t, float_t, float_t,
										float_t, float_t, float_t, unsigned int,
										MPI::Intracomm&, int = 0, int = 0, int = 0);
			float_t compute_cut_fit_error(float_t*, float_t*, float_t);

	}; // class HipGISAXS

} // namespace hig

#endif /* _HIPGISAXS_MAIN_HPP_ */
