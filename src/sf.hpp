/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf.hpp
 *  Created: Jun 18, 2012
 *  Modified: Wed 18 Sep 2013 11:33:52 AM PDT
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

#ifndef _SF_HPP_
#define _SF_HPP_

#ifdef USE_MPI
#include <mpi.h>
#include "woo/comm/multi_node_comm.hpp"
#endif // USE_MPI

#include "typedefs.hpp"
#include "globals.hpp"
#include "structure.hpp"

namespace hig {
	
	class StructureFactor {
		private:
			complex_t *sf_;
			unsigned int nx_;
			unsigned int ny_;
			unsigned int nz_;

		public:
			StructureFactor();
			~StructureFactor();

			void clear(void);
			bool compute_structure_factor(std::string, vector3_t, Lattice*, vector3_t,
											vector3_t, vector3_t, vector3_t
											#ifdef USE_MPI
												, woo::MultiNode&, const char*
											#endif
											);
			bool compute_structure_factor_gpu(std::string, vector3_t, Lattice*, vector3_t,
											vector3_t, vector3_t, vector3_t
											#ifdef USE_MPI
												, woo::MultiNode&, const char*
											#endif
											);
			complex_t operator[](unsigned int i) const { return sf_[i]; }
			void save_sf(unsigned int nqx, unsigned int nqy, unsigned int nqz, const char* filename);

			// only for testing - remove it ...
			complex_t* sf() { return sf_; }
			void printsf() {
				std::cout << "sf:" << std::endl;
				for(unsigned int i = 0; i < nx_ * ny_ * nz_; ++ i) {
					std::cout << sf_[i].real() << "," << sf_[i].imag() << "\t";
				} // for
				std::cout << std::endl;
			} // printsf()
	}; // class StructureFactor

} // namespace hig


#endif /* _SF_HPP_ */
