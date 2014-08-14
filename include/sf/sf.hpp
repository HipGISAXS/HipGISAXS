/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: sf.hpp
 *  Created: Jun 18, 2012
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

#ifndef __HIG_SF_HPP__
#define __HIG_SF_HPP__

#ifdef USE_MPI
#include <mpi.h>
#include <woo/comm/multi_node_comm.hpp>
#endif // USE_MPI

#include <common/typedefs.hpp>
#include <common/globals.hpp>
#include <model/structure.hpp>

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


#endif /* __HIG_SF_HPP__ */
