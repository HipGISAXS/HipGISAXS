/***
  *  $Id: sf.hpp 46 2012-08-23 02:01:21Z asarje $
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: sf.hpp
  *  Created: Jun 18, 2012
  *  Modified: Mon 01 Oct 2012 11:17:28 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#ifndef _SF_HPP_
#define _SF_HPP_

#include <mpi.h>

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

			bool compute_structure_factor(std::string, vector3_t, Lattice*, vector3_t,
											vector3_t, vector3_t, vector3_t, MPI::Intracomm&);
			complex_t operator[](unsigned int i) const { return sf_[i]; }

			// only for testing - remove it ...
			complex_t* sf() { return sf_; }
			void printsf() {
				std::cout << "sf:" << std::endl;
				for(unsigned int i = 0; i < nx_ * ny_ * nz_; ++ i) {
					std::cout << sf_[i].real() << "," << sf_[i].imag() << "\t";
				} // for
				std::cout << std::endl;
			} // printsf()
			void save_sf(unsigned int nqx, unsigned int nqy, unsigned int nqz, const char* filename);
	}; // class StructureFactor

} // namespace hig


#endif /* _SF_HPP_ */
