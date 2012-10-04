/***
  *  $Id$
  *
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: plot_ff.hpp
  *  Created: Aug 22, 2012
  *  Modified: Mon 01 Oct 2012 11:16:18 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  */

#include <iostream>
#include <fstream>
#include <cmath>

#include "image.hpp"

namespace hig {

	class Plot {
		public:
			Plot(unsigned int nx, unsigned int ny, unsigned int nz):
					nx_(nx), ny_(ny), nz_(nz), plot_(nx, ny, nz) {
				raw_data_ = new (std::nothrow) complex_t[nx * ny * nz];
				mag_data_ = new (std::nothrow) float_t[nx * ny * nz];
			} // Plot()

			~Plot() {
				delete[] mag_data_;
				delete[] raw_data_;
			} // ~Plot()


			bool plot(const char* infile, const char* outfile) {
				read_in(infile);
				convert_to_mag();
				plot_.construct_image(mag_data_);
				plot_.save(std::string(outfile));
				return true;
			} // plot()

		private:
			unsigned int nx_;
			unsigned int ny_;
			unsigned int nz_;
			complex_t *raw_data_;	// complex numbers
			float_t *mag_data_;		// magnitudes of corresponding complex numbers
			Image plot_;


			bool read_in(const char* infile) {
				std::ifstream fin(infile);
				for(unsigned int z = 0; z < nz_; ++ z) {
					for(unsigned int y = 0; y < ny_; ++ y) {
						for(unsigned int x = 0; x < nx_; ++ x) {
							float_t real = 0.0, imag = 0.0;
							unsigned int index = nx_ * ny_ * z + nx_ * y + x;
							fin >> real;
							fin >> imag;
							raw_data_[index] = complex_t(real, imag);
						} // for x
					} // for y
				} // for z
				fin.close();
				return true;
			} // read_in()


			bool convert_to_mag() {
				for(unsigned int i = 0; i < nx_ * ny_ * nz_; ++ i) mag_data_[i] = 0.0;
				for(int p = 0; p < 4; ++ p) {
					for(unsigned int i = 0; i < nx_ * ny_ * (nz_ / 4); ++ i) {
						complex_t temp = raw_data_[nx_ * ny_ * (nz_ / 4) * p + i];
						mag_data_[nx_ * ny_ * (nz_ / 4) * p + i] =
								sqrt(pow(temp.real(), 2.0) + pow(temp.imag(), 2.0));
					} // for
				} // for
				return true;
			} // convert_to_mag()

	}; // Plot

} // namespace
