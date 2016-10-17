/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: plot_gisaxs.hpp
 *  Created: Aug 22, 2012
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 *
 *  Licensing: The HipGISAXS software is only available to be downloaded and
 *  used by employees of academic research institutions, not-for-profit
 *  research laboratories, or governmental research facilities. Please read the
 *  accompanying LICENSE file before downloading the software. By downloading
 *  the software, you are agreeing to be bound by the terms of this
 *  NON-COMMERCIAL END USER LICENSE AGREEMENT.
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include <image/image.hpp>

namespace hig {

  class PlotGISAXS {
    public:
//      PlotGISAXS(unsigned int nx, unsigned int ny, unsigned int nz):
//          nx_(nx), ny_(ny), nz_(nz), plot_(nx, ny, nz) {
//        raw_data_ = new (std::nothrow) real_t[nx * ny * nz];
//      } // Plot()

      PlotGISAXS(unsigned int nx, unsigned int ny, unsigned int nz,
          unsigned int r, unsigned int g, unsigned int b):
          nx_(nx), ny_(ny), nz_(nz), plot_(nx, ny, nz, r, g, b) {
        raw_data_ = new (std::nothrow) real_t[nx * ny * nz];
      } // Plot()

      ~PlotGISAXS() {
        delete[] raw_data_;
      } // ~Plot()


      bool plot(const char* infile, const char* outfile) {
        read_in(infile);
        plot_.construct_image(raw_data_);
        plot_.save(std::string(outfile));
        return true;
      } // plot()

    private:
      unsigned int nx_;
      unsigned int ny_;
      unsigned int nz_;
      real_t* raw_data_;
      Image plot_;


      bool read_in(const char* infile) {
        std::ifstream fin(infile);
        for(unsigned int z = 0; z < nz_; ++ z) {
          for(unsigned int y = 0; y < ny_; ++ y) {
            for(unsigned int x = 0; x < nx_; ++ x) {
              unsigned int index = nx_ * ny_ * z + nx_ * y + x;
              fin >> raw_data_[index];
            } // for x
          } // for y
        } // for z
        fin.close();
        return true;
      } // read_in()

  }; // PlotGISAXS

} // namespace
