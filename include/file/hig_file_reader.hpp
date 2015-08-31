/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: hig_file_reader.hpp
 *  Created: Jul 11, 2012
 *  Modified: Wed 08 Oct 2014 12:13:01 PM PDT
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

#ifndef _HIG_FILE_READER_
#define _HIG_FILE_READER_

#include <fstream>

#include <common/typedefs.hpp>
#include <file/objectshape_reader.hpp>

#ifdef USE_PARALLEL_HDF5
#include <file/hdf5shape_reader.h>
#endif


namespace hig {

#ifdef USE_PARALLEL_HDF5

#ifdef __cplusplus
extern "C" {
#endif

  unsigned int c_hdf5_shape_reader(const char* filename, double** shape_def, unsigned int* num_triangles) {
    // improve this later ...
    // for now just call the old function
    h5_shape_reader(filename, shape_def, num_triangles);
    return *num_triangles;
  } // c_hdf5_shape_reader()

#ifdef __cplusplus
} // extern "C"
#endif

#endif // USE_PARALLEL_HDF5

  // this is a singleton stateless class
  class HiGFileReader {
    private:
      HiGFileReader() { }
      HiGFileReader(const HiGFileReader&);
      HiGFileReader& operator=(const HiGFileReader&);

    public:
      static HiGFileReader& instance() {
        static HiGFileReader hig_file_reader;
        return hig_file_reader;
      } // instance()

      #ifdef USE_PARALLEL_HDF5
      unsigned int hdf5_shape_reader(const char* filename,
                    double* &shape_def, unsigned int &num_triangles) {
        return c_hdf5_shape_reader(filename, &shape_def, &num_triangles);
      } // hdf5_shape_reader()
      #endif

      unsigned int shape_shape_reader(const char* filename,
                    std::vector<real_t> &shape_def, unsigned int &num_triangles) {
        std::ifstream f(filename);
        if(!f.is_open()) {
          std::cout << "Cannot open file " << filename << std::endl;
          return 0;
        } // if
        real_t s = 0.0, cx = 0.0, cy = 0.0, cz = 0.0, nx = 0.0, ny = 0.0, nz = 0.0;

        while(true) {
          f >> s;
          if(f.eof() || !f.good()) break;
          f >> nx; f >> ny; f >> nz;
          f >> cx; f >> cy; f >> cz;
          shape_def.push_back(s);
          shape_def.push_back(nx);
          shape_def.push_back(ny);
          shape_def.push_back(nz);
          shape_def.push_back(cx);
          shape_def.push_back(cy);
          shape_def.push_back(cz);
        } // while

        f.close();
        num_triangles = shape_def.size() / 7;
        return num_triangles;
      } // shape_shape_reader()

      unsigned int object_shape_reader(const char* filename,
                    double* &shape_def, unsigned int &num_triangles) {
        ObjectShapeReader obj(filename, shape_def, num_triangles);
        return num_triangles;
      } // obj_shape_reader()

  }; // class HiGFileReader

} // namespace hig


#endif /* _HIG_FILE_READER_ */
