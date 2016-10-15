/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: rawshape_reader.hpp
 *  Created: Aug 25, 2013
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

#ifndef __RAWSHAPE_READER_HPP__
#define __RAWSHAPE_READER_HPP__

#include <iostream>
#include <fstream>
#include <vector>

#include <common/typedefs.hpp>

namespace hig {

  class RawShapeReader {
    public:
      RawShapeReader(const char*, double*&, unsigned int&);
      ~RawShapeReader() { }
    private:
      bool load_raw(const char*, std::vector<real_t>&);
  }; // class RawShapeReader

} // namespace hig

#endif // __RAWSHAPE_READER_HPP__
