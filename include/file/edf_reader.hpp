/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: edf_reader.hpp
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

#ifndef __EDF_READER_HPP__
#define __EDF_READER_HPP__

#include <fstream>
#include <string>
#include <map>

#include <common/typedefs.hpp>

namespace hig {

  const size_t EDF_CHUNK_SIZE = 512;

  class EDFReader {
    public:
      EDFReader(const char*);
      ~EDFReader() { }

      bool get_data(real_t *&, unsigned int&, unsigned int&);

    private:
      bool is_white_space(char);
      bool get_next_item(char*, int&, const int, char*, int&, bool&, char*, int&, bool&);
      bool header_done(char*, const int, const int);
      bool extract_header(std::ifstream&, char*, size_t);
      bool extract_data(std::ifstream&, char*, size_t);
      unsigned int convert_to_unsigned(const std::string&);
      void print_header(void);
      void print_data(void);

    private:
      std::map <std::string, std::string> header_;  /* header key -> value map */
      std::vector <real_t> data_;          /* the data as real_ts */
      unsigned int rows_;
      unsigned int cols_;

  }; // class EDFReader


  class EDFWriter {
    public:
      EDFWriter(const char * name) :filename_(name), pixel_(172.E-06) {}
      void setCenter(real_t x, real_t y){ center_x_ = x; center_y_ = y; }
      void setEnergy(real_t e){ energy_ = e; }
      void setPixelSize(real_t px){ pixel_ = px; }
      void setSize(int r, int c){ nrow_ = r; ncol_ = c; }
      void sdd(real_t ); // calculate SDD for Xi-cam (HipIES) compatibility
      void Write(real_t *);

    private:
      const char * filename_;
      int nrow_, ncol_;
      real_t center_x_, center_y_;
      real_t energy_;
      real_t pixel_;
      real_t sdd_;
    
      EDFWriter(){} // forbid default constructor

  }; // EDFWriter

} // namespace hig

#endif // __EDF_READER_HPP__
