/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: edf_reader.cpp
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

#include <iostream>
#include <cstdlib>

#include <file/edf_reader.hpp>

namespace hig {

  EDFReader::EDFReader(const char* filename) {
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if(!infile.is_open()) {
      std::cerr << "error: could not open the EDF file " << filename << std::endl;
      exit(1);
    } // if
    char* dataseg = new (std::nothrow) char[EDF_CHUNK_SIZE];  // 512 byte blocks
    // read block of header
    if(!extract_header(infile, dataseg, EDF_CHUNK_SIZE)) exit(1);
    if(!extract_data(infile, dataseg, EDF_CHUNK_SIZE)) exit(1);
    infile.close();
    delete[] dataseg;
    //print_header();
    //print_data();
  } // EDFReader::EDFReader()


  bool EDFReader::is_white_space(char c) {
    if(c == ' '
      || c == '\n'
      || c == '\t'
      || c == '\r') return true;
    return false;
  } // EDFReader::is_white_space()


  bool EDFReader::get_next_item(char* chunk, int& counter, const int max_size,
                  char* key, int& i_k, bool& end_k,
                  char* value, int& i_v, bool& end_v) {
    if(!end_k) {
      // read key
      while(1) {
        char curr = chunk[counter ++];
        //std::cout << counter << " ==> " << curr << std::endl;
        if(curr == '=') { key[i_k] = '\0'; end_k = true; break; }    // end of keyword
        if(!is_white_space(curr)) // skip white spaces
          key[i_k ++] = curr;
        if(counter == max_size) break;
        if(curr == '}') { -- counter; return false; }
      } // while
    } // if
    if(counter == max_size) return false;
    if(!end_v) {
      // read value
      while(1) {
        char curr = chunk[counter ++];
        //std::cout << curr << std::endl;
        if(curr == ';') { value[i_v] = '\0'; end_v = true; break; }    // end of keyword
        if(!is_white_space(curr))  // skip white spaces
          value[i_v ++] = curr;
        if(counter == max_size) break;
      } // while
    } // if
    if(!end_v) return false;
    return true;
  } // EDFReader::get_next_item()


  inline bool EDFReader::header_done(char* chunk, const int counter, const int max_size) {
    int mycount = counter;
    while(1) {
      char curr = chunk[mycount ++];
      if(curr == '}') return true;
      if(mycount == max_size) break;    // reached the end of the chunk
      if(is_white_space(curr)) continue;
      break;  // curr is some other character => more stuff
    } // while
    return false;
  } // EDFReader::header_done()


  bool EDFReader::extract_header(std::ifstream& infile, char* chunk, size_t size) {
    char key[size], value[size];
    int counter = 0, i_k = 0, i_v = 0;
    bool end_k = false, end_v = false;
    infile.read(chunk, size);
    while(1) {
      if(!get_next_item(chunk, counter, size, key, i_k, end_k, value, i_v, end_v)) {
        if(header_done(chunk, counter, size)) break;
        counter = 0;
        infile.read(chunk, size);
        continue;
      } // if
      // we have a new key value pair
      header_[std::string(key)] = std::string(value);
      //std::cout << "****** [" << counter << "] " << key << " -> " << value << std::endl;
      i_k = 0; i_v = 0;
      end_k = false; end_v = false;
      if(header_done(chunk, counter, size)) break;
    } // while

    cols_ = convert_to_unsigned(std::string("Dim_1"));
    rows_ = convert_to_unsigned(std::string("Dim_2"));

    //std::cout << "*** Header reading done. " << std::endl;

    return true;  /* to signal that header was complete */
  } // EDFReader::extract_header()


  inline unsigned int EDFReader::convert_to_unsigned(const std::string& key) {
    return atoi(header_[key].c_str());
  } // EDFReader::convert_to_unsigned()


  bool EDFReader::extract_data(std::ifstream& infile, char* chunk, size_t size) {
    unsigned long int num_bytes = atol(header_[std::string("EDF_BinarySize")].c_str());
    typedef float real_t;
    real_t* data;
    unsigned long int count = 0;
    data_.clear();
    while(count != num_bytes) {
      infile.read(chunk, size);
      data = reinterpret_cast<real_t*>(chunk);
      int len = std::min(size, (num_bytes - count));
      for(int i = 0; i < len / sizeof(real_t); ++ i) data_.push_back((real_t) data[i]);
      count += len;
      //std::cout << "****** " << count << std::endl;
    } // while
    if(data_.size() != cols_ * rows_) {
      std::cerr << "error: mismatch in reference data size" << std::endl;
      std::cerr << "error: cols_ = " << cols_ << ", rows_ = " << rows_ << std::endl;
      std::cerr << "error: data_.size() = " << data_.size() << std::endl;
      return false;
    } // if
    //std::cout << "*** Data reading done." << std::endl;
    return true;
  } // EDFReader::extract_data()


  bool EDFReader::get_data(real_t*& data, unsigned int& ny, unsigned int& nz) {
    if(data_.size() < 1) return false;
    data = &(data_[0]);
    ny = cols_;
    nz = rows_;
    return true;
  } // EDFReader::get_data()
  

  void EDFReader::print_header(void) {
    for(std::map <std::string, std::string> ::iterator i = header_.begin(); i != header_.end(); ++ i)
      std::cout << "-- " << (*i).first << " --> " << (*i).second << std::endl;
  } // EDFReader::print_header()

  void EDFReader::print_data(void) {
    for(unsigned int i = 0; i < rows_; ++ i) {
      for(unsigned int j = 0; j < cols_; ++ j) {
        unsigned int idx = i * cols_ + j;
        std::cout << data_[idx] << " ";
      } // for
      std::cout << std::endl;
    } // for
  } // EDFReader::preint_data()


  /******* EDF Writer ********/
  void EDFWriter::Write(real_t * data){
    std::fstream edf(filename_, std::ios::out | std::ios::binary);
    if (!edf.is_open()){
      std::cerr << "Error: failed to open EDF file for writing" << std::endl;
      return;
    }

    std::string type("float");
    int size = nrow_ * ncol_;
    /* filename */
    edf << "Pixel Size        = " << pixel_    << std::endl;
    edf << "Center X          = " << center_x_ << std::endl;
    edf << "Center Y          = " << center_y_ << std::endl;
    edf << "Detector Distance = " << sdd_      << std::endl;
    edf << "Energy            = " << energy_   << std::endl;
    edf << "Dim_1             = " << nrow_     << std::endl;
    edf << "Dim_2             = " << ncol_     << std::endl;
    edf << "Size              = " << size * sizeof(float)  << std::endl;
    edf << "DataType          = " << type      << std::endl;
    
    // First 1024 bytes are reserved for header
    edf.seekp(1024);
    for (int i = 0; i < size; i++ ){
      float d = (float) data[i];
      edf.write(reinterpret_cast<char *>(&d), sizeof(float));
    }
    edf.close();
  }

  void EDFWriter::sdd(real_t alpha){
    real_t tan_a = std::tan(alpha);
    sdd_ = (ncol_ - center_y_) * pixel_ / tan_a;
  }

} // namespace hig
