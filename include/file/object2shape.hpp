/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: object2shape.hpp
 *  Created: Aug 25, 2012
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

#ifndef _OBJECT2SHAPE_H_
#define _OBJECT2SHAPE_H_

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <boost/tokenizer.hpp>

extern "C" {
void s2h_converter(real_t** shape_def, unsigned int num_triangles, char* hdf5_filename, MPI_Comm comm);
}

enum token_t {
  COMMENT,
  VERTEX,
  TEXTURE,
  SUB_MESH,
  MATERIAL_LIBRARY,
  MATERIAL_NAME,
  LINE,
  SMOOTH_SHADING,
  NORMAL,
  FACE,
  UNKNOWN
}; //enum

typedef struct {
  real_t x;
  real_t y;
  real_t z;
  real_t w;   // for format's completeness
} vertex_t;

typedef struct {
  int a;
  int b;
  int c;
  int d;
} poly_index_t;

typedef boost::char_separator<char> token_separator_t;
typedef boost::tokenizer<token_separator_t> tokenizer_t;

class o2s_converter {
  public:
    o2s_converter(char* filename, char* outfilename, MPI_Comm comm, bool hdf5 = true);
    ~o2s_converter() {
      if(filename_ != NULL) delete filename_;
      if(outfilename_ != NULL) delete outfilename_;
      if(shape_def_ != NULL) delete[] shape_def_;
    } // ~converter()

  private:
    void load_object(char* filename, std::vector<vertex_t> &vertices,
        std::vector<std::vector<int> > &face_list_3v,
        std::vector<std::vector<int> > &face_list_4v);
    //void load_object(char* filename, std::vector<vertex_t> &vertices,
    //    std::vector<poly_index_t> &Vn3, std::vector<poly_index_t> &Vn4,
    //    std::vector<poly_index_t> &F3, std::vector<poly_index_t> &F4);
    //real_t* convert(char* outfilename, std::vector<poly_index_t> F3,
    //    std::vector<vertex_t> vertices, bool hdf5);
    real_t* convert(char* outfilename,
        std::vector<std::vector<int> > face_list_3v,
        std::vector<vertex_t> vertices, bool hdf5);

    bool get_triangle_params(vertex_t v1, vertex_t v2, vertex_t v3,
        real_t &s_area, vertex_t &normal, vertex_t &center);
    token_t token_hash(std::string const &str);
    void findall(std::string str, char c, std::vector<int> &pos_list);
    void display_vertices(std::vector<vertex_t> &vertices);
    void display_poly_index(std::vector<poly_index_t> &indices);

    std::string *filename_;
    std::string *outfilename_;
    real_t* shape_def_;
    MPI_Comm comm_;
}; // class o2s_converter

#endif // _OBJECT2SHAPE_H_
