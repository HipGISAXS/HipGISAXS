/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: objectshape_reader.hpp
 *  Created: Aug 25, 2013
 *  Modified: Sun 29 Sep 2013 11:26:59 AM PDT
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

#ifndef __OBJECTSHAPE_READER_HPP__
#define __OBJECTSHAPE_READER_HPP__

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <boost/tokenizer.hpp>

namespace hig {

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
		float_t x;
		float_t y;
		float_t z;
		float_t w; 	// for format's completeness
	} vertex_t;

	typedef struct {
		int a;
		int b;
		int c;
		int d;
	} poly_index_t;

	typedef boost::char_separator<char> token_separator_t;
	typedef boost::tokenizer<token_separator_t> tokenizer_t;

	class ObjectShapeReader {
		public:
			ObjectShapeReader(const char* filename, std::vector<float_t>& shape_def);
			~ObjectShapeReader() { }

		private:
			bool load_object(const char* filename, std::vector<vertex_t> &vertices,
								std::vector<std::vector<int> > &face_list_3v,
								std::vector<std::vector<int> > &face_list_4v);
			float_t* convert_to_shape(std::vector<std::vector<int> > face_list_3v,
								std::vector<vertex_t> vertices);
	
			bool get_triangle_params(vertex_t v1, vertex_t v2, vertex_t v3,
								float_t &s_area, vertex_t &normal, vertex_t &center);
			token_t token_hash(std::string const &str);
			void findall(std::string str, char c, std::vector<int> &pos_list);
			void display_vertices(std::vector<vertex_t> &vertices);
			void display_poly_index(std::vector<poly_index_t> &indices);
	}; // class ObjectShapeReader

} // namespace hig

#endif // __OBJECTSHAPE_READER_HPP__
