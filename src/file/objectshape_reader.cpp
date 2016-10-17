/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: objectshape_reader.cpp
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

#include <iomanip>
#include <boost/lexical_cast.hpp>

#include <file/objectshape_reader.hpp>

namespace hig {

  /**
   * constructor: read input object file, convert to shape
   */
  ObjectShapeReader::ObjectShapeReader(const char* filename, double* &shape_def,
                      unsigned int& num_triangles) {
    std::vector<vertex_t> vertices;
    //std::vector<poly_index_t> Vn3, Vn4, F3, F4;
    std::vector<std::vector<int> > face_list_3v, face_list_4v;
    load_object(filename, vertices, face_list_3v, face_list_4v);

    //load_object(filename, vertices, Vn3, Vn4, F3, F4);
    //std::cout << "Vertices: " << std::endl;  display_vertices(vertices);
    //std::cout << "F3: " << std::endl;    display_poly_index(F3);
    //std::cout << "F4: " << std::endl;    display_poly_index(F4);

    if(face_list_4v.size() > 0)
      std::cerr << "warning: ignoring faces with four vertices. may be fix it later" << std::endl;
    std::vector<real_t> temp_shape_def;
    convert_to_shape(face_list_3v, vertices, temp_shape_def);
    num_triangles = temp_shape_def.size() / 7;
    shape_def = new (std::nothrow) double[temp_shape_def.size()];
    for(int i = 0; i < temp_shape_def.size(); ++ i) shape_def[i] = temp_shape_def[i];
  } // o2s_converter()


  /**
   * load input object file into local data structures in memory
   */
  bool ObjectShapeReader::load_object(const char* filename,
      std::vector<vertex_t> &vertices,
      /*std::vector<poly_index_t> &Vn3, std::vector<poly_index_t> &Vn4,
      std::vector<poly_index_t> &F3, std::vector<poly_index_t> &F4,*/
      std::vector<std::vector<int> > &face_list_3v,
      std::vector<std::vector<int> > &face_list_4v) {

    std::ifstream input(filename);
    if(!input.is_open()) {
      std::cout << "Unable to open file " << filename << std::endl;
      return false;
    } // if

    int num_texture = 0, num_normal = 0;
    //std::vector<int> g3num, g4num;
    //int f3num = 0, f4num = 0;
    //char discard;
    //real_t data[12], tc1[4], vn1[4];
    //poly_index_t f1;

    std::vector<int> face_vi, face_vti, face_vni;
    int face_num_3v, face_num_4v;
    tokenizer_t::iterator iter1, iter2;
    char* temp_line = new char[257];  // will not work if line has more than 256 characters

    while(!input.eof()) {
      input.getline(temp_line, 256);
      std::string temp_string(temp_line);

      temp_string.erase(std::remove(temp_string.begin(), temp_string.end(), '\r'),
                temp_string.end());
      if(temp_string.empty()) continue;

      std::stringstream temp_stream(temp_string);
      std::string token;
      temp_stream >> token;

      std::vector<int> pos_face_vertices, pos_slash;
      //int num_face_vertices = 0, num_slash = 0, dbslash = 0;

      token_separator_t sep1(" ");
      tokenizer_t tokens1(temp_string, sep1);

      switch(token_hash(token)) {
        case VERTEX:
          vertex_t v;
          temp_stream >> v.x;
          temp_stream >> v.y;
          temp_stream >> v.z;
          vertices.push_back(v);
          break;

        case TEXTURE:    // not using
          ++ num_texture;
          break;

        case SUB_MESH:    // not using
          //g3num.push_back(face_num_3v);
          //g4num.push_back(face_num_4v);
          break;

        case NORMAL:    // not using
          ++ num_normal;
          break;

        /*case FACE:
          pos_face_vertices.clear();
          pos_slash.clear();
          findall(temp_string, ' ', pos_face_vertices);
          num_face_vertices = pos_face_vertices.size();
          findall(temp_string, '/', pos_slash);
          num_slash = pos_slash.size();

          if(num_slash > 1)
            dbslash = pos_slash[1] - pos_slash[0];
          else dbslash = 0;

          if(num_slash == 0) {
            temp_stream >> f1.a;
          } else if(num_slash == num_face_vertices && dbslash > 1) {
            temp_stream >> data[0] >> discard >> data[1];
            temp_stream >> data[2] >> discard >> data[3];
            temp_stream >> data[4] >> discard >> data[5];
            if(num_face_vertices == 3) {
              f1.a = data[0]; f1.b = data[2]; f1.c = data[4];
              tc1[0] = data[1]; tc1[1] = data[3]; tc1[2] = data[5];
            } // if
            if(num_face_vertices == 4) {
              temp_stream >> data[6] >> discard >> data[7];
              f1.a = data[0]; f1.b = data[2]; f1.c = data[4]; f1.d = data[6];
              tc1[0] = data[1]; tc1[1] = data[3]; tc1[2] = data[5]; tc1[3] = data[7];
            } // if
          } else if(num_slash == 2 * num_face_vertices && dbslash == 1) {
            temp_stream >> data[0] >> discard >> discard >> data[1];
            temp_stream >> data[2] >> discard >> discard >> data[3];
            temp_stream >> data[4] >> discard >> discard >> data[5];
            if(num_face_vertices == 3) {
              f1.a = data[0]; f1.b = data[2]; f1.c = data[4]; f1.d = -1;
              vn1[0] = data[1]; vn1[1] = data[3]; vn1[2] = data[5]; vn1[3] = -1;
              Vn3.push_back(f1);
            }
            if(num_face_vertices == 4) {
              temp_stream >> data[6] >> discard >> discard >> data[7];
              f1.a = data[0]; f1.b = data[2]; f1.c = data[4]; f1.d = data[6];
              vn1[0] = data[1]; vn1[1] = data[3]; vn1[2] = data[5]; vn1[3] = data[7];
              Vn4.push_back(f1);
            }
          } else if(num_slash == 2 * num_face_vertices && dbslash > 1) {
            temp_stream >> data[0] >> discard >> data[1] >> discard >> data[2];
            temp_stream >> data[3] >> discard >> data[4] >> discard >> data[5];
            temp_stream >> data[6] >> discard >> data[7] >> discard >> data[8];
            if(num_face_vertices == 3) {
              f1.a = data[0]; f1.b = data[3]; f1.c = data[6]; f1.d = -1;
              tc1[0] = data[1]; tc1[1] = data[4]; tc1[2] = data[7]; tc1[3] = -1;
              vn1[0] = data[2]; vn1[1] = data[5]; vn1[2] = data[8]; vn1[3] = -1;
              Vn3.push_back(f1);
            } // if
            if(num_face_vertices == 4) {
              temp_stream >> data[9] >> discard >> data[10] >> discard >> data[11];
              f1.a = data[0]; f1.b = data[3]; f1.c = data[6]; f1.d = data[9];
              tc1[0] = data[1]; tc1[1] = data[4]; tc1[2] = data[7]; tc1[3] = data[10];
              vn1[0] = data[2]; vn1[1] = data[5]; vn1[2] = data[8]; vn1[3] = data[11];
              Vn4.push_back(f1);
            } // if
          } // if-else
          if(num_face_vertices == 3) {
            F3.push_back(f1);
            ++ f3num;
          } else if(num_face_vertices == 4) {
            F4.push_back(f1);
            ++ f4num;
          } else {
            std::cout << "Warning: something is not right!" << std::endl;
          } // if-else
          break;*/

        case FACE:  // new processing
          face_vi.clear();  // vertex indices
          face_vti.clear();  // vertex texture indices
          face_vni.clear();  // vertex normal indices
          iter1 = tokens1.begin();
          ++ iter1;  // ignore the first token
          for(; iter1 != tokens1.end(); ++ iter1) {
            token_separator_t sep2("/", " ", boost::keep_empty_tokens);
            tokenizer_t tokens2(*iter1, sep2);
            // cases: v v/vt v/vt/vn v//vn
            iter2 = tokens2.begin();
            face_vi.push_back(boost::lexical_cast<int>(*iter2));
            ++ iter2;
            if(iter2 != tokens2.end()) {
              if(!(*iter2).empty()) face_vti.push_back(boost::lexical_cast<int>(*iter2));
              ++ iter2;
            } // if
            if(iter2 != tokens2.end()) {
              face_vni.push_back(boost::lexical_cast<int>(*iter2));
              ++ iter2;
            } // if
            if(iter2 != tokens2.end()) {
              std::cerr << "error: could not understand the input line '"
                  << temp_string << "'" << std::endl;
              return false;
            } // if-else
          } // for
          if(face_vi.size() == 3) {
            face_list_3v.push_back(face_vi);
            ++ face_num_3v;
          } else if(face_vi.size() == 4) {
            face_list_4v.push_back(face_vi);
            ++ face_num_4v;
          } else {
            std::cerr << "warning: ignoring face with more than four vertices" << std::endl;
          } // if-else
          break;

        case MATERIAL_LIBRARY:
        case MATERIAL_NAME:
        case LINE:
        case SMOOTH_SHADING:
        case COMMENT:
          break;

        default:
          std::cout << "unprocessed: " << temp_string << std::endl;
          break;
      } // switch
    } // while

    input.close();

    return true;
  } // ObjectShapeReader::load_object()


  /**
   * convert data from object file into shape definition
   */
  bool ObjectShapeReader::convert_to_shape(std::vector<std::vector<int> > face_list_3v,
                          std::vector<vertex_t> vertices,
                          std::vector<real_t> &shape_def) {
    int num_triangles = face_list_3v.size();
    int count = 0;
    for(std::vector<std::vector<int> >::iterator i = face_list_3v.begin();
        i != face_list_3v.end(); ++ i) {
      real_t s_area = 0.0;
      vertex_t norm, center;
      if(!get_triangle_params(vertices[(*i)[0] - 1], vertices[(*i)[1] - 1], vertices[(*i)[2] - 1],
                s_area, norm, center)) continue;
      shape_def.push_back(s_area);
      shape_def.push_back(norm.x);
      shape_def.push_back(norm.y);
      shape_def.push_back(norm.z);
      shape_def.push_back(center.x);
      shape_def.push_back(center.y);
      shape_def.push_back(center.z);
    } // for

    return true;
  } // ObjectShapeReader::convert_to_shape()


  /**
   * given the vertices of a triangle in order, compute surface area, normal and center
   */
  bool ObjectShapeReader::get_triangle_params(vertex_t v1, vertex_t v2, vertex_t v3,
      real_t &s_area, vertex_t &normal, vertex_t &center) {
    center.x = (v1.x + v2.x + v3.x) / 3.0;
    center.y = (v1.y + v2.y + v3.y) / 3.0;
    center.z = (v1.z + v2.z + v3.z) / 3.0;

    vertex_t a, b;
    a.x = (v2.x - v1.x); a.y = (v2.y - v1.y); a.z = (v2.z - v1.z);
    b.x = (v3.x - v1.x); b.y = (v3.y - v1.y); b.z = (v3.z - v1.z);

    real_t norm_a = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    real_t norm_b = sqrt(b.x * b.x + b.y * b.y + b.z * b.z);
    real_t dot = a.x * b.x + a.y * b.y + a.z * b.z;
    vertex_t cross;
    cross.x = a.y * b.z - a.z * b.y; cross.y = a.z * b.x - a.x * b.z; cross.z = a.x * b.y - a.y * b.x;

    real_t norm_cross = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);

    if(norm_cross == 0 || norm_a == 0 || norm_b == 0) return false;

    normal.x = cross.x / norm_cross; normal.y = cross.y / norm_cross; normal.z = cross.z / norm_cross;

    real_t sintheta = sqrt(1 - (dot / (norm_a * norm_b)) * (dot / (norm_a * norm_b)));
    s_area = norm_a * norm_b * sintheta / 2;
    return true;
  } // ObjectShapeReader::get_triangle_params()


  /**
   * hashing of tokens. a token is the first character of a line.
   */
  token_t ObjectShapeReader::token_hash(std::string const &str) {
    if(str == "#") return COMMENT;
    if(str == "v") return VERTEX;
    if(str == "vt") return TEXTURE;
    if(str == "g") return SUB_MESH;
    if(str == "mtllib") return MATERIAL_LIBRARY;
    if(str == "usemtl") return MATERIAL_NAME;
    if(str == "l") return LINE;
    if(str == "s") return SMOOTH_SHADING;
    if(str == "vn") return NORMAL;
    if(str == "f") return FACE;
    return UNKNOWN;
  } // ObjectShapeReader::token_hash()


  /**
   * find positions of all occurances of character 'c' in string 'str'
   */
  void ObjectShapeReader::findall(std::string str, char c, std::vector<int> &pos_list) {
    int pos = 0;
    pos = str.find(c, pos);
    while(pos != (int)std::string::npos) {
      pos_list.push_back(pos);
      pos = str.find(c, pos + 1);
    } // while
  } // ObjectShapeReader::findall()


  void ObjectShapeReader::display_vertices(std::vector<vertex_t> &vertices) {
    for(std::vector<vertex_t>::iterator i = vertices.begin(); i != vertices.end(); ++ i) {
      std::cout << (*i).x << "\t" << (*i).y << "\t" << (*i).z << std::endl;
    } // for
  } // ObjectShapeReader::display_vertices()


  void ObjectShapeReader::display_poly_index(std::vector<poly_index_t> &Vn3) {
    for(std::vector<poly_index_t>::iterator i = Vn3.begin(); i != Vn3.end(); ++ i) {
      std::cout << (*i).a << "\t" << (*i).b << "\t" << (*i).c << "\t" << (*i).d << std::endl;
    } // for
  } // ObjectShapeReader::display_poly_index()

} // namespace hig
