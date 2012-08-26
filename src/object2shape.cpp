/**
 * $Id: object2shape.cpp 38 2012-08-09 23:01:20Z asarje $
 */

#include <iomanip>
#include "object2shape.hpp"


o2s_converter::o2s_converter(char* filename, char* outfilename, MPI_Comm comm, bool hdf5) {
	filename_ = NULL; outfilename_ = NULL; shape_def_ = NULL;

	filename_ = new std::string(filename);
	outfilename_ = new std::string(outfilename);
	comm_ = comm;

	std::vector<vertex_t> vertices;
	std::vector<poly_index_t> Vn3, Vn4, F3, F4;

	load_object(filename, vertices, Vn3, Vn4, F3, F4);
	//std::cout << "Vertices: " << std::endl;	display_vertices(vertices);
	//std::cout << "F3: " << std::endl;		display_poly_index(F3);
	//std::cout << "F4: " << std::endl;		display_poly_index(F4);
	convert(outfilename, F3, vertices, hdf5);
} // o2s_converter()


void o2s_converter::load_object(char* filename,
		std::vector<vertex_t> &vertices, std::vector<poly_index_t> &Vn3, std::vector<poly_index_t> &Vn4,
		std::vector<poly_index_t> &F3, std::vector<poly_index_t> &F4) {
	std::ifstream input(filename);
	if(!input.is_open()) {
		std::cout << "Unable to open file " << filename << std::endl;
		return;
	} // if

	//std::cout << "Opened file " << filename << std::endl;

	int num_texture = 0, num_normal = 0;
	std::vector<int> g3num, g4num;
	int f3num = 0, f4num = 0;

	char discard;
	float_t data[12], tc1[4], vn1[4];
	poly_index_t f1;

	int q = 0;

	char* temp_line = new char[257];
	while(!input.eof()) {
		input.getline(temp_line, 256);
		std::string temp_string(temp_line);
		if(temp_string.empty()) continue;

		std::stringstream temp_stream(temp_string);
		std::string token;
		temp_stream >> token;

		std::vector<int> pos_face_vertices, pos_slash;
		int num_face_vertices = 0, num_slash = 0, dbslash = 0;

		switch(token_hash(token)) {
			case VERTEX:
				vertex_t v;
				temp_stream >> v.x;
				temp_stream >> v.y;
				temp_stream >> v.z;
				vertices.push_back(v);
				break;

			case TEXTURE:
				++ num_texture;
				break;

			case SUB_MESH:
				g3num.push_back(f3num);
				g4num.push_back(f4num);
				break;

			case NORMAL:
				++ num_normal;
				break;

			case FACE:
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

				break;

			case COMMENT:
			case MATERIAL_LIBRARY:
			case MATERIAL_NAME:
			case LINE:
			case SMOOTH_SHADING:
				break;

			default:
				std::cout << "unprocessed: " << temp_string << std::endl;
				break;
		} // switch
	} // while

	input.close();
} // load_object()


float_t* o2s_converter::convert(char* outfilename, std::vector<poly_index_t> F3,
			std::vector<vertex_t> vertices, bool hdf5) {
	int num_triangles = F3.size();
	shape_def_ = new float_t[num_triangles * 7];

	std::ofstream output;
	if(!hdf5)
		output.open(outfilename);

	int count = 0;
	for(std::vector<poly_index_t>::iterator i = F3.begin(); i != F3.end(); ++ i) {
		float_t s_area = 0.0;
		vertex_t norm, center;
		get_triangle_params(vertices[(*i).a - 1], vertices[(*i).b - 1], vertices[(*i).c - 1],
							s_area, norm, center);
		if(!hdf5) {
			output << std::setiosflags(std::ios::fixed) << std::setprecision(6)
				<< s_area << "\t" << norm.x << "\t" << norm.y << "\t" << norm.z << "\t"
				<< center.x << "\t" << center.y << "\t" << center.z << std::endl;
		} // if
			/*	<< vertices[(*i).a - 1].x << "\t" << vertices[(*i).a - 1].y << "\t"
				<<  vertices[(*i).a - 1].z << "\t"
				<< vertices[(*i).b - 1].x << "\t" << vertices[(*i).b - 1].y << "\t"
				<<  vertices[(*i).b - 1].z << "\t"
				<< vertices[(*i).c - 1].x << "\t" << vertices[(*i).c - 1].y << "\t"
				<<  vertices[(*i).c - 1].z << "\t"
				<< std::endl; */
		shape_def_[count ++] = s_area;
		shape_def_[count ++] = norm.x;
		shape_def_[count ++] = norm.y;
		shape_def_[count ++] = norm.z;
		shape_def_[count ++] = center.x;
		shape_def_[count ++] = center.y;
		shape_def_[count ++] = center.z;
	} // for

	// call the C HDF5 function
	if(hdf5)
		s2h_converter(&shape_def_, num_triangles, outfilename, comm_);
	else
		output.close();

	return shape_def_;
} // convert()


void o2s_converter::get_triangle_params(vertex_t v1, vertex_t v2, vertex_t v3,
		float_t &s_area, vertex_t &normal, vertex_t &center) {
	center.x = (v1.x + v2.x + v3.x) / 3.0;
	center.y = (v1.y + v2.y + v3.y) / 3.0;
	center.z = (v1.z + v2.z + v3.z) / 3.0;

	vertex_t a, b;
	a.x = (v2.x - v1.x); a.y = (v2.y - v1.y); a.z = (v2.z - v1.z);
	b.x = (v3.x - v1.x); b.y = (v3.y - v1.y); b.z = (v3.z - v1.z);

	float_t norm_a = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
	float_t norm_b = sqrt(b.x * b.x + b.y * b.y + b.z * b.z);
	float_t dot = a.x * b.x + a.y * b.y + a.z * b.z;
	vertex_t cross;
	cross.x = a.y * b.z - a.z * b.y; cross.y = a.z * b.x - a.x * b.z; cross.z = a.x * b.y - a.y * b.x;

	float_t norm_cross = sqrt(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z);
	normal.x = cross.x / norm_cross; normal.y = cross.y / norm_cross; normal.z = cross.z / norm_cross;

	float_t sintheta = sqrt(1 - (dot / (norm_a * norm_b)) * (dot / (norm_a * norm_b)));
	s_area = norm_a * norm_b * sintheta / 2;
} // get_triangle_params()


token_t o2s_converter::token_hash(std::string const &str) {
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
} // token_hash


void o2s_converter::findall(std::string str, char c, std::vector<int> &pos_list) {
	int pos = 0;
	pos = str.find(c, pos);
	while(pos != std::string::npos) {
		pos_list.push_back(pos);
		pos = str.find(c, pos + 1);
	} // while
} // findall()


void o2s_converter::display_vertices(std::vector<vertex_t> &vertices) {
	for(std::vector<vertex_t>::iterator i = vertices.begin(); i != vertices.end(); ++ i) {
		std::cout << (*i).x << "\t" << (*i).y << "\t" << (*i).z << std::endl;
	} // for
} // display_vertices()


void o2s_converter::display_poly_index(std::vector<poly_index_t> &Vn3) {
	for(std::vector<poly_index_t>::iterator i = Vn3.begin(); i != Vn3.end(); ++ i) {
		std::cout << (*i).a << "\t" << (*i).b << "\t" << (*i).c << "\t" << (*i).d << std::endl;
	} // for
} // display_poly_index()
