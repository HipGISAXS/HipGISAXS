/**
 * $Id: shape2hdf5.hpp 38 2012-08-09 23:01:20Z asarje $
 *
 * Project: HipGISAXS (High-Performance GISAXS)
 */

#ifndef _ORIENT_SAMPLE_HPP_
#define _ORIENT_SAMPLE_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

namespace hig {

	class OrientSample {
		public:
			OrientSample(std::string infile, std::string outfile, int axis);
			~OrientSample() { }

		private:
			bool load_shape_def();
			bool orient_shape_def(int axis);
			bool get_rotation_matrix(int axis, std::vector<int> &rot_matrix);
			bool mat_vec_mul3x1(std::vector<int> mat, std::vector<double> vec, std::vector<double> &out);

			std::string infile_;
			std::string outfile_;
			std::vector<double> shape_def_;
			unsigned int num_triangles_;
	}; // class shape2hdf5_converter

} // namespace hig

#endif // _ORIENT_SAMPLE_HPP_
