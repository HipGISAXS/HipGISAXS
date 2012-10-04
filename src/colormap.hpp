/***
  *  Project: HipGISAXS (High-Performance GISAXS)
  *
  *  File: colormap.hpp
  *  Created: Jul 02, 2012
  *  Modified: Mon 01 Oct 2012 11:10:58 AM PDT
  *
  *  Author: Abhinav Sarje <asarje@lbl.gov>
  *
  *  Developers: Slim Chourou, Abhinav Sarje, Elaine Chan, Alexander Hexemer, Xiaoye Li
  *  Affiliation: Lawrence Berkeley National Laboratory, Berkeley, CA, USA
  */

#ifndef _COLORMAP_HPP_
#define _COLORMAP_HPP_

#include <boost/array.hpp>

namespace hig {

	typedef boost::array <unsigned char, 3> color_t;
	
	class ColorMap8 {
		private:
			unsigned char* color_palette_;	// make this a vector ... ?
			int palette_num_;
			unsigned int palette_size_;
	
		public:
			ColorMap8() { color_palette_ = NULL; init("jet"); } 	// default is jet
			ColorMap8(const char* palette) { color_palette_ = NULL; init(palette); }
			ColorMap8(const std::string palette) { color_palette_ = NULL; init(palette.c_str()); }
			~ColorMap8() { if(color_palette_ != NULL) delete[] color_palette_; }

			int palette_size() const { return palette_size_; }

			color_t operator[](unsigned int index) {
				color_t colors;
				colors[0] = colors[1] = colors[2] = 0;
				if(index >= palette_size_) {
					std::cerr << "error: color map index is out of range: "
							<< index << " / " << palette_size_ << std::endl;
					return colors;
				} // if

				colors[0] = color_palette_[index * 3];
				colors[1] = color_palette_[index * 3 + 1];
				colors[2] = color_palette_[index * 3 + 2];

				return colors;
			} // operator[]() */
	
			bool color_map(unsigned int index,
							unsigned char& red, unsigned char& green, unsigned char& blue) {
				if(index >= palette_size_) {
					std::cerr << "error: color map index is out of range: " << index << std::endl;
					return false;
				} // if
	
				red = color_palette_[index * 3];
				green = color_palette_[index * 3 + 1];
				blue = color_palette_[index * 3 + 2];

				return true;
			} // color_map()

		private:
			bool init(const char* palette_name) {
				if(palette_name == "jet") palette_num_ = 0;
				else {
					std::cerr << "error: palette '" << palette_name << "' is not defined" << std::endl;
					return false;
				} // if-else
			
				switch(palette_num_) {
					case 1:
					case 2:
			
					case 0:		// the default case - jet color map
					default:
						palette_size_ = 72;		// 72 colors (RGB)
						unsigned char color_palette[] = {
											  0,   0, 127,
											  0,   0, 141,
											  0,   0, 155,
											  0,   0, 169,
											  0,   0, 183,
											  0,   0, 198,
											  0,   0, 212,
											  0,   0, 226,
											  0,   0, 240,
											  0,   0, 255,
											  0,  14, 255,
											  0,  28, 255,
											  0,  42, 255,
											  0,  56, 255,
											  0,  70, 255,
											  0,  84, 255,
											  0,  98, 255,
											  0, 112, 255,
											  0, 127, 255,
											  0, 141, 255,
											  0, 155, 255,
											  0, 169, 255,
											  0, 183, 255,
											  0, 198, 255,
											  0, 212, 255,
											  0, 226, 255,
											  0, 240, 255,
											  0, 255, 255,
											 14, 255, 240,
											 28, 255, 226,
											 42, 255, 212,
											 56, 255, 198,
											 70, 255, 183,
											 84, 255, 169,
											 98, 255, 155,
											112, 255, 141,
											127, 255, 127,
											141, 255, 112,
											155, 255,  98,
											169, 255,  84,
											183, 255,  70,
											198, 255,  56,
											212, 255,  42,
											226, 255,  28,
											240, 255,  14,
											255, 255,   0,
											255, 240,   0,
											255, 226,   0,
											255, 212,   0,
											255, 198,   0,
											255, 183,   0,
											255, 169,   0,
											255, 155,   0,
											255, 141,   0,
											255, 127,   0,
											255, 112,   0,
											255,  98,   0,
											255,  84,   0,
											255,  70,   0,
											255,  56,   0,
											255,  42,   0,
											255,  28,   0,
											255,  14,   0,
											255,   0,   0,
											240,   0,   0,
											226,   0,   0,
											212,   0,   0,
											198,   0,   0,
											183,   0,   0,
											169,   0,   0,
											155,   0,   0,
											141,   0,   0
										};
						color_palette_ = new (std::nothrow) unsigned char[palette_size_ * 3];
						for(unsigned int i = 0; i < palette_size_ * 3; ++ i)
							color_palette_[i] = color_palette[i];
				} // switch
			
				return true;
			} // init()
	
	}; // class ColorMap8

} // namespace hig
	
#endif /* _COLORMAP_HPP_ */
