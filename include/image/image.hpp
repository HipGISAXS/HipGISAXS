/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: image.hpp
 *  Created: Jun 18, 2012
 *  Modified: Wed 08 Oct 2014 12:13:02 PM PDT
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

#ifndef _IMAGE_HPP_
#define _IMAGE_HPP_

#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/numeric/affine.hpp>

#include <common/globals.hpp>
#include <image/colormap.hpp>
#include <common/typedefs.hpp>

namespace hig {

  class Image {
    //template <typename ChannelValue, typename Layout> struct pixel;
    //typedef pixel<bits8, rgb_layout_t> rgb8_pixel_t;

    private:
      unsigned int nx_;        /* x dimension - used in 3D image construction */
      unsigned int ny_;        /* y dimension */
      unsigned int nz_;        /* z dimension */
      boost::gil::rgb8_pixel_t* image_buffer_;  /* this will hold the final rgb values */
      ColorMap8 color_map_;      /* defines mapping to colors in the defined palette */
      ColorMap new_color_map_;    /* new color mapping */

      bool scale_image(unsigned int, unsigned int, unsigned int, unsigned int,
              real_t*, real_t*&);
      bool resample_pixels(unsigned int, unsigned int, real_t*, unsigned int, unsigned int,
              real_t*&, const boost::gil::matrix3x2<real_t>&);
      bool convert_to_rgb_pixels(unsigned int, unsigned int, real_t*);
      bool convert_to_rgb_palette(unsigned int, unsigned int, real_t*);
      bool slice(Image* &img, unsigned int xval = 0);  /* obtain a slice at given x in case of 3D data */

      bool translate_pixels_to_positive(unsigned int nx, unsigned int ny, real_t * & data);
      bool normalize_pixels(unsigned int nx, unsigned int ny, real_t * & data);
      vector2_t minmax(unsigned int n, real_t* data);

      // temporary workaround ...
      void remove_nans_infs(unsigned int nx, unsigned int ny, real_t* data);

    public:
      Image(unsigned int ny, unsigned int nz);          /* initialize a 2D image object */
      Image(unsigned int ny, unsigned int nz, char* palette);
      Image(unsigned int ny, unsigned int nz, std::string palette);
      Image(unsigned int nx, unsigned int ny, unsigned int nz);  /* initialize a 3D image object */
      Image(unsigned int nx, unsigned int ny, unsigned int nz, char* palette);
      Image(unsigned int nx, unsigned int ny, unsigned int nz, std::string palette);
      Image(unsigned int nx, unsigned int ny, unsigned int nz,
          unsigned int r, unsigned int g, unsigned int b);
      ~Image();

      bool construct_image(const real_t* data, int slice);
      bool construct_image(real_t* data);
      bool construct_palette(real_t* data);
      bool save(std::string filename);      /* save the current image buffer */
      bool save(char* filename);          /* if buffer has 3D data, it will
                               save all slices */
      bool save(std::string filename, int xval);  /* save slice xval */
      bool save(char* filename, int xval);
      bool save(std::string filename, int xbegin, int xend);  /* save slices from xbegin to xend */
      bool save(char* filename, int xbegin, int xend);

  }; // class Image

} // namespace hig

#endif /* _IMAGE_HPP_ */
