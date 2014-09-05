/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: edf_reader.hpp
 *  Created: Aug 25, 2013
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

			bool get_data(float*&, unsigned int&, unsigned int&);

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
			std::map <std::string, std::string> header_;	/* header key -> value map */
			std::vector <float_t> data_;					/* the data as float_ts */
			unsigned int rows_;
			unsigned int cols_;

	}; // class EDFReader

} // namespace hig

#endif // __EDF_READER_HPP__
