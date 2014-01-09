/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: string_utils.cpp
 *  Created: Jan 08, 2014
 *  Modified: Wed 08 Jan 2014 05:13:06 PM PST
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

#include <iostream>

#include "string_utils.hpp"

namespace hig {

	bool extract_first_keyword(const std::string& str, std::string& keyword, std::string& rem_str) {
		std::size_t pos = str.find_first_of(":");
		if(pos == std::string::npos) {
			std::cerr << "error: ill-formed parameter variable name '" << str << "'" << std::endl;
			return false;
		} // if
		keyword = str.substr(0, pos);
		rem_str = str.substr(pos + 1);  // remainder string
		return true;
	} // extract_first_keyword()


	bool extract_keyword_name_and_key(const std::string& keyword, std::string& name, std::string& key) {
		std::size_t pos = keyword.find_first_of("[");
		if(pos == std::string::npos) {
			name = keyword;
			key = "";
		} else {
			name = keyword.substr(0, pos);
			std::string rem = keyword.substr(pos);
			if(rem.back() != ']') {
				std::cerr << "error: ending bracket missing in keyword '" << keyword << "'" << std::endl;
				return false;
			} // if
			rem.pop_back();
			if(rem.front() != '\'' || rem.back() != '\'') {
				std::cerr << "error: keyword key not enclosed in single-quotes in '"
							<< keyword << "'" << std::endl;
				return false;
			} // if
			rem.pop_back();
			key = rem.substr(1);
		} // if-else
		return true;
	} // extract_keyword_name_and_key()

} // namespace hig
