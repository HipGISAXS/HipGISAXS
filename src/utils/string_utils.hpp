/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: string_utils.hpp
 *  Created: Jan 08, 2014
 *  Modified: Wed 08 Jan 2014 12:00:35 PM PST
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

#ifndef __STRING_UTILS_HPP__
#define __STRING_UTILS_HPP__

#include <string>

namespace hig {

	extern bool extract_first_keyword(const std::string&, std::string&, std::string&);
	extern bool extract_keyword_name_and_key(const std::string&, std::string&, std::string&);

} // namespace hig

#endif /* __STRING_UTILS_HPP__ */
