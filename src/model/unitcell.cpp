/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: unitcell.cpp
 *  Created: Feb 18, 2015
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

#include <iostream>

#include <model/unitcell.hpp>
#include <utils/string_utils.hpp>
#include <config/token_mapper.hpp>


namespace hig {

  Unitcell::Unitcell(void) {
    init();
    key_ = string_t("");
  } // Unitcell::Unitcell()

  Unitcell::Unitcell(const string_t& key) {
    init();
    key_ = key;
  } // Unitcell::Unitcell()

  Unitcell::~Unitcell(void) {
    clear();
  } // Unitcell::~Unitcell()


  void Unitcell::insert_element(const string_t& shape_key, const vector3_t& location) {
    element_iterator_t i = elements_.find(shape_key);
    if(i == elements_.end()) {
      // insert new key
      location_list_t loc; loc.push_back(location); // copy location
      elements_[shape_key] = loc;
    } else {
      (*i).second.push_back(location);
    } // if-else
  } // Unitcell::insert_element()


  /* printing */

  void Unitcell::print() {
    std::cout << " key_ = " << key_ << std::endl;
    for(element_iterator_t e = elements_.begin(); e != elements_.end(); ++ e) {
      std::cout << "  shape:key = " << (*e).first << std::endl;
      for(location_iterator_t l = (*e).second.begin(); l != (*e).second.end(); ++ l) {
        std::cout << "    [ " << (*l)[0] << " " << (*l)[1] << " " << (*l)[2] << " ]" << std::endl;
      } // for
    } // for
    std::cout << std::endl;
  } // Unitcell::print()

} // namespace hig
