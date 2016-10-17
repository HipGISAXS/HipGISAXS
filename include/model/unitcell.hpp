/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: unitcell.hpp
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

#ifndef __UNITCELL_HPP__
#define __UNITCELL_HPP__

#include <string>
#include <unordered_map>

#include <common/enums.hpp>
#include <common/globals.hpp>

namespace hig {

  class Unitcell {

    public:
      typedef std::vector <vector3_t> location_list_t;
      typedef std::unordered_map <string_t, location_list_t> element_list_t;

      typedef location_list_t::iterator location_iterator_t;
      typedef element_list_t::iterator element_iterator_t;

    private:
      std::string key_;         /* unique key for the unit cell */
      element_list_t elements_; /* list of all elements in the unit cell */

    public:
      Unitcell(void);
      Unitcell(const string_t& key);
      ~Unitcell(void);

      inline void init(void) { clear(); }

      inline void clear(void) { key_.clear(); elements_.clear(); }

      /* setters */

      inline void key(const string_t& s) { key_ = s; }
      inline void element_list(const element_list_t& el) { elements_ = el; }
      void insert_element(const string_t& shape_key, const vector3_t& location);

      /* getters */

      inline string_t key() const { return key_; }
      inline element_iterator_t element_begin() { return elements_.begin(); }
      inline element_iterator_t element_end() { return elements_.end(); }

      /* modifiers ... TODO */
      //bool update_element(const string_t&, const vector3_t&);

      void print();

  }; // class Unitcell

  typedef std::unordered_map <string_t, Unitcell> unitcell_list_t;
  typedef unitcell_list_t::iterator unitcell_iterator_t;

} // namespace hig

#endif /* __UNITCELL_HPP__ */
