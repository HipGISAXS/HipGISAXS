/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: read_oo_output.hpp
 *  Created: Jun 09, 2012
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

#ifndef _READ_OO_INPUT_HPP_
#define _READ_OO_INPUT_HPP_

#include <istream>
#include <sstream>
#include <string>
#include <stack>

#include <common/globals.hpp>
#include <config/token_mapper.hpp>

namespace hig {

  class InputReader {
    private:
      /* singleton */
      InputReader();
      InputReader(const InputReader&);
      InputReader& operator=(const InputReader&);

      std::istringstream input_stream_;
      std::stack<TokenType> structure_stack_;  // for checking that all "{", "[", """ etc match

      Token current_token_;
      Token previous_token_;
      Token parent_token_;

      //TokenMapper mapper_;

      TokenType raw_token_lookup(char c);
      bool get_raw_token(Token& token);
      TokenType process_keyword_token(std::string& keyword);
      bool read_keyword(std::string& str);
      bool read_quoted_string(std::string& str);
      bool read_number(real_t& val);
      bool skip_white_spaces(void);
      bool skip_comments(void);

    public:
      static InputReader& instance() {
        static InputReader reader;
        return reader;
      } // instance()

      bool read_input(const char* filename);
      bool read_include_file(const string_t& filename);
      Token get_next_token(void);
      void rewind(void);

  }; // class InputReader

} // namespace hig

#endif /* _READ_OO_INPUT_HPP_ */
