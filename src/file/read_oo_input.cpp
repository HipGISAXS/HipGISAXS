/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: read_oo_input.cpp
 *  Created: Jul 04, 2012
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
#include <fstream>
#include <cstring>

#include <file/read_oo_input.hpp>

namespace hig {

  InputReader::InputReader() {
    current_token_.type_ = null_token;
    previous_token_.type_ = null_token;
  } // InputReader::InputReader()


  bool InputReader::read_input(const char* filename) {
    std::ifstream f(filename);
    if(!f.is_open()) {
      std::cerr << "fatal error: could not open input file "
            << filename << ". aborting " << std::endl;
      return false;
    } // if

    // read the whole file into a string
    std::string input_string;
    std::getline(f, input_string, (char)EOF);
    f.close();

    // create and assign the string stream
    std::istringstream temp_stream(input_string.c_str());
    input_stream_.str(temp_stream.str());

    return true;
  } // InputReader::init()


  void InputReader::rewind(void) {
    //input_stream_.seekg(0);
    input_stream_.clear();
    std::istringstream temp_stream(input_stream_.str());
    //std::cout << temp_stream.str() << std::endl;
    input_stream_.str(temp_stream.str());
    while(!structure_stack_.empty()) structure_stack_.pop();
    current_token_.type_ = null_token;
    previous_token_.type_ = null_token;
    parent_token_.type_ = null_token;
  } // InputReader::rewind()


  bool InputReader::read_include_file(const string_t& filename) {
    std::ifstream f(filename.c_str());
    if(!f.is_open()) {
      std::cerr << "fatal error: could not read include file "
                << filename << ". aborting " << std::endl;
      return false;
    } // if

    std::string include_string;
    std::getline(f, include_string, (char)EOF);
    f.close();

    int pos = input_stream_.tellg();
    std::streampos spos = input_stream_.cur;
    //std::cout << "-----" << filename << "-------" << include_string.size() << "---------- " << include_string << std::endl;
    //std::cout << "---------------------- " << pos << std::endl;

    // find the total length
    input_stream_.seekg(0, input_stream_.end);
    int len = input_stream_.tellg();
    //std::cout << "================= " << len << std::endl;
    // assuming the keyword is "include"
    int ignore = filename.size() + 12;

    // read until pos
    char* temp_buf = new char[len + include_string.size() + 1];
    input_stream_.seekg(0, input_stream_.beg);
    input_stream_.read(temp_buf, pos-ignore);
    input_stream_.ignore(ignore);
    // copy include string
    strncpy(temp_buf+pos-ignore, include_string.c_str(), include_string.size());

    //std::cout << "++++++++++++++++++++++ " << temp_buf << " ++++++++++++++++++++++" << std::endl;
    input_stream_.read(temp_buf+pos-ignore+include_string.size(), len - pos);
    temp_buf[include_string.size()+len-ignore] = '\0';
    //std::cout << "++++++++++++++++++++++ " << temp_buf << " ++++++++++++++++++++++" << std::endl;

    std::istringstream temp_stream(temp_buf);
    //std::cout << temp_stream.str() << std::endl;
    input_stream_.str(temp_stream.str());

    // set position in input_stream_ to correct place
    input_stream_.seekg(pos-ignore, input_stream_.beg);
    //std::cout << input_stream_.str() << std::endl;

    return true;
  } // InputReader::read_include_file()


  TokenType InputReader::raw_token_lookup(char c) {
    Token token;

    // first check for character
    if((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == ':') {
      return character_token;
    } else {
      switch(c) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          token.type_ = digit_token;
          break;
        case '-':
          token.type_ = negative_token;
          break;
        case ' ':  /* space */
        case '\t':  /* horizontal tab */
        case '\n':  /* newline */
        case '\v':  /* vertical tab */
        case '\f':  /* feed */
        case '\r':  /* carriage return */
          token.type_ = white_space_token;
          break;
        case '{':
          token.type_ = object_begin_token;
          break;
        case '}':
          token.type_ = object_end_token;
          break;
        case '[':
          token.type_ = array_begin_token;
          break;
        case ']':
          token.type_ = array_end_token;
          break;
        case '"':
          token.type_ = string_begin_end_token;
          break;
        case '=':
          token.type_ = assignment_token;
          break;
        case ',':
          token.type_ = separator_token;
          break;
        case '#':
          token.type_ = comment_token;
          break;
        default:
          token.type_ = error_token;
          std::cerr << "Fatal Error: Illegal character '"
                << c << "' found in input." << std::endl;
          return token.type_;
      } // switch
      return token.type_;
    } // if-else
  } // InputReader::read_token()


  bool InputReader::get_raw_token(Token& token) {
    char c;
    input_stream_.get(c);

    if(!input_stream_.good()) {
      token.type_ = null_token;
      return false;
    } // if

    token.type_ = raw_token_lookup(c);
    return true;
  } // InputReader::get_raw_token()

    
  Token InputReader::get_next_token() {
    skip_white_spaces();
    
    Token token, raw_token;
    std::string qstr;
    std::string value;

    get_raw_token(raw_token);
    switch(raw_token.type_) {

      case null_token:
        token.type_ = null_token;    // this means no more tokens left
        break;

      case digit_token:
      case negative_token:
        real_t n_value;
        input_stream_.unget();
        if(!read_number(n_value)) {
          std::cerr << "fatal error: failed while reading a number" << std::endl;
          token.type_ = error_token;
        } else {
          token.type_ = number_token;
          token.dvalue_ = n_value;
        } // if-else
        break;

      case object_begin_token:
        token.type_ = object_begin_token;
        structure_stack_.push(object_begin_token);
        break;

      case object_end_token:
        token.type_ = object_end_token;
        if(structure_stack_.top() != object_begin_token) {
          std::cerr << "fatal error: mismatched object encapsulators" << std::endl;
          token.type_ = error_token;
        } else {
          structure_stack_.pop();
        } // if-else
        break;

      case array_begin_token:
        token.type_ = array_begin_token;
        structure_stack_.push(array_begin_token);
        break;

      case array_end_token:
        token.type_ = array_end_token;
        if(structure_stack_.top() != array_begin_token) {
          std::cerr << "fatal error: mismatched array encapsulators" << std::endl;
          token.type_ = error_token;
        } else {
          structure_stack_.pop();
        } // if-else
        break;

      case string_begin_end_token:  // will always be begin since
                      // end will be removed while reading
                      //the whole string earlier
        if(!read_quoted_string(qstr)) {
          std::cerr << "fatal error: premature EOF reached while reading string" << std::endl;
          token.type_ = error_token;
        } else {
          token.type_ = string_token;
          token.svalue_ = qstr;
        } // if-else
        break;

      case character_token:
        input_stream_.unget();
        read_keyword(value);
        token.type_ = process_keyword_token(value);
        if(token.type_ == error_token) {
          std::cerr << "fatal error: unknown keyword '" << value << "'" << std::endl;
        } // if
        token.svalue_ = value;
        break;

      case assignment_token:
        token.type_ = assignment_token;
        break;

      case separator_token:
        token.type_ = separator_token;
        break;

      case comment_token:
        skip_comments();
        token.type_ = comment_token;
        break;

      default:
        std::cerr << "fatal error: unknown token" << std::endl;
        token.type_ = error_token;
    } // switch

    return token;
  } // InputReader::get_next_token()
  

  TokenType InputReader::process_keyword_token(std::string& keyword) {
    return TokenMapper::instance().get_keyword_token(keyword);
  } // InputReader::process_keyword_token()


  bool InputReader::read_keyword(std::string& str) {
    //input_stream_ >> str;    // not sure if this works
    char curr_c;
    input_stream_.get(curr_c);
    if(!input_stream_.good()) return false;    // end of input reached
    TokenType tok = raw_token_lookup(curr_c);
    while(tok == character_token || tok == digit_token || tok == negative_token) {
      str.push_back(curr_c);
      input_stream_.get(curr_c);
      if(!input_stream_.good()) return false;    // end of input reached
      tok = raw_token_lookup(curr_c);
    } // while
    input_stream_.unget();
    return true;
  } // InputReader::read_keyword()


  bool InputReader::read_quoted_string(std::string& str) {
    //input_stream_ >> str;    // not sure if this works
    char curr_c;
    input_stream_.get(curr_c);
    if(!input_stream_.good()) return false;  // premature end of file reached
    while(curr_c != '"') {
      str.push_back(curr_c);
      input_stream_.get(curr_c);
    } // while
    // not ungetting because last read character is '"'
    return true;
  } // InputReader::read_quoted_string()


  bool InputReader::read_number(real_t& val) {
    input_stream_ >> val;
    return true;
  } // InputReader::read_number()


  bool InputReader::skip_white_spaces(void) {
    while(true) {
      char c;
      input_stream_.get(c);
      if(!input_stream_.good()) return false;  // end of input reached
      if(c != ' ' && c != '\t' && c != '\v'
          && c != '\r' && c != '\n' && c != '\f') {
        input_stream_.unget();
        break;
      } // if
    } // while
    return true;
  } // InputReader::skip_spaces()


  bool InputReader::skip_comments(void) {
    //input_stream_.ignore(1 << 31, '\n');  // what if its \r or \f or \v ...
    char c;
    input_stream_.get(c);
    while(c != '\n') input_stream_.get(c);
    return true;
  } // InputReader::skip_comments()

} // namespace hig
