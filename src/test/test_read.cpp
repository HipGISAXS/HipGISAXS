/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: test_read.cpp
 *  Created: Aug 25, 2012
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

#include "hig_input.hpp"

using namespace std;
using namespace hig;

int main(int narg, char** args) {

  if(narg != 2) {
    cout << "Please specify input filename" << endl;
    exit(1);
  } // if

  if(HiGInput::instance().construct_input_config(args[1])) {
    HiGInput::instance().print_all();
    return 0;
  } // if
  return 1;
} // main()
