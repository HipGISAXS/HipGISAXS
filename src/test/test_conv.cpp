/**
 *  Project: HipGISAXS (High-Performance GISAXS)
 *
 *  File: test_conv.cpp
 *  Created: Jul 26, 2012
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
#include <cstdlib>
#include <ctime>

#include "convolutions.hpp"

using namespace std;


bool make(real_t* &m, unsigned int mx, unsigned int my) {
  m = new (nothrow) real_t[mx * my];
  for(int i = 0; i < mx * my; ++ i) m[i] = ((real_t)rand() / RAND_MAX);
} // make()

void printm(real_t* m, unsigned int mx, unsigned int my) {
  cout << "[";
  for(int y = 0; y < my; ++ y) {
    for(int x = 0; x < mx; ++ x) {
      cout << m[mx * y + x] << " ";
    } // for
    cout << ";\n";
  } // for
  cout << "]" << endl;
} // printm


int main(int narg, char** args) {
  srand(time(NULL));

  int asize = 0, bsize = 0;
  if(narg == 3) {
    asize = atoi(args[1]);
    bsize = atoi(args[2]);
  } else {
    asize = 3;
    bsize = 2;
  } // if

  real_t *a, *b;
  make(a, asize, asize);
  make(b, bsize, bsize);
//  real_t a[] = {  0.7790, 0.1520, 0.6449,
//          0.6142, 0.5742, 0.0352,
//          0.0535, 0.2173, 0.7655};
//  real_t b[] = {  0.4603, 0.8580,
//          0.1763, 0.5797};

  real_t *c;
  unsigned int cx, cy;
  //c[] = {  0.3586, 0.7383, 0.4273, 0.5534,
  //      0.4201, 1.2697, 0.7106, 0.4041,
  //      0.1329, 0.6032, 0.8779, 0.6772,
  //      0.0094, 0.0693, 0.2609, 0.4437};
  //c[] = {  1.2697, 0.7106,
  //      0.6032, 0.8779};

  hig::Convolutions::instance().convolution_2d(hig::Convolutions::conv_valid,
                        asize, asize, a, bsize, bsize, b,
                        cx, cy, c);
  cout << "a: " << endl;
  printm(a, asize, asize);
  cout << "b: " << endl;
  printm(b, bsize, bsize);
  cout << "computed c: " << cx << "x" << cy << endl;
  printm(c, cx, cy);
//  cout << "actual c full: " << endl
//      << "0.3586  0.7383  0.4273  0.5534" << endl
//      << "0.4201  1.2697  0.7106  0.4041" << endl
//      << "0.1329  0.6032  0.8779  0.6772" << endl
//      << "0.0094  0.0693  0.2609  0.4437" << endl;

  return 0;
} // main()
 main()
