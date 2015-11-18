/**
 *  Project:
 *
 *  File: cbesselj.cpp
 *  Created: Nov 16, 2015
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */


#include <iostream>
#include <fstream>

#include <numerics/numeric_utils.hpp>
#include <woo/timer/woo_boostchronotimers.hpp>


typedef std::vector<hig::complex_t> cvec_t;

bool read_input(char* fname, cvec_t& data) {
  std::ifstream input(fname);
  while(!input.eof()) {
    hig::real_t r, i;
    input >> r;
    input >> i;
    data.push_back(hig::complex_t(r, i));
  } // while
  input.close();
  return true;
} // read_input()


bool write_output(char* fname, const cvec_t& data) {
  std::ofstream output(fname);
  for(cvec_t::const_iterator i = data.begin(); i != data.end(); ++ i)
    output << (*i).real() << "\t" << (*i).imag() << std::endl;
  output.close();
  return true;
} // write_output()


bool compute_vals(const cvec_t& data, cvec_t& vals, hig::complex_t (*func)(hig::complex_t, int)) {
  vals.clear();
  vals.resize(data.size());
  #pragma omp parallel for
  for(int i = 0; i < data.size(); ++ i) {
    vals[i] = ((*func)(data[i], 1));
  } // for
  return true;
} // compute_vals()


hig::real_t error_function(const hig::complex_t a, const hig::complex_t b) {
  hig::complex_t d = a - b;
  return std::abs(d);
} // error_function()


hig::real_t compute_error(const cvec_t& d1, const cvec_t& d2) {
  hig::real_t e = 0.;
  if(d1.size() != d2.size()) {
    std::cerr << "error: invalid comparison due to unequal data lengths" << std::endl;
    return -1.;
  } // if
  #pragma omp parallel for reduction(+:e)
  for(int i = 0; i < d1.size(); ++ i) {
    e += error_function(d1[i], d2[i]);
  } // for
  return e;
} // compute_error()


int main(int narg, char** args) {

  if(narg != 3) {
    std::cout << "usage: cbesselj <inputfile> <write_output>" << std::endl;
    return 0;
  } // if

  cvec_t data, vals_ref, vals_new;
  woo::BoostChronoTimer timer;

  // read input data
  std::cout << "** reading input data ..." << std::flush;
  timer.start();
  read_input(args[1], data);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ data size = " << data.size() << std::endl;

  // compute cbesselj using hipgisaxs implemented version
  std::cout << "** computing values ..." << std::flush;
  timer.start();
  compute_vals(data, vals_new, &hig::cbesselj);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  // compute cbesselj using GNU math (reference)
  std::cout << "** computing reference ..." << std::flush;
  timer.start();
  compute_vals(data, vals_ref, &hig::cbessj);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  // compute error
  std::cout << "** computing error ..." << std::flush;
  timer.start();
  hig::real_t e = compute_error(vals_new, vals_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e << " [average = " << e / data.size() << "]" << std::endl;

  if(atoi(args[2])) {
    // save the computed values
    std::cout << "** saving values ... " << std::flush;
    timer.start();
    write_output("ref.vals", vals_ref);
    write_output("new.vals", vals_new);
    timer.stop();
    std::cout << "done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  } // if

  return 0;
} // main()
