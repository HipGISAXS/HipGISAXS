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
#include <woo/timer/woo_boostchronotimers.hpp>

#include <numerics/numeric_utils.hpp>


typedef std::vector<hig::complex_t> cvec_t;
typedef std::vector<hig::avx_m256c_t> avx_m256cvec_t;

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

bool convert_to_avx(const cvec_t& data, avx_m256cvec_t& avx_data) {
  avx_data.clear();
  __attribute__((aligned(32))) hig::real_t real[hig::AVX_VEC_LEN_], imag[hig::AVX_VEC_LEN_];
  for(int i = 0; i < data.size(); i += hig::AVX_VEC_LEN_) {
    for(int j = 0; j < hig::AVX_VEC_LEN_; ++ j) {
      real[j] = data[i + j].real();
      imag[j] = data[i + j].imag();
    } // for j
    hig::avx_m256c_t temp;
    temp.real = _mm256_load_pd(real);
    temp.imag = _mm256_load_pd(imag);
    avx_data.push_back(temp);
  } // for i
  return true;
} // convert_to_avx()

bool convert_from_avx(const avx_m256cvec_t& avx_data, cvec_t& data) {
  data.clear();
  __attribute__((aligned(32))) hig::real_t real[hig::AVX_VEC_LEN_], imag[hig::AVX_VEC_LEN_];
  for(int i = 0; i < avx_data.size(); ++ i) {
    _mm256_store_pd(real, avx_data[i].real);
    _mm256_store_pd(imag, avx_data[i].imag);
    for(int j = 0; j < hig::AVX_VEC_LEN_; ++ j) {
      hig::complex_t temp(real[j], imag[j]);
      data.push_back(temp);
    } // for j
  } // for i
  return true;
} // convert_to_avx()


bool compute_vals(const cvec_t& data, cvec_t& vals, hig::complex_t (*func)(hig::complex_t, int)) {
  vals.clear();
  vals.resize(data.size());
  #pragma omp parallel for
  for(int i = 0; i < data.size(); ++ i) {
    vals[i] = ((*func)(data[i], 1));
  } // for
  return true;
} // compute_vals()


bool compute_vals_vec(const avx_m256cvec_t& data, avx_m256cvec_t& vals,
                      hig::avx_m256c_t (*func)(hig::avx_m256c_t, int)) {
  vals.clear();
  vals.resize(data.size());
  #pragma omp parallel for
  for(int i = 0; i < data.size(); ++ i) {
    vals[i] = ((*func)(data[i], 1));
  } // for
  return true;
} // compute_vals_vec()


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

  std::cout << "** scalar tests **" << std::endl;

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
  hig::real_t e1 = compute_error(vals_new, vals_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e1 << " [average = " << e1 / data.size() << "]" << std::endl;

  if(atoi(args[2])) {
    // save the computed values
    std::cout << "** saving values ... " << std::flush;
    timer.start();
    write_output("ref.vals", vals_ref);
    write_output("new.vals", vals_new);
    timer.stop();
    std::cout << "done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  } // if

  std::cout << "** vector tests **" << std::endl;

  // convert data to vectors
  avx_m256cvec_t avx_data, avx_vals_new, avx_vals_ref;
  convert_to_avx(data, avx_data);

  // using the vectorized version from hipgisaxs
  std::cout << "** computing values [vectorized] ..." << std::flush;
  timer.start();
  compute_vals_vec(avx_data, avx_vals_new, &hig::cbesselj_vec);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;

  // using vectorized GNU math (ref)
  std::cout << "** computing values [vectorized] ..." << std::flush;
  timer.start();
  compute_vals_vec(avx_data, avx_vals_new, &hig::cbessj_vec);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  
  // convert to scalars
  convert_from_avx(avx_vals_new, vals_new);
  // compute error
  std::cout << "** computing error (vector vs. scalar) ..." << std::flush;
  timer.start();
  hig::real_t e2 = compute_error(vals_new, vals_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e2 << " [average = " << e2 / data.size() << "]" << std::endl;

  convert_from_avx(avx_vals_ref, vals_ref);
  // compute error
  std::cout << "** computing error (vector vs. vector) ..." << std::flush;
  timer.start();
  hig::real_t e3 = compute_error(vals_new, vals_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e3 << " [average = " << e3 / data.size() << "]" << std::endl;

  return 0;
} // main()
