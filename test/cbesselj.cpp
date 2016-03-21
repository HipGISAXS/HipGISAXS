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
#include <malloc.h>
#include <woo/timer/woo_boostchronotimers.hpp>

#include <numerics/numeric_utils.hpp>
#include <numerics/cpu/avx_numerics.hpp>


typedef std::vector<hig::complex_t> cvec_t;
typedef hig::avx_m256c_t m256c_t;
typedef std::vector<m256c_t> avx_m256cvec_t;

template <std::size_t N>
struct MyAllocator {
  char data[N];
  void* p;
  std::size_t sz;
  MyAllocator() : p(data), sz(N) {}
  template <typename T>
  T* aligned_alloc(std::size_t a = alignof(T)) {
    if (std::align(a, sizeof(T), p, sz)) {
      T* result = reinterpret_cast<T*>(p);
      p = (char*)p + sizeof(T);
      sz -= sizeof(T);
      return result;
    } // if
    return nullptr;
  } // aligned_alloc()
}; // struct MyAllocator


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

bool convert_to_avx(const cvec_t& data, size_t sz, size_t pad, hig::avx_m256c_t* avx_data, size_t& vsz) {
  __attribute__((aligned(32))) hig::real_t real[hig::AVX_VEC_LEN_], imag[hig::AVX_VEC_LEN_];
  int k = 0;
  for(int i = 0; i < sz - hig::AVX_VEC_LEN_; i += hig::AVX_VEC_LEN_, ++ k) {
    for(int j = 0; j < hig::AVX_VEC_LEN_; ++ j) {
      real[j] = data[i + j].real();
      imag[j] = data[i + j].imag();
    } // for j
    avx_data[k].real = _mm256_load_pd(real);
    avx_data[k].imag = _mm256_load_pd(imag);
  } // for i
  // last iteration, to take care of padding if there
  for(int j = 0; j < hig::AVX_VEC_LEN_ - pad; ++ j) {
    real[j] = data[sz - hig::AVX_VEC_LEN_ + j].real();
    imag[j] = data[sz - hig::AVX_VEC_LEN_ + j].imag();
  } // for j
  for(int j = hig::AVX_VEC_LEN_ - 1; j >= pad; -- j) {
    real[j] = 0.0;
    imag[j] = 0.0;
  } // for j
  avx_data[k].real = _mm256_load_pd(real);
  avx_data[k].imag = _mm256_load_pd(imag);
  vsz = k + 1;
  return true;
} // convert_to_avx()

bool convert_to_scalar(cvec_t& data, const m256c_t* avx_data, size_t sz, size_t pad) {
  data.clear();
  __attribute__((aligned(32))) hig::real_t real[hig::AVX_VEC_LEN_], imag[hig::AVX_VEC_LEN_];
  for(int i = 0; i < sz - 1; ++ i) {
    _mm256_store_pd(real, avx_data[i].real);
    _mm256_store_pd(imag, avx_data[i].imag);
    for(int j = 0; j < hig::AVX_VEC_LEN_; ++ j) {
      hig::complex_t temp(real[j], imag[j]);
      data.push_back(temp);
    } // for j
  } // for i
  _mm256_store_pd(real, avx_data[sz - 1].real);
  _mm256_store_pd(imag, avx_data[sz - 1].imag);
  for(int j = 0; j < hig::AVX_VEC_LEN_ - pad; ++ j) {
    hig::complex_t temp(real[j], imag[j]);
    data.push_back(temp);
  } // for j
  return true;
} // convert_to_scalar()


bool compute_vals(const cvec_t& data, cvec_t& vals, hig::complex_t (*func)(hig::complex_t, int)) {
  vals.clear();
  vals.resize(data.size());
  #pragma omp parallel for
  for(int i = 0; i < data.size(); ++ i) {
    vals[i] = ((*func)(data[i], 1));
  } // for
  return true;
} // compute_vals()


bool compute_vals_vec(m256c_t* &data, size_t sz, m256c_t* vals, m256c_t (*func)(m256c_t, int)) {
  #pragma omp parallel for
  for(int i = 0; i < sz; ++ i) {
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
    std::cerr << std::endl << "error: invalid comparison due to unequal data lengths" << std::endl;
    std::cerr << "       " << d1.size() << " != " << d2.size() << std::endl;
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

  cvec_t data, vals_ref, vals_new, vals_vec_new, vals_vec_ref;
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
  size_t sz = ceil((float) data.size() / hig::AVX_VEC_LEN_);
  size_t rem = data.size() % hig::AVX_VEC_LEN_;
  size_t pad = (rem == 0) ? 0 : hig::AVX_VEC_LEN_ - rem;
  std::cout << "** converting scalars into avx vectors [padding: " << pad << "] ..." << std::endl;
  // need aligned buffers
  hig::avx_m256c_t *avx_data = (hig::avx_m256c_t*)_mm_malloc(sz * sizeof(hig::avx_m256c_t), 64);
  hig::avx_m256c_t *avx_vals_new = (hig::avx_m256c_t*)_mm_malloc(sz * sizeof(hig::avx_m256c_t), 64);
  hig::avx_m256c_t *avx_vals_ref = (hig::avx_m256c_t*)_mm_malloc(sz * sizeof(hig::avx_m256c_t), 64);
  size_t vsz;
  convert_to_avx(data, data.size(), pad, avx_data, vsz);

  std::cout << "data.size(): " << data.size() << " sz: " << sz << " pad: " << pad
            << " vsz: " << vsz << std::endl;
  assert(sz == vsz);

  // using the vectorized version from hipgisaxs
  std::cout << "** computing values [vectorized] ..." << std::flush;
  timer.start();
  compute_vals_vec(avx_data, sz, avx_vals_new, &hig::avx_cbesselj_cp);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;

  // using vectorized GNU math (ref)
  std::cout << "** computing reference [vectorized] ..." << std::flush;
  timer.start();
  compute_vals_vec(avx_data, sz, avx_vals_ref, &hig::avx_cbessj_cp);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  
  // convert to scalars
  convert_to_scalar(vals_vec_new, avx_vals_new, sz, pad);
  // compute error
  std::cout << "** computing error (vector new vs. scalar ref) ..." << std::flush;
  timer.start();
  hig::real_t e2 = compute_error(vals_vec_new, vals_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e2 << " [average = " << e2 / data.size() << "]" << std::endl;

  // convert to scalars
  convert_to_scalar(vals_vec_ref, avx_vals_ref, sz, pad);
  // compute error
  std::cout << "** computing error (vector new vs. vector ref) ..." << std::flush;
  timer.start();
  hig::real_t e3 = compute_error(vals_vec_new, vals_vec_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e3 << " [average = " << e3 / data.size() << "]" << std::endl;

  // for sanity of the reference
  std::cout << "** computing sanity (vector ref vs. scalar ref) ..." << std::flush;
  timer.start();
  hig::real_t e4 = compute_error(vals_vec_ref, vals_ref);
  timer.stop();
  std::cout << " done. [" << timer.elapsed_msec() << " ms.]" << std::endl;
  std::cout << "++ total error = " << e4 << " [average = " << e4 / data.size() << "]" << std::endl;

  return 0;
} // main()
