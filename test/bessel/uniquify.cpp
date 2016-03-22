/**
 *  Project:
 *
 *  File: uniquify.cpp
 *  Created: Nov 16, 2015
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <fstream>
#include <complex>
#include <unordered_set>
#include <functional>
#include <cmath>


typedef std::complex <double> complex_t;

class hash_t {
  public:
  size_t operator()(const complex_t n) const {
    std::hash<double> hash;
    return hash(std::norm(n));
  } // operator()
}; // class hash_t

class equal_t {
  bool operator()(const complex_t a, const complex_t b) {
    return (a.real() == b.real() && a.imag() == b.imag());
  } // operator()
}; // class equal_t

inline bool is_int(double n) {
  return (rint(n) == n);
} // is_int()


inline bool is_int(complex_t n) {
  return (rint(n.real()) == n.real() && rint(n.imag()) == n.imag());
} // is_int()


int main(int narg, char** args) {
  if(narg != 4) {
    std::cout << "usage: uniquify <inputfile> <outputfile> <is_complex>" << std::endl;
    return 0;
  } // if

  int iscomplex = atoi(args[3]);

  if(!iscomplex) {
    std::ifstream in(args[1]);
    double n;
    in >> n;
    bool isint = is_int(n);
    std::unordered_set <double> data;
    data.insert(n);
    while(!in.eof()) {
      in >> n;
      data.insert(n);
    } // while
    in.close();

    std::ofstream out(args[2]);
    for(std::unordered_set <double>::const_iterator i = data.begin(); i != data.end(); ++ i) {
      if(isint) out << (int) *i << std::endl;
      else out << *i << std::endl;
    } // for
    out.close();
  } else {
    // for complex data
    std::ifstream in(args[1]);
    double nr, ni;
    in >> nr; in >> ni;
    bool isint = is_int(nr);
    std::unordered_set<complex_t, hash_t> data;
    data.insert(complex_t(nr, ni));
    while(!in.eof()) {
      in >> nr; in >> ni;
      data.insert(complex_t(nr, ni));
    } // while
    in.close();

    std::ofstream out(args[2]);
    for(std::unordered_set<complex_t, hash_t>::const_iterator i = data.begin(); i != data.end(); ++ i) {
      if(isint) out << (int) (*i).real() << "\t" << (int) (*i).imag() << std::endl;
      else out << (*i).real() << "\t" << (*i).imag() << std::endl;
    } // for
    out.close();
  } // if-else

  return 0;
} // main()
