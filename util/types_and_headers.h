#ifndef BASIC_DATA_TYPES_AND_HEADER
#define BASIC_DATA_TYPES_AND_HEADER

//------------------------------------------------------------------------------
// Standard headers
//------------------------------------------------------------------------------
#include <iostream>
#include <complex>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <initializer_list>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <sys/time.h>
#include <utility>
#include <numeric>
#include "omp.h"
//------------------------------------------------------------------------------
using std::vector;
using std::pair;
using std::string;
using std::to_string;
using std::set;
using std::map;
using std::unordered_set;
using std::unordered_map;
using std::initializer_list;
using std::complex;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// External library 1 -- TBLIS -- necessary
//------------------------------------------------------------------------------
// By Devin Matthews, et al.
// https://github.com/devinamatthews/tblis
// Intro:
// TBLIS is a library and framework for performing tensor operations,
// especially tensor contraction, using efficient native algorithms.
// In other words, TBLIS does not need tensor transposition to perform tensor
// contraction.
// In comparison, a popular but not optimal tensor contraction strategy consists
// of 3 steps: transposition, d/s/zgemm, transposition.
// However, TBLIS does not support sparse tensor (at least for now).
#include "tblis.h"
typedef vector<tblis::len_type>    len_vec;
typedef vector<tblis::label_type>  lab_vec;
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// External library 2 -- Eigen -- needed, but only for tests and part of "EZH5"
//------------------------------------------------------------------------------
// http://eigen.tuxfamily.org
#include <Eigen/Dense>
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// External library 3 -- HDF5 -- optional, but required when saving/loading data
//------------------------------------------------------------------------------
// Check https://support.hdfgroup.org/HDF5/ for instructions on
// how to install HDF5
//------------------------------------------------------------------------------
// For easier user interface of HDF5 functions, a wrapper called "EZH5" is used.
// By M. Chen
// https://github.com/mileschen360/ezh5
// Notice:
// If not enabled, save/load functionality for many things will not work
#ifdef USE_EZH5
#include "hdf5.h"
#include "../util/ezh5.h"
#endif
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// External library 4 -- HPTT -- optional
//------------------------------------------------------------------------------
// By Paul Springer
// https://github.com/springer13/hptt
// Intro:
// HPTT is a high-performance C++ library for out-of-place tensor transpositions.
// Notice:
// Tensor transposition is used to reshape tensors before sending them to QR/SVD.
// If enabled, tensor transposition calls will be sent to HPTT driver.
// If not enabled, tensor transposition will be sent to a not optimal function.
#ifdef USE_HPTT
#include <hptt.h>
#endif
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Single Abelian quantum number
// first: qn;    second: dimension
//------------------------------------------------------------------------------
typedef pair<int, unsigned> quantum_number;
//------------------------------------------------------------------------------
// For sorting quantum numbers by qn
bool qnCompare(const quantum_number& qn1, const quantum_number& qn2) {
  return qn1.first < qn2.first;
}
// For sorting singular values of different quantum blocks
bool pairCompare(const pair<double, unsigned>& p1, const pair<double, unsigned>& p2) {
  return p1.first > p2.first;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Enum
//------------------------------------------------------------------------------
// Index type
//------------------------------------------------------------------------------
enum index_type {Link, Site};
//------------------------------------------------------------------------------
// Arrow type, for qtensor
//------------------------------------------------------------------------------
enum arrow_type {Inward, Outward};
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// For svd, and sweeping procedures
//------------------------------------------------------------------------------
#define MoveFromLeft 0
#define MoveFromRight 1
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Deal with complex conjugation in class template
inline double cconj(double val) {return val;}
inline std::complex<double> cconj(std::complex<double> val) {return std::conj(val);}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Thread safe random numbers
//------------------------------------------------------------------------------
double thread_safe_random_double() {
  // static thread_local std::mt19937 generator(std::random_device{}());
  static thread_local std::mt19937 generator;
  std::uniform_real_distribution<double> distribution(0, 1);
  return distribution(generator);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Deal with complex random numbers in class template
//------------------------------------------------------------------------------
inline void random_array(double* val, unsigned size){
  for (size_t i = 0; i < size; i++) {
    val[i] = 2*(thread_safe_random_double() - 0.5);
  }
}
inline void random_array(std::complex<double>* val, unsigned size){
  for (size_t i = 0; i < size; i++) {
    val[i] = std::complex<double> (2*(thread_safe_random_double() - 0.5), 2*(thread_safe_random_double() - 0.5));
  }
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// short hand
//------------------------------------------------------------------------------
typedef initializer_list<unsigned>    uint_list;
typedef initializer_list<int>         int_list;
typedef initializer_list<string>      str_list;
typedef initializer_list<arrow_type>  arr_list;
typedef initializer_list<index_type>  typ_list;
typedef vector<unsigned>              uint_vec;
typedef vector<int>                   int_vec;
typedef vector<string>                str_vec;
typedef vector<arrow_type>            arr_vec;
typedef vector<index_type>            typ_vec;
//------------------------------------------------------------------------------


#endif
