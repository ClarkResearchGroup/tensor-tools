#ifndef BASIC_DATA_TYPES_AND_HEADER
#define BASIC_DATA_TYPES_AND_HEADER

// Standard headers
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
#include <omp.h>

using std::vector;
using std::string;
using std::to_string;
using std::set;
using std::map;
using std::unordered_set;
using std::unordered_map;
using std::initializer_list;

// Eigen3
#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif
#ifdef USE_OPENBLAS
#define EIGEN_USE_BLAS
#endif
#include <Eigen/Dense>
typedef Eigen::MatrixXd Mxd;
typedef Eigen::MatrixXcd Mxc;

// EZH5 (by M. Chen)
// The code included is adpated from the original library
// from https://github.com/mileschen360/ezh5
#ifdef USE_EZH5
#include "../util/ezh5.h"
#endif

// TBLIS (by Devin Matthews, et al.)
// TBLIS is a library and framework for performing tensor operations,
// especially tensor contraction, using efficient native algorithms.
// https://github.com/devinamatthews/tblis
#include "tblis.h"
using tblis::tensor;
using tblis::label_type;
using tblis::len_type;

// Tensor index type
enum index_type {Link, Site};

#endif
