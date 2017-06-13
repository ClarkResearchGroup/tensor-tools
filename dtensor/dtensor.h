#ifndef DENSE_TENSOR_CLASS_HEADER
#define DENSE_TENSOR_CLASS_HEADER

#include "../util/types_and_headers.h"
#include "tensor_index.h"

//---------------------------------------------------------------------------
// Type definitions
typedef initializer_list<unsigned>    int_list;
typedef initializer_list<string>      str_list;
typedef initializer_list<index_type>  typ_list;

typedef vector<unsigned>    int_vec;
typedef vector<string>      str_vec;
typedef vector<index_type>  typ_vec;
typedef vector<len_type>    len_vec;
typedef vector<label_type>  lab_vec;
//---------------------------------------------------------------------------


/*
Design of the templated dtensor class:

(1) dtensor's underlying storage and algorithms are provided
    by TBLIS.
    TBLIS is a library and framework for performing tensor operations,
    especially tensor contraction, using efficient native algorithms.
    https://github.com/devinamatthews/tblis

(2) Tensor ordering scheme
    dtensor's elements are basically stored in an (hidden) array contiguous in memory.
    Assume that we label the elements of a rank-4 dtensor as T(i1,i2,i3,i4),
    where i1 \in [0,s1), i2 \in [0,s2), i3 \in [0,s3), i4 \in [0,s4).
    When first index i1 increases by 1, the tensor element moves to the next 1
    in the underlying array, so i1 has stride 1.
    Similarly, i2 has stride s1, i3 has stride s1*s2, i4 has stride s1*s2*s3.

(3) Tensor multiplication
    operator * assumes Einstein summation rule, with one caveat.
    Normally, we would only sum over an index if it is repeated 2 times.
    And we would frown at an expression where an index is repeated more than 2 times,
    treating it as a mistake.
    However, under the current implementation, summation is applied even if an index
    is repeated more than 2 times.

(3) When performing dtensor contraction, Einstein summation rule is assumed.
Repeated indices will be summed over.
*/

template <typename T>
class dtensor{
public:
  //---------------------------------------------------------------------------
  // Constructors
  dtensor();                               // default constructor
  dtensor(int_list idx_sizes);                 // given idx_sizes, random names, default to type Link
  dtensor(int_vec idx_sizes);                  // given idx_sizes, random names, default to type Link
  dtensor(int_list idx_sizes, str_list names); // given idx_sizes and names, default to type Link
  dtensor(int_vec idx_sizes, str_vec names);   // given idx_sizes and names, default to type Link
  dtensor(int_list idx_sizes, str_list names, typ_list types); // given idx_sizes, names, and types
  dtensor(int_vec idx_sizes, str_vec names, typ_vec types);    // given idx_sizes, names, and types
  dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels); // given idx_sizes, names, types, levels
  dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels);    // given idx_sizes, names, types, levels
  dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels, T* data_array); // given idx_sizes, names, types, levels, data
  dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels, T* data_array);    // given idx_sizes, names, types, levels, data
  dtensor(vector<tensor_index>& idx_vec);
  dtensor(const dtensor<T>& other);  // copy constructor
  dtensor(dtensor<T>&& other);       // move constructor
  //---------------------------------------------------------------------------
  // Destructor
  ~dtensor(){}

  //---------------------------------------------------------------------------
  // Storage
  unsigned size;        // total size (number of elements) of the tensor
  unsigned rank;        // number of indices
  vector<tensor_index> idx_set; // full set of indices
  tensor<T> _T;         // tblis::tensor<T>
  bool _initted;        // initilization flag

  //---------------------------------------------------------------------------
  // Initializer
  void setRandom();
  void setZero();

  //---------------------------------------------------------------------------
  // Overloaded operator
  dtensor& operator = (const dtensor<T>& other); // copy assignment
  dtensor& operator = (dtensor<T>&& other);      // move assignment
  dtensor operator * (dtensor<T>& A);      // repeated indices are summed over
  dtensor& operator *= (const T c);              // scaling
  dtensor& operator /= (const T c);              // scaling
  dtensor operator * (const T c);                // scaling
  dtensor operator / (const T c);                // scaling

  //---------------------------------------------------------------------------
  // Full contraction (ends in a scalar)
  T contract(dtensor<T>& A);

  //---------------------------------------------------------------------------
  // Get diagonal subtensor
  // only possible when tensor indices must com in "pairs",
  // meaning same name string but different prime level
  dtensor diagonal();
  dtensor diagonal(index_type type);

  //---------------------------------------------------------------------------
  // Prime level manipulation
  void prime(int inc=1);
  void primeLink();
  void primeSite();
  void mapPrime(unsigned from, unsigned to);
  void mapPrime(unsigned from, unsigned to, index_type type);

  //---------------------------------------------------------------------------
  // Save/Load
#ifdef USE_EZH5
  void save(std::string fn);
  void load(std::string fn);
#endif

  //---------------------------------------------------------------------------
  // Get norm
  double norm();
  void normalize();

  //---------------------------------------------------------------------------
  // Print
  void print(unsigned print_level=0);

};

#endif
