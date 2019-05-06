#ifndef DENSE_TENSOR_CLASS_HEADER
#define DENSE_TENSOR_CLASS_HEADER

#include "../util/types_and_headers.h"
#include "../util/ezh5.h"
#include "dtensor_index.h"
#include "dtensor_index_op.h"
//#include "dtensor_view.h"
#include <ctf.hpp>

string indToStr(vector<dtensor_index> &indices,unordered_map<string,char> &charMap);
string indicesToChar(vector<dtensor_index> &indices, unordered_map<string,char> &charMap);
unsigned indicesToSize(vector<dtensor_index> &indices);

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

template <typename T> class dtensor_view;


template <typename T>
class dtensor{
public:
  //---------------------------------------------------------------------------
  // Constructors
  dtensor();                               // default constructor
  dtensor(uint_list idx_sizes);                 // given idx_sizes, random names, default to type Link
  dtensor(uint_vec& idx_sizes);                  // given idx_sizes, random names, default to type Link
  dtensor(uint_list idx_sizes, str_list names); // given idx_sizes and names, default to type Link
  dtensor(uint_vec& idx_sizes, str_vec& names);   // given idx_sizes and names, default to type Link
  dtensor(uint_list idx_sizes, str_list names, typ_list types); // given idx_sizes, names, and types
  dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types);    // given idx_sizes, names, and types
  dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels); // given idx_sizes, names, types, levels
  dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels);    // given idx_sizes, names, types, levels
  dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, T* data_array); // given idx_sizes, names, types, levels, data
  dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, T* data_array);    // given idx_sizes, names, types, levels, data
  dtensor(vector<dtensor_index>& idx_vec);
  dtensor(initializer_list<dtensor_index> idx_list);
  dtensor(vector<dtensor_index>& idx_vec, T* data_array);
  dtensor(vector<dtensor_index>& idx_vec, CTF::Tensor<T>& data_array);
  dtensor(const dtensor<T>& other);       // copy constructor
  //dtensor(const dtensor_view<T>& other);  // copy constructor
  dtensor(dtensor<T>&& other);            // move constructor
  //---------------------------------------------------------------------------
  // Destructor
  ~dtensor(){}
  //---------------------------------------------------------------------------
  // Reset
  void reset(vector<dtensor_index>& idx_vec, bool makeZero=true);
  //---------------------------------------------------------------------------
  // Resize indices
  // (data preserved when dimension of indices lowered, filled with val when enlarged)
  void resize(uint_vec& new_sizes, T val=T());
  void resize(uint_list new_sizes, T val=T());
  void resize(vector<dtensor_index>& new_idx_set, T val=T());

  //---------------------------------------------------------------------------
  // Storage
  unsigned size;                  // total size (number of elements) of the tensor
  unsigned rank;                  // number of indices
  vector<dtensor_index> idx_set;   // full set of tensor indices (dtensor_index.h)
  //tblis::tensor<T> _T;            // tblis::tensor_view<T>, provide tensor functionality (does not own data)
  CTF::Tensor<T> __T; 
  bool _initted;                  // initilization flag

  //---------------------------------------------------------------------------
  // Initializer
  void setRandom();
  void setZero();
  void setOne();

  //---------------------------------------------------------------------------
  // Permutate
  void permute(uint_vec& perm);
  void permute(uint_list perm);

  //---------------------------------------------------------------------------
  // Overloaded operator
  dtensor<T>& operator = (const dtensor<T>& other);      // copy assignment
  //dtensor<T>& operator = (const dtensor_view<T>& other); // copy assignment
  dtensor<T>& operator = (dtensor<T>&& other);           // move assignment
  dtensor<T>  operator * (dtensor<T>& A);                // repeated indices are summed over
  //dtensor<T>  operator * (dtensor_view<T>& A);           // repeated indices are summed over
  dtensor<T>& operator += (dtensor<T>& A);               //
  //dtensor<T>& operator += (dtensor_view<T>& A);          //
  dtensor<T>& operator -= (dtensor<T>& A);               //
  //dtensor<T>& operator -= (dtensor_view<T>& A);          //
  dtensor<T>  operator + (dtensor<T>& A);                //
  //dtensor<T>  operator + (dtensor_view<T>& A);           //
  dtensor<T>  operator - (dtensor<T>& A);                //
  //dtensor<T>  operator - (dtensor_view<T>& A);           //
  dtensor<T>& operator *= (const T c);                   // scaling
  dtensor<T>& operator /= (const T c);                   // scaling
  dtensor<T>  operator * (const T c);                    // scaling
  dtensor<T>  operator / (const T c);                    // scaling

  //---------------------------------------------------------------------------
  // Full contraction (ends in a scalar)
  T contract(dtensor<T>& A);
  //T contract(dtensor_view<T>& A);

  //---------------------------------------------------------------------------
  // Get diagonal subtensor
  // only possible when tensor indices come in "pairs",
  // meaning same name string but different prime level
  //dtensor<T> diagonal();
  dtensor<T> diagonal(index_type type=All);

  //---------------------------------------------------------------------------
  // Prime level manipulation
  void prime(int inc=1);
  void primeLink(int inc=1);
  void primeSite(int inc=1);
  void mapPrime(unsigned from, unsigned to);
  void mapPrime(unsigned from, unsigned to, index_type type);

  void dag(){};
  void conj();

  //---------------------------------------------------------------------------
  // Save/Load
  void save(string fn, int64_t offset=0);
  void load(string fn);
  void save(ezh5::Node& fW, int64_t offset=0);
  void load(ezh5::Node& fR);

  //---------------------------------------------------------------------------
  // Get norm
  double norm();
  double normalize();

  //---------------------------------------------------------------------------
  // special arithmetic operations with another tensor in the same format/pattern
  void add(dtensor<T>& A, T c = T(1));
  T inner_product(dtensor<T>& A);

  //---------------------------------------------------------------------------
  // Print
  void print(unsigned print_level=0);

  //---------------------------------------------------------------------------
  // index functions for CTF
  string getIndices(); //prints string based on rank
  string getIndices(unordered_map<string,char> &charMap);

private:
  string _indices;

};

//Calculate Entanglement Entropy from a rank 1 CTF tensor

template<typename T>
double calcEntropy(dtensor<T>& S, double cutoff=1e-24){
  assert(S.rank==1); assert(S._initted==true);
  CTF::Scalar<double> vNEE = 0.0;
  vNEE[""] += CTF::Function<T,double>([cutoff](T sg){ if(real(sg)>cutoff) return -norm(sg)*std::log(norm(sg)); else return 0.0; })(S.__T["i"]);
  return vNEE;
}
#endif
