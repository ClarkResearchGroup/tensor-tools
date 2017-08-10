#ifndef DENSE_TENSOR_VIEW_CLASS_HEADER
#define DENSE_TENSOR_VIEW_CLASS_HEADER

#include "../util/types_and_headers.h"
#include "../util/ezh5.h"
#include "dtensor_index.h"
#include "dtensor_index_op.h"
#include "dtensor.h"

template <typename T> class dtensor;

/*
dtensor_view class is the same as dtensor class,
except that it does not own the data container.
*/

template <typename T>
class dtensor_view{
public:
  //---------------------------------------------------------------------------
  // Constructors
  dtensor_view(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, T* data_array); // given idx_sizes, names, types, levels, data
  dtensor_view(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, T* data_array);    // given idx_sizes, names, types, levels, data
  dtensor_view(vector<dtensor_index>& idx_vec, T* data_array);
  dtensor_view(const dtensor<T>& other);       // copy constructor
  dtensor_view(const dtensor_view<T>& other);  // copy constructor
  dtensor_view(dtensor_view<T>&& other);       // move constructor
  //---------------------------------------------------------------------------
  // Destructor
  ~dtensor_view(){}
  //---------------------------------------------------------------------------
  // Update data pointer
  void update(T* data_array);

  //---------------------------------------------------------------------------
  // Storage
  unsigned size;                  // total size (number of elements) of the tensor
  unsigned rank;                  // number of indices
  vector<dtensor_index> idx_set;   // full set of tensor indices (dtensor_index.h)
  tblis::tensor_view<T> _T;            // tblis::tensor_view<T>, provide tensor functionality (does not own data)
  bool _initted;                  // initilization flag

  //---------------------------------------------------------------------------
  // Initializer
  void setRandom();
  void setZero();

  //---------------------------------------------------------------------------
  // Permutate
  void permute(uint_vec& perm);
  void permute(uint_list perm);

  //---------------------------------------------------------------------------
  // Overloaded operator
  dtensor_view<T>& operator = (const dtensor_view<T>& other); // copy assignment
  dtensor_view<T>& operator = (dtensor_view<T>&& other);      // move assignment
  dtensor<T>       operator * (dtensor_view<T>& A);           // repeated indices are summed over
  dtensor<T>       operator * (dtensor<T>& A);                // repeated indices are summed over
  dtensor_view<T>& operator += (dtensor_view<T>& A);          //
  dtensor_view<T>& operator += (dtensor<T>& A);               //
  dtensor_view<T>& operator -= (dtensor_view<T>& A);          //
  dtensor_view<T>& operator -= (dtensor<T>& A);               //
  dtensor<T>       operator + (dtensor_view<T>& A);           //
  dtensor<T>       operator + (dtensor<T>& A);                //
  dtensor<T>       operator - (dtensor_view<T>& A);           //
  dtensor<T>       operator - (dtensor<T>& A);                //
  dtensor_view<T>& operator *= (const T c);                   // scaling
  dtensor_view<T>& operator /= (const T c);                   // scaling
  dtensor<T>       operator * (const T c);                    // scaling
  dtensor<T>       operator / (const T c);                    // scaling

  //---------------------------------------------------------------------------
  // Full contraction (ends in a scalar)
  T contract(dtensor_view<T>& A);
  T contract(dtensor<T>& A);

  //---------------------------------------------------------------------------
  // Get diagonal subtensor
  // only possible when tensor indices come in "pairs",
  // meaning same name string but different prime level
  dtensor<T> diagonal();
  dtensor<T> diagonal(index_type type);

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
  void save(std::string fn);
  // void load(std::string fn);

  //---------------------------------------------------------------------------
  // Get norm
  double norm();
  double normalize();

  //---------------------------------------------------------------------------
  // Print
  void print(unsigned print_level=0);

};

#endif
