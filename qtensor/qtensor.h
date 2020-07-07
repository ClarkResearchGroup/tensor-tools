/*
 * Copyright 2020 Ryan Levy, Xiongjie Yu, and Bryan K. Clark
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef QUANTUM_NUMBERED_TENSOR_CLASS_HEADER
#define QUANTUM_NUMBERED_TENSOR_CLASS_HEADER

#include "../util/types_and_headers.h"
#include "../util/ezh5.h"
#include "../linalg/lapack_wrapper.h"
#include "qtensor_index.h"
#include "qtensor_index_op.h"
#include <ctf.hpp>

/*
Design of the templated qtensor class:

(1) qtensor's underlying storage and algorithms are provided
    by CTF.
    CTF is a library and framework for performing tensor operations,
    especially tensor contraction, using efficient native algorithms.
    https://github.com/cyclops-community/ctf/

(2) When performing qtensor contraction, Einstein summation rule is assumed.
    Repeated indices will be summed over.
*/


string indToStr(vector<qtensor_index> &indices,unordered_map<string,char> &charMap);
string indicesToChar(vector<qtensor_index> &indices, unordered_map<string,char> &charMap);
string indToStrNP(vector<qtensor_index> &indices,unordered_map<string,char> &charMap);
string indicesToCharNP(vector<qtensor_index> &indices, unordered_map<string,char> &charMap);

template <typename T>
class qtensor{
public:
  //---------------------------------------------------------------------------
  // Constructors -- sets up the qtensor_index
  qtensor();                               // default constructor
  qtensor(arr_list arrows);                 // given idx_sizes, random names, default to type Link
  qtensor(arr_vec& arrows);                  // given idx_sizes, random names, default to type Link
  qtensor(arr_list arrows, str_list names); // given idx_sizes and names, default to type Link
  qtensor(arr_vec& arrows, str_vec& names);   // given idx_sizes and names, default to type Link
  qtensor(arr_list arrows, str_list names, typ_list types); // given idx_sizes, names, and types
  qtensor(arr_vec& arrows, str_vec& names, typ_vec& types);    // given idx_sizes, names, and types
  qtensor(arr_list arrows, str_list names, typ_list types, uint_list levels); // given idx_sizes, names, types, levels
  qtensor(arr_vec& arrows, str_vec& names, typ_vec& types, uint_vec& levels);    // given idx_sizes, names, types, levels
  qtensor(vector<qtensor_index>& idx_vec);
  qtensor(vector<qtensor_index>&& idx_vec);
  qtensor(initializer_list<qtensor_index> idx_list);

  qtensor(const qtensor<T>& other);       // copy constructor
  qtensor(qtensor<T>&& other);            // move constructor
  //---------------------------------------------------------------------------
  // Destructor
  ~qtensor(){}
  //---------------------------------------------------------------------------
  // set up quantum number information {(qn, qdim)} for a selected qtensor_index
  void addQNtoIndex(unsigned idx, quantum_number qn);
  //---------------------------------------------------------------------------
  // Initialize legal tensor block
  void initBlock();
  void clearBlock();

  //---------------------------------------------------------------------------
  // Reset
  void reset(vector<qtensor_index>& qidx_vec);

  //---------------------------------------------------------------------------
  // Storage
  unsigned rank;
  vector<qtensor_index> idx_set;
  vector<vector<T>> block; //TODO: remove!
  vector<CTF::Tensor<T> > _block;
  vector< qn_vec >  block_index_qn;
  vector< int_vec > block_index_qd;
  vector< uint_vec > block_index_qi;
  unordered_map< string, unsigned > block_id_by_qn_str;
  bool _initted;
  //---------------------------------------------------------------------------
  // Index help
  string _indices;
  string getIndices();
  string getIndices(unordered_map<string,char> &charMap);
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
  qtensor<T>& operator = (const qtensor<T>& other);      // copy assignment
  qtensor<T>& operator = (qtensor<T>&& other);           // move assignment
  qtensor<T>  operator * (qtensor<T>& A);                // repeated indices are summed over
  qtensor<T>& operator += (qtensor<T>& A);               //
  qtensor<T>& operator -= (qtensor<T>& A);               //
  qtensor<T>  operator + (qtensor<T>& A);                //
  qtensor<T>  operator - (qtensor<T>& A);                //
  qtensor<T>& operator *= (const T c);                   // scaling
  qtensor<T>& operator /= (const T c);                   // scaling
  qtensor<T>  operator * (const T c);                    // scaling
  qtensor<T>  operator / (const T c);                    // scaling

  //---------------------------------------------------------------------------
  // Full contraction (ends in a scalar)
  T contract(qtensor<T>& A);

  //---------------------------------------------------------------------------
  // Get diagonal subtensor
  // only possible when tensor indices come in "pairs",
  // meaning same name string but different prime level
  qtensor<T> diagonal();

  //---------------------------------------------------------------------------
  // special arithmetic operations with another tensor in the same format/pattern
  void add(qtensor<T>& A, T c = T(1));
  T inner_product(qtensor<T>& A);

  //---------------------------------------------------------------------------
  // Prime level manipulation
  void prime(int inc=1);
  void primeLink(int inc=1);
  void primeSite(int inc=1);
  void mapPrime(unsigned from, unsigned to);
  void mapPrime(unsigned from, unsigned to, index_type type);
  void dag();
  void conj();

  //---------------------------------------------------------------------------
  // Save/Load
  void save(string fn);
  void load(string fn);
  void save(ezh5::Node& fW);
  void load(ezh5::Node& fR);

  //---------------------------------------------------------------------------
  // Get norm
  double norm();
  double normalize();

  //---------------------------------------------------------------------------
  // Print
  void print(unsigned print_level=0);

};

template<typename T>
double calcEntropy(qtensor<T>& S, double cutoff=1e-24){
  assert(S.rank==1); assert(S._initted==true);
  CTF::Scalar<double> vNEE = 0.0;
  for(auto block : S._block)
    vNEE[""] += CTF::Function<T,double>([cutoff](T sg){ if(real(sg)>cutoff) return -norm(sg)*std::log(norm(sg)); else return 0.0; })(block["i"]);
  return vNEE;
}

#endif
