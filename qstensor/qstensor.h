#ifndef QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS_HEADER
#define QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS_HEADER

#include "../util/types_and_headers.h"
#include "../util/ezh5.h"
#include "../linalg/lapack_wrapper.h"
#include "../qtensor/qtensor_index.h"
#include "../qtensor/qtensor_index_op.h"
#include "../qtensor/qtensor.h"
#include <ctf.hpp>

/*
Design of the templated qstensor class:

(1) qstensor's underlying storage and algorithms are provided
    by TBLIS.
    TBLIS is a library and framework for performing tensor operations,
    especially tensor contraction, using efficient native algorithms.
    https://github.com/devinamatthews/tblis

(2) When performing qstensor contraction, Einstein summation rule is assumed.
    Repeated indices will be summed over.
*/



template <typename T>
class qstensor{
public:
  //---------------------------------------------------------------------------
  // Constructors -- sets up the qtensor_index
  qstensor();                               // default constructor
  qstensor(arr_list arrows);                 // given idx_sizes, random names, default to type Link
  qstensor(arr_vec& arrows);                  // given idx_sizes, random names, default to type Link
  qstensor(arr_list arrows, str_list names); // given idx_sizes and names, default to type Link
  qstensor(arr_vec& arrows, str_vec& names);   // given idx_sizes and names, default to type Link
  qstensor(arr_list arrows, str_list names, typ_list types); // given idx_sizes, names, and types
  qstensor(arr_vec& arrows, str_vec& names, typ_vec& types);    // given idx_sizes, names, and types
  qstensor(arr_list arrows, str_list names, typ_list types, uint_list levels); // given idx_sizes, names, types, levels
  qstensor(arr_vec& arrows, str_vec& names, typ_vec& types, uint_vec& levels);    // given idx_sizes, names, types, levels
  qstensor(vector<qtensor_index>& idx_vec);
  qstensor(vector<qtensor_index>&& idx_vec);
  qstensor(initializer_list<qtensor_index> idx_list);

  qstensor(const qstensor<T>& other);       // copy constructor
  qstensor(qstensor<T>&& other);            // move constructor
  
  qstensor(const qtensor<T>& other);       // copy convert constructor
  qstensor(qtensor<T>&& other);            // move convert constructor
  //---------------------------------------------------------------------------
  // Destructor
  ~qstensor(){}
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
  CTF::Tensor<T> _T;
  vector< int_vec >  block_index_qn;
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
  qstensor<T>& operator = (const qstensor<T>& other);      // copy assignment
  qstensor<T>& operator = (qstensor<T>&& other);           // move assignment
  qstensor<T>& operator = (const qtensor<T>& other);      // convert copy assignment
  qstensor<T>& operator = (qtensor<T>&& other);           // convert move assignment
  qstensor<T>  operator * (qstensor<T>& A);                // repeated indices are summed over
  qstensor<T>& operator += (qstensor<T>& A);               //
  qstensor<T>& operator -= (qstensor<T>& A);               //
  qstensor<T>  operator + (qstensor<T>& A);                //
  qstensor<T>  operator - (qstensor<T>& A);                //
  qstensor<T>& operator *= (const T c);                   // scaling
  qstensor<T>& operator /= (const T c);                   // scaling
  qstensor<T>  operator * (const T c);                    // scaling
  qstensor<T>  operator / (const T c);                    // scaling

  //---------------------------------------------------------------------------
  // Full contraction (ends in a scalar)
  T contract(qstensor<T>& A);

  //---------------------------------------------------------------------------
  // Get diagonal subtensor
  // only possible when tensor indices come in "pairs",
  // meaning same name string but different prime level
  qstensor<T> diagonal();

  //---------------------------------------------------------------------------
  // special arithmetic operations with another tensor in the same format/pattern
  void add(qstensor<T>& A, T c = T(1));
  T inner_product(qstensor<T>& A);

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
double calcEntropy(qstensor<T>& S, double cutoff=1e-24){
  assert(S.rank==1); assert(S._initted==true);
  CTF::Scalar<double> vNEE = 0.0;
  for(auto block : S._block)
    vNEE[""] += CTF::Function<T,double>([cutoff](T sg){ if(real(sg)>cutoff) return -norm(sg)*std::log(norm(sg)); else return 0.0; })(block["i"]);
  return vNEE;
}

#endif
