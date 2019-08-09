#ifndef QUANTUM_NUMBERED_SPARSE_TENSORTRAIN_CLASS_H
#define QUANTUM_NUMBERED_SPARSE_TENSORTRAIN_CLASS_H

#include "../util/ezh5.h"
#include "../linalg/lapack_wrapper.h"
#include "../qstensor/qstensor_all.h"
#include "qtt.h"
#include "../models/sites/sites.h"


/*
This class is the base class of MPS and MPO.
*/

// N is the number of physical indices per site (1:MPS, 2:MPO)
template <typename T, unsigned N>
class qsTensorTrain
{
public:
  //---------------------------------------------------------------------------
  // ID used for naming the the tensors
  unsigned _id; // assigned randomly when an object is created

  //---------------------------------------------------------------------------
  // Basic properties
  unsigned length;         // length of the tensor train
  unsigned phy_dim;        // dimension of each physical bond
  uint_vec bond_dims;       // dimensions of virtual bonds
  int center;              // canonicalization center
  int totalQ;              // total abelian quantum number
  vector<int> phy_qn;      // quantum numbers of physical bonds

  //---------------------------------------------------------------------------
  // Data containers
  vector< qstensor<T> > A;

  //---------------------------------------------------------------------------
  // Flag -- whether the tensors are allocated
  bool tensors_allocated;

  //---------------------------------------------------------------------------
  // Constructors
  qsTensorTrain();
  qsTensorTrain(abstract_sites* s, int Q = 0);
  qsTensorTrain(abstract_sites* s, str_vec product_string);
  qsTensorTrain(unsigned L, unsigned pD, vector<int>& phyQN, int Q);
  qsTensorTrain(unsigned L, unsigned pD, vector<int>& phyQN, uint_vec product_state);
  qsTensorTrain(const qsTensorTrain<T, N>& other);
  qsTensorTrain(qsTensorTrain<T, N>&& other);
  qsTensorTrain(const qTensorTrain<T, N>& other);
  qsTensorTrain(qTensorTrain<T, N>&& other);
  // Destructor
  ~qsTensorTrain(){};

  //---------------------------------------------------------------------------
  // Set shapes
  void setLength(int L);
  void setPhysicalDim(int s);

  //---------------------------------------------------------------------------
  // Tensor data management
  void allocateTensors(unsigned* product_state=nullptr);
  void freeTensors();

  //---------------------------------------------------------------------------
  // Set tensors values
  void setZero();
  void setRandom();

  //---------------------------------------------------------------------------
  // Basic arithmetic operations
  qsTensorTrain<T, N>& operator = (const qsTensorTrain<T, N>& other);
  qsTensorTrain<T, N>& operator = (qsTensorTrain<T, N>&& other);
  qsTensorTrain<T, N>& operator = (const qTensorTrain<T, N>& other); //convert
  qsTensorTrain<T, N>& operator = (qTensorTrain<T, N>&& other);      //convert
  qsTensorTrain<T, N>& operator *=(const T c);
  qsTensorTrain<T, N>& operator /=(const T c);
  qsTensorTrain<T, N>  operator * (const T c) const;
  qsTensorTrain<T, N>  operator / (const T c) const;

  //---------------------------------------------------------------------------
  // Print qTensorTrain information
  void print(int level=0);

  //---------------------------------------------------------------------------
  // txt or HDF5 storage
  void save(std::string fn, std::string wfn="");
  void load(std::string fn);
  void save(ezh5::Node& fW);
  void load(ezh5::Node& fR);

  //---------------------------------------------------------------------------
  // Canonicalization
  void rc();
  void lc();

  //---------------------------------------------------------------------------
  // Adjust canonical center
  double position(int site);

  //---------------------------------------------------------------------------
  // Norm
  void normalize();
  double norm();
};


// Definition of MPS and MPO through qTensorTrain
template <typename T>
using qsMPS =  qsTensorTrain<T, 1>;

template <typename T>
using qsMPO =  qsTensorTrain<T, 2>;


#endif
