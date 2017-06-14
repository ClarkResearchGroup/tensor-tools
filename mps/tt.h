#ifndef My_TENSORTRAIN_CLASS_H
#define My_TENSORTRAIN_CLASS_H

#include "../linalg/lapack_wrapper.h"
#include "../util/ezh5.h"

/*
This class is the base class of MPS and MPO.
(1) It handles the allocation of the data structure (Eigen3 dense matrix).
(2) It defines basic arithmetic operations (+,-,*,/,=) on the tensor trains.
(3) It provides the file input/output utility of the tensor trains.

As a side note, for the MatrixProduct (or TensorTrain) representation,
the underlying data structure is more Matrix-like than rank-3-tensor-like.
Using Eigen3's matrix as storage type also makes it to manipulate the linear
algebras.

Therefore, even there is a dtensor (dense tensor) class, it is not used here.
dtensor is more suitable when faced with a more complicated network, like
the networks in DMRG's eigen value problem, or SIMPS's linear equation problem.
For those situations, the best choice is to tensorize matrices in the MPS/MPO.
*/

template <typename T>
class TensorTrain
{
public:
  unsigned TTID;
  int length;
  int index_size; // For mps, index_size=phy_size; For mpo, index_size=phy_size^2
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>** M;
  std::vector<int> bond_dims;
  bool tensors_allocated;

  TensorTrain();
  TensorTrain(const TensorTrain<T>& other);
  TensorTrain(TensorTrain<T>&& other);
  ~TensorTrain();

  // set shapes
  void setLength(int L);
  void setIndexSize(int s);
  void setBondDim(int bd);
  void setBondDims(std::vector<int> bds);

  // tensor container management
  void allocateTensors();
  void freeTensors();

  // set contents
  void setZero();
  void setZero(int bd);
  void setRandom();
  void setRandom(int bd);

  // basic arithmetic operations
  TensorTrain& operator = (const TensorTrain<T>& other);
  TensorTrain& operator = (TensorTrain<T>&& other);
  TensorTrain& operator +=(const TensorTrain<T>& A);
  TensorTrain& operator -=(const TensorTrain<T>& A);
  TensorTrain  operator + (const TensorTrain<T>& A);
  TensorTrain  operator - (const TensorTrain<T>& A);
  TensorTrain& operator *=(const T& c);
  TensorTrain& operator /=(const T& c);
  TensorTrain  operator * (const T& c);
  TensorTrain  operator / (const T& c);

  // display
  void print(int level=0);

  // file stream
  void save(std::string fn);
  void load(std::string fn);

  // canonicalization
  void rc();
  void lc();
  void normalize();
  double norm();
  void moveLeft(int site, bool print_EE=false);
  void moveRight(int site, bool print_EE=false);

  // svd
  void svd(double tol, bool dry_run=false);

};


#endif
