#ifndef BIG_DENSE_TENSOR_CLASS_HEADER
#define BIG_DENSE_TENSOR_CLASS_HEADER

#include "dtensor.h"
#include "dtensor_view.h"

/*
Big dtensor class:
  1. Purpose:
  Assume we need to multiply a tensor structure (composed of several tensors,
  say, A*B*C*D) by some other tensors (E_i, i=1,2,3,...) for a lot of times,
  like in any iterative solver, where A*B*C*D are identified as a "matrix".
  Often it happens that contracting A*B*C*D into a single tensor is very costly,
  probably due to its large tensor rank. And often we would be better off
  doing A*B*E_i*C*D, because contracting E_i in the middle makes the tensor rank
  small throughout the process.

  So, in some situation, we do not want to contract A,B,C,D first. But in order
  to make the the solver code more coherent, it is helpful to hold them in a
  single object that we can identify as a "matrix". Then we can use it to perform
  "matrix-vector" multiplication in the old style.
  big_dtensor class is designed exactly to fullfill this purpose.

  "matrix", when contracted, must produce indices in pairs of prime levels
  of 0 and 1.

  2. Usage:
  The big_tensor class provides a uniform interface for performing the following
  operations:

       |--   --|
       |   |   |
  y =  L--mid--R
       |   |   |
       |-- x --|

  where effectively, our "matrix" A is of the structure

       |--   --|
       |   |   |
  A =  L--mid--R
       |   |   |
       |--   --|

  Special case:
  L and R are initialized to be nullptrs. When not set, they are effectively
  1 (scalar).

*/

template <typename T>
class big_dtensor{
public:
  big_dtensor(){
    L=nullptr;
    R=nullptr;
  }
  ~big_dtensor(){}

  void setLeft(dtensor<T>* Left) {L=Left;}
  void setRight(dtensor<T>* Right) {R=Right;}
  void addMid(dtensor<T>* m) {mid.push_back(m);}

  dtensor<T> product(dtensor<T>& v);
  dtensor<T> product(dtensor_view<T>& v);

  dtensor<T> diagonal(); // does not include bias term

  T expec(dtensor<T>& v);
  T expec(dtensor_view<T>& v);

private:
  // Left/Right environment
  dtensor<T>* L;
  dtensor<T>* R;
  // Ops between L and R
  vector< dtensor<T>* > mid;
};

#endif
