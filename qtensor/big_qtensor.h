#ifndef BIG_QUANTUM_NUMBERED_TENSOR_CLASS_HEADER
#define BIG_QUANTUM_NUMBERED_TENSOR_CLASS_HEADER

#include "qtensor.h"

/*
Big qtensor class:
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
  big_qtensor class is designed exactly to fullfill this purpose.

  "matrix", when contracted, must produce indices in pairs of prime levels
  of 0 and 1.
*/

template <typename T>
class big_qtensor{
public:
  big_qtensor(){
    L=nullptr;
    R=nullptr;
  }
  ~big_qtensor(){}

  void setLeft(qtensor<T>* Left) {L=Left; A=std::move(A*(*L));}
  void setRight(qtensor<T>* Right) {R=Right; }
  void addMid(qtensor<T>* m) {mid.push_back(m); A=std::move(A*(*m));}

  qtensor<T> product(qtensor<T>& v);
  qtensor<T> diagonal();
  T expec(qtensor<T>& v);

private:
  // Left/Right environment
  qtensor<T>* L;
  qtensor<T>* R;
  qtensor<T> A;
  // (1) One site
  vector< qtensor<T>* > mid;
};

#endif
