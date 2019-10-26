#ifndef BIG_QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS_HEADER
#define BIG_QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS_HEADER

#include "qstensor.h"

/*
Big qstensor class:
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
  big_qstensor class is designed exactly to fullfill this purpose.

  "matrix", when contracted, must produce indices in pairs of prime levels
  of 0 and 1.
*/

template <typename T>
class big_qstensor{
public:
  big_qstensor(){
    L=nullptr;
    R=nullptr;
    size_=0;
  }
  ~big_qstensor(){}

  void setLeft(qstensor<T>* Left) {L=Left; }
  void setRight(qstensor<T>* Right) {R=Right;}
  void addMid(qstensor<T>* m) {  mid=std::move(mid*(*m)); }

  qstensor<T> product(qstensor<T>& v);
  qstensor<T> diagonal();
  T expec(qstensor<T>& v);
  size_t size() const;

private:
  // Left/Right environment
  qstensor<T>* L;
  qstensor<T>* R;
  qstensor<T> mid;
  // (1) One site
  //vector< qstensor<T>* > mid;
  mutable size_t size_;
};

#endif
