#ifndef DENSE_TENSOR_BASED_CG_SOLVER
#define DENSE_TENSOR_BASED_CG_SOLVER

#include "tensor_cg.h"

template <typename T>
void elemWiseDivide(dtensor<T>& A, dtensor<T>& B){
  assert(A._initted && B._initted);
  assert(A.size == B.size);
  for (unsigned i = 0; i < A.size; i++) {
    if( B._T.data()[i]!=T(0) )
      A._T.data()[i] = A._T.data()[i]/B._T.data()[i];
  }
}
template void elemWiseDivide(dtensor<double>& A, dtensor<double>& B);
template void elemWiseDivide(dtensor< std::complex<double> >& A, dtensor< std::complex<double> >& B);

template <typename T>
void LRWx(dtensor<T>& LEnv, dtensor<T>& REnv, dtensor<T>& W, dtensor<T>& x, dtensor<T>& y){
  dtensor<T> t1, t2;
  if(REnv.size==0){
    t1 = std::move(LEnv * x);
    y  = std::move(t1 * W);
  }else if(LEnv.size==0){
    t1 = std::move(REnv * x);
    y  = std::move(t1 * W);
  }else{
    t1 = std::move(REnv * x);
    t2 = std::move(t1 * W);
    y  = std::move(LEnv * t2);
  }
  y.prime(-1);
}

template <typename T>
int tensor_CG(dtensor<T>& LEnv, dtensor<T>& REnv, dtensor<T>& W, dtensor<T>& x, dtensor<T>& b, int max_iter, double tol){
  dtensor<T> A, D, t1, t2;
  if(REnv.size==0){
    A = std::move(LEnv * W);
  }else if(LEnv.size==0){
    A = std::move(REnv * W);
  }else{
    t1 = std::move(REnv * W);
    A  = std::move(LEnv * t1);
  }
  D = std::move(A.diagonal());

  T resid, alpha, beta, rho, rho_1, normb;
  normb = b.norm();

  dtensor<T> r, p, z, q;
  LRWx(LEnv,REnv,W,x,t1); r = b - t1;

  if (normb == 0.0)
    normb = 1;

  if ((resid = r.norm() / normb) <= tol) {
    // tol = std::abs(resid);
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++){
    z = r; // elemWiseDivide(z,D);
    rho = r.contract(z);
    if(i==1)
      p = z;
    else{
      beta = rho / rho_1;
      t1 = p * beta;  p = z + t1;
    }

    LRWx(LEnv,REnv,W,p,q);
    alpha = rho / p.contract(q);
    t1 = p * alpha; x += t1;
    t1 = q * alpha; r -= t1;

    if ((resid = r.norm() / normb) <= tol) {
      // tol = std::abs(resid);
      max_iter = i;
      return 0;
    }

    rho_1 = rho;

    // std::cout << "resid = " << resid << " " << i << " " << alpha << " " << beta << '\n';
  }
  return 1;
}
template int tensor_CG(dtensor<double>& LEnv, dtensor<double>& REnv, dtensor<double>& W, dtensor<double>& x, dtensor<double>& b, int max_iter, double tol);
template int tensor_CG(dtensor< std::complex<double> >& LEnv, dtensor< std::complex<double> >& REnv, dtensor< std::complex<double> >& W, dtensor< std::complex<double> >& x, dtensor< std::complex<double> >& b, int max_iter, double tol);


#endif
