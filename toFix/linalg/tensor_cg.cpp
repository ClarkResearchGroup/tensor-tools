#ifndef DENSE_TENSOR_BASED_CG_SOLVER
#define DENSE_TENSOR_BASED_CG_SOLVER

#include "tensor_cg.h"

template <typename T>
void elemWiseDivide(dtensor<T>& A, dtensor<T>& B){
  assert(1==2);
  /*assert(A._initted && B._initted);
  assert(A.size == B.size);
  for (size_t i = 0; i < A.size; i++) {
    if( B._T.data()[i]!=T(0) && B._T.data()[i]>1e-16 )
      A._T.data()[i] = A._T.data()[i]/B._T.data()[i];
  }*/
}
template void elemWiseDivide(dtensor<double>& A, dtensor<double>& B);
template void elemWiseDivide(dtensor< std::complex<double> >& A, dtensor< std::complex<double> >& B);


template <typename T>
void elemWiseDivide(qtensor<T>& A, qtensor<T>& B){
  assert(A._initted && B._initted);
  assert(A.rank == B.rank);
  for (auto i = A.block_id_by_qn_str.begin(); i != A.block_id_by_qn_str.end(); ++i){
    string qn_str = i->first;
    unsigned A_id = i->second;
    if(B.block_id_by_qn_str.find(qn_str)!=B.block_id_by_qn_str.end()){
      unsigned B_id = B.block_id_by_qn_str.at(qn_str);
      assert(A.block[A_id].size() == B.block[B_id].size());
      for (size_t j = 0; j < A.block[A_id].size(); j++) {
        if( B.block[B_id][j]!=T(0) )
          A.block[A_id][j] = A.block[A_id][j]/B.block[B_id][j];
      }
    }
  }
}
template void elemWiseDivide(qtensor<double>& A, qtensor<double>& B);
template void elemWiseDivide(qtensor< std::complex<double> >& A, qtensor< std::complex<double> >& B);


template <typename T, template <typename> class BigTensorType, template <typename> class TensorType1, template <typename> class TensorType2>
int tensor_CG(BigTensorType<T>& A, TensorType1<T>& x, TensorType2<T>& b, int max_iter, double tol, bool preconditioning){
  TensorType1<T> D, t1, t2;
  if(preconditioning){
    D = std::move(A.diagonal());
    // uint_vec perm;
    // find_index_permutation(D.idx_set, b.idx_set, perm);
    // D.permute(perm);
  }
  T resid, alpha, beta, rho, rho_1, normb;
  normb = b.norm();
  TensorType1<T> r, p, z, q;
  t1 = std::move(A.product(x));
  r = std::move(b - t1);
  if (normb == 0.0)
    normb = 1;
  if ((resid = r.norm() / normb) <= tol) {
    return 0;
  }
  for (int i = 1; i <= max_iter; i++){
    z = r;
    if(preconditioning) elemWiseDivide(z,D);
    if(preconditioning){
      r.dag(); r.conj();
      rho = r.contract(z);
      r.dag(); r.conj();
    }else{
      rho = r.norm();
      rho *= rho;
    }
    if(i==1)
      p = z;
    else{
      beta = rho / rho_1;
      t1 = std::move(p * beta);  p = std::move(z + t1);
    }
    q = std::move(A.product(p)); p.conj();
    alpha = rho / p.contract(q); q.dag(); p.conj();
    t1 = std::move(p * alpha); x += t1;
    t1 = std::move(q * alpha); r -= t1;
    if ((resid = r.norm() / normb) <= tol) {
      return 0;
    }
    rho_1 = rho;
    // std::cout << "resid = " << resid << " " << i << " " << alpha << " " << beta << '\n';
  }
  return 1;
}
template int tensor_CG(big_dtensor<double>& A, dtensor<double>& x, dtensor<double>& b, int max_iter, double tol, bool preconditioning);
template int tensor_CG(big_dtensor< std::complex<double> >& A, dtensor< std::complex<double> >& x, dtensor< std::complex<double> >& b, int max_iter, double tol, bool preconditioning);
template int tensor_CG(big_dtensor<double>& A, dtensor<double>& x, dtensor_view<double>& b, int max_iter, double tol, bool preconditioning);
template int tensor_CG(big_dtensor< std::complex<double> >& A, dtensor< std::complex<double> >& x, dtensor_view< std::complex<double> >& b, int max_iter, double tol, bool preconditioning);

template int tensor_CG(big_qtensor<double>& A, qtensor<double>& x, qtensor<double>& b, int max_iter, double tol, bool preconditioning);
template int tensor_CG(big_qtensor< std::complex<double> >& A, qtensor< std::complex<double> >& x, qtensor< std::complex<double> >& b, int max_iter, double tol, bool preconditioning);

#endif
