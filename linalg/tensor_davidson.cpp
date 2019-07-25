#ifndef DENSE_TENSOR_BASED_DAVIDSON_SOLVER
#define DENSE_TENSOR_BASED_DAVIDSON_SOLVER

#include "tensor_davidson.h"

template <typename T>
void elemWiseDivide(dtensor<T>& A, double theta, dtensor<T>& B){
  assert(A._initted && B._initted);
  assert(A.size == B.size);
  /*for (unsigned i = 0; i < A.size; i++) {
    if( B._T.data()[i]!=T(theta) )
      A._T.data()[i] = A._T.data()[i]/(B._T.data()[i] - theta);
  }*/
 unordered_map<string,char> charMap;
 auto indA = indicesToChar(A.idx_set,charMap);
 auto indB = indicesToChar(B.idx_set,charMap);
 CTF::Transform<T,T>([theta](T b, T& a){ if(b!=T(theta)) a/=(b-theta); })(B.__T[indB.c_str()],A.__T[indA.c_str()]);
}
template void elemWiseDivide(dtensor<double>& A, double theta, dtensor<double>& B);
template void elemWiseDivide(dtensor< std::complex<double> >& A, double theta, dtensor< std::complex<double> >& B);


template <typename T>
void elemWiseDivide(qtensor<T>& A, double theta, qtensor<T>& B){
  assert(A._initted && B._initted);
  assert(A.rank == B.rank);
  unordered_map<string,char> charMap;
  auto indA = indicesToChar(A.idx_set,charMap);
  auto indB = indicesToChar(B.idx_set,charMap);
  for (auto i = A.block_id_by_qn_str.begin(); i != A.block_id_by_qn_str.end(); ++i){
    string qn_str = i->first;
    unsigned A_id = i->second;
    if(B.block_id_by_qn_str.find(qn_str)!=B.block_id_by_qn_str.end()){
      unsigned B_id = B.block_id_by_qn_str.at(qn_str);
      assert(A._block[A_id].get_tot_size(false) == B._block[B_id].get_tot_size(false));
      CTF::Transform<T,T>([theta](T b, T& a){ if(b!=T(theta)) a/=(b-theta); })(B._block[B_id][indB.c_str()],A._block[A_id][indA.c_str()]);
    }
  }
}
template void elemWiseDivide(qtensor<double>& A, double theta, qtensor<double>& B);
template void elemWiseDivide(qtensor< std::complex<double> >& A, double theta, qtensor< std::complex<double> >& B);

template <typename T>
void elemWiseDivide(qstensor<T>& A, double theta, qstensor<T>& B){
  assert(A._initted && B._initted);
  assert(A.rank == B.rank);
 unordered_map<string,char> charMap;
 auto indA = indicesToChar(A.idx_set,charMap);
 auto indB = indicesToChar(B.idx_set,charMap);
 CTF::Transform<T,T>([theta](T b, T& a){ if(b!=T(theta)) a/=(b-theta); })(B._T[indB.c_str()],A._T[indA.c_str()]);
}
template void elemWiseDivide(qstensor<double>& A, double theta, qstensor<double>& B);
template void elemWiseDivide(qstensor< std::complex<double> >& A, double theta, qstensor< std::complex<double> >& B);

template <typename T, template <typename> class BigTensorType, template <typename> class TensorType>
double tensor_davidson(BigTensorType<T>& A, TensorType<T>& x, int m, int max_restart, double tol, char mode){
  assert(m>=1 && max_restart>=1 && tol>0);
  TensorType<T> u, ua, r;
  TensorType<T>* v  = new TensorType<T> [m];
  TensorType<T>* va = new TensorType<T> [m];
  TensorType<T> D   = std::move(A.diagonal());
  // D.dag();
  // uint_vec perm;
  // find_index_permutation(D.idx_set, x.idx_set, perm);
  // D.permute(perm);
  T* M              = new T [m*m];
  double* evals     = new double [m];
  double eval       = 0;
  double r_norm     = 0.0;
  double u_norm     = 0.0;
  // Restarded Jacobi-Davidson
  for (int i = 0; i < max_restart; i++) {
    // Restart
    for (int j = 0; j < m; j++) {
      // Davidson method
      v[j] = x;
      // (2) Orthogonalize v
      for (int k = j-1; k >= 0; k--) {
        T alpha = v[k].inner_product(v[j]);
        v[j].add(v[k], -alpha);
      }
      v[j].normalize();
      // (3) Get va = (A*v)^{\dage}
      /*for (int k = 0; k < j+1; k++) {
        va[k] = std::move(A.product(v[k]));
      }*/
      va[j] = std::move(A.product(v[j]));
      // (4) Update Hamiltonian projected to subspace
      for (int k = 0; k < j+1; k++) {
        for (int l = 0; l <= k; l++) {
          M[l + k*(j+1)] = va[l].inner_product(v[k]);
          M[k + l*(j+1)] = cconj(M[l + k*(j+1)]);
        }
      }
      if(j>0){
        // (5) Diagonalize M
        DIAG(j+1, M, evals);
        // (6) Get the Ritz vector u with smallest/largest eigenvalue
        if(mode == 'S'){
          eval = evals[0];
          for (int k = 0; k < j+1; k++) {
            if(k==0){
              u = v[k];
              u *= M[0];
              ua = va[k];
              ua *= M[0];
            }else{
              u.add(v[k], M[k]);
              ua.add(va[k], M[k]);
            }
          }
        }else if(mode == 'L'){
          eval = evals[j];
          for (int k = 0; k < j+1; k++) {
            if(k==0){
              u = v[k];
              u *= M[j*(j+1)];
              ua = va[k];
              ua *= M[j*(j+1)];
            }else{
              u.add(v[k], M[k+j*(j+1)]);
              ua.add(va[k], M[k+j*(j+1)]);
            }
          }
        }
      }else{
        eval = std::real(M[0]);
        u = v[0];
        ua = va[0];
      }
      // (8) Get residue vector r = ua - evals[0] * u
      // u_norm = u.norm();
      ua.add(u, -eval);
      r = ua;
      r_norm = r.norm();
      // (9) Expand search space
      if(j<m-1){
        x = r;
        elemWiseDivide(x, eval, D);
        x.normalize();
      }
      // std::cout << "Residue norm = " << r_norm << '\n';
    }
    x = u;
    if(r_norm<tol) break;
  }
  delete [] v;
  delete [] va;
  delete [] evals;
  delete [] M;

  return eval;
}
template double tensor_davidson(big_dtensor<double>& A, dtensor<double>& x, int m, int max_restart, double tol, char mode);
template double tensor_davidson(big_dtensor< std::complex<double> >& A, dtensor< std::complex<double> >& x, int m, int max_restart, double tol, char mode);
template double tensor_davidson(big_qtensor<double>& A, qtensor<double>& x, int m, int max_restart, double tol, char mode);
template double tensor_davidson(big_qtensor< std::complex<double> >& A, qtensor< std::complex<double> >& x, int m, int max_restart, double tol, char mode);
template double tensor_davidson(big_qstensor<double>& A, qstensor<double>& x, int m, int max_restart, double tol, char mode);
template double tensor_davidson(big_qstensor< std::complex<double> >& A, qstensor< std::complex<double> >& x, int m, int max_restart, double tol, char mode);


#endif
