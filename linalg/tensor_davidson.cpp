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

template <typename T, template <typename> class BigTensorType, template <typename> class TensorType>
double tensor_davidsonIT(BigTensorType<T>& A, TensorType<T>& x, int m, int max_restart, double tol, char mode){
  assert(m>=1 && max_restart>=1 && tol>0);
  auto x_norm = x.norm();
  x *= 1./x_norm;
  auto maxsize = 4;//A.rank //not sure
  auto actual_maxiter = std::min(max_restart,maxsize-1);

  TensorType<T>* v  = new TensorType<T> [actual_maxiter+2];
  TensorType<T>* va = new TensorType<T> [actual_maxiter+2];
  
  T* M              = new T [(actual_maxiter+2)*(actual_maxiter+2)];
  T* Mtemp          = new T [(actual_maxiter+2)*(actual_maxiter+2)];
  double* evals     = new double [(actual_maxiter+2)];
  double* evecs     = new double [(actual_maxiter+2)*(actual_maxiter+2)];
  double eval       = 0;
  double qnorm = 0.;

  double last_lambda = 1000.;
  double Approx0 = 1e-12;

  v[0] = x;
  va[0] = std::move(A.product(v[0]));

  double initEn = std::real(va[0].inner_product(v[0])); //TODO check for cmplx
  size_t t = 0;
  size_t iter = 0;

  double lambda = 0.;
  for(size_t ii=0;ii<actual_maxiter+1;ii++){
    //diag v*A*v
    //compute residual q
    int ni = ii+1;
    auto& q = v[ni];

    //Step A (or I) of Davidson (1975)
    if(ii==0){
      lambda = initEn;
      for(int mi=0;mi<ni;mi++){ //M is (ni,ni)
        M[mi] = lambda;
        Mtemp[mi] = lambda;
      }
      q = va[0];
      q.add(v[0],-lambda); 
    }
    else{
      /*for(int mi=0;mi<ni*ni;mi++)
        M[mi] *= -1;*/
      //copy M->Mtemp since it gets destroyed
      //assume that we'll add a row/column
      for(int mi =0;mi<ni;mi++){
        for(int mj=0;mj<=mi;mj++){
          Mtemp[mj + mi*(ni+1)] = M[mj + mi*(ni)];
          Mtemp[mi + mj*(ni+1)] = M[mi + mj*(ni)];
        }
      }
      DIAG(ni,M,evals);
      /*for(int mi=0;mi<ni*ni;mi++)
        M[mi] *= -1;
      for(int k=0;k<ni;k++)
        evals[k]*=-1;*/
      lambda = evals[0];
      //perr<<"L="<<lambda<<endl;
      x = v[0]; x*=M[0];//x*=evals[0];
      q = va[0]; q*=M[0];//q*=evals[0];
      for(int k=1;k<ii;k++){
        x.add(v[k],M[k]);
        q.add(va[k],M[k]);
      }
      //Step B of Davidson (1975)
      //Calculate residual q
      q.add(x,-lambda);
      /*if(std::real(M[0]) < 0){
        x *= -1.;
        q *= -1.;
      }*/
    }
    //Step C of Davidson (1975)
    //Check convergence
    qnorm = q.norm();
    bool converged = (qnorm < tol && std::abs(lambda-last_lambda) < tol) || qnorm < std::max(Approx0,tol*1e-3);
    last_lambda = lambda;
    if((qnorm < 1e-20) || (converged && ii >= 1) || (ii==actual_maxiter))
      goto done;
    //Step D of Davidson (1975)
    //Apply Davidson preconditioner
    //TODO

    //Step E and F of Davidson (1975)
    //Do Gram-Schmidt on d (Npass times)
    //to include it in the subbasis

    int Npass = 1;
    auto Vq = std::vector<T>(ni);
    int pass = 1;
    int tot_pass = 0;
    while(pass <=Npass){
      ++tot_pass;
      for(int k=0;k<ni;k++){
        Vq[k] = q.contract(v[k]); //includes dag
      }
      for(int k=0;k<ni;k++){
        q.add(v[k],-Vq[k]);
      }
      auto qnrm = q.norm();
      if(qnrm<1e-10){
        //assert(1==2);
        //Orthogonalization failure,
        //try randomizing
        q = v[ni-1];
        q.setRandom();
        qnrm = q.norm();
        --pass;
        if(ni >= maxsize){
          //Not be possible to orthogonalize if
          //max size of q (vecSize after randomize)
          //is size of current basis
          goto done;
        }
        if(tot_pass > Npass*3){
          goto done;
        }
      }
      q *= 1./qnrm;
      ++pass; 
    }
    //Step G of Davidson (1975)
    //Expand AV and M
    //for next step
    va[ni] = std::move(A.product(v[ni]));
    //Step H of Davidson (1975)
    //Add new row and column to M
    swap(Mtemp,M);
    //for(int mi=0;mi<ni*ni;mi++) perr<<M[mi]<<" ";perr<<endl;
    int l = ni;
    for (int k = 0; k < ni+1; k++) {
      M[l + k*(ni+1)] = va[l].inner_product(v[k]);
      M[k + l*(ni+1)] = cconj(M[l + k*(ni+1)]);
    }
    //for(int mi=0;mi<(ni+1)*(ni+1);mi++) perr<<M[mi]<<" ";perr<<endl;

  } //(ii)
  
   done:
  x.normalize();
  delete [] v;
  delete [] va;
  delete [] evals;
  delete [] M;
  delete [] Mtemp;
  
  return lambda;


  


}
template double tensor_davidsonIT(big_dtensor<double>& A, dtensor<double>& x, int m, int max_restart, double tol, char mode);
template double tensor_davidsonIT(big_dtensor< std::complex<double> >& A, dtensor< std::complex<double> >& x, int m, int max_restart, double tol, char mode);
template double tensor_davidsonIT(big_qtensor<double>& A, qtensor<double>& x, int m, int max_restart, double tol, char mode);
template double tensor_davidsonIT(big_qtensor< std::complex<double> >& A, qtensor< std::complex<double> >& x, int m, int max_restart, double tol, char mode);
template double tensor_davidsonIT(big_qstensor<double>& A, qstensor<double>& x, int m, int max_restart, double tol, char mode);
template double tensor_davidsonIT(big_qstensor< std::complex<double> >& A, qstensor< std::complex<double> >& x, int m, int max_restart, double tol, char mode);

#endif
