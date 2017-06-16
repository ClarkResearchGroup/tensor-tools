#ifndef MPS_MPO_OPERATIONS
#define MPS_MPO_OPERATIONS

#include "mps_mpo_ops.h"


template <typename T>
void KroneckerProd(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& B, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& res){
  assert(A.size()>0 && B.size()>0);
  res.setZero(A.rows()*B.rows(), A.cols()*B.cols());
  int r1 = A.rows();
  int c1 = A.cols();
  int r2 = B.rows();
  int c2 = B.cols();
  for(int i = 0; i < r1; ++i)
  {
    for(int j = 0; j < c1; ++j)
    {
      res.block(i*r2, j*c2, r2, c2) = A(i,j)*B;
    }
  }
}
template void KroneckerProd(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& A, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& B, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& res);
template void KroneckerProd(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& A, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& B, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>& res);

// res = A * psi (A on top of psi)
template <typename T>
void fitApplyMPO(MPO<T>& A, MPS<T>& psi, MPS<T>& res, char transA, bool to_compress, double cutoff){
  assert(psi.tensors_allocated && A.tensors_allocated);
  assert(psi.index_size==A.phy_size);
  assert(psi.length==A.length);
  MPS<T> phi(psi.length, psi.index_size, 1);
  int pD = psi.index_size;
  int L  = psi.length;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* tp = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> [pD];
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp;
  for(int i = 0; i < L; ++i){
    for(int j = 0; j < pD; ++j){
      for(int k = 0; k < pD; ++k){
        // note that when transA = 'T', only transposition is applied
        // to achieve Hermitian conjugation, modify the code so that it includes
        // a conjugation as well.
        if(transA!='T' || transA!='t')
          KroneckerProd<T>(A.M[i][j*pD+k], psi.M[i][k], temp);
        else
          KroneckerProd<T>(A.M[i][k*pD+j], psi.M[i][k], temp);
        if(k==0)
          tp[j] = temp;
        else
          tp[j] += temp;
      }
    }
    for(int j = 0; j < pD; ++j){
      phi.M[i][j] = tp[j];
    }
  }
  for(int i = 0; i < L; ++i){
    phi.bond_dims[i] = phi.M[i][0].rows();
  }
  delete [] tp;
  if (to_compress) {
    int max_iter = 10;
    res.fit(phi, max_iter, cutoff);
  }else{
    res = phi;
  }
}
template void fitApplyMPO(MPO<double>& A, MPS<double>& psi, MPS<double>& res, char transA, bool to_compress, double cutoff);
template void fitApplyMPO(MPO< std::complex<double> >& A, MPS< std::complex<double> >& psi, MPS< std::complex<double> >& res, char transA, bool to_compress, double cutoff);

// res = A * B (A on top of B)
template <typename T>
void fitApplyMPO(MPO<T>& A, MPO<T>& B, MPO<T>& res, char transA, bool to_compress, double cutoff){
  assert(A.tensors_allocated && B.tensors_allocated);
  assert(A.index_size==B.index_size);
  assert(A.length==B.length);
  MPO<T> C(A.length, A.phy_size, 1);
  int pD = A.phy_size;
  int L  = A.length;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* tp = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> [pD*pD];
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> temp;
	for(int i = 0; i < L; ++i){
		for(int j = 0; j < pD*pD; ++j){
			int topIdx = j/pD; int botIdx = j%pD;
      // note that when transA = 'T', only transposition is applied
      // to achieve Hermitian conjugation, modify the code so that it includes
      // a conjugation as well.
			for(int k = 0; k < pD; ++k){
				if(transA!='T' || transA!='t')
					KroneckerProd(A.M[i][topIdx*pD+k], B.M[i][k*pD+botIdx], temp);
				else
					KroneckerProd(A.M[i][k*pD+topIdx], B.M[i][k*pD+botIdx], temp);
				if(k==0)
					tp[j] = temp;
				else
					tp[j] += temp;
			}
		}
		for(int j = 0; j < pD*pD; ++j){
			C.M[i][j] = tp[j];
		}
	}
  for(int i = 0; i < L; ++i){
    C.bond_dims[i] = C.M[i][0].rows();
  }
  delete [] tp;
  if (to_compress) {
    int max_iter = 10;
    res.fit(C, max_iter, cutoff);
  }else{
    res = C;
  }
}
template void fitApplyMPO(MPO<double>& A, MPO<double>& B, MPO<double>& res, char transA, bool to_compress, double cutoff);
template void fitApplyMPO(MPO< std::complex<double> >& A, MPO< std::complex<double> >& B, MPO< std::complex<double> >& res, char transA, bool to_compress, double cutoff);

template <typename T>
MPO<T> diagonal(const MPO<T>& H){
	assert(H.tensors_allocated);
	MPO<T> A = H;
	int L  = H.length;
	int pD = H.phy_size;
	for(int i = 0; i < L; i++){
		for(int j = 0; j < pD; j++){
			for(int k = j+1; k < pD; k++){
				A.M[i][j*pD+k].setZero();
			}
		}
	}
	return A;
}
template MPO<double> diagonal(const MPO<double>& H);
template MPO< std::complex<double> > diagonal(const MPO< std::complex<double> >& H);

template <typename T>
MPS<T> DiagonalMPOAsMPS(const MPO<T>& H){
  assert(H.tensors_allocated);
  MPS<T> psi(H.length, H.phy_size, 1);
  for (int i = 0; i < psi.length; i++) {
    for (int j = 0; j < psi.index_size; j++) {
      psi.M[i][j] = H.M[i][j*H.phy_size+j];
    }
  }
  for(int i = 0; i < psi.length; ++i){
    psi.bond_dims[i] = psi.M[i][0].rows();
  }
  return psi;
}
template MPS<double> DiagonalMPOAsMPS(const MPO<double>& H);
template MPS< std::complex<double> > DiagonalMPOAsMPS(const MPO< std::complex<double> >& H);

template <typename T>
MPO<T> MPSAsDiagonalMPO(const MPS<T>& psi){
  assert(psi.tensors_allocated);
  MPO<T> H(psi.length, psi.index_size, 1);
  H.setBondDims(psi.bond_dims);
  H.setZero();
  for (int i = 0; i < H.length; i++) {
    for (int j = 0; j < H.phy_size; j++) {
      H.M[i][j*H.phy_size+j] = psi.M[i][j];
    }
  }
  return H;
}
template MPO<double> MPSAsDiagonalMPO(const MPS<double>& psi);
template MPO< std::complex<double> > MPSAsDiagonalMPO(const MPS< std::complex<double> >& psi);

template <typename T>
MPO<T> offdiagonal(const MPO<T>& H){
	MPO<T> A=H, B=diagonal(H);
	A = A - B;
	return A;
}
template MPO<double> offdiagonal(const MPO<double>& H);
template MPO< std::complex<double> > offdiagonal(const MPO< std::complex<double> >& H);

template <typename T>
double var(const MPO<T>& H){
	MPO<T> A = diagonal(H);
	return (l2norm(H)-l2norm(A))/std::pow(H.phy_size,H.length);
}
template double var(const MPO<double>& H);
template double var(const MPO< std::complex<double> >& H);

template <typename T>
T trace(const MPO<T>& H){
	assert(H.tensors_allocated);
	int pD = H.phy_size;
	int L  = H.length;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> A = H.M[0][0]/2;
	for(int j = 1; j < pD; ++j) {
		A += H.M[0][j*pD+j]/2;
	}
	for(int i = 1; i < L; ++i){
		Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tp = H.M[i][0]/2;
		for(int j = 1; j < pD; ++j){
			tp += H.M[i][j*pD+j]/2;
		}
		A = A * tp;
	}
	return A(0,0);
}
template double trace(const MPO<double>& H);
template std::complex<double> trace(const MPO< std::complex<double> >& H);

template <typename T>
double l2norm(const MPO<T>& H){
	MPO<T> A = H;
	double nm = A.norm();
	return nm*nm;
}
template double l2norm(const MPO<double>& H);
template double l2norm(const MPO< std::complex<double> >& H);

#endif
