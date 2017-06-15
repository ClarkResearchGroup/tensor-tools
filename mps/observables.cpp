#ifndef My_OBSERVABLES
#define My_OBSERVABLES

#include "observables.h"

//-------------------------------------
// Expectation velues
template <typename T>
T psiHphi (MPS<T>& psi, MPO<T>& H, MPS<T>& phi) {
	assert(psi.tensors_allocated);
	assert(H.tensors_allocated);
	assert(phi.tensors_allocated);
	assert(psi.index_size==phi.index_size);
	assert(psi.length==phi.length);
	assert(psi.length==H.length);
	T res = T(0);
	dtensor<T> acc, t1, t2, t3, t4, t5;
	for (int i = 0; i < psi.length; i++) {
		t1 = std::move( psi.tensorize(i) ); t1.dag();
		t2 = std::move( H.tensorize(i) );
		t3 = std::move( phi.tensorize(i) );
		if(i==0){
			t4 = std::move(t1*t2);
			acc = std::move(t4*t3);
		}else if(i==psi.length-1){
			t4 = std::move(acc*t1);
			t5 = std::move(t4*t2);
			res = t5.contract(t3);
		}else{
			t4 = std::move(acc*t1);
			t5 = std::move(t4*t2);
			acc = std::move(t5*t3);
			// acc.print();
		}
	}
	return res;
}
template double psiHphi (MPS<double>& psi, MPO<double>& H, MPS<double>& phi);
template  std::complex<double>  psiHphi (MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, MPS< std::complex<double> >& phi);

template <typename T>
T psiphi (MPS<T>& psi, MPS<T>& phi) {
	assert(psi.tensors_allocated);
	assert(phi.tensors_allocated);
	assert(psi.index_size==phi.index_size);
	assert(psi.length==phi.length);
	T res = T(0);
	dtensor<T> acc, t1, t2, t3;
	for (int i = 0; i < psi.length; i++) {
		t1 = std::move( psi.tensorize(i) ); t1.dag(); t1.mapPrime(1,0,Site);
		t2 = std::move( phi.tensorize(i) );
		if(i==0){
			acc = (t1*t2);
		}else if(i==psi.length-1){
			t3 = std::move(acc*t1);
			res = t3.contract(t2);
		}else{
			t3 = std::move(acc*t1);
			acc = std::move(t3*t2);
		}
	}
	return res;
}
template double psiphi (MPS<double>& psi, MPS<double>& phi);
template  std::complex<double>  psiphi (MPS< std::complex<double> >& psi, MPS< std::complex<double> >& phi);


//-------------------------------------
// MPO operations
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
