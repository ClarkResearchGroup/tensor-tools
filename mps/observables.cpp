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
			acc = std::move(t1*t2);
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


// void EE(MPS& phi, std::vector<double>& vec){
// 	assert(phi.tensors_allocated);
// 	assert(phi.norm()>1e-12);
// 	MPS psi = phi;
// 	psi.normalize();
// 	int L  = psi.length;
// 	int pD = psi.index_size;
// 	int bD = *std::max_element(psi.bond_dims.begin(), psi.bond_dims.end());
// 	double * sv = new double [pD*bD]();
// 	int row, col, tDim;
// 	Mxd TM, U, V;
// 	for(int i = 0; i < L-1; i++){
// 		row=psi.M[i][0].rows();
// 		col=psi.M[i][0].cols();
// 		TM.resize(pD*row,col);
// 		tDim = std::min(TM.rows(),TM.cols());
// 		for(int j = 0; j < pD; j++){
// 			TM.block(j*row,0,row,col) = psi.M[i][j].block(0,0,row,col);
// 		}
// 		rSVD(TM,tDim,sv,U,V,'r');
// 		double EEn = 0;
// 		for(int j = 0; j < tDim; ++j){
// 			if(sv[j]>1e-15) EEn -= sv[j]*sv[j]*log(sv[j]*sv[j]);
// 		}
// 		vec.push_back(EEn);
// 		for(int tid = 0; tid < pD; tid++){
// 			psi.M[i][tid].setZero();
// 			psi.M[i][tid].block(0,0,row,U.cols()) = U.block(tid*row,0,row,U.cols());
// 			Mxd tempM;
// 			tempM.noalias() = V * psi.M[i+1][tid];
// 			psi.M[i+1][tid].setZero();
// 			psi.M[i+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
// 		}
// 	}
// 	delete [] sv;
// }
//
// //-------------------------------------
// // MPO
// MPO diagonal(const MPO& H){
// 	assert(H.tensors_allocated);
// 	MPO A = H;
// 	int L  = H.length;
// 	int pD = std::round(std::sqrt(H.index_size));
// 	for(int i = 0; i < L; i++){
// 		for(int j = 0; j < pD; j++){
// 			for(int k = j+1; k < pD; k++){
// 				A.M[i][j*pD+k].setZero();
// 			}
// 		}
// 	}
// 	return A;
// }
//
// MPO offdiagonal(const MPO& H){
// 	MPO A=H, B=diagonal(H);
// 	A = A - B;
// 	return A;
// }
//
// void EE(MPO& HH, std::vector<double>& vec){
// 	assert(HH.tensors_allocated);
// 	assert(HH.norm()>1e-15);
// 	MPO H = HH;
// 	H.normalize();
// 	int L  = H.length;
// 	int pD = std::round(std::sqrt(H.index_size));
// 	int bD = *std::max_element(H.bond_dims.begin(), H.bond_dims.end());
// 	double * sv = new double [pD*pD*bD]();
// 	int row, col, tDim;
// 	Mxd TM, U, V;
// 	for(int i = 0; i < L-1; i++){
// 		row=H.M[i][0].rows();
// 		col=H.M[i][0].cols();
// 		TM.resize(pD*pD*row,col);
// 		tDim = std::min(TM.rows(),TM.cols());
// 		for(int j = 0; j < pD*pD; j++){
// 			TM.block(j*row,0,row,col) = H.M[i][j].block(0,0,row,col);
// 		}
// 		rSVD(TM,tDim,sv,U,V,'r');
// 		double EEn = 0;
// 		for(int j = 0; j < tDim; ++j){
// 			if(sv[j]>1e-15) EEn -= sv[j]*sv[j]*log(sv[j]*sv[j]);
// 		}
// 		vec.push_back(EEn);
// 		for(int tid = 0; tid < pD*pD; tid++){
// 			H.M[i][tid].setZero();
// 			H.M[i][tid].block(0,0,row,U.cols()) = U.block(tid*row,0,row,U.cols());
// 			Mxd tempM;
// 			tempM.noalias() = V * H.M[i+1][tid];
// 			H.M[i+1][tid].setZero();
// 			H.M[i+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
// 		}
// 	}
// 	delete [] sv;
// }
//
// double var(const MPO& H){
// 	MPO A = diagonal(H);
// 	return (l2norm(H)-l2norm(A))/std::pow(H.index_size,H.length/2.0);
// }
//
// double trace(const MPO& H){
// 	assert(H.tensors_allocated);
// 	int pD = std::round(std::sqrt(H.index_size));
// 	int L  = H.length;
// 	Mxd A = H.M[0][0]/2;
// 	for(int j = 1; j < pD; ++j) {
// 		A += H.M[0][j*pD+j]/2;
// 	}
// 	for(int i = 1; i < L; ++i){
// 		Mxd tp = H.M[i][0]/2;
// 		for(int j = 1; j < pD; ++j){
// 			tp += H.M[i][j*pD+j]/2;
// 		}
// 		A = A * tp;
// 	}
// 	return A(0,0);
// }
//
// double l2norm(const MPO& H){
// 	MPO A = H;
// 	double nm = A.norm();
// 	return nm*nm;
// }

#endif
