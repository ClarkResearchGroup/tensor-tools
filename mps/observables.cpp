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

#endif
