#ifndef BIG_QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS
#define BIG_QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS

#include "big_qstensor.h"

template <typename T>
qstensor<T> big_qstensor<T>::product(qstensor<T>& v){
  // contract pattern:
  // L-v-mid-R
  // qstensor<T> res;
  // if(L!=nullptr) res = std::move((*L)*v);
  // for(auto m : mid){
  //   if(res.rank>0)
  //     res = std::move(res*(*m));
  //   else
  //     res = std::move(v*(*m));
  // }
  // if(R!=nullptr) res = std::move(res*(*R));

  qstensor<T> res = std::move(A*v);
  if(R!=nullptr) res = std::move(res*(*R));

  res.prime(-1);
  return res;
}
template qstensor<double> big_qstensor<double>::product(qstensor<double>& v);
template qstensor< std::complex<double> > big_qstensor< std::complex<double> >::product(qstensor< std::complex<double> >& v);


template <typename T>
T big_qstensor<T>::expec(qstensor<T>& v){
  // contract pattern:
  // L-v-mid-R
  // qstensor<T> res;
  // if(L!=nullptr) res = std::move((*L)*v);
  // for(auto m : mid){
  //   if(res.rank>0)
  //     res = std::move(res*(*m));
  //   else
  //     res = std::move(v*(*m));
  // }
  // if(R!=nullptr) res = std::move(res*(*R));

  qstensor<T> res = std::move(A*v);
  if(R!=nullptr) res = std::move(res*(*R));

  res.prime(-1);
  return res.contract(v);
}
template double big_qstensor<double>::expec(qstensor<double>& v);
template std::complex<double> big_qstensor< std::complex<double> >::expec(qstensor< std::complex<double> >& v);


template <typename T>
qstensor<T> big_qstensor<T>::diagonal(){
  // qstensor<T> res;
  // if(L!=nullptr) res = (*L);
  // for(auto m : mid){
  //   res = std::move(res*(*m));
  // }
  // if(R!=nullptr) res = std::move(res*(*R));
  // res = std::move(res.diagonal());

  qstensor<T> res;
  if(R==nullptr){
    res = std::move(A.diagonal());
  }else{
    res = std::move(A*(*R));
    res = std::move(res.diagonal());
    //qstensor<T> Rd = R->diagonal();
    //res = std::move(A.diagonal());
    //res = std::move(res*Rd);
  }

  // qstensor<T> res = std::move(A.diagonal());

  return res;
}
template qstensor<double> big_qstensor<double>::diagonal();
template qstensor< std::complex<double> > big_qstensor< std::complex<double> >::diagonal();

#endif
