#ifndef BIG_DENSE_TENSOR_CLASS
#define BIG_DENSE_TENSOR_CLASS

#include "big_dtensor.h"

template <typename T>
dtensor<T> big_dtensor<T>::product(dtensor<T>& v){
  //assert(1==2);
  // contract pattern:
  // L-v-mid-R
  dtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  for(auto m : mid){
    if(res.rank>0)
      res = std::move(res*(*m));
    else
      res = std::move(v*(*m));
  }
  if(R!=nullptr) res = std::move(res*(*R));
  res.prime(-1);
  return res;
}
template dtensor<double> big_dtensor<double>::product(dtensor<double>& v);
template dtensor< std::complex<double> > big_dtensor< std::complex<double> >::product(dtensor< std::complex<double> >& v);

template <typename T>
dtensor<T> big_dtensor<T>::product(dtensor_view<T>& v){
    assert(1==2);
  // contract pattern:
  // L-v-mid-R
  dtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  for(auto m : mid){
    if(res.rank>0)
      res = std::move(res*(*m));
    else
      res = std::move(v*(*m));
  }
  if(R!=nullptr) res = std::move(res*(*R));
  res.prime(-1);
  return res;
}
template dtensor<double> big_dtensor<double>::product(dtensor_view<double>& v);
template dtensor< std::complex<double> > big_dtensor< std::complex<double> >::product(dtensor_view< std::complex<double> >& v);

template <typename T>
T big_dtensor<T>::expec(dtensor<T>& v){
    assert(1==2);
  // contract pattern:
  // L-v-mid-R
  dtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  for(auto m : mid){
    if(res.rank>0)
      res = std::move(res*(*m));
    else
      res = std::move(v*(*m));
  }
  if(R!=nullptr) res = std::move(res*(*R));
  res.prime(-1);
  return res.contract(v);
}
template double big_dtensor<double>::expec(dtensor<double>& v);
template std::complex<double> big_dtensor< std::complex<double> >::expec(dtensor< std::complex<double> >& v);

template <typename T>
T big_dtensor<T>::expec(dtensor_view<T>& v){
    assert(1==2);
  // contract pattern:
  // L-v-mid-R
  dtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  for(auto m : mid){
    if(res.rank>0)
      res = std::move(res*(*m));
    else
      res = std::move(v*(*m));
  }
  if(R!=nullptr) res = std::move(res*(*R));
  res.prime(-1);
  v.conj();
  T ans = res.contract(v);
  v.conj();
  return ans;
}
template double big_dtensor<double>::expec(dtensor_view<double>& v);
template std::complex<double> big_dtensor< std::complex<double> >::expec(dtensor_view< std::complex<double> >& v);

template <typename T>
dtensor<T> big_dtensor<T>::diagonal(){
  //  assert(1==2);
  dtensor<T> res;
  // if(L!=nullptr) (*L).print();
  if(L!=nullptr) res = std::move((*L).diagonal(Link));
  for(auto m : mid){
    res = std::move(res*(*m));
  }
  res = std::move(res.diagonal(Site));
  // if(R!=nullptr) (*R).print();
  if(R!=nullptr) res = std::move(res*(*R));
  res = std::move(res.diagonal(Link));
  return res;
}
template dtensor<double> big_dtensor<double>::diagonal();
template dtensor< std::complex<double> > big_dtensor< std::complex<double> >::diagonal();

#endif
