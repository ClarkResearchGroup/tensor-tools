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

/*template <typename T>
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
template dtensor< std::complex<double> > big_dtensor< std::complex<double> >::product(dtensor_view< std::complex<double> >& v);*/

template <typename T>
T big_dtensor<T>::expec(dtensor<T>& v){
  // contract pattern:
  // L-v-mid-R
  /*dtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  for(auto m : mid){
    if(res.rank>0)
      res = std::move(res*(*m));
    else
      res = std::move(v*(*m));
  }
  if(R!=nullptr) res = std::move(res*(*R));
  res.prime(-1);*/
  dtensor<T> res = std::move(product(v));
  return res.contract(v);
}
template double big_dtensor<double>::expec(dtensor<double>& v);
template std::complex<double> big_dtensor< std::complex<double> >::expec(dtensor< std::complex<double> >& v);

/*template <typename T>
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
template std::complex<double> big_dtensor< std::complex<double> >::expec(dtensor_view< std::complex<double> >& v);*/

template <typename T>
dtensor<T> big_dtensor<T>::diagonal(){

  //  assert(1==2);
  dtensor<T> res;
  // if(L!=nullptr) (*L).print();
  if(L!=nullptr) res = std::move((*L).diagonal(Link));
  for(auto m : mid){
    dtensor<T> dm = m->diagonal(Site);
    res = std::move(res*dm);
  }
  //res = std::move(res.diagonal(Site));
  // if(R!=nullptr) (*R).print();
  if(R!=nullptr){
    dtensor<T> dR = R->diagonal(Link);
    res = std::move(res*dR);
  }
  //res = std::move(res.diagonal(Link));

  return res;
}
template dtensor<double> big_dtensor<double>::diagonal();
template dtensor< std::complex<double> > big_dtensor< std::complex<double> >::diagonal();


template <typename T>
size_t big_dtensor<T>::size() const{
  if(size_ == 0){
    size_ = 1;
    if(L!=nullptr){
      for(unsigned l=0; l< L->rank; l++){
        if(L->idx_set[l].level() > 0)
          size_ *= L->idx_set[l].size();
      }
    }
    if(R!=nullptr){
      for(unsigned l=0; l< R->rank; l++){
        if(R->idx_set[l].level() > 0)
          size_ *= R->idx_set[l].size();
      }
    }
  }
  /*for(int l=0;l<mid.rank;l++){
    if(mid.idx_set[l].type() == Site)
      size_*= mid.idx_set[l].size();
  }*/
  for(auto m : mid){
    for(unsigned l=0;l<m->rank;l++){
      if(m->idx_set[l].type() == Site)
        size_*= m->idx_set[l].size();
    }
  }
  return size_;
}

template size_t big_dtensor<double>::size() const;
template size_t big_dtensor< std::complex<double> >::size() const;
#endif
