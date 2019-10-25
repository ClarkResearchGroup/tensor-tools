#ifndef BIG_QUANTUM_NUMBERED_TENSOR_CLASS
#define BIG_QUANTUM_NUMBERED_TENSOR_CLASS

#include "big_qtensor.h"

template <typename T>
qtensor<T> big_qtensor<T>::product(qtensor<T>& v){
  // contract pattern:
  // L-v-mid-R
  qtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  if(res.rank>0)
    res = std::move(res*(mid));
  else
    res = std::move(v*(mid));
  if(R!=nullptr) res = std::move(res*(*R));

  /*qtensor<T> res = std::move(A*v);
  if(R!=nullptr) res = std::move(res*(*R));*/

  res.prime(-1);
  return res;
}
template qtensor<double> big_qtensor<double>::product(qtensor<double>& v);
template qtensor< std::complex<double> > big_qtensor< std::complex<double> >::product(qtensor< std::complex<double> >& v);


template <typename T>
T big_qtensor<T>::expec(qtensor<T>& v){
  // contract pattern:
  // L-v-mid-R
  /*qtensor<T> res;
  if(L!=nullptr) res = std::move((*L)*v);
  if(res.rank>0)
    res = std::move(res*(mid));
  else
    res = std::move(v*(mid));
  
  if(R!=nullptr) res = std::move(res*(*R));*/

  /*qtensor<T> res = std::move(A*v);
  if(R!=nullptr) res = std::move(res*(*R));*/

  //res.prime(-1);
  qtensor<T> res = std::move(product(v));
  return res.contract(v);
}
template double big_qtensor<double>::expec(qtensor<double>& v);
template std::complex<double> big_qtensor< std::complex<double> >::expec(qtensor< std::complex<double> >& v);


template <typename T>
qtensor<T> big_qtensor<T>::diagonal(){
  qtensor<T> res;
  if(L!=nullptr) res = (*L);
  res = std::move(res*(mid));
  
  if(R!=nullptr) res = std::move(res*(*R));
  res = std::move(res.diagonal());

  /*qtensor<T> res;
  if(R==nullptr){
    res = std::move(A.diagonal());
  }else{
    res = std::move(A*(*R));
    res = std::move(res.diagonal());
    //qtensor<T> Rd = R->diagonal();
    //res = std::move(A.diagonal());
    //res = std::move(res*Rd);
  }*/

  // qtensor<T> res = std::move(A.diagonal());

  return res;
}
template qtensor<double> big_qtensor<double>::diagonal();
template qtensor< std::complex<double> > big_qtensor< std::complex<double> >::diagonal();

template <typename T>
size_t big_qtensor<T>::size() const{
  if(size_ == -1){
    size_ = 1;
    if(L!=nullptr){
      for(int l=0; l< L->rank; l++){
        if(L->idx_set[l].level() > 0)
          size_ *= L->idx_set[l].size();
      }
    }
    if(R!=nullptr){
      for(int l=0; l< R->rank; l++){
        if(R->idx_set[l].level() > 0)
          size_ *= R->idx_set[l].size();
      }
    }
  }
  for(int l=0;l<mid.rank;l++){
    if(mid.idx_set[l].type() == Site)
      size_*= mid.idx_set[l].size();
  }
  return size_;
}

template size_t big_qtensor<double>::size() const;
template size_t big_qtensor< std::complex<double> >::size() const;
#endif
