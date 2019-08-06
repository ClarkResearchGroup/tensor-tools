#ifndef BIG_QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS
#define BIG_QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS

#include "big_qstensor.h"

template <typename T>
qstensor<T> big_qstensor<T>::product(qstensor<T>& v){
  // contract pattern:
  // L-v-mid-R
   qstensor<T> res;
   if(L!=nullptr) res = std::move((*L)*v);
   if(res.rank>0)
     res = std::move(res*(mid));
   else
     res = std::move(v*(mid));
  if(R!=nullptr) res = std::move(res*(*R));

  /*qstensor<T> res = std::move(A*v);
  if(R!=nullptr) res = std::move(res*(*R));*/

  res.prime(-1);
  return res;
}
template qstensor<double> big_qstensor<double>::product(qstensor<double>& v);
template qstensor< std::complex<double> > big_qstensor< std::complex<double> >::product(qstensor< std::complex<double> >& v);


template <typename T>
T big_qstensor<T>::expec(qstensor<T>& v){
  // contract pattern:
  // L-v-mid-R
   qstensor<T> res;
   if(L!=nullptr) res = std::move((*L)*v);
   if(res.rank>0)
     res = std::move(res*mid);
   else
     res = std::move(v*(mid));
   if(R!=nullptr) res = std::move(res*(*R));

  assert(1==2);
  //qstensor<T> res = std::move(A*v);
  //if(R!=nullptr) res = std::move(res*(*R));

  res.prime(-1);
  return res.contract(v);
}
template double big_qstensor<double>::expec(qstensor<double>& v);
template std::complex<double> big_qstensor< std::complex<double> >::expec(qstensor< std::complex<double> >& v);

template <typename T>
qstensor<T> big_qstensor<T>::diagonal(){
   /*qstensor<T> res;
   if(L!=nullptr) res = (*L).diagonal();
   for(auto m : mid){
     qstensor<T> dm = m->diagonal();
     res = std::move(res*(dm));
   }
   if(R!=nullptr){
     qstensor<T> dR = R->diagonal();
     res = std::move(res*dR);
   }*/

  qstensor<T> res;
  if(R==nullptr){
    res = std::move((*L)*mid);
    res = std::move(res.diagonal());
  }else{
    //fake operators without storage so we can get final index structure
    qstensor<T> midtemp(mid.idx_set); idxToSparse(mid.idx_set, midtemp._T); midtemp._initted=true;
    qstensor<T> Ltemp(L->idx_set); idxToSparse(L->idx_set,Ltemp._T); Ltemp._initted=true;
    qstensor<T> Rtemp(R->idx_set); idxToSparse(R->idx_set,Rtemp._T); Rtemp._initted=true;
    //res = std::move(A*(*R));
    //res = std::move(res.diagonal());
    res = std::move(Ltemp*midtemp);
    res = std::move(res*Rtemp);
    res = std::move(res.diagonal(All));

    unordered_map<string,char> charMap;
    auto indMid = indicesToCharNP(mid.idx_set,charMap);
    auto indL = indicesToCharNP(L->idx_set,charMap);
    auto indR = indicesToCharNP(R->idx_set,charMap);
    auto indRes = indicesToCharNP(res.idx_set,charMap);
    //perr<<indA<< " "<<indR<< " "<<indRes<<endl;
    //TODO: check effeciency of this contraction
    res._T[indRes.c_str()] = L->_T[indL.c_str()]*mid._T[indMid.c_str()]*R->_T[indR.c_str()];
  }

  // qstensor<T> res = std::move(A.diagonal());

  return res;
}
template qstensor<double> big_qstensor<double>::diagonal();
template qstensor< std::complex<double> > big_qstensor< std::complex<double> >::diagonal();

#endif
