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
  /*// L-v-mid-R
   qstensor<T> res;
   if(L!=nullptr) res = std::move((*L)*v);
   if(res.rank>0)
     res = std::move(res*mid);
   else
     res = std::move(v*(mid));
   if(R!=nullptr) res = std::move(res*(*R));

  assert(1==2);*/
  //qstensor<T> res = std::move(A*v);
  //if(R!=nullptr) res = std::move(res*(*R));
  qstensor<T> res = std::move(product(v));

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
    /*vector<qtensor_index> L_diag = {L->idx_set[1],L->idx_set[2]};
    vector<qtensor_index> R_diag = {R->idx_set[1],R->idx_set[2]};
    CTF::Tensor<T> LT(2,L->_T.lens+1); LT.sparsify();
    CTF::Tensor<T> RT(2,R->_T.lens+1); RT.sparsify();
    vector<qtensor_index> mid_diag = {mid.idx_set[0],mid.idx_set[1],
                                      mid.idx_set[3],mid.idx_set[5]};
    vector<int64_t> mid_sizes = {mid._T.lens[0],mid._T.lens[1],
                                       mid._T.lens[3],mid._T.lens[5]};
    CTF::Tensor<T> midT(4,mid_sizes.data());midT.sparsify();

    //TODO: only have LT or RT in memory at once?

    auto indMid_diag = indicesToCharNP(mid_diag,charMap);
    auto indL_diag = indicesToCharNP(L_diag,charMap);
    auto indR_diag = indicesToCharNP(R_diag,charMap);
    LT[indL_diag.c_str()] = L->_T[indL.c_str()];
    RT[indR_diag.c_str()] = R->_T[indR.c_str()];
    midT[indMid_diag.c_str()] = mid._T[indMid.c_str()];
    //perr<<indA<< " "<<indR<< " "<<indRes<<endl;
    //TODO: check effeciency of this contraction
    res._T[indRes.c_str()] = LT[indL_diag.c_str()]*midT[indMid_diag.c_str()]*RT[indR_diag.c_str()];*/
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

template <typename T>
size_t big_qstensor<T>::size() const{
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
template size_t big_qstensor<double>::size() const;
template size_t big_qstensor< std::complex<double> >::size() const;

#endif
