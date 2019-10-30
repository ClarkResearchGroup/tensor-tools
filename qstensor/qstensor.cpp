#ifndef QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS
#define QUANTUM_NUMBERED_SPARSE_TENSOR_CLASS
#include "qstensor.h"

template<typename T>
void idxToSparse(vector<qtensor_index> &idx_set, CTF::Tensor<T> &M){
 vector<int64_t> idx_sizes(idx_set.size(),0);//dense storage index sizes
 for(size_t i=0;i<idx_set.size();i++){
   for(size_t j=0;j<idx_set[i].size();j++){
     idx_sizes[i]+=idx_set[i].qdim(j);
   }
  }
 //M = CTF::Tensor<>(idx_sizes.size(),true,idx_sizes.data());
 M = CTF::Tensor<>(idx_sizes.size(),idx_sizes.data());
}
template<typename T,template <typename> class TensorType >
void blockToSparse(const TensorType<T> &A, CTF::Tensor<T> &M){
 vector<int64_t> idx_sizes(A.rank,0);//dense storage index sizes
 vector<unordered_map<int,int64_t> > offsets(A.rank); //find corner of block in dense

 for(size_t i=0;i<A.idx_set.size();i++){
   int64_t offset =0;
   for(size_t j=0;j<A.idx_set[i].size();j++){
     idx_sizes[i]+=A.idx_set[i].qdim(j);
     offsets[i][A.idx_set[i].qn(j)] = offset;
     offset+= A.idx_set[i].qdim(j);
   }
  }
 //M = CTF::Tensor<>(idx_sizes.size(),true,idx_sizes.data());
 M = CTF::Tensor<>(idx_sizes.size(),idx_sizes.data());
 //place data into M
 vector<int64_t> zeros(A.rank,0);
 for(size_t i=0;i<A._block.size();i++){
   vector<int64_t> starts(A.rank);
   vector<int64_t> ends(A.rank);
   for(unsigned j=0;j<A.rank;j++){
    starts[j] = offsets[j][A.block_index_qn[i][j]];
    ends[j]   = starts[j]+A.block_index_qd[i][j];
   }

   M.slice(starts.data(),ends.data(),0.,A._block[i],zeros.data(),A._block[i].lens,1.);
 }
 M.sparsify();
}
template void blockToSparse(const qtensor<double> &A,  CTF::Tensor<double> &M);
template void blockToSparse(const qstensor<complex<double> > &A, CTF::Tensor<complex<double> > &M);
template void blockToSparse(const qtensor<complex<double> > &A,  CTF::Tensor<complex<double> > &M);
template void blockToSparse(const qstensor<double> &A, CTF::Tensor<double> &M);

template<typename T>
void getOffsets(qstensor<T> &A, vector<unordered_map<int,int64_t> >& blockOffsets){
 blockOffsets.resize(A.rank); //find corner of block in dense
 for(size_t i=0;i<A.idx_set.size();i++){
   int64_t offset =0;
   for(size_t j=0;j<A.idx_set[i].size();j++){
     blockOffsets[i][A.idx_set[i].qn(j)] = offset;
     offset+= A.idx_set[i].qdim(j);
   }
  }
}
template void getOffsets(qstensor<double> &A, vector<unordered_map<int,int64_t> >& blockOffsets);
template void getOffsets(qstensor<std::complex<double> > &A, vector<unordered_map<int,int64_t> >& blockOffsets);
//-----------------------------------------------------------------------------
// Constructors
template <typename T>
qstensor<T>::qstensor(){
  rank = 0;
  _initted = false;
}
template qstensor<double>::qstensor();
template qstensor< std::complex<double> >::qstensor();


template <typename T>
qstensor<T>::qstensor(arr_list arrows){
  rank = 0;
  for(auto a :arrows){
    idx_set.emplace_back(a);
    ++rank;
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_list arrows);
template qstensor< std::complex<double> >::qstensor(arr_list arrows);


template <typename T>
qstensor<T>::qstensor(arr_vec& arrows){
  rank = 0;
  for(auto a :arrows){
    idx_set.emplace_back(a);
    ++rank;
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_vec& arrows);
template qstensor< std::complex<double> >::qstensor(arr_vec& arrows);


template <typename T>
qstensor<T>::qstensor(arr_list arrows, str_list names){
  rank = 0;
  arr_vec arrows_vec;
  str_vec names_vec;
  for(auto a :arrows){
    arrows_vec.push_back(a);
    ++rank;
  }
  for(auto n :names){
    names_vec.push_back(n);
  }
  for (size_t i = 0; i < arrows_vec.size(); i++) {
    idx_set.emplace_back(arrows_vec[i], names_vec[i]);
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_list arrows, str_list names);
template qstensor< std::complex<double> >::qstensor(arr_list arrows, str_list names);


template <typename T>
qstensor<T>::qstensor(arr_vec& arrows, str_vec& names){
  rank = arrows.size();
  for (size_t i = 0; i < arrows.size(); i++) {
    idx_set.emplace_back(arrows[i], names[i]);
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_vec& arrows, str_vec& names);
template qstensor< std::complex<double> >::qstensor(arr_vec& arrows, str_vec& names);


template <typename T>
qstensor<T>::qstensor(arr_list arrows, str_list names, typ_list types){
  rank = 0;
  arr_vec arrows_vec;
  str_vec names_vec;
  typ_vec types_vec;
  for(auto a :arrows){
    arrows_vec.push_back(a);
    ++rank;
  }
  for(auto n :names){
    names_vec.push_back(n);
  }
  for(auto t :types){
    types_vec.push_back(t);
  }
  for (size_t i = 0; i < arrows_vec.size(); i++) {
    idx_set.emplace_back(arrows_vec[i], names_vec[i], types_vec[i]);
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_list arrows, str_list names, typ_list types);
template qstensor< std::complex<double> >::qstensor(arr_list arrows, str_list names, typ_list types);


template <typename T>
qstensor<T>::qstensor(arr_vec& arrows, str_vec& names, typ_vec& types){
  rank = arrows.size();
  for (size_t i = 0; i < arrows.size(); i++) {
    idx_set.emplace_back(arrows[i], names[i], types[i]);
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_vec& arrows, str_vec& names, typ_vec& types);
template qstensor< std::complex<double> >::qstensor(arr_vec& arrows, str_vec& names, typ_vec& types);


template <typename T>
qstensor<T>::qstensor(arr_list arrows, str_list names, typ_list types, uint_list levels){
  rank = 0;
  arr_vec arrows_vec;
  str_vec names_vec;
  typ_vec types_vec;
  uint_vec levels_vec;
  for(auto a :arrows){
    arrows_vec.push_back(a);
    ++rank;
  }
  for(auto n :names){
    names_vec.push_back(n);
  }
  for(auto t :types){
    types_vec.push_back(t);
  }
  for(auto l :levels){
    levels_vec.push_back(l);
  }
  for (size_t i = 0; i < rank; i++) {
    idx_set.emplace_back(arrows_vec[i], names_vec[i], types_vec[i], levels_vec[i]);
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_list arrows, str_list names, typ_list types, uint_list levels);
template qstensor< std::complex<double> >::qstensor(arr_list arrows, str_list names, typ_list types, uint_list levels);


template <typename T>
qstensor<T>::qstensor(arr_vec& arrows, str_vec& names, typ_vec& types, uint_vec& levels){
  rank = arrows.size();
  for (size_t i = 0; i < arrows.size(); i++) {
    idx_set.emplace_back(arrows[i], names[i], types[i], levels[i]);
  }
  _initted = false;
}
template qstensor<double>::qstensor(arr_vec& arrows, str_vec& names, typ_vec& types, uint_vec& levels);
template qstensor< std::complex<double> >::qstensor(arr_vec& arrows, str_vec& names, typ_vec& types, uint_vec& levels);


template <typename T>
qstensor<T>::qstensor(vector<qtensor_index>& qidx_vec){
  rank = qidx_vec.size();
  idx_set = qidx_vec;
  _initted = false;
}
template qstensor<double>::qstensor(vector<qtensor_index>& idx_vec);
template qstensor< std::complex<double> >::qstensor(vector<qtensor_index>& idx_vec);


template <typename T>
qstensor<T>::qstensor(vector<qtensor_index>&& qidx_vec){
  idx_set = std::move(qidx_vec);
  rank = idx_set.size();
  _initted = false;
}
template qstensor<double>::qstensor(vector<qtensor_index>&& idx_vec);
template qstensor< std::complex<double> >::qstensor(vector<qtensor_index>&& idx_vec);


template <typename T>
qstensor<T>::qstensor(initializer_list<qtensor_index> qidx_list){
  rank = 0;
  for(auto i : qidx_list){
    idx_set.push_back(i);
    ++rank;
  }
  _initted = false;
}
template qstensor<double>::qstensor(initializer_list<qtensor_index> qidx_list);
template qstensor< std::complex<double> >::qstensor(initializer_list<qtensor_index> qidx_list);


template <typename T>
qstensor<T>::qstensor(const qstensor<T>& other){
  rank               = other.rank;
  idx_set            = other.idx_set;
  _block             = other._block;
  _T                 = other._T;
  block_index_qn     = other.block_index_qn;
  block_index_qd     = other.block_index_qd;
  block_index_qi     = other.block_index_qi;
  block_id_by_qn_str = other.block_id_by_qn_str;
  _initted           = other._initted;
}
template qstensor<double>::qstensor(const qstensor<double>& other);
template qstensor< std::complex<double> >::qstensor(const qstensor< std::complex<double> >& other);


template <typename T>
qstensor<T>::qstensor(qstensor<T>&& other){
  rank               = other.rank;
  idx_set            = std::move(other.idx_set);
  _block             = std::move(other._block);
  _T                 = std::move(other._T);
  block_index_qn     = std::move(other.block_index_qn);
  block_index_qd     = std::move(other.block_index_qd);
  block_index_qi     = std::move(other.block_index_qi);
  block_id_by_qn_str = std::move(other.block_id_by_qn_str);
  _initted           = other._initted;
}
template qstensor<double>::qstensor(qstensor<double>&& other);
template qstensor< std::complex<double> >::qstensor(qstensor< std::complex<double> >&& other);

template <typename T>
qstensor<T>::qstensor(const qtensor<T>& other){
  rank               = other.rank;
  idx_set            = other.idx_set;
  //_block             = other._block;
  block_index_qn     = other.block_index_qn;
  block_index_qd     = other.block_index_qd;
  block_index_qi     = other.block_index_qi;
  block_id_by_qn_str = other.block_id_by_qn_str;
  _initted           = other._initted;
  //make _T
  blockToSparse(other,_T);
  _block.clear();
}
template qstensor<double>::qstensor(const qtensor<double>& other);
template qstensor< std::complex<double> >::qstensor(const qtensor< std::complex<double> >& other);


template <typename T>
qstensor<T>::qstensor(qtensor<T>&& other){
  rank               = other.rank;
  idx_set            = std::move(other.idx_set);
  //_block             = std::move(other._block);
  block_index_qn     = std::move(other.block_index_qn);
  block_index_qd     = std::move(other.block_index_qd);
  block_index_qi     = std::move(other.block_index_qi);
  block_id_by_qn_str = std::move(other.block_id_by_qn_str);
  _initted           = other._initted;
  //make _T
  blockToSparse(other,_T);
  _block.clear();
}
template qstensor<double>::qstensor(qtensor<double>&& other);
template qstensor< std::complex<double> >::qstensor(qtensor< std::complex<double> >&& other);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// set up quantum number information (qn, qdim) for selected qtensor_index
template <typename T>
void qstensor<T>::addQNtoIndex(unsigned idx, quantum_number qn){
  assert(idx<rank);
  idx_set[idx].addQN(qn);
}
template void qstensor<double>::addQNtoIndex(unsigned idx, quantum_number qn);
template void qstensor< std::complex<double> >::addQNtoIndex(unsigned idx, quantum_number qn);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Initialize legal block
template <typename T>
void qstensor<T>::initBlock(){
  if(rank>0 && !_initted){
    // get strides
    uint_vec stride;
    stride.push_back(1);
    for (size_t i = 0; i < rank; i++) {
      stride.push_back(idx_set[i].size() * stride[i]);
    }
    unsigned n = stride.back();
    // iterate over all combinations of qns
    for (size_t i = 0; i < n; i++) {
      uint_vec pos(rank);
      int_vec qds(rank);
      int_vec  qns(rank);
      unsigned size = 1;
      int totalQN = 0;
      for (size_t j = 0; j < rank; j++) {
        pos[j] = unsigned(i/stride[j])%idx_set[j].size();
        qns[j] = idx_set[j].qn(pos[j]);
        qds[j] = idx_set[j].qdim(pos[j]);
        size *= qds[j];
        if(idx_set[j].arrow()==Inward){
          totalQN += qns[j];
        }else{
          totalQN -= qns[j];
        }
      }
      // legal combination
      if(totalQN==0){
        string qn_str;
        for (size_t j = 0; j < rank; j++) {
          qn_str += (to_string(qns[j])+" ");
        }
        //_block.emplace_back(rank,qds.data());
        //block.emplace_back(size);
        block_index_qn.push_back(std::move(qns));
        block_index_qd.push_back(std::move(qds));
        block_index_qi.push_back(std::move(pos));
        block_id_by_qn_str[qn_str] = block_index_qn.size()-1;
      }
    }
  }
  idxToSparse(idx_set,_T);
  _initted = true;
}
template void qstensor<double>::initBlock();
template void qstensor< std::complex<double> >::initBlock();


template <typename T>
void qstensor<T>::clearBlock(){
  rank = idx_set.size();
  //block.clear();
  _block.clear();
  if(_initted) _T.set_zero();
  block_index_qn.clear();
  block_index_qd.clear();
  block_index_qi.clear();
  block_id_by_qn_str.clear();
  _initted = false;
}
template void qstensor<double>::clearBlock();
template void qstensor< std::complex<double> >::clearBlock();
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Reset
template <typename T>
void qstensor<T>::reset(vector<qtensor_index>& idx_vec){
  clearBlock();
  idx_set = idx_vec;
  rank = idx_set.size();
}
template void qstensor<double>::reset(vector<qtensor_index>& idx_vec);
template void qstensor< std::complex<double> >::reset(vector<qtensor_index>& idx_vec);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Set values
template <typename T>
void qstensor<T>::setRandom(){
  assert(_initted);
  vector<unordered_map<int,int64_t> >offsets;
  getOffsets(*this,offsets);
  vector<int64_t> zeros(rank,0);
  for(size_t i=0;i< block_index_qd.size();i++){
    CTF::Tensor<T> ran(rank,block_index_qd[i].data());
    ran.fill_random(0,1);
   vector<int64_t> starts(rank);
   vector<int64_t> ends(rank);
   for(unsigned j=0;j<rank;j++){
    starts[j] = offsets[j][block_index_qn[i][j]];
    ends[j]   = starts[j]+block_index_qd[i][j];
   }
   _T.slice(starts.data(),ends.data(),0.,ran,zeros.data(),ran.lens,1.);

  }
}
template void qstensor<double>::setRandom();
template void qstensor< std::complex<double> >::setRandom();

template <typename T>
void qstensor<T>::setZero(){
  assert(_initted);
  /*for (size_t i = 0; i < _block.size(); i++) {
    _block[i].set_zero();
  }*/
  _T.set_zero();
}
template void qstensor<double>::setZero();
template void qstensor< std::complex<double> >::setZero();

template <typename T>
void qstensor<T>::setOne(){
  assert(_initted);
  string ind = getIndices();
  /*for (size_t i = 0; i < _block.size(); i++) {
    _block[i][ind.c_str()] = 1.;
  }*/
  _T[ind.c_str()] = 1.;
}
template void qstensor<double>::setOne();
template void qstensor< std::complex<double> >::setOne();
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Permute
template <typename T>
void qstensor<T>::permute(uint_vec& perm){
  assert(_initted);
  bool perm_needed = false;
  for (size_t i = 0; i < perm.size(); i++) {
    if(i!=perm[i]){
      perm_needed = true;
      break;
    }
  }
  if (perm_needed){
    qstensor<T> A(*this);
    vector<int64_t> idx_sizes(rank);
    for (size_t i = 0; i < rank; i++) {
      A.idx_set[i] = idx_set[perm[i]];
      idx_sizes[i] = _T.lens[perm[i]];
    }
    A.block_id_by_qn_str.clear();
    T alpha = 1;
    T beta  = 0;
    unordered_map<string,char> charMap;
    string indOld = getIndices(charMap);
    string indNew = A.getIndices(charMap);
    for (size_t i = 0; i < block_index_qn.size(); i++) {
      string A_qn_str;
      for (size_t j = 0; j < rank; j++) {
        A.block_index_qn[i][j] = block_index_qn[i][perm[j]];
        A.block_index_qd[i][j] = block_index_qd[i][perm[j]];
        A.block_index_qi[i][j] = block_index_qi[i][perm[j]];
        A_qn_str += (to_string(A.block_index_qn[i][j])+" ");
      }
      A.block_id_by_qn_str[A_qn_str] = i;
    }
    //permute by "hand"
    CTF::Tensor<> B(A.rank,true,idx_sizes.data());
    B[indNew.c_str()] = A._T[indOld.c_str()];
    A._T = std::move(B);
    *this = A;
  }
}
template void qstensor<double>::permute(uint_vec& perm);
template void qstensor< std::complex<double> >::permute(uint_vec& perm);

template <typename T>
void qstensor<T>::permute(uint_list perm){
  uint_vec perm_vec;
  for(auto s : perm){
    perm_vec.push_back(s);
  }
  permute(perm_vec);
}
template void qstensor<double>::permute(uint_list perm);
template void qstensor< std::complex<double> >::permute(uint_list perm);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Operator overloading
template <typename T>
qstensor<T>& qstensor<T>::operator=(const qstensor<T>& other){
  if(this!=&other){
    clearBlock();
    rank               = other.rank;
    idx_set            = other.idx_set;
    _block              = other._block;
    _T                 = other._T;
    block_index_qn     = other.block_index_qn;
    block_index_qd     = other.block_index_qd;
    block_index_qi     = other.block_index_qi;
    block_id_by_qn_str = other.block_id_by_qn_str;
    _initted           = other._initted;
  }
  return *this;
}
template qstensor<double>& qstensor<double>::operator=(const qstensor<double> &other);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator=(const qstensor< std::complex<double> > &other);


template <typename T>
qstensor<T>& qstensor<T>::operator=(qstensor<T>&& other){
  if(this!=&other){
    clearBlock();
    rank               = other.rank;
    idx_set            = std::move(other.idx_set);
    _block              = std::move(other._block);
    _T                 = std::move(other._T);
    block_index_qn     = std::move(other.block_index_qn);
    block_index_qd     = std::move(other.block_index_qd);
    block_index_qi     = std::move(other.block_index_qi);
    block_id_by_qn_str = std::move(other.block_id_by_qn_str);
    _initted           = other._initted;
  }
  return *this;
}
template qstensor<double>& qstensor<double>::operator=(qstensor<double>&& other);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator=(qstensor< std::complex<double> >&& other);

template <typename T>
qstensor<T>& qstensor<T>::operator=(const qtensor<T>& other){//convert
  if(true){
    clearBlock();
    rank               = other.rank;
    idx_set            = other.idx_set;
    //_block              = other._block;
    //_T                 = other._T;
    blockToSparse(other,_T);
    block_index_qn     = other.block_index_qn;
    block_index_qd     = other.block_index_qd;
    block_index_qi     = other.block_index_qi;
    block_id_by_qn_str = other.block_id_by_qn_str;
    _initted           = other._initted;
    _block.clear();
  }
  return *this;
}
template qstensor<double>& qstensor<double>::operator=(const qtensor<double> &other);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator=(const qtensor< std::complex<double> > &other);


template <typename T>
qstensor<T>& qstensor<T>::operator=(qtensor<T>&& other){//convert
  if(true){
    clearBlock();
    rank               = other.rank;
    idx_set            = std::move(other.idx_set);
    //_block              = std::move(other._block);
    //_T                 = std::move(other._T);
    blockToSparse(other,_T);
    block_index_qn     = std::move(other.block_index_qn);
    block_index_qd     = std::move(other.block_index_qd);
    block_index_qi     = std::move(other.block_index_qi);
    block_id_by_qn_str = std::move(other.block_id_by_qn_str);
    _initted           = other._initted;
    _block.clear();
  }
  return *this;
}
template qstensor<double>& qstensor<double>::operator=(qtensor<double>&& other);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator=(qtensor< std::complex<double> >&& other);

template <typename T>
qstensor<T> qstensor<T>::operator * (qstensor<T>& other){
  assert(_initted || other._initted);
  if( _initted && !other._initted ){
    qstensor<T> res(*this);
    return res;
  }
  if( other._initted && !_initted ){
    qstensor<T> res(other);
    return res;
  }
  //copy
  //qstensor<T> A(*this);
  //qstensor<T> B(other);
  vector<qtensor_index> A_idx_set = idx_set;
  vector<qtensor_index> B_idx_set = other.idx_set;
  // mark repeated indices
  unordered_map<string,int>  index_marker;
  vector< pair<unsigned,unsigned> > idx_pair;
  for (size_t i = 0; i < rank; i++) {
    index_marker[A_idx_set[i].tag()] = i;
  }
  for (size_t i = 0; i < other.rank; i++) {
    qtensor_index tp = B_idx_set[i];
    tp.dag();
    if(index_marker.find(tp.tag()) == index_marker.end()){
      index_marker[B_idx_set[i].tag()] = i;
    }else{
      idx_pair.push_back( std::make_pair(unsigned(index_marker[tp.tag()]), i) );
      index_marker[tp.tag()] = -1;
      index_marker[B_idx_set[i].tag()] = -1;
    }
  }
  // permute
  uint_vec A_perm;
  uint_vec B_perm;
  for (size_t i = 0; i < rank; i++) {
    if (index_marker[A_idx_set[i].tag()] != -1){
      A_perm.push_back(i);
    }
  }
  for (size_t i = 0; i < idx_pair.size(); i++) {
    A_perm.push_back(idx_pair[i].first);
    B_perm.push_back(idx_pair[i].second);
  }
  for (size_t i = 0; i < other.rank; i++) {
    if (index_marker[B_idx_set[i].tag()] != -1){
      B_perm.push_back(i);
    }
  }
  //A.permute(A_perm);
  //B.permute(B_perm);
  for (size_t i = 0; i < rank; i++) {
    A_idx_set[i] = idx_set[A_perm[i]];
  }
  for (size_t i = 0; i < other.rank; i++) {
    B_idx_set[i] = other.idx_set[B_perm[i]];
  }
  // number of repeated indices
  unsigned num_rep = idx_pair.size();
  // set up new qtensor_index
  vector<qtensor_index> res_index_set;
  for (size_t i = 0; i < rank-num_rep; i++) {
    res_index_set.push_back(A_idx_set[i]);
  }
  for (size_t i = num_rep; i < other.rank; i++) {
    res_index_set.push_back(B_idx_set[i]);
  }
  qstensor<T> res(res_index_set);
  res._initted = true;
  // get blocks info
  set<int> mid_QN_set;
  unordered_map< int, set<uint_vec> > left_index_qi;
  unordered_map< int, set<uint_vec> > mid_index_qi;
  unordered_map< int, set<uint_vec> > right_index_qi;
  for (size_t i = 0; i < block_index_qn.size(); i++) {
    int mid_QN = 0;
    uint_vec left_qi;
    for (size_t j = 0; j < rank-num_rep; j++) {
      left_qi.push_back(block_index_qi[i][A_perm[j]]);
      if(A_idx_set[j].arrow()==Inward){
        mid_QN += block_index_qn[i][A_perm[j]];
      }else{
        mid_QN -= block_index_qn[i][A_perm[j]];
      }
    }
    mid_QN_set.insert(mid_QN);
    left_index_qi[mid_QN].insert(left_qi);
  }
  for (size_t i = 0; i < other.block_index_qn.size(); i++) {
    int mid_QN = 0;
    uint_vec right_qi;
    for (size_t j = num_rep; j < other.rank; j++) {
      right_qi.push_back(other.block_index_qi[i][B_perm[j]]);
      if(B_idx_set[j].arrow()==Inward){
        mid_QN -= other.block_index_qn[i][B_perm[j]];
      }else{
        mid_QN += other.block_index_qn[i][B_perm[j]];
      }
    }
    uint_vec mid_qi;
    for (size_t j = 0; j < num_rep; j++) {
      mid_qi.push_back(other.block_index_qi[i][B_perm[j]]);
    }
    right_index_qi[mid_QN].insert(right_qi);
    mid_index_qi[mid_QN].insert(mid_qi);
  }
  unordered_map<string,char> charMap;
  string indA_L = getIndices(charMap);
  string indB_R = other.getIndices(charMap);
  string indC   = res.getIndices(charMap);
  idxToSparse(res_index_set, res._T);
  res._T[indC.c_str()] = _T[indA_L.c_str()]*other._T[indB_R.c_str()];
  // merge blocks
  for (auto i1 = mid_QN_set.begin(); i1 != mid_QN_set.end(); ++i1){
    int q = *i1;
    const set<uint_vec>& left_qi_set  = left_index_qi[q];
    const set<uint_vec>& right_qi_set = right_index_qi[q];
    const set<uint_vec>& mid_qi_set   = mid_index_qi[q];
    for (auto i2 = right_qi_set.begin(); i2 != right_qi_set.end(); ++i2){
      const uint_vec& right_qi = *i2;
      for (auto i3 = left_qi_set.begin(); i3 != left_qi_set.end(); ++i3){
        const uint_vec& left_qi = *i3;
        uint_vec res_block_index_qi;
        int_vec res_block_index_qd;
        int_vec  res_block_index_qn;
        unsigned res_block_size = 1;
        string   res_qn_str;
        int_vec  A_block_qn(rank);
        int_vec  B_block_qn(other.rank);
        int M=1,N=1;
        for (size_t i = 0; i < left_qi.size(); i++) {
          res_block_index_qi.push_back(left_qi[i]);
          res_block_index_qd.push_back(A_idx_set[i].qdim(left_qi[i]));
          res_block_index_qn.push_back(A_idx_set[i].qn(left_qi[i]));
          res_block_size *= res_block_index_qd.back();
          res_qn_str += (to_string(res_block_index_qn.back())+" ");
          A_block_qn[i] = res_block_index_qn.back();
          M *= res_block_index_qd.back();
        }
        for (size_t i = 0; i < right_qi.size(); i++) {
          res_block_index_qi.push_back(right_qi[i]);
          res_block_index_qd.push_back(B_idx_set[num_rep+i].qdim(right_qi[i]));
          res_block_index_qn.push_back(B_idx_set[num_rep+i].qn(right_qi[i]));
          res_block_size *= res_block_index_qd.back();
          res_qn_str += (to_string(res_block_index_qn.back())+" ");
          B_block_qn[num_rep+i] = res_block_index_qn.back();
          N *= res_block_index_qd.back();
        }
        // std::cout << "\n" << '\n';
        //res._block.emplace_back(res_index_set.size(),res_block_index_qd.data());
        res.block_index_qn.push_back(res_block_index_qn);
        res.block_index_qd.push_back(res_block_index_qd);
        res.block_index_qi.push_back(res_block_index_qi);
        res.block_id_by_qn_str[res_qn_str] = res.block_index_qn.size()-1;
        // sum over all blocks
      }
    }
  }

  /*for(int ii=0;ii<res._block.size();ii++){
    for(int l=0;l<res.rank;l++){
      assert(res._block[ii].lens[l]==res.block_index_qd[ii][l]);
    }
  }*/
  return res;
}
template qstensor<double> qstensor<double>::operator * (qstensor<double>& other);
template qstensor< std::complex<double> > qstensor< std::complex<double> >::operator * (qstensor< std::complex<double> >& other);


#if !defined(USE_HPTT)
template <typename T>
qstensor<T> qstensor<T>::operator + (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // copy
  qstensor<T> res = A;
  // permute
  uint_vec perm;
  find_index_permutation(res.idx_set, idx_set, perm);
  res.permute(perm);
  // add
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < block.size(); i++) {
    string res_qn_str;
    for (size_t j = 0; j < rank; j++) {
      res_qn_str += (to_string(res.block_index_qn[i][j])+" ");
    }
    unsigned idx = res.block_id_by_qn_str[res_qn_str];
    // iterate over all elements
    for (size_t k = 0; k < block[i].size(); k++) {
      res.block[idx][k] += block[i][k];
    }
  }
  return res;
}
#else
template <typename T>
qstensor<T> qstensor<T>::operator + (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // copy
  qstensor<T> res(*this);
  // find permutation
  uint_vec perm;
  find_index_permutation(A.idx_set, res.idx_set, perm);
  // add
  T alpha = 1;
  T beta  = 1;
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  for (size_t i = 0; i < A.block.size(); i++) {
    string res_qn_str;
    int_vec idx_sizes(rank);
    for (size_t j = 0; j < rank; j++) {
      res_qn_str += (to_string(A.block_index_qn[i][perm[j]])+" ");
      idx_sizes[j] = A.block_index_qd[i][j];
    }
    unsigned res_idx = res.block_id_by_qn_str[res_qn_str];
    /*auto plan = hptt::create_plan(
      (int *)perm.data(), rank,
      alpha, A.block[i].data(), idx_sizes.data(), NULL,
      beta, res.block[res_idx].data(), NULL,
      hptt::ESTIMATE,numThreads);
    plan->execute();*/
  }
  return res;
}
#endif
template qstensor<double> qstensor<double>::operator + (qstensor<double>& A);
template qstensor< std::complex<double> > qstensor< std::complex<double> >::operator + (qstensor< std::complex<double> >& A);


#if !defined(USE_HPTT)
template <typename T>
qstensor<T> qstensor<T>::operator - (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // copy
  qstensor<T> res = A;
  // permute
  uint_vec perm;
  find_index_permutation(res.idx_set, idx_set, perm);
  res.permute(perm);
  // add
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < block.size(); i++) {
    string res_qn_str;
    for (size_t j = 0; j < rank; j++) {
      res_qn_str += (to_string(res.block_index_qn[i][j])+" ");
    }
    unsigned idx = res.block_id_by_qn_str[res_qn_str];
    // iterate over all elements
    for (size_t k = 0; k < block[i].size(); k++) {
      res.block[idx][k] = block[i][k] - res.block[idx][k];
    }
  }
  return res;
}
#else
template <typename T>
qstensor<T> qstensor<T>::operator - (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // copy
  qstensor<T> res(*this);
  // find permutation
  uint_vec perm;
  find_index_permutation(A.idx_set, res.idx_set, perm);
  // add
  T alpha = -1;
  T beta  = 1;
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  for (size_t i = 0; i < A.block.size(); i++) {
    string res_qn_str;
    int_vec idx_sizes(rank);
    for (size_t j = 0; j < rank; j++) {
      res_qn_str += (to_string(A.block_index_qn[i][perm[j]])+" ");
      idx_sizes[j] = A.block_index_qd[i][j];
    }
    unsigned res_idx = res.block_id_by_qn_str[res_qn_str];
    /*auto plan = hptt::create_plan(
      (int *)perm.data(), rank,
      alpha, A.block[i].data(), idx_sizes.data(), NULL,
      beta, res.block[res_idx].data(), NULL,
      hptt::ESTIMATE,numThreads);
    plan->execute();*/
  }
  return res;
}
#endif
template qstensor<double> qstensor<double>::operator - (qstensor<double>& A);
template qstensor< std::complex<double> > qstensor< std::complex<double> >::operator - (qstensor< std::complex<double> >& A);


#if !defined(USE_HPTT)
template <typename T>
qstensor<T>& qstensor<T>::operator += (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // permute
  uint_vec perm;
  find_index_permutation(A.idx_set, idx_set, perm);
  A.permute(perm);
  // add
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < block.size(); i++) {
    string A_qn_str;
    for (size_t j = 0; j < rank; j++) {
      A_qn_str += (to_string(A.block_index_qn[i][j])+" ");
    }
    unsigned idx = A.block_id_by_qn_str[A_qn_str];
    // iterate over all elements
    for (size_t k = 0; k < block[i].size(); k++) {
      block[i][k] += A.block[idx][k];
    }
  }
  return *this;
}
#else
template <typename T>
qstensor<T>& qstensor<T>::operator += (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // find permutation
  uint_vec perm;
  find_index_permutation(A.idx_set, idx_set, perm);
  // add
  T alpha = 1;
  T beta  = 1;
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  for (size_t i = 0; i < A.block.size(); i++) {
    string qn_str;
    int_vec idx_sizes(rank);
    for (size_t j = 0; j < rank; j++) {
      qn_str += (to_string(A.block_index_qn[i][perm[j]])+" ");
      idx_sizes[j] = A.block_index_qd[i][j];
    }
    unsigned idx = block_id_by_qn_str[qn_str];
    /*auto plan = hptt::create_plan(
      (int *)perm.data(), rank,
      alpha, A.block[i].data(), idx_sizes.data(), NULL,
      beta, block[idx].data(), NULL,
      hptt::ESTIMATE,numThreads);
    plan->execute();*/
  }
  return *this;
}
#endif
template qstensor<double>& qstensor<double>::operator += (qstensor<double>& A);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator += (qstensor< std::complex<double> >& A);


#if !defined(USE_HPTT)
template <typename T>
qstensor<T>& qstensor<T>::operator -= (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // permute
  uint_vec perm;
  find_index_permutation(A.idx_set, idx_set, perm);
  A.permute(perm);
  // add
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < block.size(); i++) {
    string A_qn_str;
    for (size_t j = 0; j < rank; j++) {
      A_qn_str += (to_string(A.block_index_qn[i][j])+" ");
    }
    unsigned idx = A.block_id_by_qn_str[A_qn_str];
    // iterate over all elements
    for (size_t k = 0; k < block[i].size(); k++) {
      block[i][k] -= A.block[idx][k];
    }
  }
  return *this;
}
#else
template <typename T>
qstensor<T>& qstensor<T>::operator -= (qstensor<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  // find permutation
  uint_vec perm;
  find_index_permutation(A.idx_set, idx_set, perm);
  // add
  T alpha = -1;
  T beta  = 1;
  char* p = std::getenv("OMP_NUM_THREADS");
  int numThreads = 1;
  if(p){
    numThreads = atoi(p);
  }
  omp_set_num_threads(numThreads);
  for (size_t i = 0; i < A.block.size(); i++) {
    string qn_str;
    int_vec idx_sizes(rank);
    for (size_t j = 0; j < rank; j++) {
      qn_str += (to_string(A.block_index_qn[i][perm[j]])+" ");
      idx_sizes[j] = A.block_index_qd[i][j];
    }
    unsigned idx = block_id_by_qn_str[qn_str];
    /*auto plan = hptt::create_plan(
      (int *)perm.data(), rank,
      alpha, A.block[i].data(), idx_sizes.data(), NULL,
      beta, block[idx].data(), NULL,
      hptt::ESTIMATE,numThreads);
    plan->execute();*/
  }
  return *this;
}
#endif
template qstensor<double>& qstensor<double>::operator -= (qstensor<double>& A);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator -= (qstensor< std::complex<double> >& A);


template <typename T>
qstensor<T>& qstensor<T>::operator *= (const T c){
  assert(_initted);
  CTF::Scalar<T> cs(c);
  auto ind = getIndices();
  /*for (size_t i = 0; i < _block.size(); i++) {
    _block[i][ind.c_str()]*=cs[""];
  }*/
  _T[ind.c_str()]*=cs[""];
  return *this;
}
template qstensor<double>& qstensor<double>::operator*=(const double c);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator*=(const std::complex<double> c);


template <typename T>
qstensor<T>& qstensor<T>::operator /= (const T c){
  assert(_initted);
  auto invC = 1./c;
  CTF::Scalar<T> cs(invC);
  string indA = getIndices();
  /*for (size_t i = 0; i < _block.size(); i++) {
    _block[i][indA.c_str()]*=cs[""];
  }*/
  _T[indA.c_str()]*=cs[""];
  return *this;
}
template qstensor<double>& qstensor<double>::operator/=(const double c);
template qstensor< std::complex<double> >& qstensor< std::complex<double> >::operator/=(const std::complex<double> c);


template <typename T>
qstensor<T> qstensor<T>::operator*(const T c){
  assert(1==2);
  assert(_initted);
  qstensor A(*this);
  for (size_t i = 0; i < A.block.size(); i++) {
    for (size_t j = 0; j < A.block[i].size(); j++) {
      A.block[i][j] *= c;
    }
  }
  return A;
}
template qstensor<double> qstensor<double>::operator*(const double c);
template qstensor< std::complex<double> > qstensor< std::complex<double> >::operator*(const std::complex<double> c);


template <typename T>
qstensor<T> qstensor<T>::operator/(const T c){
  assert(1==2);
  assert(_initted);
  qstensor A(*this);
  for (size_t i = 0; i < A.block.size(); i++) {
    for (size_t j = 0; j < A.block[i].size(); j++) {
      A.block[i][j] *= c;
    }
  }
  return A;
}
template qstensor<double> qstensor<double>::operator/(const double c);
template qstensor< std::complex<double> > qstensor< std::complex<double> >::operator/(const std::complex<double> c);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Contract to scalar
template <typename T>
T qstensor<T>::contract(qstensor<T>& A){
  assert(_initted && A._initted);
  assert(rank>0 && A.rank>0);
  uint_vec perm;
  //TODO: more effecient dag (copy vectors only)
  //qstensor<T> B(A); B.dag();
  vector<qtensor_index> A_idx_set = A.idx_set;
  //dag the set
  for(unsigned i=0;i<A.rank;i++) A_idx_set[i].dag();

  //find_index_permutation(B.idx_set, idx_set, perm);
  //B.permute(perm);
  T res = 0;
  unordered_map<string,char> charMap;
  string indA = indicesToChar(A_idx_set,charMap);//B.getIndices(charMap);
  string ind = getIndices(charMap);
  res=_T[ind.c_str()]*A._T[indA.c_str()];
  return res;
}
template double qstensor<double>::contract(qstensor<double>& A);
template std::complex<double> qstensor< std::complex<double> >::contract(qstensor< std::complex<double> >& A);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// special arithmetic operations with another tensor in the format/pattern
template <typename T>
void qstensor<T>::add(qstensor<T>& A, T c){
  assert(A._initted && _initted);
  assert(A.rank == rank);
  //assert(A.block_index_qn.size() == block_index_qn.size()); //assumed but not necessary
  unordered_map<string,char> charMap;
  auto ind = getIndices(charMap);
  auto indA    = A.getIndices(charMap);
  CTF::Scalar<T> cs(c);
  _T[ind.c_str()] += cs[""]*A._T[indA.c_str()];
  /*for (auto i = block_id_by_qn_str.begin(); i != block_id_by_qn_str.end(); ++i){
    string qn_str = i->first;
    unsigned t_id = i->second;
    if(A.block_id_by_qn_str.find(qn_str)!=A.block_id_by_qn_str.end()){
      unsigned A_id = A.block_id_by_qn_str.at(qn_str);
      assert(A._block[A_id].get_tot_size(false) == _block[t_id].get_tot_size(false));
      assert(A._block[A_id].order == _block[t_id].order);
      _block[t_id][ind.c_str()] += cs[""]*A._block[A_id][indA.c_str()];
    }
  }*/
  for (auto i = A.block_id_by_qn_str.begin(); i != A.block_id_by_qn_str.end(); ++i){
    string qn_str = i->first;
    unsigned A_id = i->second;
    //TODO this can be made slightly more effecient
    if(block_id_by_qn_str.find(qn_str)==block_id_by_qn_str.end()){
      block_index_qn.push_back(A.block_index_qn[A_id]);
      block_index_qd.push_back(A.block_index_qd[A_id]);
      block_index_qi.push_back(A.block_index_qi[A_id]);
      block_id_by_qn_str[qn_str] = block_index_qn.size();
    }
  }
}
template void qstensor<double>::add(qstensor<double>& A, double c);
template void qstensor< std::complex<double> >::add(qstensor< std::complex<double> >& A, std::complex<double> c);


template <>
double qstensor<double>::inner_product(qstensor<double>& A){
  assert(A._initted && _initted);
  assert(A.rank == rank);
  unordered_map<string,char> charMap;
  auto ind  =   getIndices(charMap);
  auto indA = A.getIndices(charMap);
  double res = 0;
  res += _T[ind.c_str()]*A._T[indA.c_str()];
  return res;
}
template<> 
std::complex<double> qstensor< std::complex<double> >::inner_product(qstensor< std::complex<double> >& A){

  assert(A._initted && _initted);
  assert(A.rank == rank);
  unordered_map<string,char> charMap;
  auto ind  =   getIndices(charMap);
  auto indA = A.getIndices(charMap);
  double res = 0;
  using cmplx = std::complex<double>;
  res += CTF::Function<double,cmplx,cmplx>([](cmplx l, cmplx r){ return std::real(cconj(l)*r);})(_T[ind.c_str()],A._T[indA.c_str()]);
  return res;
}
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Get diagonal subtensor
// e.g. For a qstensor with index (i1, i1p, i2, i2p), assume
// (i1, i1p) cancel each other's quantum number, and
// (i2, i2p) cancel each other's quantum number.
// Then the diagonal is then made by only two index (i1, i2).
// But this only makes sense if i1 is able to cancel i2 as well.
// So not every qstensor with paired indices can have a sensible
// diagonal subtensor.
template <typename T>
qstensor<T> qstensor<T>::diagonal(){
  assert(_initted && rank>1);
  assert(rank%2==0);
  vector<qtensor_index> new_idx_set;
  vector< std::pair<int,int> > idx_pair;
  uint_vec left_index;
  uint_vec right_index;
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = i+1; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j])){
        if(idx_set[i].level() < idx_set[j].level()){
          new_idx_set.push_back(idx_set[i]);
          left_index.push_back(i);
          right_index.push_back(j);
          idx_pair.push_back(std::make_pair(i,j));
        }else{
          new_idx_set.push_back(idx_set[j]);
          left_index.push_back(j);
          right_index.push_back(i);
          idx_pair.push_back(std::make_pair(j,i));
        }
        break;
      }
    }
  }
  assert(2*new_idx_set.size()==rank);
  // Set up new qstensor
  qstensor<T> res(new_idx_set);
  unordered_map<string,char> charMap;
  auto indNew  = indicesToCharNP(new_idx_set,charMap);
  auto indOrig = indicesToCharNP(idx_set,charMap);
  // extract diagonal blocks
  for (size_t i = 0; i < block_index_qn.size(); i++) {
    uint_vec t_stride;
    t_stride.push_back(1);
    for (size_t j = 0; j < rank; j++) {
      t_stride.push_back(block_index_qd[i][j] * t_stride[j]);
    }
    uint_vec left_qi;
    bool is_diag_block = true;
    for (size_t j = 0; j < left_index.size(); j++) {
      if(block_index_qi[i][left_index[j]] != block_index_qi[i][right_index[j]]){
        is_diag_block = false;
        break;
      }else{
        left_qi.push_back(block_index_qi[i][left_index[j]]);
      }
    }
    if(is_diag_block){
      const uint_vec& res_block_index_qi = left_qi;
      int_vec  res_block_index_qd;
      int_vec  res_block_index_qn;
      unsigned res_block_size = 1;
      string   res_qn_str;
      for (size_t j = 0; j < left_index.size(); j++) {
        res_block_index_qn.push_back(block_index_qn[i][left_index[j]]);
        res_block_index_qd.push_back(block_index_qd[i][left_index[j]]);
        res_block_size *= res_block_index_qd.back();
        res_qn_str += (to_string(res_block_index_qn.back())+" ");
      }
      //res._block.emplace_back(res_block_index_qd.size(),res_block_index_qd.data());
      res.block_index_qn.push_back(res_block_index_qn);
      res.block_index_qd.push_back(res_block_index_qd);
      res.block_index_qi.push_back(res_block_index_qi);
      res.block_id_by_qn_str[res_qn_str] = res.block_index_qn.size()-1;
      // copy values
      //res._block.back()[indNew.c_str()] = _block[i][indOrig.c_str()];
    }
  }
  idxToSparse(new_idx_set,res._T);
  res._T[indNew.c_str()] = _T[indOrig.c_str()];
  res._initted = true;
  return res;
}
template qstensor<double> qstensor<double>::diagonal();
template qstensor< std::complex<double> > qstensor< std::complex<double> >::diagonal();

//this is a dummy function that doesn't allocate new memory,
//so that the proper indices can be setup for the final call 
//as CTF cannot figure out the proper indices along the way
template <typename T>
qstensor<T> qstensor<T>::diagonal(index_type type){
  assert(_initted && rank>1);
  assert(rank%2==0);
  vector<qtensor_index> new_idx_set;
  vector< std::pair<int,int> > idx_pair;
  uint_vec left_index;
  uint_vec right_index;
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = i+1; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j])){
        if(idx_set[i].level() < idx_set[j].level()){
          new_idx_set.push_back(idx_set[i]);
          left_index.push_back(i);
          right_index.push_back(j);
          idx_pair.push_back(std::make_pair(i,j));
        }else{
          new_idx_set.push_back(idx_set[j]);
          left_index.push_back(j);
          right_index.push_back(i);
          idx_pair.push_back(std::make_pair(j,i));
        }
        break;
      }
    }
  }
  assert(2*new_idx_set.size()==rank);
  // Set up new qstensor
  qstensor<T> res(new_idx_set);
  unordered_map<string,char> charMap;
  auto indNew  = indicesToCharNP(new_idx_set,charMap);
  auto indOrig = indicesToCharNP(idx_set,charMap);
  // extract diagonal blocks
  for (size_t i = 0; i < block_index_qn.size(); i++) {
    uint_vec t_stride;
    t_stride.push_back(1);
    for (size_t j = 0; j < rank; j++) {
      t_stride.push_back(block_index_qd[i][j] * t_stride[j]);
    }
    uint_vec left_qi;
    bool is_diag_block = true;
    for (size_t j = 0; j < left_index.size(); j++) {
      if(block_index_qi[i][left_index[j]] != block_index_qi[i][right_index[j]]){
        is_diag_block = false;
        break;
      }else{
        left_qi.push_back(block_index_qi[i][left_index[j]]);
      }
    }
    if(is_diag_block){
      const uint_vec& res_block_index_qi = left_qi;
      int_vec  res_block_index_qd;
      int_vec  res_block_index_qn;
      unsigned res_block_size = 1;
      string   res_qn_str;
      for (size_t j = 0; j < left_index.size(); j++) {
        res_block_index_qn.push_back(block_index_qn[i][left_index[j]]);
        res_block_index_qd.push_back(block_index_qd[i][left_index[j]]);
        res_block_size *= res_block_index_qd.back();
        res_qn_str += (to_string(res_block_index_qn.back())+" ");
      }
      //res._block.emplace_back(res_block_index_qd.size(),res_block_index_qd.data());
      res.block_index_qn.push_back(res_block_index_qn);
      res.block_index_qd.push_back(res_block_index_qd);
      res.block_index_qi.push_back(res_block_index_qi);
      res.block_id_by_qn_str[res_qn_str] = res.block_index_qn.size()-1;
      // copy values
      //res._block.back()[indNew.c_str()] = _block[i][indOrig.c_str()];
    }
  }
  idxToSparse(new_idx_set,res._T);
  //res._T[indNew.c_str()] = _T[indOrig.c_str()];
  res._initted = true;
  return res;
}
template qstensor<double> qstensor<double>::diagonal(index_type type);
template qstensor< std::complex<double> > qstensor< std::complex<double> >::diagonal(index_type type);
//---------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Print Info
template <typename T>
void qstensor<T>::print(unsigned print_level){
  pout<<"-------------------------------------"<<'\n';
  pout<<"(1) Tensor's rank = "<<rank<<'\n';
  pout<<"(2) Tensor's (# non-zero, % sparsity) = ("<<_T.nnz_tot<<","<<(double)_T.nnz_tot/(_T.get_tot_size(false))<<")\n";
  pout<<"(3) Tensor's real edge lengths = ("; for(unsigned l=0;l<rank;l++) pout<<_T.lens[l]<<","; pout<<")\n";
  pout<<"(4) Tensor's index (arrow, name, type, prime level), {(qn, qdim)}"<<'\n';
  for (size_t i = 0; i < rank; i++) {
    pout << "    ";
    if(idx_set[i].arrow()==Inward){
      pout<<"(Inward"<<", ";
    }else{
      pout<<"(Outward"<<", ";
    }
    pout<<idx_set[i].name()<<", ";
    if(idx_set[i].type()==Link){
      pout<<"Link"<<", ";
    }else{
      pout<<"Site"<<", ";
    }
    pout<<idx_set[i].level()<<")"<<", {";
    for (size_t j = 0; j < idx_set[i].size(); j++) {
      pout << "(" << idx_set[i].qn(j) << ", " << idx_set[i].qdim(j) << ") ";
    }
    pout << "}" << '\n';
  }
  if (print_level>0) {
    pout<<"(5) Number of legal QN block = "<<block_index_qn.size()<<'\n';
    pout<<"(6) QN block"<<'\n';
    for (size_t i = 0; i < block_index_qi.size(); i++){
      uint_vec v = block_index_qi[i];
      pout<<"Block "<<i<<" indices: ";
      for (size_t j = 0; j < rank; j++) {
        pout<<"(";
        if(idx_set[j].arrow()==Inward){
          pout<<"Inward"<<", qn=";
        }else{
          pout<<"Outward"<<", qn=";
        }
        pout<<block_index_qn[i][j]<<", qdim="<<block_index_qd[i][j]<<")"<<" ";
      }
      pout << '\n';
    }
    if(print_level>1){
      _T.print();
      pout<<'\n';
    }
  
  }
  pout << "-------------------------------------" << '\n';
}
template void qstensor<double>::print(unsigned print_level);
template void qstensor< std::complex<double> >::print(unsigned print_level);


//-----------------------------------------------------------------------------
// Save/Load
template <typename T>
void qstensor<T>::save(string fn){
  assert(1==2);
  assert(_initted);
  uint_vec idx_arrows;
  str_vec  idx_names;
  uint_vec idx_types;
  uint_vec idx_levels;
  vector<int_vec>  idx_qn(rank);
  vector<uint_vec> idx_qdim(rank);
  for (size_t i = 0; i < rank; i++) {
    idx_arrows.push_back(idx_set[i].arrow());
    idx_names.push_back(idx_set[i].name());
    idx_types.push_back(idx_set[i].type());
    idx_levels.push_back(idx_set[i].level());
    for (size_t j = 0; j < idx_set[i].size(); j++) {
      idx_qn[i].push_back(idx_set[i].qn(j));
      idx_qdim[i].push_back(idx_set[i].qdim(j));
    }
  }
  ezh5::File fh5W (fn, H5F_ACC_TRUNC);
  fh5W["num_blocks"] = block.size();
  fh5W["rank"] = rank;
  fh5W["idx_arrows"] = idx_arrows;
  fh5W["idx_types"] = idx_types;
  fh5W["idx_levels"] = idx_levels;
  for (size_t i = 0; i < rank; i++) {
    fh5W["idx_qn_"+to_string(i)] = idx_qn[i];
    fh5W["idx_qdim_"+to_string(i)] = idx_qdim[i];
  }
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec(idx_names[i].begin(),idx_names[i].end());
    fh5W["idx_name_"+std::to_string(i)] = vec;
  }
  for (size_t i = 0; i < block.size(); i++) {
    fh5W["block_"+to_string(i)] = block[i];
    fh5W["block_"+to_string(i)+"_qi"] = block_index_qi[i];
  }
}
template void qstensor<double>::save(string fn);
template void qstensor< std::complex<double> >::save(string fn);


template <typename T>
void qstensor<T>::save(ezh5::Node& fh5W){
  assert(_initted);
  uint_vec idx_arrows;
  str_vec  idx_names;
  uint_vec idx_types;
  uint_vec idx_levels;
  vector<int_vec>  idx_qn(rank);
  vector<uint_vec> idx_qdim(rank);
  for (size_t i = 0; i < rank; i++) {
    idx_arrows.push_back(idx_set[i].arrow());
    idx_names.push_back(idx_set[i].name());
    idx_types.push_back(idx_set[i].type());
    idx_levels.push_back(idx_set[i].level());
    for (size_t j = 0; j < idx_set[i].size(); j++) {
      idx_qn[i].push_back(idx_set[i].qn(j));
      idx_qdim[i].push_back(idx_set[i].qdim(j));
    }
  }
  fh5W["num_blocks"] = block_index_qn.size();
  fh5W["rank"] = rank;
  fh5W["idx_arrows"] = idx_arrows;
  fh5W["idx_types"] = idx_types;
  fh5W["idx_levels"] = idx_levels;
  for (size_t i = 0; i < rank; i++) {
    fh5W["idx_qn_"+to_string(i)] = idx_qn[i];
    fh5W["idx_qdim_"+to_string(i)] = idx_qdim[i];
  }
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec(idx_names[i].begin(),idx_names[i].end());
    fh5W["idx_name_"+std::to_string(i)] = vec;
  }
  for (size_t i = 0; i < block_index_qn.size(); i++) {
    fh5W["block_"+to_string(i)+"_qi"] = block_index_qi[i];
  }
}
template void qstensor<double>::save(ezh5::Node& fW);
template void qstensor< std::complex<double> >::save(ezh5::Node& fW);


template <typename T>
void qstensor<T>::load(string fn){
  assert(1==2);
  uint_vec idx_arrows_int;
  arr_vec  idx_arrows;
  str_vec  idx_names;
  uint_vec idx_types_int;
  typ_vec  idx_types;
  uint_vec idx_levels;
  unsigned num_blocks;
  ezh5::File fh5R (fn, H5F_ACC_RDONLY);
  fh5R["num_blocks"] >> num_blocks;
  fh5R["rank"] >> rank;
  fh5R["idx_arrows"] >> idx_arrows_int;
  fh5R["idx_types"] >> idx_types_int;
  fh5R["idx_levels"] >> idx_levels;
  idx_set.clear();
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec;
    fh5R["idx_name_"+to_string(i)] >> vec;
    string s = string(vec.begin(),vec.end());
    idx_names.push_back(s);
    idx_arrows.push_back(arrow_type(idx_arrows_int[i]));
    idx_types.push_back(index_type(idx_types_int[i]));
    idx_set.push_back(qtensor_index(idx_arrows[i],idx_names[i],idx_types[i],idx_levels[i]));
    int_vec idx_qn;
    uint_vec idx_qdim;
    fh5R["idx_qn_"+to_string(i)] >> idx_qn;
    fh5R["idx_qdim_"+to_string(i)] >> idx_qdim;
    idx_set[i].addQN(idx_qn, idx_qdim);
  }
  qstensor<T> A(idx_set);
  A._initted = true;
  for (size_t i = 0; i < num_blocks; i++) {
    uint_vec qi_vec;
    int_vec qd_vec;
    int_vec  qn_vec;
    string   qn_str;
    fh5R["block_"+to_string(i)+"_qi"] >> qi_vec;
    for (size_t j = 0; j < rank; j++) {
      qn_vec.push_back(idx_set[j].qn(qi_vec[j]));
      qd_vec.push_back(idx_set[j].qdim(qi_vec[j]));
      qn_str += (to_string(qn_vec.back())+" ");
    }
    A.block.push_back(vector<T>(1));
    A.block_index_qi.push_back(qi_vec);
    A.block_index_qd.push_back(qd_vec);
    A.block_index_qn.push_back(qn_vec);
    A.block_id_by_qn_str[qn_str] = i;
    fh5R["block_"+to_string(i)] >> A.block.back();
  }
  (*this) = A;
}
template void qstensor<double>::load(string fn);
template void qstensor< std::complex<double> >::load(string fn);


template <typename T>
void qstensor<T>::load(ezh5::Node& fh5R){
  uint_vec idx_arrows_int;
  arr_vec  idx_arrows;
  str_vec  idx_names;
  uint_vec idx_types_int;
  typ_vec  idx_types;
  uint_vec idx_levels;
  unsigned num_blocks;
  fh5R["num_blocks"] >> num_blocks;
  fh5R["rank"] >> rank;
  fh5R["idx_arrows"] >> idx_arrows_int;
  fh5R["idx_types"] >> idx_types_int;
  fh5R["idx_levels"] >> idx_levels;
  idx_set.clear();
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec;
    fh5R["idx_name_"+to_string(i)] >> vec;
    string s = string(vec.begin(),vec.end());
    idx_names.push_back(s);
    idx_arrows.push_back(arrow_type(idx_arrows_int[i]));
    idx_types.push_back(index_type(idx_types_int[i]));
    idx_set.push_back(qtensor_index(idx_arrows[i],idx_names[i],idx_types[i],idx_levels[i]));
    int_vec idx_qn;
    uint_vec idx_qdim;
    fh5R["idx_qn_"+to_string(i)] >> idx_qn;
    fh5R["idx_qdim_"+to_string(i)] >> idx_qdim;
    idx_set[i].addQN(idx_qn, idx_qdim);
  }
  qstensor<T> A(idx_set);
  A._initted = true;
  for (size_t i = 0; i < num_blocks; i++) {
    uint_vec qi_vec;
    int_vec qd_vec;
    int_vec  qn_vec;
    string   qn_str;
    fh5R["block_"+to_string(i)+"_qi"] >> qi_vec;
    for (size_t j = 0; j < rank; j++) {
      qn_vec.push_back(idx_set[j].qn(qi_vec[j]));
      qd_vec.push_back(idx_set[j].qdim(qi_vec[j]));
      qn_str += (to_string(qn_vec.back())+" ");
    }
    A.block_index_qi.push_back(qi_vec);
    A.block_index_qd.push_back(qd_vec);
    A.block_index_qn.push_back(qn_vec);
    A.block_id_by_qn_str[qn_str] = i;
  }
  assert(A.block_index_qi.size()==num_blocks);
  assert(A.block_index_qd.size()==num_blocks);
  assert(A.block_index_qn.size()==num_blocks);
  assert(A.block_id_by_qn_str.size()==num_blocks);
  A.initBlock();
  (*this) = A;
}
template void qstensor<double>::load(ezh5::Node& fR);
template void qstensor< std::complex<double> >::load(ezh5::Node& fR);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Prime level manipulation
template <typename T>
void qstensor<T>::prime(int inc){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].prime(inc);
  }
}
template void qstensor<double>::prime(int inc);
template void qstensor< std::complex<double> >::prime(int inc);

template <typename T>
void qstensor<T>::primeLink(int inc){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeLink(inc);
  }
}
template void qstensor<double>::primeLink(int inc);
template void qstensor< std::complex<double> >::primeLink(int inc);

template <typename T>
void qstensor<T>::primeSite(int inc){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeSite(inc);
  }
}
template void qstensor<double>::primeSite(int inc);
template void qstensor< std::complex<double> >::primeSite(int inc);

template <typename T>
void qstensor<T>::mapPrime(unsigned from, unsigned to){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].mapPrime(from, to);
  }
}
template void qstensor<double>::mapPrime(unsigned from, unsigned to);
template void qstensor< std::complex<double> >::mapPrime(unsigned from, unsigned to);

template <typename T>
void qstensor<T>::mapPrime(unsigned from, unsigned to, index_type type){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].mapPrime(from, to, type);
  }
}
template void qstensor<double>::mapPrime(unsigned from, unsigned to, index_type type);
template void qstensor< std::complex<double> >::mapPrime(unsigned from, unsigned to, index_type type);

template <typename T>
void qstensor<T>::dag(){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].dag();
  }
}
template void qstensor<double>::dag();
template void qstensor< std::complex<double> >::dag();

template <>
void qstensor<double>::conj(){}

template <>
void qstensor<std::complex<double> >::conj(){
  auto ind = getIndices();
  using cmplx = std::complex<double>;
  CTF::Transform<cmplx>([ind](cmplx & d){ d= std::conj(d); })(_T[ind.c_str()]);
}
//---------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Norm
template<> double qstensor<double>::norm(){
  assert(_initted);
  double res = 0.0;

  return std::real(_T.norm2());
}
//template double qstensor<double>::norm();
//template double qstensor< std::complex<double> >::norm();

template<> double qstensor<std::complex<double> >::norm(){
  assert(_initted);
  double res = 0.0;
  string indthis = getIndices();
  using cmplx = std::complex<double>;
  CTF::Scalar<double> tot;
  auto A  = _T[indthis.c_str()];
  tot[""] += CTF::Function<double,cmplx,cmplx>([](cmplx l, cmplx r){ return std::real(cconj(l)*r);})(A,A);
  res+= tot;
  return std::sqrt(res);
}

template <typename T>
double qstensor<T>::normalize(){
  assert(_initted);
  double res = norm();
  auto ind = getIndices();
  if( res!=0.) _T[ind.c_str()]*=1./res;
  //else         assert(1==2);
  return res;
}
template double qstensor<double>::normalize();
template double qstensor< std::complex<double> >::normalize();
//-----------------------------------------------------------------------------
// Index help
template <typename T>
string qstensor<T>::getIndices(){
  char ch='a';
  _indices = string(ch,rank);
  for (unsigned i=0;i<rank;i++){ _indices[i]=ch++; }
  return _indices;
}
template string qstensor<double>::getIndices();
template string qstensor<std::complex<double> >::getIndices();

template <typename T>
string qstensor<T>::getIndices(unordered_map<string,char> &charMap){
  _indices = indicesToChar(idx_set,charMap);
  return _indices;
}
template string qstensor<double>::getIndices(unordered_map<string,char> &charMap);
template string qstensor<std::complex<double> >::getIndices(unordered_map<string,char> &charMap);
//-----------------------------------------------------------------------------


#endif
