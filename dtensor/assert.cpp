#include "dtensor.h"
#include "dtensor_op.h"

typedef vector<unsigned int>    len_vec;
typedef vector<char>  lab_vec;

template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto s : idx_sizes) {
    idx_set.push_back(dtensor_index(s));
    //idx_lens.push_back(tblis::len_type(s));
    ++rank;
    size *= s;
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto s : idx_sizes) {
    idx_set.push_back(dtensor_index(s));
    //idx_lens.push_back(tblis::len_type(s));
    ++rank;
    size *= s;
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  uint_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(v_sizes[i],v_names[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names, typ_list types){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  uint_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  typ_vec v_types(types.begin(), types.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(v_sizes[i],v_names[i],v_types[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names, typ_list types);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names, typ_list types);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i],types[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types);



template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i],types[i],levels[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  //_T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, T* data_array){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  uint_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  typ_vec v_types(types.begin(), types.end());
  uint_vec v_levels(levels.begin(), levels.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(v_sizes[i],v_names[i],v_types[i],v_levels[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  //_T.reset(idx_lens);
  //std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, double* data_array);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, std::complex<double>* data_array);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, T* data_array){
  assert(1==2);
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i],types[i],levels[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  //_T.reset(idx_lens);
  //std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, double* data_array);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, std::complex<double>* data_array);



template <typename T>
dtensor<T>::dtensor(vector<dtensor_index>& idx_vec, T* data_array){
  assert(1==2);
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    idx_lens.push_back(idx_vec[i].size());
    ++rank;
    size *= idx_vec[i].size();
  }
  //_T.reset(idx_lens);
  //std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(vector<dtensor_index>& idx_vec, double* data_array);
template dtensor< std::complex<double> >::dtensor(vector<dtensor_index>& idx_vec, std::complex<double>* data_array);


/*template <typename T>
dtensor<T>::dtensor(const dtensor_view<T>& other){
  assert(1==2);
  rank = other.rank;
  size = other.size;
  idx_set = other.idx_set;
  //_T.reset(other._T);
  _initted = other._initted;
}
template dtensor<double>::dtensor(const dtensor_view<double>& other);
template dtensor< std::complex<double> >::dtensor(const dtensor_view< std::complex<double> >& other);*/


//---------------------------------------------------------------------------
// Resize indices
// (data preserved when dimension of indices lowered, filled with val when enlarged)
template <typename T>
void dtensor<T>::resize(uint_vec& new_sizes, T val){
  assert(1==2);
  assert(_initted);
  assert(new_sizes.size()==rank);
  std::vector<long int> v;
  size = 1;
  for (size_t i = 0; i < rank; i++) {
    if(idx_set[i].size()!=new_sizes[i]) idx_set[i].resize(new_sizes[i]);
    v.push_back(new_sizes[i]);
    size *= new_sizes[i];
  }
  //_T.resize(v, val);
}
template void dtensor<double>::resize(uint_vec& new_sizes, double val);
template void dtensor< std::complex<double> >::resize(uint_vec& new_sizes, std::complex<double> val);


template <typename T>
void dtensor<T>::resize(uint_list new_sizes, T val){
  assert(1==2);
  uint_vec new_sizes_vec;
  for(auto v : new_sizes){
    new_sizes_vec.push_back(v);
  }
  resize(new_sizes_vec, val);
}
template void dtensor<double>::resize(uint_list new_sizes, double val);
template void dtensor< std::complex<double> >::resize(uint_list new_sizes, std::complex<double> val);



template <typename T>
void dtensor<T>::resize(vector<dtensor_index>& new_idx_set, T val){
  assert(1==2);
  assert(_initted);
  assert(new_idx_set.size()==rank);
  idx_set = new_idx_set;
  std::vector<long int> v;
  size = 1;
  for (size_t i = 0; i < rank; i++) {
    v.push_back(idx_set[i].size());
    size *= idx_set[i].size();
  }
  //_T.resize(v, val);
}
template void dtensor<double>::resize(vector<dtensor_index>& new_idx_set, double val);
template void dtensor< std::complex<double> >::resize(vector<dtensor_index>& new_idx_set, std::complex<double> val);



//-----------------------------------------------------------------------------
// Permute
#if !defined(USE_HPTT)
template <typename T>
void dtensor<T>::permute(uint_vec& perm)
{
    assert(1==2);
  assert(_initted);
  bool perm_needed = false;
  for (size_t i = 0; i < perm.size(); i++) {
    if(i!=perm[i]){
      perm_needed = true;
      break;
    }
  }
  if (perm_needed){
    vector<dtensor_index> idx_set_p(idx_set);
    for (size_t i = 0; i < rank; i++) {
      idx_set[i] = idx_set_p[perm[i]];
    }
    uint_vec stride1, stride2;
    int s1 = 1, s2 = 1;
    stride1.push_back(s1);
    stride2.push_back(s2);
    for (size_t i = 0; i < rank-1; i++) {
      s1 *= idx_set_p[i].size();
      s2 *= idx_set[i].size();
      stride1.push_back(s1);
      stride2.push_back(s2);
    }
    T* A = new T [size]();
    //std::copy(_T.data(), _T.data()+size, A);
    ////////////////////////////////
    char* p = std::getenv("OMP_NUM_THREADS");
    int numThreads = 1;
    if(p){
      numThreads = atoi(p);
    }
    omp_set_num_threads(numThreads);
    #pragma omp parallel for default(shared)
    for (size_t i = 0; i < size; i++) {
      int old_idx[rank];
      int new_idx[rank];
      for (size_t j = 0; j < rank; j++) {
        old_idx[j] = int(i/stride1[j])%idx_set_p[j].size();
      }
      for (size_t j = 0; j < rank; j++) {
        new_idx[j] = old_idx[perm[j]];
      }
      int idx = 0;
      for (size_t j = 0; j < rank; j++) {
        idx += new_idx[j] * stride2[j];
      }
      //_T.data()[idx] = A[i];
    }
    delete [] A;
  }
}
#else
template <typename T>
void dtensor<T>::permute(uint_vec& perm){
    assert(1==2);
  assert(_initted);
  bool perm_needed = false;
  for (size_t i = 0; i < perm.size(); i++) {
    if(i!=perm[i]){
      perm_needed = true;
      break;
    }
  }
  if (perm_needed){
    vector<dtensor_index> idx_set_p(idx_set);
    T* A = new T [size];
    //std::copy(_T.data(), _T.data()+size, A);
    T* B; //= _T.data();
    T alpha = 1;
    T beta  = 0;
    int idx_sizes[rank];
    for (size_t i = 0; i < rank; i++) {
      idx_sizes[i] = idx_set[i].size();
      idx_set[i] = idx_set_p[perm[i]];
    }
    char* p = std::getenv("OMP_NUM_THREADS");
    int numThreads = 1;
    if(p){
      numThreads = atoi(p);
    }
    omp_set_num_threads(numThreads);
    // auto plan = hptt::create_plan((int *)perm.data(),rank,alpha,A,idx_sizes,NULL,beta,B,NULL,hptt::ESTIMATE,numThreads);
    auto plan = hptt::create_plan((int *)perm.data(),rank,alpha,A,idx_sizes,NULL,beta,B,NULL,hptt::PATIENT,numThreads);
    // auto plan = hptt::create_plan(perm.data(),rank,alpha,A,idx_sizes,NULL,beta,B,NULL,hptt::MEASURE,numThreads);
    // auto plan = hptt::create_plan(perm.data(),rank,alpha,A,idx_sizes,NULL,beta,B,NULL,hptt::PATIENT,numThreads);
    plan->execute();
    delete [] A;
  }
}
#endif
template void dtensor<double>::permute(uint_vec& perm);
template void dtensor< std::complex<double> >::permute(uint_vec& perm);

template <typename T>
void dtensor<T>::permute(uint_list perm){
    assert(1==2);
  uint_vec perm_vec;
  for(auto s : perm){
    perm_vec.push_back(s);
  }
  permute(perm_vec);
}
template void dtensor<double>::permute(uint_list perm);
template void dtensor< std::complex<double> >::permute(uint_list perm);
//-----------------------------------------------------------------------------


/*template <typename T>
dtensor<T>& dtensor<T>::operator=(const dtensor_view<T>& other){
    assert(1==2);
  rank = other.rank;
  size = other.size;
  idx_set = other.idx_set;
  //_T.reset(other._T);
  
  _initted = other._initted;
  return *this;
}
template dtensor<double>& dtensor<double>::operator=(const dtensor_view<double> &other);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator=(const dtensor_view< std::complex<double> > &other);*/


template <typename T>
dtensor<T> dtensor<T>::operator*(const T c){
    assert(1==2);
  assert(_initted);
  dtensor A(*this);
  /*for (size_t i = 0; i < A.size; i++) {
    A._T.data()[i] *= c;
  }*/
  return A;
}



template dtensor<double> dtensor<double>::operator*(const double c);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator*(const std::complex<double> c);



// template <typename T>
// dtensor<T> dtensor<T>::operator * (dtensor<T>& A){
//   assert(_initted || A._initted);
//   if( _initted && !A._initted ){
//     dtensor<T> res(*this);
//     return res;
//   }
//   if( A._initted && !_initted ){
//     dtensor<T> res(A);
//     return res;
//   }
//   vector<dtensor_index> res_index_set;
//   index_sets_difference(idx_set, A.idx_set, res_index_set);
//   assert(res_index_set.size()>0); // result cannnot be a scalar
//   lab_vec this_labels;
//   lab_vec A_labels;
//   lab_vec res_labels;
//   char ch = 'a';
//   unordered_map<string,char> labels_map;
//   for (size_t i = 0; i < rank; i++) {
//     if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
//       this_labels.push_back(ch);
//       labels_map[idx_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       this_labels.push_back(labels_map.at(idx_set[i].tag()));
//     }
//   }
//   for (size_t i = 0; i < A.rank; i++) {
//     if(labels_map.find(A.idx_set[i].tag()) == labels_map.end()){
//       A_labels.push_back(ch);
//       labels_map[A.idx_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
//     }
//   }
//   for (size_t i = 0; i < res_index_set.size(); i++) {
//     if(labels_map.find(res_index_set[i].tag()) == labels_map.end()){
//       res_labels.push_back(ch);
//       labels_map[res_index_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       res_labels.push_back(labels_map.at(res_index_set[i].tag()));
//     }
//   }
//   dtensor<T> res(res_index_set);
//   tblis::mult(T(1),_T,this_labels.data(),A._T,A_labels.data(),T(0),res._T,res_labels.data());
//   return res;
// }
// template dtensor<double> dtensor<double>::operator * (dtensor<double>& A);
// template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor< std::complex<double> >& A);

//
// template <typename T>
// dtensor<T> dtensor<T>::operator * (dtensor<T>& A){
//   assert(_initted || A._initted);
//   if( _initted && !A._initted ){
//     dtensor<T> res(*this);
//     return res;
//   }
//   if( A._initted && !_initted ){
//     dtensor<T> res(A);
//     return res;
//   }
//   unordered_map<string,int>  labels_num_map;
//   unordered_map<string,char> labels_char_map;
//   char ch = 'a';
//   // Get the number of times a index appears
//   for (size_t i = 0; i < rank; i++) {
//     if(labels_num_map.find(idx_set[i].tag()) == labels_num_map.end()){
//       labels_num_map[idx_set[i].tag()] = 1;
//       labels_char_map[idx_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       labels_num_map[idx_set[i].tag()] += 1;
//     }
//   }
//   for (size_t i = 0; i < A.rank; i++) {
//     if(labels_num_map.find(A.idx_set[i].tag()) == labels_num_map.end()){
//       labels_num_map[A.idx_set[i].tag()] = 1;
//       labels_char_map[A.idx_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       labels_num_map[A.idx_set[i].tag()] += 1;
//     }
//   }
//   // Set up new dtensor_index
//   vector<dtensor_index> res_index_set;
//   lab_vec this_labels;
//   lab_vec A_labels;
//   lab_vec res_labels;
//   for (size_t i = 0; i < rank; i++) {
//     this_labels.push_back(labels_char_map[idx_set[i].tag()]);
//     if (labels_num_map[idx_set[i].tag()] == 1){
//       res_index_set.push_back(idx_set[i]);
//       res_labels.push_back(labels_char_map[idx_set[i].tag()]);
//     }
//   }
//   for (size_t i = 0; i < A.rank; i++) {
//     A_labels.push_back(labels_char_map[A.idx_set[i].tag()]);
//     if (labels_num_map[A.idx_set[i].tag()] == 1){
//       res_index_set.push_back(A.idx_set[i]);
//       res_labels.push_back(labels_char_map[A.idx_set[i].tag()]);
//     }
//   }
//   dtensor<T> res(res_index_set);
//   tblis::mult(T(1),_T,this_labels.data(),A._T,A_labels.data(),T(0),res._T,res_labels.data());
//   return res;
// }
// template dtensor<double> dtensor<double>::operator * (dtensor<double>& A);
// template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor< std::complex<double> >& A);



/*template <typename T>
dtensor<T> dtensor<T>::operator * (dtensor_view<T>& A){
    assert(1==2);
  assert(_initted || A._initted);
  if( _initted && !A._initted ){
    dtensor<T> res(*this);
    return res;
  }
  if( A._initted && !_initted ){
    dtensor<T> res(A);
    return res;
  }
  unordered_map<string,int>  labels_num_map;
  unordered_map<string,char> labels_char_map;
  char ch = 'a';
  // Get the number of times a index appears
  for (size_t i = 0; i < rank; i++) {
    labels_num_map[idx_set[i].tag()] = 1;
    labels_char_map[idx_set[i].tag()] = ch;
    ++ch;
  }
  for (size_t i = 0; i < A.rank; i++) {
    if(labels_num_map.find(A.idx_set[i].tag()) == labels_num_map.end()){
      labels_num_map[A.idx_set[i].tag()] = 1;
      labels_char_map[A.idx_set[i].tag()] = ch;
      ++ch;
    }else{
      labels_num_map[A.idx_set[i].tag()] += 1;
    }
  }
  // Set up new dtensor_index
  vector<dtensor_index> res_index_set;
  lab_vec this_labels;
  lab_vec A_labels;
  lab_vec res_labels;
  for (size_t i = 0; i < rank; i++) {
    this_labels.push_back(labels_char_map[idx_set[i].tag()]);
    if (labels_num_map[idx_set[i].tag()] == 1){
      res_index_set.push_back(idx_set[i]);
      res_labels.push_back(labels_char_map[idx_set[i].tag()]);
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_char_map[A.idx_set[i].tag()]);
    if (labels_num_map[A.idx_set[i].tag()] == 1){
      res_index_set.push_back(A.idx_set[i]);
      res_labels.push_back(labels_char_map[A.idx_set[i].tag()]);
    }
  }
  dtensor<T> res(res_index_set);
  //tblis::mult(T(1),_T,this_labels.data(),A._T,A_labels.data(),T(0),res._T,res_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator * (dtensor_view<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor_view< std::complex<double> >& A);*/


template <typename T>
dtensor<T> dtensor<T>::operator + (dtensor<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  dtensor<T> res = *this;
  //tblis::add(T(1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator + (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator + (dtensor< std::complex<double> >& A);


/*template <typename T>
dtensor<T> dtensor<T>::operator + (dtensor_view<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  dtensor<T> res = *this;
  //tblis::add(T(1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator + (dtensor_view<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator + (dtensor_view< std::complex<double> >& A);*/


template <typename T>
dtensor<T> dtensor<T>::operator - (dtensor<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  dtensor<T> res = *this;
  //tblis::add(T(-1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator - (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator - (dtensor< std::complex<double> >& A);


/*template <typename T>
dtensor<T> dtensor<T>::operator - (dtensor_view<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  dtensor<T> res = *this;
  //tblis::add(T(-1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator - (dtensor_view<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator - (dtensor_view< std::complex<double> >& A);*/


template <typename T>
dtensor<T>& dtensor<T>::operator += (dtensor<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  //tblis::add(T(1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator += (dtensor<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator += (dtensor< std::complex<double> >& A);


/*template <typename T>
dtensor<T>& dtensor<T>::operator += (dtensor_view<T>& A){
  assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  //tblis::add(T(1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator += (dtensor_view<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator += (dtensor_view< std::complex<double> >& A);*/



template <typename T>
dtensor<T>& dtensor<T>::operator -= (dtensor<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  //tblis::add(T(-1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator -= (dtensor<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator -= (dtensor< std::complex<double> >& A);


/*template <typename T>
dtensor<T>& dtensor<T>::operator -= (dtensor_view<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
      this_labels.push_back(ch);
      labels_map[idx_set[i].tag()] = ch;
      ++ch;
    }else{
      this_labels.push_back(labels_map.at(idx_set[i].tag()));
    }
  }
  for (size_t i = 0; i < A.rank; i++) {
    A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
  }
  //tblis::add(T(-1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator -= (dtensor_view<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator -= (dtensor_view< std::complex<double> >& A);*/



template <typename T>
dtensor<T> dtensor<T>::operator/(const T c){
    assert(1==2);
  assert(_initted);
  dtensor A(*this);
  /*for (size_t i = 0; i < A.size; i++) {
    A._T.data()[i] /= c;
  }*/
  return A;
}
template dtensor<double> dtensor<double>::operator/(const double c);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator/(const std::complex<double> c);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Contract to scalar
// template <typename T>
// T dtensor<T>::contract(dtensor<T>& A){
//   assert(_initted && A._initted);
//   assert(rank>0 && A.rank>0);
//   vector<dtensor_index> res_index_set;
//   index_sets_difference(idx_set, A.idx_set, res_index_set);
//   assert(res_index_set.size()==0); // contract to a scalar
//   lab_vec this_labels;
//   lab_vec A_labels;
//   char ch = 'a';
//   unordered_map<string,char> labels_map;
//   for (size_t i = 0; i < rank; i++) {
//     if(labels_map.find(idx_set[i].tag()) == labels_map.end()){
//       this_labels.push_back(ch);
//       labels_map[idx_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       this_labels.push_back(labels_map.at(idx_set[i].tag()));
//     }
//   }
//   for (size_t i = 0; i < A.rank; i++) {
//     if(labels_map.find(A.idx_set[i].tag()) == labels_map.end()){
//       A_labels.push_back(ch);
//       labels_map[A.idx_set[i].tag()] = ch;
//       ++ch;
//     }else{
//       A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
//     }
//   }
//   T res = 0;
//   tblis::dot(_T,this_labels.data(),A._T,A_labels.data(),res);
//   return res;
// }
// template double dtensor<double>::contract(dtensor<double>& A);
// template std::complex<double> dtensor< std::complex<double> >::contract(dtensor< std::complex<double> >& A);



/*template <typename T>
T dtensor<T>::contract(dtensor_view<T>& A){
    assert(1==2);
  assert(_initted && A._initted);
  assert(rank>0 && A.rank>0);
  vector<dtensor_index> res_index_set;
  index_sets_difference(idx_set, A.idx_set, res_index_set);
  assert(res_index_set.size()==0); // contract to a scalar
  lab_vec this_labels;
  lab_vec A_labels;
  char ch = 'a';
  unordered_map<string,char> labels_map;
  for (size_t i = 0; i < rank; i++) {
    this_labels.push_back(ch);
    labels_map[idx_set[i].tag()] = ch;
    ++ch;
  }
  for (size_t i = 0; i < A.rank; i++) {
    if(labels_map.find(A.idx_set[i].tag()) == labels_map.end()){
      A_labels.push_back(ch);
      labels_map[A.idx_set[i].tag()] = ch;
      ++ch;
    }else{
      A_labels.push_back(labels_map.at(A.idx_set[i].tag()));
    }
  }
  T res = 0;
  //tblis::dot(_T,this_labels.data(),A._T,A_labels.data(),res);
  return res;
}
template double dtensor<double>::contract(dtensor_view<double>& A);
template std::complex<double> dtensor< std::complex<double> >::contract(dtensor_view< std::complex<double> >& A);*/
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Save/Load
template <typename T>
void dtensor<T>::save(string fn){
    assert(1==2);
  assert(_initted);
  uint_vec idx_sizes;
  str_vec idx_names;
  uint_vec idx_types;
  uint_vec idx_levels;
  vector<T> d;//(_T.data(),_T.data()+size);
  for (size_t i = 0; i < rank; i++) {
    idx_sizes.push_back(idx_set[i].size());
    idx_names.push_back(idx_set[i].name());
    idx_types.push_back(idx_set[i].type());
    idx_levels.push_back(idx_set[i].level());
  }
  std::string idx_name_pref = "idx_name_";
  ezh5::File fh5W (fn, H5F_ACC_TRUNC);
  fh5W["rank"] = rank;
  fh5W["size"] = size;
  fh5W["idx_sizes"] = idx_sizes;
  fh5W["idx_types"] = idx_types;
  fh5W["idx_levels"] = idx_levels;
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec(idx_names[i].begin(),idx_names[i].end());
    fh5W[idx_name_pref+std::to_string(i)] = vec;
  }
  fh5W["T"] = d;
}
template void dtensor<double>::save(string fn);
template void dtensor< std::complex<double> >::save(string fn);


template <typename T>
void dtensor<T>::save(ezh5::Node& fW){
    assert(1==2);
  assert(_initted);
  uint_vec idx_sizes;
  str_vec idx_names;
  uint_vec idx_types;
  uint_vec idx_levels;
  vector<T> d;//(_T.data(),_T.data()+size);
  for (size_t i = 0; i < rank; i++) {
    idx_sizes.push_back(idx_set[i].size());
    idx_names.push_back(idx_set[i].name());
    idx_types.push_back(idx_set[i].type());
    idx_levels.push_back(idx_set[i].level());
  }
  std::string idx_name_pref = "idx_name_";
  fW["rank"] = rank;
  fW["size"] = size;
  fW["idx_sizes"] = idx_sizes;
  fW["idx_types"] = idx_types;
  fW["idx_levels"] = idx_levels;
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec(idx_names[i].begin(),idx_names[i].end());
    fW[idx_name_pref+std::to_string(i)] = vec;
  }
  fW["T"] = d;
}
template void dtensor<double>::save(ezh5::Node& fW);
template void dtensor< std::complex<double> >::save(ezh5::Node& fW);


template <typename T>
void dtensor<T>::load(string fn){
    assert(1==2);
  uint_vec idx_sizes;
  str_vec idx_names;
  uint_vec idx_types_int;
  typ_vec idx_types;
  uint_vec idx_levels;
  len_vec idx_lens;
  vector<T> d;
  std::string idx_name_pref = "idx_name_";
  ezh5::File fh5R (fn, H5F_ACC_RDONLY);
  fh5R["rank"] >> rank;
  fh5R["size"] >> size;
  fh5R["idx_sizes"] >> idx_sizes;
  fh5R["idx_types"] >> idx_types_int;
  fh5R["idx_levels"] >> idx_levels;
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec;
    fh5R[idx_name_pref+std::to_string(i)] >> vec;
    std::string s = std::string(vec.begin(),vec.end());
    idx_names.push_back(s);
  }
  idx_set.clear();
  for (size_t i = 0; i < rank; i++) {
    idx_types.push_back(index_type(idx_types_int[i]));
    idx_set.push_back(dtensor_index(idx_sizes[i],idx_names[i],idx_types[i],idx_levels[i]));
    //idx_lens.push_back(tblis::len_type(idx_sizes[i]));
  }
  //_T.reset(idx_lens);
  fh5R["T"] >> d;
  //std::copy(d.begin(),d.end(),_T.data());
  _initted = true;
}
template void dtensor<double>::load(string fn);
template void dtensor< std::complex<double> >::load(string fn);


template <typename T>
void dtensor<T>::load(ezh5::Node& fR){
    assert(1==2);
  uint_vec idx_sizes;
  str_vec idx_names;
  uint_vec idx_types_int;
  typ_vec idx_types;
  uint_vec idx_levels;
  len_vec idx_lens;
  vector<T> d;
  std::string idx_name_pref = "idx_name_";
  fR["rank"] >> rank;
  fR["size"] >> size;
  fR["idx_sizes"] >> idx_sizes;
  fR["idx_types"] >> idx_types_int;
  fR["idx_levels"] >> idx_levels;
  for (size_t i = 0; i < rank; i++) {
    std::vector<char> vec;
    fR[idx_name_pref+std::to_string(i)] >> vec;
    std::string s = std::string(vec.begin(),vec.end());
    idx_names.push_back(s);
  }
  idx_set.clear();
  for (size_t i = 0; i < rank; i++) {
    idx_types.push_back(index_type(idx_types_int[i]));
    idx_set.push_back(dtensor_index(idx_sizes[i],idx_names[i],idx_types[i],idx_levels[i]));
    //idx_lens.push_back(tblis::len_type(idx_sizes[i]));
  }
  //_T.reset(idx_lens);
  fR["T"] >> d;
  //std::copy(d.begin(),d.end(),_T.data());
  _initted = true;
}
template void dtensor<double>::load(ezh5::Node& fR);
template void dtensor< std::complex<double> >::load(ezh5::Node& fR);
//-----------------------------------------------------------------------------



template <typename T>
void dtensor<T>::primeSite(int inc){
    assert(1==2);
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeSite(inc);
  }
}
template void dtensor<double>::primeSite(int inc);
template void dtensor< std::complex<double> >::primeSite(int inc);

template <typename T>
void dtensor<T>::mapPrime(unsigned from, unsigned to){
    assert(1==2);
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].mapPrime(from, to);
  }
}
template void dtensor<double>::mapPrime(unsigned from, unsigned to);
template void dtensor< std::complex<double> >::mapPrime(unsigned from, unsigned to);

template <typename T>
void dtensor<T>::mapPrime(unsigned from, unsigned to, index_type type){
    assert(1==2);
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].mapPrime(from, to, type);
  }
}
template void dtensor<double>::mapPrime(unsigned from, unsigned to, index_type type);
template void dtensor< std::complex<double> >::mapPrime(unsigned from, unsigned to, index_type type);




template <typename T>
void qr(dtensor<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& Q, dtensor<T>& R)
{
  assert(1==2);
  // Permute dtensor
  unsigned r=1, c=1;
  vector<dtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
    r *= v.size();
  }
  for(auto v : right){
    new_idx_set.push_back(v);
    c *= v.size();
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Set up mid index
  dtensor_index mid(std::min(r,c));
  // Set up Q and R
  vector<dtensor_index> Q_idx_set(left);
  Q_idx_set.push_back(mid);
  vector<dtensor_index> R_idx_set;
  R_idx_set.push_back(mid);
  for(auto v : right){
    R_idx_set.push_back(v);
  }
  Q.reset(Q_idx_set);
  R.reset(R_idx_set);
  // perform QR
  //QR(r, c, A._T.data(), Q._T.data(), R._T.data());
}
template void qr(dtensor<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& Q, dtensor<double>& R);
template void qr(dtensor< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& Q, dtensor< std::complex<double> >& R);

/*template <typename T>
void qr(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& Q, dtensor<T>& R)
{
    assert(1==2);
  // Permute dtensor
  unsigned r=1, c=1;
  vector<dtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
    r *= v.size();
  }
  for(auto v : right){
    new_idx_set.push_back(v);
    c *= v.size();
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Set up mid index
  dtensor_index mid(std::min(r,c));
  // Set up Q and R
  vector<dtensor_index> Q_idx_set(left);
  Q_idx_set.push_back(mid);
  vector<dtensor_index> R_idx_set;
  R_idx_set.push_back(mid);
  for(auto v : right){
    R_idx_set.push_back(v);
  }
  Q.reset(Q_idx_set);
  R.reset(R_idx_set);
  // perform QR
  //QR(r, c, A._T.data(), Q._T.data(), R._T.data());
}
template void qr(dtensor_view<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& Q, dtensor<double>& R);
template void qr(dtensor_view< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& Q, dtensor< std::complex<double> >& R);*/


/*template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction)
{
    assert(1==2);
  // Permute dtensor
  unsigned r=1, c=1;
  vector<dtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
    r *= v.size();
  }
  for(auto v : right){
    new_idx_set.push_back(v);
    c *= v.size();
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Set up mid index
  dtensor_index mid(std::min(r,c));
  // Set up U and V
  vector<dtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<dtensor_index> V_idx_set;
  V_idx_set.push_back(mid);
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  // perform SVD
  if(direction==MoveFromLeft){
    U.reset(U_idx_set);
    //SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L');
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    //SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R');
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }
}
template void svd(dtensor_view<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction);
template void svd(dtensor_view< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction);*/




template <typename T>
void svd(dtensor<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff)
{
    assert(1==2);
  // Permute dtensor
  unsigned r=1, c=1;
  vector<dtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
    r *= v.size();
  }
  for(auto v : right){
    new_idx_set.push_back(v);
    c *= v.size();
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Set up mid index
  dtensor_index mid(std::min(r,c));
  // Set up U and V
  vector<dtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<dtensor_index> V_idx_set;
  V_idx_set.push_back(mid);
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  // perform SVD
  if(direction==MoveFromLeft){
    U.reset(U_idx_set);
    //SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L', cutoff);
    if(S.size() != mid.size()){
      uint_vec U_sizes;
      for(auto v : left){
        U_sizes.push_back(v.size());
      }
      U_sizes.push_back(S.size());
      U.resize(U_sizes);
    }
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    //SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R', cutoff);
    if(S.size() != mid.size()){
      uint_vec V_sizes;
      V_sizes.push_back(S.size());
      for(auto v : right){
        V_sizes.push_back(v.size());
      }
      V.resize(V_sizes);
    }
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }
}
template void svd(dtensor<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction, double cutoff);
template void svd(dtensor< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff);

/*template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff)
{
    assert(1==2);
  // Permute dtensor
  unsigned r=1, c=1;
  vector<dtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
    r *= v.size();
  }
  for(auto v : right){
    new_idx_set.push_back(v);
    c *= v.size();
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Set up mid index
  dtensor_index mid(std::min(r,c));
  // Set up U and V
  vector<dtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<dtensor_index> V_idx_set;
  V_idx_set.push_back(mid);
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  // perform SVD
  if(direction==MoveFromLeft){
    U.reset(U_idx_set);
    //SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L', cutoff);
    if(S.size() != mid.size()){
      uint_vec U_sizes;
      for(auto v : left){
        U_sizes.push_back(v.size());
      }
      U_sizes.push_back(S.size());
      U.resize(U_sizes);
    }
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    //SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R', cutoff);
    if(S.size() != mid.size()){
      uint_vec V_sizes;
      V_sizes.push_back(S.size());
      for(auto v : right){
        V_sizes.push_back(v.size());
      }
      V.resize(V_sizes);
    }
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }
}
template void svd(dtensor_view<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction, double cutoff);
template void svd(dtensor_view< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff);*/




/*template <typename T>
void svd(dtensor<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  svd(A,left,right,U,V,S,direction);*/
  //HERE
  /*  assert(1==2);
  // Permute dtensor
  unsigned r=1, c=1;
  vector<dtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
    r *= v.size();
  }
  for(auto v : right){
    new_idx_set.push_back(v);
    c *= v.size();
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Set up mid index
  dtensor_index mid(std::min(r,c));
  // Set up U and V
  vector<dtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<dtensor_index> V_idx_set;
  V_idx_set.push_back(mid);
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  // perform SVD
  if(direction==MoveFromLeft){
    U.reset(U_idx_set);
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L', cutoff);
    if(S.size() > K) S.resize(K);
    if(S.size() != mid.size()){
      uint_vec U_sizes;
      for(auto v : left){
        U_sizes.push_back(v.size());
      }
      U_sizes.push_back( S.size() );
      U.resize(U_sizes);
    }
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R', cutoff);
    if(S.size() > K) S.resize(K);
    if(S.size() != mid.size()){
      uint_vec V_sizes;
      V_sizes.push_back( S.size() );
      for(auto v : right){
        V_sizes.push_back(v.size());
      }
      V.resize(V_sizes);
    }
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }*/


// template <typename T>
// void svd(dtensor_view<T>& A,
//         vector<dtensor_index>& left,
//         vector<dtensor_index>& right,
//         dtensor<T>& U, dtensor<T>& V, vector<double>& S,
//         int direction, double cutoff, long unsigned K)
// {
//     assert(1==2);
//   // Permute dtensor
//   unsigned r=1, c=1;
//   vector<dtensor_index> new_idx_set;
//   for(auto v : left){
//     new_idx_set.push_back(v);
//     r *= v.size();
//   }
//   for(auto v : right){
//     new_idx_set.push_back(v);
//     c *= v.size();
//   }
//   uint_vec perm;
//   find_index_permutation(A.idx_set, new_idx_set, perm);
//   A.permute(perm);
//   // Set up mid index
//   dtensor_index mid(std::min(r,c));
//   // Set up U and V
//   vector<dtensor_index> U_idx_set(left);
//   U_idx_set.push_back(mid);
//   vector<dtensor_index> V_idx_set;
//   V_idx_set.push_back(mid);
//   for(auto v : right){
//     V_idx_set.push_back(v);
//   }
//   // perform SVD
//   if(direction==MoveFromLeft){
//     U.reset(U_idx_set);
//     SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L', cutoff);
//     if(S.size() > K) S.resize(K);
//     if(S.size() != mid.size()){
//       uint_vec U_sizes;
//       for(auto v : left){
//         U_sizes.push_back(v.size());
//       }
//       U_sizes.push_back( S.size() );
//       U.resize(U_sizes);
//     }
//     V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
//   }
//   else if(direction==MoveFromRight){
//     V.reset(V_idx_set);
//     SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R', cutoff);
//     if(S.size() > K) S.resize(K);
//     if(S.size() != mid.size()){
//       uint_vec V_sizes;
//       V_sizes.push_back( S.size() );
//       for(auto v : right){
//         V_sizes.push_back(v.size());
//       }
//       V.resize(V_sizes);
//     }
//     U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
//   }
// }


/*template <typename T>
void svd_bond(dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction)
{
    assert(1==2);
  dtensor<T> U,V;
  dtensor<T> combined = std::move(A_left * A_right);
  vector<dtensor_index> left;
  vector<dtensor_index> right;
  string tag = mid.tag();
  // Separate dtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != tag) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_left.rank; j++) {
    string idx_tag = A_left.idx_set[j].tag();
    if (idx_tag != tag) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction);
  mid.resize(S.size());
  A_right = V;
  A_right.idx_set[0] = mid;
  A_left = U;
  A_left.idx_set.back() = mid;
}
template void svd_bond(dtensor<double>& A_left, dtensor<double>& A_right, dtensor_index& mid, vector<double>& S, int direction);
template void svd_bond(dtensor< std::complex<double> >& A_left, dtensor< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction);*/


/*template <typename T>
void svd_bond(dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction)
{
    assert(1==2);
  dtensor<T> U,V;
  dtensor<T> combined = std::move(A_left * A_right);
  vector<dtensor_index> left;
  vector<dtensor_index> right;
  string tag = mid.tag();
  // Separate dtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != tag) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_left.rank; j++) {
    string idx_tag = A_left.idx_set[j].tag();
    if (idx_tag != tag) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction);
  mid.resize(S.size());
  A_right = V;
  A_right.idx_set[0] = mid;
  A_left = U;
  A_left.idx_set.back() = mid;
}
template void svd_bond(dtensor_view<double>& A_left, dtensor_view<double>& A_right, dtensor_index& mid, vector<double>& S, int direction);
template void svd_bond(dtensor_view< std::complex<double> >& A_left, dtensor_view< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction);*/



template <typename T>
void svd_bond(dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
    assert(1==2);

  dtensor<T> U,V;
  dtensor<T> combined = std::move(A_left * A_right);
  vector<dtensor_index> left;
  vector<dtensor_index> right;
  string tag = mid.tag();
  // Separate dtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != tag) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_left.rank; j++) {
    string idx_tag = A_left.idx_set[j].tag();
    if (idx_tag != tag) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction,cutoff,K);
  mid.resize(S.size());
  A_right = V;
  A_right.idx_set[0] = mid;
  A_left = U;
  A_left.idx_set.back() = mid;
}
template void svd_bond(dtensor<double>& A_left, dtensor<double>& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(dtensor< std::complex<double> >& A_left, dtensor< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);


/*template <typename T>
void svd_bond(dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  assert(1==2);
  dtensor<T> U,V;
  dtensor<T> combined = std::move(A_left * A_right);
  vector<dtensor_index> left;
  vector<dtensor_index> right;
  string tag = mid.tag();
  // Separate dtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != tag) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_left.rank; j++) {
    string idx_tag = A_left.idx_set[j].tag();
    if (idx_tag != tag) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction,cutoff,K);
  mid.resize(S.size());
  A_right = V;
  A_right.idx_set[0] = mid;
  A_left = U;
  A_left.idx_set.back() = mid;
}
template void svd_bond(dtensor_view<double>& A_left, dtensor_view<double>& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(dtensor_view< std::complex<double> >& A_left, dtensor_view< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);*/


/*template <typename T>
void svd_bond(dtensor<T>& combined, dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  assert(1==2);
  dtensor<T> U,V;
  vector<dtensor_index> left;
  vector<dtensor_index> right;
  string tag = mid.tag();
  // Separate dtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != tag) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_left.rank; j++) {
    string idx_tag = A_left.idx_set[j].tag();
    if (idx_tag != tag) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction,cutoff,K);
  mid.resize(S.size());
  A_right = V;
  A_right.idx_set[0] = mid;
  A_left = U;
  A_left.idx_set.back() = mid;
}
template void svd_bond(dtensor<double>& combined, dtensor_view<double>& A_left, dtensor_view<double>& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(dtensor< std::complex<double> >& combined, dtensor_view< std::complex<double> >& A_left, dtensor_view< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);*/

/*template <typename T>
void svd(dtensor<T>& A,
         vector<dtensor_index>& left,
         vector<dtensor_index>& right,
         dtensor<T>& U, dtensor<T>& V, vector<double>& S,
         int direction)
{
  assert(1==2);
  unsigned r=indicesToSize(left);
  unsigned c=indicesToSize(right);
// Set up mid index
  dtensor_index mid(std::min(r,c));
  
  // Set up U and V  
  vector<dtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<dtensor_index> V_idx_set;
  V_idx_set.push_back(mid);
  for(auto v : right){
    V_idx_set.push_back(v);
  }

  unordered_map<string,char> charMap;
  auto indU = indicesToChar(U_idx_set,charMap);
  auto indV = indicesToChar(V_idx_set,charMap);
  
  //CTF::Tensor<T> _S;
  CTF::Vector<T> _S;
  auto indS = string(1,charMap[mid.tag()]);
  auto indA = indToStr(A.idx_set,charMap);
  A.__T[indA.c_str()].svd(U.__T[indU.c_str()],_S[indS.c_str()],V.__T[indV.c_str()]); //,R+3); //Does this work independent of what's in U and V

  //convert _S into a vector
  int64_t np;
  int64_t * inds;
  T * data;
  _S.get_local_data(&np, &inds, &data);
  free(inds);
  S.resize(_S.len);
  //std::copy(S.begin(),S.end(),data);
  std::fill(S.begin(),S.end(),0);

  if(direction==MoveFromLeft){
      U.reset(U_idx_set,false);
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set,false);
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }
}
template void svd(dtensor<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction);
template void svd(dtensor< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction);*/

template <typename T>
void svd_bond(dtensor<T>& combined, dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  assert(1==2);
  dtensor<T> U,V;
  vector<dtensor_index> left;
  vector<dtensor_index> right;
  string tag = mid.tag();
  // Separate dtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != tag) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_left.rank; j++) {
    string idx_tag = A_left.idx_set[j].tag();
    if (idx_tag != tag) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction,cutoff,K);
  mid.resize(S.size());
  A_right = V;
  A_right.idx_set[0] = mid;
  A_left = U;
  A_left.idx_set.back() = mid;
}
template void svd_bond(dtensor<double>& combined, dtensor<double>& A_left, dtensor<double>& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(dtensor< std::complex<double> >& combined, dtensor< std::complex<double> >& A_left, dtensor< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);

template <typename T>
void svd(dtensor<T>& A,
         vector<dtensor_index>& left,
         vector<dtensor_index>& right,
         dtensor<T>& U, dtensor<T>& V, vector<double>& S,
         int direction,double cutoff, long unsigned K)
{
  assert(1==2);
  unsigned r=indicesToSize(left);
  unsigned c=indicesToSize(right);
  dtensor_index mid(std::min(r,c));
  
  // Set up U and V  
  vector<dtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<dtensor_index> V_idx_set;
  V_idx_set.push_back(mid);
  for(auto v : right){
    V_idx_set.push_back(v);
  }

  unordered_map<string,char> charMap;
  auto indU = indicesToChar(U_idx_set,charMap);
  auto indV = indicesToChar(V_idx_set,charMap);
  
  CTF::Tensor<T> _S;
  auto indS = string(1,charMap[mid.tag()]);
  auto indA = indToStr(A.idx_set,charMap);
  //auto indU = indToStr(U_idx_set,charMap);
  //auto indV = indToStr(V_idx_set,charMap);
  //cout<<indA<<" "<<indU<<" "<<indS<<" "<<indV<<endl;
  //A.__T.print();
  //A.__T[indA.c_str()].svd(_U[indU.c_str()],_S[indS.c_str()],_V[indV.c_str()]); //,R+3);
  //U.__T = CTF::Tensor<T>;  //Is this ok?
  //V.__T = CTF::
  //  A.__T[indA.c_str()].svd(U.__T[indU.c_str()],_S[indS.c_str()],V.__T[indV.c_str()]); //,R+3); //Does this work independent of what's in U and V
  
  A.__T[indA.c_str()].svd(U.__T[indU.c_str()],_S[indS.c_str()],V.__T[indV.c_str()],K); //,R+3); //Does this work independent of what's in U and V
  //V.__T.print();
  //cerr<<"U*U:"<< U.__T.reduce(CTF::OP_SUMSQ) << endl;
  //cerr<<"V*V:"<< V.__T.reduce(CTF::OP_SUMSQ) << endl;
  //convert _S into a vector
  int64_t np;
  int64_t * inds;
  T * data;
  _S.get_local_data(&np, &inds, &data);
  free(inds);
  //cerr<<"S SIZE:"<<np<<endl;
  S.resize(_S.lens[0]);
  //S = std::move(std::vector<double>(np,reinterpret_cast<double*>(data)));
  //std::copy(S.begin(),S.end(),data);
  std::fill(S.begin(),S.end(),0);
  /// ,


  
  if(direction==MoveFromLeft){
    /*for (int j=0;j<U.__T.order;j++)
    cerr<<"Length U: "<<j<<" is  "<<U.__T.lens[j]<<" "<<U.idx_set[j].tag()<<endl;
    cerr<<endl;

    for (int j=0;j<A.__T.order;j++)
      cerr<<"Length: "<<j<<" is  "<<A.__T.lens[j]<<" "<<A.idx_set[j].tag()<<endl;*/

    
    U.reset(U_idx_set,false);
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  //cerr<<"U*U:"<< U.__T.reduce(CTF::OP_SUMSQ) << " "<< U.contract(U) << endl;
  //cerr<<"V*V:"<< V.__T.reduce(CTF::OP_SUMSQ) << " "<< V.contract(V) << endl;
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set,false);
    
    /*for (int j=0;j<A.__T.order;j++)
      cerr<<"Length: "<<j<<" is  "<<A.__T.lens[j]<<" "<<A.idx_set[j].tag()<<endl;

    for (int j=0;j<V.__T.order;j++)
    cerr<<"Length V: "<<j<<" is  "<<V.__T.lens[j]<<" "<<V.idx_set[j].tag()<<endl;
    cerr<<endl;


    cerr<<"Length U: "<<i<<" is  "<<U.__T.lens[j]<<endl;
    cerr<<endl;
    for (int j=0;j<A.order;j++)
      cerr<<"Length: "<<i<<" is  "<<A.__T.lens[j]<<endl;*/

    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
    //cerr<<V.__T.order<<endl;
  }
  //cerr<<"Post SVD: "<<endl;
  //cerr<<"U*U:"<< U.__T.reduce(CTF::OP_SUMSQ) << " "<< U.contract(U) << endl;
  //cerr<<"V*V:"<< V.__T.reduce(CTF::OP_SUMSQ) << " "<< V.contract(V) << endl;
  return; 
  exit(1);
  
  //  cout<<"Currently my tensor is "<<A.__T<<endl;
  /*for (auto a : A.idx_set)
    cerr<<a.tag()<<endl;
  cerr<<endl;
  for (auto l : left)
    cerr<<l.tag()<<endl;
  cerr<<endl;
  for (auto r : right)
    cerr<<r.tag()<<endl;*/
  //Tensor<T> U, S, V;
}


template void svd(dtensor<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd(dtensor< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff, long unsigned K);
