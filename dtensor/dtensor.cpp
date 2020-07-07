/*
 * Copyright 2020 Xiongjie Yu, Ryan Levy, and Bryan K. Clark
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef DENSE_TENSOR_CLASS
#define DENSE_TENSOR_CLASS
#include "dtensor.h"

//-----------------------------------------------------------------------------
// Constructors
template <typename T>
dtensor<T>::dtensor(){
  rank = 0;
  size = 0;
  _initted = false;
}
template dtensor<double>::dtensor();
template dtensor< std::complex<double> >::dtensor();


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto s : idx_sizes) {
    idx_set.push_back(dtensor_index(s));
    idx_lens.push_back(tblis::len_type(s));
    ++rank;
    size *= s;
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto s : idx_sizes) {
    idx_set.push_back(dtensor_index(s));
    idx_lens.push_back(tblis::len_type(s));
    ++rank;
    size *= s;
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names){
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
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names, typ_list types){
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
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names, typ_list types);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names, typ_list types);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i],types[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels){
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
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i],types[i],levels[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels);


template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, T* data_array){
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
  _T.reset(idx_lens,tblis::uninitialized);
  std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, double* data_array);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels, std::complex<double>* data_array);


template <typename T>
dtensor<T>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, T* data_array){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(idx_sizes[i],names[i],types[i],levels[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens,tblis::uninitialized);
  std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, double* data_array);
template dtensor< std::complex<double> >::dtensor(uint_vec& idx_sizes, str_vec& names, typ_vec& types, uint_vec& levels, std::complex<double>* data_array);


template <typename T>
dtensor<T>::dtensor(vector<dtensor_index>& idx_vec){
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    idx_lens.push_back(idx_vec[i].size());
    ++rank;
    size *= idx_vec[i].size();
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(vector<dtensor_index>& idx_vec);
template dtensor< std::complex<double> >::dtensor(vector<dtensor_index>& idx_vec);


template <typename T>
dtensor<T>::dtensor(initializer_list<dtensor_index> idx_list){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto i : idx_list){
    idx_set.push_back(i);
    idx_lens.push_back(i.size());
    ++rank;
    size *= i.size();
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(initializer_list<dtensor_index> idx_list);
template dtensor< std::complex<double> >::dtensor(initializer_list<dtensor_index> idx_list);


template <typename T>
dtensor<T>::dtensor(vector<dtensor_index>& idx_vec, T* data_array){
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    idx_lens.push_back(idx_vec[i].size());
    ++rank;
    size *= idx_vec[i].size();
  }
  _T.reset(idx_lens,tblis::uninitialized);
  std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(vector<dtensor_index>& idx_vec, double* data_array);
template dtensor< std::complex<double> >::dtensor(vector<dtensor_index>& idx_vec, std::complex<double>* data_array);


template <typename T>
dtensor<T>::dtensor(const dtensor<T>& other){
  rank = other.rank;
  size = other.size;
  idx_set = other.idx_set;
  _T.reset(other._T);
  _initted = other._initted;
}
template dtensor<double>::dtensor(const dtensor<double>& other);
template dtensor< std::complex<double> >::dtensor(const dtensor< std::complex<double> >& other);


template <typename T>
dtensor<T>::dtensor(const dtensor_view<T>& other){
  rank = other.rank;
  size = other.size;
  idx_set = other.idx_set;
  _T.reset(other._T);
  _initted = other._initted;
}
template dtensor<double>::dtensor(const dtensor_view<double>& other);
template dtensor< std::complex<double> >::dtensor(const dtensor_view< std::complex<double> >& other);


template <typename T>
dtensor<T>::dtensor(dtensor<T>&& other){
  rank = other.rank;
  size = other.size;
  idx_set = std::move(other.idx_set);
  _T.reset(std::move(other._T));
  _initted = other._initted;
}
template dtensor<double>::dtensor(dtensor<double>&& other);
template dtensor< std::complex<double> >::dtensor(dtensor< std::complex<double> >&& other);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Reset
template <typename T>
void dtensor<T>::reset(vector<dtensor_index>& idx_vec){
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    idx_lens.push_back(idx_vec[i].size());
    ++rank;
    size *= idx_vec[i].size();
  }
  _T.reset(idx_lens);
  _initted = true;
  setZero();
}
template void dtensor<double>::reset(vector<dtensor_index>& idx_vec);
template void dtensor< std::complex<double> >::reset(vector<dtensor_index>& idx_vec);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Resize indices
// (data preserved when dimension of indices lowered, filled with val when enlarged)
template <typename T>
void dtensor<T>::resize(uint_vec& new_sizes, T val){
  assert(_initted);
  assert(new_sizes.size()==rank);
  std::vector<long int> v;
  size = 1;
  for (size_t i = 0; i < rank; i++) {
    if(idx_set[i].size()!=new_sizes[i]) idx_set[i].resize(new_sizes[i]);
    v.push_back(new_sizes[i]);
    size *= new_sizes[i];
  }
  _T.resize(v, val);
}
template void dtensor<double>::resize(uint_vec& new_sizes, double val);
template void dtensor< std::complex<double> >::resize(uint_vec& new_sizes, std::complex<double> val);

template <typename T>
void dtensor<T>::resize(uint_list new_sizes, T val){
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
  assert(_initted);
  assert(new_idx_set.size()==rank);
  idx_set = new_idx_set;
  std::vector<long int> v;
  size = 1;
  for (size_t i = 0; i < rank; i++) {
    v.push_back(idx_set[i].size());
    size *= idx_set[i].size();
  }
  _T.resize(v, val);
}
template void dtensor<double>::resize(vector<dtensor_index>& new_idx_set, double val);
template void dtensor< std::complex<double> >::resize(vector<dtensor_index>& new_idx_set, std::complex<double> val);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Set values
template <typename T>
void dtensor<T>::setRandom(){
  assert(_initted);
  random_array(_T.data(), size);
}
template void dtensor<double>::setRandom();
template void dtensor< std::complex<double> >::setRandom();

template <typename T>
void dtensor<T>::setZero(){
  assert(_initted);
  for (size_t i = 0; i < size; i++) {
    _T.data()[i] = T(0.0);
  }
}
template void dtensor<double>::setZero();
template void dtensor< std::complex<double> >::setZero();

template <typename T>
void dtensor<T>::setOne(){
  assert(_initted);
  for (size_t i = 0; i < size; i++) {
    _T.data()[i] = T(1);
  }
}
template void dtensor<double>::setOne();
template void dtensor< std::complex<double> >::setOne();
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Permute
#if !defined(USE_HPTT)
template <typename T>
void dtensor<T>::permute(uint_vec& perm)
{
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
    len_vec idx_lens(rank);
    for (size_t i = 0; i < rank; i++) {
      idx_set[i]  = idx_set_p[perm[i]];
      idx_lens[i] = idx_set[i].size();
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
    std::copy(_T.data(), _T.data()+size, A);
    _T.reset(idx_lens);
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
      _T.data()[idx] = A[i];
    }
    delete [] A;
  }
}
#else
template <typename T>
void dtensor<T>::permute(uint_vec& perm){
  assert(_initted);
  bool perm_needed = false;
  for (size_t i = 0; i < perm.size(); i++) {
    if(i!=perm[i]){
      perm_needed = true;
      break;
    }
  }
  if (perm_needed){
    //std::cerr<<"Need perm!"<<std::endl;
    vector<dtensor_index> idx_set_p(idx_set);
    T* A = new T [size];
    std::copy(_T.data(), _T.data()+size, A);
    T alpha = 1;
    T beta  = 0;
    len_vec idx_lens(rank);
    int idx_sizes[rank];
    for (size_t i = 0; i < rank; i++) {
      idx_sizes[i] = idx_set[i].size();
      idx_set[i]   = idx_set_p[perm[i]];
      idx_lens[i]  = idx_set[i].size();
    }
    _T.reset(idx_lens,tblis::uninitialized);
    T* B = _T.data(); //ensure B points to new data
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
  uint_vec perm_vec;
  for(auto s : perm){
    perm_vec.push_back(s);
  }
  permute(perm_vec);
}
template void dtensor<double>::permute(uint_list perm);
template void dtensor< std::complex<double> >::permute(uint_list perm);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Operator overloading
template <typename T>
dtensor<T>& dtensor<T>::operator=(const dtensor<T>& other){
  if(this!=&other){
    rank = other.rank;
    size = other.size;
    idx_set = other.idx_set;
    _T.reset(other._T);
    _initted = other._initted;
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator=(const dtensor<double> &other);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator=(const dtensor< std::complex<double> > &other);


template <typename T>
dtensor<T>& dtensor<T>::operator=(const dtensor_view<T>& other){
  rank = other.rank;
  size = other.size;
  idx_set = other.idx_set;
  _T.reset(other._T);
  _initted = other._initted;
  return *this;
}
template dtensor<double>& dtensor<double>::operator=(const dtensor_view<double> &other);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator=(const dtensor_view< std::complex<double> > &other);


template <typename T>
dtensor<T>& dtensor<T>::operator=(dtensor<T>&& other){
  if(this!=&other){
    rank = other.rank;
    size = other.size;
    idx_set = std::move(other.idx_set);
    _T.reset(std::move(other._T));
    _initted = other._initted;
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator=(dtensor<double>&& other);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator=(dtensor< std::complex<double> >&& other);


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


template <typename T>
dtensor<T> dtensor<T>::operator * (dtensor<T>& A){
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
  tblis::mult<T>(T(1),_T,this_labels.data(),A._T,A_labels.data(),T(0),res._T,res_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator * (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor< std::complex<double> >& A);


template <typename T>
dtensor<T> dtensor<T>::operator * (dtensor_view<T>& A){
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
  tblis::mult<T>(T(1),_T,this_labels.data(),A._T,A_labels.data(),T(0),res._T,res_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator * (dtensor_view<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor_view< std::complex<double> >& A);


template <typename T>
dtensor<T> dtensor<T>::operator + (dtensor<T>& A){
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
  tblis::add<T>(T(1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator + (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator + (dtensor< std::complex<double> >& A);


template <typename T>
dtensor<T> dtensor<T>::operator + (dtensor_view<T>& A){
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
  tblis::add<T>(T(1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator + (dtensor_view<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator + (dtensor_view< std::complex<double> >& A);


template <typename T>
dtensor<T> dtensor<T>::operator - (dtensor<T>& A){
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
  tblis::add<T>(T(-1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator - (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator - (dtensor< std::complex<double> >& A);


template <typename T>
dtensor<T> dtensor<T>::operator - (dtensor_view<T>& A){
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
  tblis::add<T>(T(-1),A._T,A_labels.data(),T(1),res._T,this_labels.data());
  return res;
}
template dtensor<double> dtensor<double>::operator - (dtensor_view<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator - (dtensor_view< std::complex<double> >& A);


template <typename T>
dtensor<T>& dtensor<T>::operator += (dtensor<T>& A){
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
  tblis::add<T>(T(1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator += (dtensor<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator += (dtensor< std::complex<double> >& A);


template <typename T>
dtensor<T>& dtensor<T>::operator += (dtensor_view<T>& A){
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
  tblis::add<T>(T(1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator += (dtensor_view<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator += (dtensor_view< std::complex<double> >& A);


template <typename T>
dtensor<T>& dtensor<T>::operator -= (dtensor<T>& A){
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
  tblis::add<T>(T(-1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator -= (dtensor<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator -= (dtensor< std::complex<double> >& A);


template <typename T>
dtensor<T>& dtensor<T>::operator -= (dtensor_view<T>& A){
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
  tblis::add<T>(T(-1),A._T,A_labels.data(),T(1),_T,this_labels.data());
  return *this;
}
template dtensor<double>& dtensor<double>::operator -= (dtensor_view<double>& A);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator -= (dtensor_view< std::complex<double> >& A);


template <typename T>
dtensor<T>& dtensor<T>::operator *= (const T c){
  assert(_initted);
  for (size_t i = 0; i < size; i++) {
    _T.data()[i] *= c;
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator*=(const double c);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator*=(const std::complex<double> c);

template <typename T>
dtensor<T>& dtensor<T>::operator /= (const T c){
  assert(_initted);
  for (size_t i = 0; i < size; i++) {
    _T.data()[i] /= c;
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator/=(const double c);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator/=(const std::complex<double> c);

template <typename T>
dtensor<T> dtensor<T>::operator*(const T c){
  assert(_initted);
  dtensor A(*this);
  for (size_t i = 0; i < A.size; i++) {
    A._T.data()[i] *= c;
  }
  return A;
}
template dtensor<double> dtensor<double>::operator*(const double c);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator*(const std::complex<double> c);

template <typename T>
dtensor<T> dtensor<T>::operator/(const T c){
  assert(_initted);
  dtensor A(*this);
  for (size_t i = 0; i < A.size; i++) {
    A._T.data()[i] /= c;
  }
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

template <typename T>
T dtensor<T>::contract(dtensor<T>& A){
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
  tblis::dot<T>(_T,this_labels.data(),A._T,A_labels.data(),res);
  return res;
}
template double dtensor<double>::contract(dtensor<double>& A);
template std::complex<double> dtensor< std::complex<double> >::contract(dtensor< std::complex<double> >& A);


template <typename T>
T dtensor<T>::contract(dtensor_view<T>& A){
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
  tblis::dot<T>(_T,this_labels.data(),A._T,A_labels.data(),res);
  return res;
}
template double dtensor<double>::contract(dtensor_view<double>& A);
template std::complex<double> dtensor< std::complex<double> >::contract(dtensor_view< std::complex<double> >& A);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Get diagonal subtensor
// only possible when some tensor indices com in "pairs",
// meaning same name but different prime level
template <typename T>
dtensor<T> dtensor<T>::diagonal(){
  assert(_initted && rank>1);
  vector<dtensor_index> new_idx_set;
  vector< std::pair<int,int> > idx_pairs;
  for (size_t i = 0; i < rank; i++) {
    bool is_shared = false;
    for (size_t j = 0; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j])){
        is_shared = true;
        break;
      }
    }
    if(!is_shared){
      new_idx_set.push_back(idx_set[i]);
      idx_pairs.push_back(std::make_pair(i,-1));
    }
  }
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = i+1; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j])){
        if(idx_set[i].level() < idx_set[j].level()){
          new_idx_set.push_back(idx_set[i]);
        }else{
          new_idx_set.push_back(idx_set[j]);
        }
        idx_pairs.push_back(std::make_pair(i,j));
        break;
      }
    }
  }
  if(new_idx_set.size()==rank){
    dtensor<T> A(*this);
    return A;
  }
  dtensor<T> A(new_idx_set);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < A.size; i++) {
    int A_idx[A.rank];
    int this_idx[rank];
    for (size_t j = 0; j < A.rank; j++) {
      A_idx[j] = int(i/A._T.stride(j))%new_idx_set[j].size();
      int i1 = idx_pairs[j].first;
      int i2 = idx_pairs[j].second;
      if(i1>=0) this_idx[i1] = A_idx[j];
      if(i2>=0) this_idx[i2] = A_idx[j];
    }
    int idx = 0;
    for (size_t j = 0; j < rank; j++) {
      idx += this_idx[j] * _T.stride(j);
    }
    A._T.data()[i] = _T.data()[idx];
  }
  return A;
}
template dtensor<double> dtensor<double>::diagonal();
template dtensor< std::complex<double> > dtensor< std::complex<double> >::diagonal();

template <typename T>
dtensor<T> dtensor<T>::diagonal(index_type type){
  assert(_initted && rank>1);
  vector<dtensor_index> new_idx_set;
  vector< std::pair<int,int> > idx_pairs;
  for (size_t i = 0; i < rank; i++) {
    bool is_shared = false;
    for (size_t j = 0; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j], type)){
        is_shared = true;
        break;
      }
    }
    if(!is_shared){
      new_idx_set.push_back(idx_set[i]);
      idx_pairs.push_back(std::make_pair(i,-1));
    }
  }
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = i+1; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j], type)){
        if(idx_set[i].level() < idx_set[j].level()){
          new_idx_set.push_back(idx_set[i]);
        }else{
          new_idx_set.push_back(idx_set[j]);
        }
        idx_pairs.push_back(std::make_pair(i,j));
        break;
      }
    }
  }
  if(new_idx_set.size()==rank){
    dtensor<T> A(*this);
    return A;
  }
  dtensor<T> A(new_idx_set);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < A.size; i++) {
    int A_idx[A.rank];
    int this_idx[rank];
    for (size_t j = 0; j < A.rank; j++) {
      A_idx[j] = int(i/A._T.stride(j))%new_idx_set[j].size();
      int i1 = idx_pairs[j].first;
      int i2 = idx_pairs[j].second;
      if(i1>=0) this_idx[i1] = A_idx[j];
      if(i2>=0) this_idx[i2] = A_idx[j];
    }
    int idx = 0;
    for (size_t j = 0; j < rank; j++) {
      idx += this_idx[j] * _T.stride(j);
    }
    A._T.data()[i] = _T.data()[idx];
  }
  return A;
}
template dtensor<double> dtensor<double>::diagonal(index_type type);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::diagonal(index_type type);
//---------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Print Info
template <typename T>
void dtensor<T>::print(unsigned print_level){
  std::cout<<"-------------------------------------"<<'\n';
  std::cout<<"(1) Tensor's rank = "<<rank<<'\n';
  std::cout<<"(2) Tensor's size = "<<size<<'\n';
  std::cout<<"(3) Tensor's index (size, name, type, prime level)"<<'\n';
  for (size_t i = 0; i < rank; i++) {
    std::cout<<"                   ("<<idx_set[i].size()<<", "<<idx_set[i].name()<<", ";
    if(idx_set[i].type()==Link){
      std::cout<<"Link"<<", ";
    }else{
      std::cout<<"Site"<<", ";
    }
    std::cout<<idx_set[i].level()<<")"<<'\n';
  }
  if (print_level>0) {
    std::cout<<"(4) Tensor data:"<<'\n';
    for (size_t i = 0; i < size; i++) {
      std::cout<<_T.data()[i]<<" ";
    }
    std::cout<<'\n';
  }
  std::cout << "-------------------------------------" << '\n';
}
template void dtensor<double>::print(unsigned print_level);
template void dtensor< std::complex<double> >::print(unsigned print_level);


//-----------------------------------------------------------------------------
// Save/Load
template <typename T>
void dtensor<T>::save(string fn){
  assert(_initted);
  uint_vec idx_sizes;
  str_vec idx_names;
  uint_vec idx_types;
  uint_vec idx_levels;
  vector<T> d(_T.data(),_T.data()+size);
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
  assert(_initted);
  uint_vec idx_sizes;
  str_vec idx_names;
  uint_vec idx_types;
  uint_vec idx_levels;
  vector<T> d(_T.data(),_T.data()+size);
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
    idx_lens.push_back(tblis::len_type(idx_sizes[i]));
  }
  _T.reset(idx_lens);
  fh5R["T"] >> d;
  std::copy(d.begin(),d.end(),_T.data());
  _initted = true;
}
template void dtensor<double>::load(string fn);
template void dtensor< std::complex<double> >::load(string fn);


template <typename T>
void dtensor<T>::load(ezh5::Node& fR){
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
    idx_lens.push_back(tblis::len_type(idx_sizes[i]));
  }
  _T.reset(idx_lens);
  fR["T"] >> d;
  std::copy(d.begin(),d.end(),_T.data());
  _initted = true;
}
template void dtensor<double>::load(ezh5::Node& fR);
template void dtensor< std::complex<double> >::load(ezh5::Node& fR);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Prime level manipulation
template <typename T>
void dtensor<T>::prime(int inc){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].prime(inc);
  }
}
template void dtensor<double>::prime(int inc);
template void dtensor< std::complex<double> >::prime(int inc);

template <typename T>
void dtensor<T>::primeLink(int inc){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeLink(inc);
  }
}
template void dtensor<double>::primeLink(int inc);
template void dtensor< std::complex<double> >::primeLink(int inc);

template <typename T>
void dtensor<T>::primeSite(int inc){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeSite(inc);
  }
}
template void dtensor<double>::primeSite(int inc);
template void dtensor< std::complex<double> >::primeSite(int inc);

template <typename T>
void dtensor<T>::mapPrime(unsigned from, unsigned to){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].mapPrime(from, to);
  }
}
template void dtensor<double>::mapPrime(unsigned from, unsigned to);
template void dtensor< std::complex<double> >::mapPrime(unsigned from, unsigned to);

template <typename T>
void dtensor<T>::mapPrime(unsigned from, unsigned to, index_type type){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].mapPrime(from, to, type);
  }
}
template void dtensor<double>::mapPrime(unsigned from, unsigned to, index_type type);
template void dtensor< std::complex<double> >::mapPrime(unsigned from, unsigned to, index_type type);

template <typename T>
void dtensor<T>::mapPrime(dtensor_index& in, unsigned from, unsigned to){
  for (size_t i = 0; i < rank; i++) {
    if(idx_set[i] == in){
      idx_set[i].mapPrime(from, to);
      break;
    }
  }
}
template void dtensor<double>::mapPrime(dtensor_index& in, unsigned from, unsigned to);
template void dtensor< std::complex<double> >::mapPrime(dtensor_index& in, unsigned from, unsigned to);

template <typename T>
void dtensor<T>::mapPrime(std::vector<dtensor_index>& vec_in, unsigned from, unsigned to){
  for(auto in : vec_in){
    for (size_t i = 0; i < rank; i++) {
      if(idx_set[i] == in){
        idx_set[i].mapPrime(from, to);
        break;
      }
    }
  }
}
template void dtensor<double>::mapPrime(std::vector<dtensor_index>& vec_in, unsigned from, unsigned to);
template void dtensor< std::complex<double> >::mapPrime(std::vector<dtensor_index>& vec_in, unsigned from, unsigned to);

template <typename T>
void dtensor<T>::noPrime(index_type type){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].noPrime(type);
  }
}
template void dtensor<double>::noPrime(index_type type);
template void dtensor< std::complex<double> >::noPrime(index_type type);

template <typename T>
void dtensor<T>::conj(){
  if (std::is_same<T, std::complex<double>>::value) {
    for (size_t i = 0; i < size; i++) _T.data()[i] = cconj(_T.data()[i]);
  }
}
template void dtensor<double>::conj();
template void dtensor< std::complex<double> >::conj();
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
// special arithmetic operations with another tensor in the same format/pattern
template <typename T>
void dtensor<T>::add(dtensor<T>& A, T c){
  assert(_initted && A._initted);
  assert(rank==A.rank);
  for (size_t i = 0; i < size; i++) _T.data()[i] += c*A._T.data()[i];
}
template void dtensor<double>::add(dtensor<double>& A, double c);
template void dtensor< std::complex<double> >::add(dtensor< std::complex<double> >& A, std::complex<double> c);


template <typename T>
T dtensor<T>::inner_product(dtensor<T>& A){
  assert(_initted && A._initted);
  assert(rank==A.rank);
  T res = 0;
  for (size_t i = 0; i < size; i++) res += cconj(_T.data()[i])*A._T.data()[i];
  return res;
}
template double dtensor<double>::inner_product(dtensor<double>& A);
template std::complex<double> dtensor< std::complex<double> >::inner_product(dtensor< std::complex<double> >& A);

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Norm
template <typename T>
double dtensor<T>::norm(){
  double res = 0.0;
  for (size_t i = 0; i < size; i++) {
    res += std::real(_T.data()[i]*std::conj(_T.data()[i]));
  }
  return std::sqrt(res);
}
template double dtensor<double>::norm();
template double dtensor< std::complex<double> >::norm();


template <typename T>
double dtensor<T>::normalize(){
  double res = 0.0;
  for (size_t i = 0; i < size; i++) {
    res += std::real(_T.data()[i]*std::conj(_T.data()[i]));
  }
  res = sqrt(res);
  for (size_t i = 0; i < size; i++) {
      _T.data()[i] /= res;
  }
  return res;
}
template double dtensor<double>::normalize();
template double dtensor< std::complex<double> >::normalize();
//-----------------------------------------------------------------------------


#endif
