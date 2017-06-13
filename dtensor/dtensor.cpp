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
dtensor<T>::dtensor(int_list idx_sizes){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto s : idx_sizes) {
    idx_set.push_back(tensor_index(s));
    idx_lens.push_back(len_type(s));
    ++rank;
    size *= s;
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_list idx_sizes);
template dtensor< std::complex<double> >::dtensor(int_list idx_sizes);


template <typename T>
dtensor<T>::dtensor(int_vec idx_sizes){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (auto s : idx_sizes) {
    idx_set.push_back(tensor_index(s));
    idx_lens.push_back(len_type(s));
    ++rank;
    size *= s;
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_vec idx_sizes);
template dtensor< std::complex<double> >::dtensor(int_vec idx_sizes);


template <typename T>
dtensor<T>::dtensor(int_list idx_sizes, str_list names){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  int_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(tensor_index(v_sizes[i],v_names[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_list idx_sizes, str_list names);
template dtensor< std::complex<double> >::dtensor(int_list idx_sizes, str_list names);


template <typename T>
dtensor<T>::dtensor(int_vec idx_sizes, str_vec names){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(tensor_index(idx_sizes[i],names[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_vec idx_sizes, str_vec names);
template dtensor< std::complex<double> >::dtensor(int_vec idx_sizes, str_vec names);


template <typename T>
dtensor<T>::dtensor(int_list idx_sizes, str_list names, typ_list types){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  int_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  typ_vec v_types(types.begin(), types.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(tensor_index(v_sizes[i],v_names[i],v_types[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_list idx_sizes, str_list names, typ_list types);
template dtensor< std::complex<double> >::dtensor(int_list idx_sizes, str_list names, typ_list types);


template <typename T>
dtensor<T>::dtensor(int_vec idx_sizes, str_vec names, typ_vec types){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(tensor_index(idx_sizes[i],names[i],types[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_vec idx_sizes, str_vec names, typ_vec types);
template dtensor< std::complex<double> >::dtensor(int_vec idx_sizes, str_vec names, typ_vec types);


template <typename T>
dtensor<T>::dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  int_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  typ_vec v_types(types.begin(), types.end());
  int_vec v_levels(levels.begin(), levels.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(tensor_index(v_sizes[i],v_names[i],v_types[i],v_levels[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels);
template dtensor< std::complex<double> >::dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels);


template <typename T>
dtensor<T>::dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(tensor_index(idx_sizes[i],names[i],types[i],levels[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels);
template dtensor< std::complex<double> >::dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels);


template <typename T>
dtensor<T>::dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels, T* data_array){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  int_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  typ_vec v_types(types.begin(), types.end());
  int_vec v_levels(levels.begin(), levels.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(tensor_index(v_sizes[i],v_names[i],v_types[i],v_levels[i]));
    idx_lens.push_back(v_sizes[i]);
    ++rank;
    size *= v_sizes[i];
  }
  _T.reset(idx_lens);
  std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels, double* data_array);
template dtensor< std::complex<double> >::dtensor(int_list idx_sizes, str_list names, typ_list types, int_list levels, std::complex<double>* data_array);


template <typename T>
dtensor<T>::dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels, T* data_array){
  rank = 0;
  size = 1;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_sizes.size(); i++) {
    idx_set.push_back(tensor_index(idx_sizes[i],names[i],types[i],levels[i]));
    idx_lens.push_back(idx_sizes[i]);
    ++rank;
    size *= idx_sizes[i];
  }
  _T.reset(idx_lens);
  std::copy(data_array,data_array+size,_T.data());
  _initted = true;
}
template dtensor<double>::dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels, double* data_array);
template dtensor< std::complex<double> >::dtensor(int_vec idx_sizes, str_vec names, typ_vec types, int_vec levels, std::complex<double>* data_array);


template <typename T>
dtensor<T>::dtensor(vector<tensor_index>& idx_vec){
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  len_vec idx_lens;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    idx_lens.push_back(idx_vec[i]._size);
    ++rank;
    size *= idx_vec[i]._size;
  }
  _T.reset(idx_lens);
  _initted = true;
}
template dtensor<double>::dtensor(vector<tensor_index>& idx_vec);
template dtensor< std::complex<double> >::dtensor(vector<tensor_index>& idx_vec);


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


//-----------------------------------------------------------------------------
// Set values
template <typename T>
void dtensor<T>::setRandom(){
  assert(_initted);
  for (size_t i = 0; i < size; i++) {
    _T.data()[i] = T(2*(drand48()-0.5));
  }
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


template <typename T>
dtensor<T> dtensor<T>::operator * (dtensor<T>& A){
  assert(_initted && A._initted);
  assert(rank>0 && A.rank>0);
  int_vec new_idx_sizes(rank+A.rank);
  str_vec new_idx_names(rank+A.rank);
  typ_vec new_idx_types(rank+A.rank);
  int_vec new_idx_levels(rank+A.rank);
  lab_vec this_labels;
  lab_vec A_labels;
  lab_vec new_labels(rank+A.rank);
  int_vec to_remove;
  unsigned num_shared_idx = 0;
  char ch='a';
  for (size_t i = 0; i < rank; i++) {
    this_labels.push_back(ch);
    ch += 1;
  }
  for (size_t i = 0; i < A.rank; i++) {
    bool is_shared = false;
    for (size_t j = 0; j < rank; j++) {
      if(A.idx_set[i] == idx_set[j]){
        is_shared = true;
        A_labels.push_back(this_labels[j]);
        ++num_shared_idx;
        to_remove.push_back(j);
        to_remove.push_back(i+rank);
        break;
      }
    }
    if(!is_shared){
      A_labels.push_back(ch);
      ch += 1;
    }
  }
  std::sort(to_remove.begin(),to_remove.end());
  for (size_t i = 0; i < rank; i++) {
    new_idx_sizes[i] = idx_set[i]._size;
    new_idx_names[i] = idx_set[i]._name;
    new_idx_types[i] = idx_set[i]._type;
    new_idx_levels[i] = idx_set[i]._level;
    new_labels[i] = this_labels[i];
  }
  for (size_t i = 0; i < A.rank; i++) {
    new_idx_sizes[i+rank] = A.idx_set[i]._size;
    new_idx_names[i+rank] = A.idx_set[i]._name;
    new_idx_types[i+rank] = A.idx_set[i]._type;
    new_idx_levels[i+rank] = A.idx_set[i]._level;
    new_labels[i+rank] = A_labels[i];
  }
  while(to_remove.size()>0){
    int idx = to_remove.back();
    to_remove.pop_back();
    new_idx_sizes.erase(new_idx_sizes.begin()+idx);
    new_idx_names.erase(new_idx_names.begin()+idx);
    new_idx_types.erase(new_idx_types.begin()+idx);
    new_idx_levels.erase(new_idx_levels.begin()+idx);
    new_labels.erase(new_labels.begin()+idx);
  }
  // does not deal with full contraction
  // because it ends in a scalar rather than
  // a tensor
  assert(rank+A.rank-2*num_shared_idx>0);
  dtensor<T> B(new_idx_sizes,new_idx_names,new_idx_types,new_idx_levels);
  tblis::mult(T(1),_T,this_labels.data(),A._T,A_labels.data(),T(0),B._T,new_labels.data());
  return B;
}
template dtensor<double> dtensor<double>::operator * (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor< std::complex<double> >& A);

template <typename T>
dtensor<T>& dtensor<T>::operator *= (const T c){
  assert(_initted);
  lab_vec labels;
  char ch = 'a';
  for (size_t i = 0; i < size; i++) {
    labels.push_back(ch);
    ++ch;
    tblis::scale(c, _T, labels.data());
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator*=(const double c);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator*=(const std::complex<double> c);

template <typename T>
dtensor<T>& dtensor<T>::operator /= (const T c){
  assert(_initted);
  lab_vec labels;
  char ch = 'a';
  for (size_t i = 0; i < size; i++) {
    labels.push_back(ch);
    ++ch;
    tblis::scale(1.0/c, _T, labels.data());
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator/=(const double c);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator/=(const std::complex<double> c);

template <typename T>
dtensor<T> dtensor<T>::operator*(const T c){
  assert(_initted);
  dtensor A(*this);
  lab_vec labels;
  char ch = 'a';
  for (size_t i = 0; i < size; i++) {
    labels.push_back(ch);
    ++ch;
    tblis::scale(c, A._T, labels.data());
  }
  return A;
}
template dtensor<double> dtensor<double>::operator*(const double c);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator*(const std::complex<double> c);

template <typename T>
dtensor<T> dtensor<T>::operator/(const T c){
  assert(_initted);
  dtensor A(*this);
  lab_vec labels;
  char ch = 'a';
  for (size_t i = 0; i < size; i++) {
    labels.push_back(ch);
    ++ch;
    tblis::scale(1.0/c, A._T, labels.data());
  }
  return A;
}
template dtensor<double> dtensor<double>::operator/(const double c);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator/(const std::complex<double> c);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Contract to scalar
template <typename T>
T dtensor<T>::contract(dtensor<T>& A){
  assert(_initted && A._initted);
  assert(rank>0 && A.rank>0);
  lab_vec this_labels;
  lab_vec A_labels;
  unsigned num_shared_idx = 0;
  char ch='a';
  for (size_t i = 0; i < rank; i++) {
    this_labels.push_back(ch);
    ch += 1;
  }
  for (size_t i = 0; i < A.rank; i++) {
    bool is_shared = false;
    for (size_t j = 0; j < rank; j++) {
      if(A.idx_set[i] == idx_set[j]){
        is_shared = true;
        A_labels.push_back(this_labels[j]);
        ++num_shared_idx;
        break;
      }
    }
    if(!is_shared){
      A_labels.push_back(ch);
      ch += 1;
    }
  }
  // does not deal with full contraction
  // because it ends in a scalar rather than
  // a tensor
  assert(rank+A.rank-2*num_shared_idx==0);
  T res;
  tblis::dot(_T,this_labels.data(),A._T,A_labels.data(),res);
  return res;
}
template double dtensor<double>::contract(dtensor<double>& A);
template std::complex<double> dtensor< std::complex<double> >::contract(dtensor< std::complex<double> >& A);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Get diagonal subtensor
// only possible when some tensor indices com in "pairs",
// meaning same name string but different prime level
template <typename T>
dtensor<T> dtensor<T>::diagonal(){
  assert(_initted && rank>1);
  int_vec new_idx_sizes;
  str_vec new_idx_names;
  typ_vec new_idx_types;
  int_vec new_idx_levels;
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
      new_idx_sizes.push_back(idx_set[i]._size);
      new_idx_names.push_back(idx_set[i]._name);
      new_idx_types.push_back(idx_set[i]._type);
      new_idx_levels.push_back(idx_set[i]._level);
      idx_pairs.push_back(std::make_pair(i,-1));
    }
  }
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = i+1; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j])){
        new_idx_sizes.push_back(idx_set[i]._size);
        new_idx_names.push_back(idx_set[i]._name);
        new_idx_types.push_back(idx_set[i]._type);
        new_idx_levels.push_back(0);
        idx_pairs.push_back(std::make_pair(i,j));
        break;
      }
    }
  }
  dtensor<T> A(new_idx_sizes,new_idx_names,new_idx_types,new_idx_levels);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < A.size; i++) {
    int A_idx[A.rank];
    int this_idx[rank];
    for (size_t j = 0; j < A.rank; j++) {
      A_idx[j] = int(i/A._T.stride(j))%new_idx_sizes[j];
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
  int_vec new_idx_sizes;
  str_vec new_idx_names;
  typ_vec new_idx_types;
  int_vec new_idx_levels;
  vector< std::pair<int,int> > idx_pairs;
  for (size_t i = 0; i < rank; i++) {
    bool is_shared = false;
    for (size_t j = 0; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j],type)){
        is_shared = true;
        break;
      }
    }
    if(!is_shared){
      new_idx_sizes.push_back(idx_set[i]._size);
      new_idx_names.push_back(idx_set[i]._name);
      new_idx_types.push_back(idx_set[i]._type);
      new_idx_levels.push_back(idx_set[i]._level);
      idx_pairs.push_back(std::make_pair(i,-1));
    }
  }
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = i+1; j < rank; j++) {
      if(idx_set[i].similar(idx_set[j],type)){
        new_idx_sizes.push_back(idx_set[i]._size);
        new_idx_names.push_back(idx_set[i]._name);
        new_idx_types.push_back(idx_set[i]._type);
        new_idx_levels.push_back(0);
        idx_pairs.push_back(std::make_pair(i,j));
        break;
      }
    }
  }
  dtensor<T> A(new_idx_sizes,new_idx_names,new_idx_types,new_idx_levels);
  #pragma omp parallel for default(shared)
  for (size_t i = 0; i < A.size; i++) {
    int A_idx[A.rank];
    int this_idx[rank];
    for (size_t j = 0; j < A.rank; j++) {
      A_idx[j] = int(i/A._T.stride(j))%new_idx_sizes[j];
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
  std::cout<<"(2) Tensor's index (size, name, type, prime level)"<<'\n';
  for (int i = 0; i < rank; i++) {
    std::cout<<"                   ("<<idx_set[i]._size<<", "<<idx_set[i]._name<<", ";
    if(idx_set[i]._type==Link){
      std::cout<<"Link"<<", ";
    }else{
      std::cout<<"Site"<<", ";
    }
    std::cout<<idx_set[i]._level<<")"<<'\n';
  }
  if (print_level>0) {
    std::cout<<"(3) Matricized tensor:"<<'\n';
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
#ifdef USE_EZH5
template <typename T>
void dtensor<T>::save(string fn){
  assert(_initted);
  int_vec idx_sizes;
  str_vec idx_names;
  int_vec idx_types;
  int_vec idx_levels;
  vector<T> d(_T.data(),_T.data()+size);
  for (size_t i = 0; i < rank; i++) {
    idx_sizes.push_back(idx_set[i]._size);
    idx_names.push_back(idx_set[i]._name);
    idx_types.push_back(idx_set[i]._type);
    idx_levels.push_back(idx_set[i]._level);
  }
  std::string idx_name_pref = "idx_name_";
  ezh5::File fh5W (fn, H5F_ACC_TRUNC);
  fh5W["rank"] = rank;
  fh5W["size"] = size;
  fh5W["idx_sizes"] = idx_sizes;
  fh5W["idx_types"] = idx_types;
  fh5W["idx_levels"] = idx_levels;
  for (int i = 0; i < rank; i++) {
    std::vector<char> vec(idx_names[i].begin(),idx_names[i].end());
    fh5W[idx_name_pref+std::to_string(i)] = vec;
  }
  fh5W["T"] = d;
}
template void dtensor<double>::save(string fn);
template void dtensor< std::complex<double> >::save(string fn);

template <typename T>
void dtensor<T>::load(string fn){
  int_vec idx_sizes;
  str_vec idx_names;
  int_vec idx_types_int;
  typ_vec idx_types;
  int_vec idx_levels;
  len_vec idx_lens;
  vector<T> d;
  std::string idx_name_pref = "idx_name_";
  ezh5::File fh5R (fn, H5F_ACC_RDONLY);
  fh5R["rank"] >> rank;
  fh5R["size"] >> size;
  fh5R["idx_sizes"] >> idx_sizes;
  fh5R["idx_types"] >> idx_types_int;
  fh5R["idx_levels"] >> idx_levels;
  for (int i = 0; i < rank; i++) {
    std::vector<char> vec;
    fh5R[idx_name_pref+std::to_string(i)] >> vec;
    std::string s = std::string(vec.begin(),vec.end());
    idx_names.push_back(s);
  }
  idx_set.clear();
  for (size_t i = 0; i < rank; i++) {
    idx_types.push_back(index_type(idx_types_int[i]));
    idx_set.push_back(tensor_index(idx_sizes[i],idx_names[i],idx_types[i],idx_levels[i]));
    idx_lens.push_back(len_type(idx_sizes[i]));
  }
  _T.reset(idx_lens);
  fh5R["T"] >> d;
  std::copy(d.begin(),d.end(),_T.data());
  _initted = true;
}
template void dtensor<double>::load(string fn);
template void dtensor< std::complex<double> >::load(string fn);
#endif
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
void dtensor<T>::primeLink(){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeLink();
  }
}
template void dtensor<double>::primeLink();
template void dtensor< std::complex<double> >::primeLink();

template <typename T>
void dtensor<T>::primeSite(){
  for (size_t i = 0; i < rank; i++) {
    idx_set[i].primeSite();
  }
}
template void dtensor<double>::primeSite();
template void dtensor< std::complex<double> >::primeSite();

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
//---------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Norm
template <typename T>
double dtensor<T>::norm(){
  double res = 0.0;
  if(std::is_same< T, std::complex<double> >::value){
    for (size_t i = 0; i < size; i++) {
        res += std::real(_T.data()[i]*std::conj(_T.data()[i]));
    }
  }else{
    for (size_t i = 0; i < size; i++) {
        res += std::real(_T.data()[i]*std::conj(_T.data()[i]));
    }
  }
  return std::sqrt(res);
}
template double dtensor<double>::norm();
template double dtensor< std::complex<double> >::norm();


template <typename T>
void dtensor<T>::normalize(){
  double res = 0.0;
  if(std::is_same< T, std::complex<double> >::value){
    for (size_t i = 0; i < size; i++) {
        res += std::real(_T.data()[i]*std::conj(_T.data()[i]));
    }
  }else{
    for (size_t i = 0; i < size; i++) {
        res += std::real(_T.data()[i]*std::conj(_T.data()[i]));
    }
  }
  res = sqrt(res);
  for (size_t i = 0; i < size; i++) {
      _T.data()[i] /= res;
  }
}
template void dtensor<double>::normalize();
template void dtensor< std::complex<double> >::normalize();
//-----------------------------------------------------------------------------


#endif
