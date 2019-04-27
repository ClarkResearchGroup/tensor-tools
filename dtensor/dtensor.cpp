#ifndef DENSE_TENSOR_CLASS
#define DENSE_TENSOR_CLASS
#include "dtensor.h"

string indToStr(vector<dtensor_index> &indices,unordered_map<string,char> &charMap)
{
  assert(charMap.size() >= indices.size());
  string myString="";
  for (auto i: indices){
    const string thisTag = i.tag();
    myString+=charMap[thisTag];
  }
  return myString;
}
string indicesToChar(vector<dtensor_index> &indices, unordered_map<string,char> &charMap)
{
  char ch='a';
  for (auto i : charMap) //find highest char
    if (i.second > ch)
      ch=i.second;
  if(ch!='a') ++ch; //increment from latest
  for (auto i : indices){
    auto it= charMap.find(i.tag());
    if (it==charMap.end()){ //new tag, add it to map
      charMap[i.tag()]=ch;
      ++ch;
    }
  }
  return indToStr(indices,charMap);
}

unsigned indicesToSize(vector<dtensor_index> &indices)
{
  unsigned r=1;
  for (auto v : indices)
    r*=v.size();
  return r;
}


vector<char> indToVec(vector<dtensor_index> &indices,unordered_map<string,char> &charMap)
{
  vector<char> myVec;
  for (auto i : indices){
    const string thisTag=i.tag();
    myVec.push_back(charMap[thisTag]);
  }
  return myVec;

}

string indToStrNP(vector<dtensor_index> &indices,unordered_map<string,char> &charMap)
{
  string myString="";
  for (auto i: indices){
    const string thisTag = noPrime(i).tag();
    //cerr<<thisTag<<endl;
    myString+=charMap[thisTag];
  }
  return myString;
}
string indicesToCharNP(vector<dtensor_index> &indices, unordered_map<string,char> &charMap)
{
  char ch='a';
  for (auto i : charMap) //find highest char
    if (i.second > ch)
      ch=i.second;
  if(ch!='a') ++ch; //increment from latest
  for (auto i : indices){
    auto it= charMap.find(noPrime(i).tag());
    if (it==charMap.end()){ //new tag, add it to map
      charMap[noPrime(i).tag()]=ch;
      ++ch;
    }
  }
  return indToStrNP(indices,charMap);
}

//-----------------------------------------------------------------------------
// Constructors
template <typename T>
dtensor<T>::dtensor(){
  rank = 0;
  size = 0;
  _initted = false;
  //__T=CTF::Tensor<>();
}
template dtensor<double>::dtensor();
template dtensor< std::complex<double> >::dtensor();




template <typename T>
dtensor<T>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels){
  rank = 0;
  size = 1;

  uint_vec v_sizes(idx_sizes.begin(), idx_sizes.end());
  str_vec v_names(names.begin(), names.end());
  typ_vec v_types(types.begin(), types.end());
  uint_vec v_levels(levels.begin(), levels.end());
  for (size_t i = 0; i < v_sizes.size(); i++) {
    idx_set.push_back(dtensor_index(v_sizes[i],v_names[i],v_types[i],v_levels[i]));
    ++rank;
    size *= v_sizes[i];
  }

  vector<int> idx_sizes_int(begin(idx_sizes),end(idx_sizes));
  __T=CTF::Tensor<>(idx_sizes.size(),idx_sizes_int.data());
  _initted = true;
}
template dtensor<double>::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels);
template dtensor< std::complex<double> >::dtensor(uint_list idx_sizes, str_list names, typ_list types, uint_list levels);

template <typename T>
dtensor<T>::dtensor(vector<dtensor_index>& idx_vec){
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  vector<int> idx_sizes;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    ++rank;
    size *= idx_vec[i].size();
    idx_sizes.push_back(idx_vec[i].size());
  }
  __T=CTF::Tensor<>(idx_vec.size(),idx_sizes.data());
  _initted = true;
}
template dtensor<double>::dtensor(vector<dtensor_index>& idx_vec);
template dtensor< std::complex<double> >::dtensor(vector<dtensor_index>& idx_vec);


template <typename T>
dtensor<T>::dtensor(initializer_list<dtensor_index> idx_list){
  rank = 0;
  size = 1;
  std::vector<int> idx_lens;
  for (auto i : idx_list){
    idx_set.push_back(i);
    idx_lens.push_back(i.size());
    ++rank;
    size *= i.size();
  }
  __T=CTF::Tensor<>(idx_lens.size(),idx_lens.data());
  _initted = true;
}
template dtensor<double>::dtensor(initializer_list<dtensor_index> idx_list);
template dtensor< std::complex<double> >::dtensor(initializer_list<dtensor_index> idx_list);


//Will those actually work?
template <typename T>
dtensor<T>::dtensor(vector<dtensor_index>& idx_vec, CTF::Tensor<T>& data_array){
  rank = 0;
  size = 1;
  idx_set = idx_vec;
  for (size_t i = 0; i < idx_vec.size(); i++) {
    ++rank;
    size *= idx_vec[i].size();
  }
  __T = std::move(data_array);
  _initted = true;
}
template dtensor<double>::dtensor(vector<dtensor_index>& idx_vec, CTF::Tensor<double>& data_array);
//:template dtensor< std::complex<double> >::dtensor(vector<dtensor_index>& idx_vec, CTF::Tensor<std::complex<double> >& data_array);

template <typename T>
dtensor<T>::dtensor(const dtensor<T>& other){
  rank = other.rank;
  size = other.size;
  idx_set = other.idx_set;
  __T=other.__T;
  _initted = other._initted;
}
template dtensor<double>::dtensor(const dtensor<double>& other);
template dtensor< std::complex<double> >::dtensor(const dtensor< std::complex<double> >& other);



template <typename T>
dtensor<T>::dtensor(dtensor<T>&& other){
  rank = other.rank;
  size = other.size;
  idx_set = std::move(other.idx_set);
  __T = std::move(other.__T);
  _initted = other._initted;
}
template dtensor<double>::dtensor(dtensor<double>&& other);
template dtensor< std::complex<double> >::dtensor(dtensor< std::complex<double> >&& other);
//-----------------------------------------------------------------------------


//---------------------------------------------------------------------------
// Reset

//What's this being used for.  Should we be resetting the CTF tensor sizes??
//What does a tblis reset have to do?
template <typename T>
void dtensor<T>::reset(vector<dtensor_index>& idx_vec, bool makeZero){
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
  _initted = true;
  if(makeZero) setZero();
}
template void dtensor<double>::reset(vector<dtensor_index>& idx_vec, bool makeZero);
template void dtensor< std::complex<double> >::reset(vector<dtensor_index>& idx_vec, bool makeZero);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Set values
template <typename T>
void dtensor<T>::setRandom(){
  assert(_initted);
  __T.fill_random(0,1);
}
template void dtensor<double>::setRandom();
template void dtensor< std::complex<double> >::setRandom();

template <typename T>
void dtensor<T>::setZero(){
  assert(_initted);
  __T.set_zero();
}
template void dtensor<double>::setZero();
template void dtensor< std::complex<double> >::setZero();

//This should be cleaner
template <typename T>
void dtensor<T>::setOne(){
  assert(_initted);
  //build generic set of indices
  string idxs;
  char i = 'a';
  for(auto idx: idx_set) idxs+= i++;
  __T[idxs.c_str()] = 1.0;
}
template void dtensor<double>::setOne();
template void dtensor< std::complex<double> >::setOne();
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Operator overloading
template <typename T>
dtensor<T>& dtensor<T>::operator=(const dtensor<T>& other){
  if(this!=&other){
    rank = other.rank;
    size = other.size;
    idx_set = other.idx_set;
    __T=other.__T;
    _initted = other._initted;
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator=(const dtensor<double> &other);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator=(const dtensor< std::complex<double> > &other);


template <typename T>
dtensor<T>& dtensor<T>::operator=(dtensor<T>&& other){ //this is when you have A = std::move(A*V)
  if(this!=&other){
    rank = other.rank;
    size = other.size;
    idx_set = std::move(other.idx_set);
    _initted = other._initted;
    __T=std::move(other.__T);
  }
  return *this;
}
template dtensor<double>& dtensor<double>::operator=(dtensor<double>&& other);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator=(dtensor< std::complex<double> >&& other);


//This should be cleaned up
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
  // cerr<<"START"<<endl;
  // for (int i=0;i<res.__T.order;i++){
  //   cerr<<res.__T.lens[i]<<" "<<res.idx_set[i].tag()<<endl;
  // }
  // cerr<<endl;
  // for (int i=0;i<__T.order;i++){
  //   cerr<<__T.lens[i]<<" "<<idx_set[i].tag()<<endl;
  // }
  // cerr<<endl;
  // for (int i=0;i<A.__T.order;i++){
  //   cerr<<A.__T.lens[i]<<" "<<A.idx_set[i].tag()<<endl;
  // }
  // cerr<<"END"<<endl;
  // cerr<<endl;

  auto a=  res.__T[res_labels.data()];
   auto b=  __T[this_labels.data()];
   auto c= A.__T[A_labels.data()];
  T one=T(1);
  T zero=T(0);
  c.mult_scl((char*)&one);
  a.mult_scl((char*)&zero);
  //    exit(1);
  //cerr<<"Pre contraction"<<endl;
  a+=b*c;
  //cerr<<"post contraction"<<endl;
  for (int i=0;i<__T.order;i++)
    assert(__T.lens[i]==idx_set[i].size());
   

  return res;
}
template dtensor<double> dtensor<double>::operator * (dtensor<double>& A);
template dtensor< std::complex<double> > dtensor< std::complex<double> >::operator * (dtensor< std::complex<double> >& A);


template <typename T>
dtensor<T>& dtensor<T>::operator *= (const T c){
  assert(_initted);
  unordered_map<string,char> charMap;
  CTF::Scalar<T> cs(c);
  __T[getIndices(charMap).c_str()]*=cs[""];
  return *this;
}
template dtensor<double>& dtensor<double>::operator*=(const double c);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator*=(const std::complex<double> c);

//This should be cleaned up
template <typename T>
dtensor<T>& dtensor<T>::operator /= (const T c){
  assert(_initted);
  string letters="";
  char ch='a';
  for (unsigned i=0;i<rank;i++){
    letters.push_back(ch);
    ++ch;
  }
  CTF::Transform<T>([c](T & d){ d= d/c; })(__T[letters.c_str()]);
  return *this;
}
template dtensor<double>& dtensor<double>::operator/=(const double c);
template dtensor< std::complex<double> >& dtensor< std::complex<double> >::operator/=(const std::complex<double> c);

//This shoud be cleaned up.
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
  res=__T[this_labels.data()]*A.__T[A_labels.data()];
 
  //  tblis::dot(_T,this_labels.data(),A._T,A_labels.data(),res);
  return res;
}
template double dtensor<double>::contract(dtensor<double>& A);
template std::complex<double> dtensor< std::complex<double> >::contract(dtensor< std::complex<double> >& A);

//---------------------------------------------------------------------------


// Get diagonal subtensor
// only possible when some tensor indices com in "pairs",
// meaning same name but different prime level
//This should be cleaned up.
template <typename T>
dtensor<T> dtensor<T>::diagonal(){
  //assert(1==2);
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
  //dtensor<T> A(new_idx_set);
  /*#pragma omp parallel for default(shared)
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
  }*/
  unordered_map<string,char> charMap;

  auto indOrig = indicesToCharNP(idx_set,charMap);
  auto indNew  = indicesToCharNP(new_idx_set,charMap);
  CTF::Tensor<T> tempT;
  tempT[indNew.c_str()] = __T[indOrig.c_str()];

  dtensor<T> A(new_idx_set,tempT);

  return A;
}
template dtensor<double> dtensor<double>::diagonal();
template dtensor< std::complex<double> > dtensor< std::complex<double> >::diagonal();


//This should be cleaned up.
template <typename T>
dtensor<T> dtensor<T>::diagonal(index_type type){
  //assert(1==2);
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
  /*dtensor<T> A(new_idx_set);
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
  }*/
  unordered_map<string,char> charMap;

  dtensor<T> A(new_idx_set);
  auto indNew  = indicesToCharNP(new_idx_set,charMap);
  auto indOrig = indicesToCharNP(idx_set,charMap);
  //cerr<<indNew << " " << indOrig<<endl;
  //print();
  //__T.print();
  A.__T[indNew.c_str()] = __T[indOrig.c_str()];

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
    /*for (size_t i = 0; i < size; i++) {
      std::cout<<_T.data()[i]<<" ";
    }*/
    __T.print();
    std::cout<<'\n';
  }
  std::cout << "-------------------------------------" << '\n';
}
template void dtensor<double>::print(unsigned print_level);
template void dtensor< std::complex<double> >::print(unsigned print_level);




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


//This shoudl be cleaned up
template <typename T>
void dtensor<T>::conj(){
  string letters="";
  char ch='a';
  for (unsigned i=0;i<rank;i++){
    letters.push_back(ch);
    ++ch;
  }
  if (std::is_same<T, std::complex<double>>::value) {
    ////C++ sdf    
    CTF::Transform<T>([](T & d){ d= std::conj(d); })(__T[letters.c_str()]);
  }
  
}
template void dtensor<double>::conj();
template void dtensor< std::complex<double> >::conj();
//---------------------------------------------------------------------------

///What does this do and this should be cleaned up
//---------------------------------------------------------------------------
// special arithmetic operations with another tensor in the same format/pattern
template <typename T>
void dtensor<T>::add(dtensor<T>& A, T c){
   // assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  //for (size_t i = 0; i < size; i++) _T.data()[i] += c*A._T.data()[i];
  unordered_map<string,char> charMap;
  auto indThis = getIndices(charMap);
  auto indA    = A.getIndices(charMap);
  //HACK
  CTF::Scalar<T> cs(c);
  __T[indThis.c_str()] += cs[""]*(A.__T[indA.c_str()]);
  //__T[indThis.c_str()] += std::real(c)*(A.__T[indA.c_str()]);
}
template void dtensor<double>::add(dtensor<double>& A, double c);
template void dtensor< std::complex<double> >::add(dtensor< std::complex<double> >& A, std::complex<double> c);


///Is this the best way to do this!  Hack with conj stuff running around here
template <typename T>
T dtensor<T>::inner_product(dtensor<T>& A){
  //assert(1==2);
  assert(_initted && A._initted);
  assert(rank==A.rank);
  /*T res = 0;
  for (size_t i = 0; i < size; i++) res += cconj(_T.data()[i])*A._T.data()[i];*/
  unordered_map<string,char> charMap;
  CTF::Scalar<T> tot;
  auto left  =   __T[getIndices(charMap).c_str()];
  auto right = A.__T[A.getIndices(charMap).c_str()];
  tot[""] += CTF::Function<T,T,T>([](T l, T r){ return cconj(l)*r;})(left,right);

  return tot;
}
template double dtensor<double>::inner_product(dtensor<double>& A);
template std::complex<double> dtensor< std::complex<double> >::inner_product(dtensor< std::complex<double> >& A);

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Norm
template <typename T>
double dtensor<T>::norm(){
  double res = 0.0;
  return std::real(__T.norm2()); 
}
template double dtensor<double>::norm();
template double dtensor< std::complex<double> >::norm();


template <typename T>
double dtensor<T>::normalize(){
  double res = 0.0;
  res = norm();
  (*this)/=res;
  return res;
}
template double dtensor<double>::normalize();
template double dtensor< std::complex<double> >::normalize();
//-----------------------------------------------------------------------------

template <typename T>
string dtensor<T>::getIndices(unordered_map<string,char> &charMap){
  _indices = indicesToChar(idx_set,charMap);
  return _indices;
}
template string dtensor<double>::getIndices(unordered_map<string,char> &charMap);
template string dtensor<std::complex<double> >::getIndices(unordered_map<string,char> &charMap);

#endif
