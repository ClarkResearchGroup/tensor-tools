#ifndef DENSE_TENSORTRAIN_CLASS
#define DENSE_TENSORTRAIN_CLASS

#include "tt.h"

//---------------------------------------------------------------------------
// Constructors
template <typename T, unsigned N>
dTensorTrain<T, N>::dTensorTrain() {
  tensors_allocated = false;
  length  = 0;
  phy_dim = 0;
  center  = -1;
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());

}
template dTensorTrain<double, 1>::dTensorTrain();
template dTensorTrain<double, 2>::dTensorTrain();
template dTensorTrain<std::complex<double>, 1>::dTensorTrain();
template dTensorTrain<std::complex<double>, 2>::dTensorTrain();


template <typename T, unsigned N>
dTensorTrain<T, N>::dTensorTrain(abstract_sites* s, unsigned bd){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength((*s).N());
  setPhysicalDim((*s).phy_dim());
  setBondDim(bd);
  allocateTensors();
  setZero();
  center  = -1;
}
template dTensorTrain<double, 1>::dTensorTrain(abstract_sites* s, unsigned bd);
template dTensorTrain<double, 2>::dTensorTrain(abstract_sites* s, unsigned bd);
template dTensorTrain<std::complex<double>, 1>::dTensorTrain(abstract_sites* s, unsigned bd);
template dTensorTrain<std::complex<double>, 2>::dTensorTrain(abstract_sites* s, unsigned bd);


template <typename T, unsigned N>
dTensorTrain<T, N>::dTensorTrain(abstract_sites* s, str_vec product_string){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength((*s).N());
  setPhysicalDim((*s).phy_dim());
  setBondDim(1);
  uint_vec product_state = (*s).product_state(product_string);
  allocateTensors(product_state.data());
  center  = -1;
}
template dTensorTrain<double, 1>::dTensorTrain(abstract_sites* s, str_vec product_string);
template dTensorTrain<double, 2>::dTensorTrain(abstract_sites* s, str_vec product_string);
template dTensorTrain<std::complex<double>, 1>::dTensorTrain(abstract_sites* s, str_vec product_string);
template dTensorTrain<std::complex<double>, 2>::dTensorTrain(abstract_sites* s, str_vec product_string);


template <typename T, unsigned N>
dTensorTrain<T, N>::dTensorTrain(unsigned l, unsigned pd, unsigned bd){
  _id = unsigned(1e9*thread_safe_random_double()+1e7*thread_safe_random_double()+1e5*thread_safe_random_double()+1e3*thread_safe_random_double());
  tensors_allocated = false;
  setLength(l);
  setPhysicalDim(pd);
  setBondDim(bd);
  allocateTensors();
  setZero();
  center  = -1;
}
template dTensorTrain<double, 1>::dTensorTrain(unsigned l, unsigned pd, unsigned bd);
template dTensorTrain<double, 2>::dTensorTrain(unsigned l, unsigned pd, unsigned bd);
template dTensorTrain<std::complex<double>, 1>::dTensorTrain(unsigned l, unsigned pd, unsigned bd);
template dTensorTrain<std::complex<double>, 2>::dTensorTrain(unsigned l, unsigned pd, unsigned bd);


template <typename T, unsigned N>
dTensorTrain<T, N>::dTensorTrain(const dTensorTrain<T, N>& other){
  tensors_allocated = false;
  if(other.tensors_allocated){
    _id = other._id;
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    setBondDims(other.bond_dims);
    center = other.center;
    A.clear();
    A = other.A;
    tensors_allocated = true;
  }
}
template dTensorTrain<double, 1>::dTensorTrain(const dTensorTrain<double, 1>& other);
template dTensorTrain<double, 2>::dTensorTrain(const dTensorTrain<double, 2>& other);
template dTensorTrain<std::complex<double>, 1>::dTensorTrain(const dTensorTrain< std::complex<double>, 1>& other);
template dTensorTrain<std::complex<double>, 2>::dTensorTrain(const dTensorTrain< std::complex<double>, 2>& other);


template <typename T, unsigned N>
dTensorTrain<T, N>::dTensorTrain(dTensorTrain<T, N>&& other){
  tensors_allocated = false;
  if(other.tensors_allocated){
    _id = other._id;
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    setBondDims(other.bond_dims);
    center = other.center;
    A.clear();
    A = std::move(other.A);
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
}
template dTensorTrain<double, 1>::dTensorTrain(dTensorTrain<double, 1>&& other);
template dTensorTrain<double, 2>::dTensorTrain(dTensorTrain<double, 2>&& other);
template dTensorTrain<std::complex<double>, 1>::dTensorTrain(dTensorTrain< std::complex<double>, 1>&& other);
template dTensorTrain<std::complex<double>, 2>::dTensorTrain(dTensorTrain< std::complex<double>, 2>&& other);


//---------------------------------------------------------------------------
// Set shapes
template <typename T, unsigned N>
void dTensorTrain<T, N>::setLength(int L) {
  assert(L>2);
  assert(!tensors_allocated);
  length = L;
  bond_dims.resize(length+1);
}
template void dTensorTrain<double, 1>::setLength(int L);
template void dTensorTrain<double, 2>::setLength(int L);
template void dTensorTrain<std::complex<double>, 1>::setLength(int L);
template void dTensorTrain<std::complex<double>, 2>::setLength(int L);


template <typename T, unsigned N>
void dTensorTrain<T, N>::setPhysicalDim(int s){
  assert(s>0);
  assert(!tensors_allocated);
  phy_dim = s;
}
template void dTensorTrain<double, 1>::setPhysicalDim(int s);
template void dTensorTrain<double, 2>::setPhysicalDim(int s);
template void dTensorTrain<std::complex<double>, 1>::setPhysicalDim(int s);
template void dTensorTrain<std::complex<double>, 2>::setPhysicalDim(int s);


template <typename T, unsigned N>
void dTensorTrain<T, N>::setBondDim(int bd){
  assert(bd>0);
  bond_dims.at(0) = 1;
  for (size_t i = 1; i < length; i++) bond_dims.at(i) = bd;
  bond_dims.at(length) = 1;
}
template void dTensorTrain<double, 1>::setBondDim(int bd);
template void dTensorTrain<double, 2>::setBondDim(int bd);
template void dTensorTrain<std::complex<double>, 1>::setBondDim(int bd);
template void dTensorTrain<std::complex<double>, 2>::setBondDim(int bd);


template <typename T, unsigned N>
void dTensorTrain<T, N>::setBondDims(const uint_vec& bds){
  for (size_t i = 0; i < length+1; i++) {
    bond_dims[i] = bds.at(i);
    assert(bond_dims[i] > 0);
  }
  assert(bond_dims[0]==1 && bond_dims[length]==1);
}
template void dTensorTrain<double, 1>::setBondDims(const uint_vec& bds);
template void dTensorTrain<double, 2>::setBondDims(const uint_vec& bds);
template void dTensorTrain<std::complex<double>, 1>::setBondDims(const uint_vec& bds);
template void dTensorTrain<std::complex<double>, 2>::setBondDims(const uint_vec& bds);


//---------------------------------------------------------------------------
// Tensor data management
template <typename T, unsigned N>
void dTensorTrain<T, N>::allocateTensors(unsigned* product_state){
  if(!tensors_allocated) {
    string Link_name_pref = "ID"+to_string(_id)+"Link";
    string Site_name_pref = "Site";
    if(N==1){
      // MPS
      for (size_t i = 0; i < length; i++) {
        string left_link_name  = Link_name_pref+to_string(i);
        string right_link_name = Link_name_pref+to_string(i+1);
        string site_name       = Site_name_pref+to_string(i);
        A.push_back(
          std::move(
            dtensor<T>(
              {bond_dims[i], phy_dim, bond_dims[i+1]},
              {left_link_name, site_name, right_link_name},
              {Link, Site, Link},
              {0,0,0})
          )
        );
      }
      // Set product state
      if(product_state!=nullptr){
        for (size_t i = 0; i < length; i++) {
          A[i].setZero();
          unsigned s0 = A[i]._T.stride(0);
          unsigned s1 = A[i]._T.stride(1);
          unsigned s2 = A[i]._T.stride(2);
          for (size_t j = 0; j < A[i].size; j++) {
            unsigned phy_idx = unsigned(j/s1)%phy_dim;
            if(phy_idx==product_state[i]){
              A[i]._T.data()[j] = 1.0;
            }
          }
        }
      }
    }
    if(N==2){
      // MPO
      for (size_t i = 0; i < length; i++) {
        string left_link_name  = Link_name_pref+to_string(i);
        string right_link_name = Link_name_pref+to_string(i+1);
        string site_name       = Site_name_pref+to_string(i);
        A.push_back(
          std::move(
            dtensor<T>(
              {bond_dims[i], phy_dim, phy_dim, bond_dims[i+1]},
              {left_link_name, site_name, site_name, right_link_name},
              {Link, Site, Site, Link},
              {0,0,1,0})
          )
        );
      }
    }
    tensors_allocated = true;
  }
}
template void dTensorTrain<double, 1>::allocateTensors(unsigned* product_state);
template void dTensorTrain<double, 2>::allocateTensors(unsigned* product_state);
template void dTensorTrain<std::complex<double>, 1>::allocateTensors(unsigned* product_state);
template void dTensorTrain<std::complex<double>, 2>::allocateTensors(unsigned* product_state);


template <typename T, unsigned N>
void dTensorTrain<T, N>::freeTensors(){
  if(tensors_allocated) {
    A.clear();
    tensors_allocated = false;
  }
}
template void dTensorTrain<double, 1>::freeTensors();
template void dTensorTrain<double, 2>::freeTensors();
template void dTensorTrain<std::complex<double>, 1>::freeTensors();
template void dTensorTrain<std::complex<double>, 2>::freeTensors();


//---------------------------------------------------------------------------
// Set tensors values
template <typename T, unsigned N>
void dTensorTrain<T, N>::setZero(){
  if(tensors_allocated) {
    for (size_t i = 0; i < length; i++) {
      A[i].setZero();
    }
  }
}
template void dTensorTrain<double, 1>::setZero();
template void dTensorTrain<double, 2>::setZero();
template void dTensorTrain<std::complex<double>, 1>::setZero();
template void dTensorTrain<std::complex<double>, 2>::setZero();


template <typename T, unsigned N>
void dTensorTrain<T, N>::setRandom(){
  if(tensors_allocated) {
    for (size_t i = 0; i < length; i++) {
      A[i].setRandom();
    }
  }
}
template void dTensorTrain<double, 1>::setRandom();
template void dTensorTrain<double, 2>::setRandom();
template void dTensorTrain<std::complex<double>, 1>::setRandom();
template void dTensorTrain<std::complex<double>, 2>::setRandom();


//---------------------------------------------------------------------------
// Print dTensorTrain information
template <typename T, unsigned N>
void dTensorTrain<T, N>::print(int level){
  std::cout<<"---------------------------------------------------"<<'\n';
  std::cout << "length = " << length << ",  phy_dim = " << phy_dim << ", center = " << center << ", id = " << _id << '\n';
  std::cout << "bond_dims vector = " << " ";
  for(auto v : bond_dims) std::cout << v << " ";
  std::cout<<std::endl;
  if(level>=1){
    std::cout << "Information of individual tensors:" << '\n';
    if(tensors_allocated){
      for (size_t i = 0; i < length; i++) {
        A[i].print(std::max(level-1,0));
      }
    }
  }
  std::cout<<"---------------------------------------------------"<<'\n';
}
template void dTensorTrain<double, 1>::print(int level);
template void dTensorTrain<double, 2>::print(int level);
template void dTensorTrain<std::complex<double>, 1>::print(int level);
template void dTensorTrain<std::complex<double>, 2>::print(int level);

//---------------------------------------------------------------------------
// Basic arithmetic operations
template <typename T, unsigned N>
dTensorTrain<T, N>& dTensorTrain<T, N>::operator = (const dTensorTrain<T, N>& other){
  if(this!=&other && other.tensors_allocated){
    _id = other._id;
    freeTensors();
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    setBondDims(other.bond_dims);
    A = other.A;
    tensors_allocated = true;
  }
  return *this;
}
template dTensorTrain<double, 1>& dTensorTrain<double, 1>::operator = (const dTensorTrain<double, 1>& other);
template dTensorTrain<double, 2>& dTensorTrain<double, 2>::operator = (const dTensorTrain<double, 2>& other);
template dTensorTrain<std::complex<double>, 1>& dTensorTrain<std::complex<double>, 1>::operator = (const dTensorTrain<std::complex<double>, 1>& other);
template dTensorTrain<std::complex<double>, 2>& dTensorTrain<std::complex<double>, 2>::operator = (const dTensorTrain<std::complex<double>, 2>& other);


template <typename T, unsigned N>
dTensorTrain<T, N>& dTensorTrain<T, N>::operator = (dTensorTrain<T, N>&& other){
  if(this!=&other && other.tensors_allocated){
    _id = other._id;
    freeTensors();
    setLength(other.length);
    setPhysicalDim(other.phy_dim);
    setBondDims(other.bond_dims);
    A = std::move(other.A);
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
  return *this;
}
template dTensorTrain<double, 1>& dTensorTrain<double, 1>::operator = (dTensorTrain<double, 1>&& other);
template dTensorTrain<double, 2>& dTensorTrain<double, 2>::operator = (dTensorTrain<double, 2>&& other);
template dTensorTrain<std::complex<double>, 1>& dTensorTrain<std::complex<double>, 1>::operator = (dTensorTrain<std::complex<double>, 1>&& other);
template dTensorTrain<std::complex<double>, 2>& dTensorTrain<std::complex<double>, 2>::operator = (dTensorTrain<std::complex<double>, 2>&& other);


template <typename T, unsigned N>
dTensorTrain<T, N>& dTensorTrain<T, N>::operator *= (const T c){
  assert(tensors_allocated);
  A[0] *= c;
  return *this;
}
template dTensorTrain<double, 1>& dTensorTrain<double, 1>::operator *= (const double c);
template dTensorTrain<double, 2>& dTensorTrain<double, 2>::operator *= (const double c);
template dTensorTrain<std::complex<double>, 1>& dTensorTrain<std::complex<double>, 1>::operator *= (const  std::complex<double> c);
template dTensorTrain<std::complex<double>, 2>& dTensorTrain<std::complex<double>, 2>::operator *= (const  std::complex<double> c);


template <typename T, unsigned N>
dTensorTrain<T, N>& dTensorTrain<T, N>::operator /= (const T c){
  assert(tensors_allocated);
  A[0] /= c;
  return *this;
}
template dTensorTrain<double, 1>& dTensorTrain<double, 1>::operator /= (const double c);
template dTensorTrain<double, 2>& dTensorTrain<double, 2>::operator /= (const double c);
template dTensorTrain< std::complex<double>, 1>& dTensorTrain<std::complex<double>, 1>::operator /= (const  std::complex<double> c);
template dTensorTrain< std::complex<double>, 2>& dTensorTrain<std::complex<double>, 2>::operator /= (const  std::complex<double> c);


template <typename T, unsigned N>
dTensorTrain<T, N> dTensorTrain<T, N>::operator * (const T c) const{
  assert(tensors_allocated);
  dTensorTrain<T, N> t;
  t = *this;
  t.A[0] *= c;
  return t;
}
template dTensorTrain<double, 1> dTensorTrain<double, 1>::operator * (const double c) const;
template dTensorTrain<double, 2> dTensorTrain<double, 2>::operator * (const double c) const;
template dTensorTrain< std::complex<double>, 1> dTensorTrain<std::complex<double>, 1>::operator * (const  std::complex<double> c) const;
template dTensorTrain< std::complex<double>, 2 > dTensorTrain<std::complex<double>, 2>::operator * (const  std::complex<double> c) const;


template <typename T, unsigned N>
dTensorTrain<T, N> dTensorTrain<T, N>::operator / (const T c) const{
  assert(tensors_allocated);
  dTensorTrain<T, N> t;
  t = *this;
  t.A[0] /= c;
  return t;
}
template dTensorTrain<double, 1> dTensorTrain<double, 1>::operator / (const double c) const;
template dTensorTrain<double, 2> dTensorTrain<double, 2>::operator / (const double c) const;
template dTensorTrain< std::complex<double>, 1> dTensorTrain<std::complex<double>, 1>::operator / (const  std::complex<double> c) const;
template dTensorTrain< std::complex<double>, 2> dTensorTrain<std::complex<double>, 2>::operator / (const  std::complex<double> c) const;


template <typename T, unsigned N>
void dTensorTrain<T, N>::save(std::string fn){
  assert(tensors_allocated);
  ezh5::File fh5 (fn, H5F_ACC_TRUNC);
  fh5["length"] = length;
  fh5["phy_dim"] = phy_dim;
  fh5["bond_dims"] = bond_dims;
  fh5["center"] = center;
  fh5["id"] = _id;
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].save(nd);
  }
}
template void dTensorTrain<double, 1>::save(std::string fn);
template void dTensorTrain<double, 2>::save(std::string fn);
template void dTensorTrain<std::complex<double>, 1>::save(std::string fn);
template void dTensorTrain<std::complex<double>, 2>::save(std::string fn);


template <typename T, unsigned N>
void dTensorTrain<T, N>::load(std::string fn){
  freeTensors();
  ezh5::File fh5 (fn, H5F_ACC_RDONLY);
  fh5["length"] >> length;
  fh5["phy_dim"] >> phy_dim;
  fh5["bond_dims"] >> bond_dims;
  fh5["center"] >> center;
  fh5["id"] >> _id;
  allocateTensors();
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].load(nd);
  }
}
template void dTensorTrain<double, 1>::load(std::string fn);
template void dTensorTrain<double, 2>::load(std::string fn);
template void dTensorTrain<std::complex<double>, 1>::load(std::string fn);
template void dTensorTrain<std::complex<double>, 2>::load(std::string fn);


template <typename T, unsigned N>
void dTensorTrain<T, N>::save(ezh5::Node& fh5){
  assert(tensors_allocated);
  fh5["length"] = length;
  fh5["phy_dim"] = phy_dim;
  fh5["bond_dims"] = bond_dims;
  fh5["center"] = center;
  fh5["id"] = _id;
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].save(nd);
  }
}
template void dTensorTrain<double, 1>::save(ezh5::Node& fh5);
template void dTensorTrain<double, 2>::save(ezh5::Node& fh5);
template void dTensorTrain<std::complex<double>, 1>::save(ezh5::Node& fh5);
template void dTensorTrain<std::complex<double>, 2>::save(ezh5::Node& fh5);


template <typename T, unsigned N>
void dTensorTrain<T, N>::load(ezh5::Node& fh5){
  freeTensors();
  fh5["length"] >> length;
  fh5["phy_dim"] >> phy_dim;
  fh5["bond_dims"] >> bond_dims;
  fh5["center"] >> center;
  fh5["id"] >> _id;
  allocateTensors();
  for (size_t i = 0; i < length; i++){
    ezh5::Node nd = fh5["Tensor"+to_string(i)];
    A[i].load(nd);
  }
}
template void dTensorTrain<double, 1>::load(ezh5::Node& fh5);
template void dTensorTrain<double, 2>::load(ezh5::Node& fh5);
template void dTensorTrain<std::complex<double>, 1>::load(ezh5::Node& fh5);
template void dTensorTrain<std::complex<double>, 2>::load(ezh5::Node& fh5);


template <typename T, unsigned N>
void dTensorTrain<T, N>::rc(){
  assert(tensors_allocated);
  dtensor<T> U,V;
  vector<double> S;
  for (size_t i = length-1; i > 0; i--) {
    vector<dtensor_index> left;
    vector<dtensor_index> right;
    dtensor_index mid;
    string LinkName = "ID"+to_string(_id)+"Link"+to_string(i);
    // Separate dtensor_index
    for (size_t j = 0; j < A[i].rank; j++) {
      string idx_name = A[i].idx_set[j].name();
      if (idx_name == LinkName) {
        left.push_back(A[i].idx_set[j]);
        mid = A[i].idx_set[j];
      }else{
        right.push_back(A[i].idx_set[j]);
      }
    }
    // SVD
    svd(A[i],left,right,U,V,S,MoveFromRight);
    mid.resize(S.size());
    bond_dims[i] = S.size();
    A[i] = V;
    A[i].idx_set[0] = mid;
    A[i-1] = std::move(A[i-1]*U);
    A[i-1].idx_set.back() = mid;
  }
  center = 0;
}
template void dTensorTrain<double, 1>::rc();
template void dTensorTrain<double, 2>::rc();
template void dTensorTrain<std::complex<double>, 1>::rc();
template void dTensorTrain<std::complex<double>, 2>::rc();


template <typename T, unsigned N>
void dTensorTrain<T, N>::lc(){
  assert(tensors_allocated);
  dtensor<T> U,V;
  vector<double> S;
  for (size_t i = 0; i < length-1; i++) {
    vector<dtensor_index> left;
    vector<dtensor_index> right;
    dtensor_index mid;
    string LinkName = "ID"+to_string(_id)+"Link"+to_string(i+1);
    // Separate dtensor_index
    for (size_t j = 0; j < A[i].rank; j++) {
      string idx_name = A[i].idx_set[j].name();
      if (idx_name == LinkName) {
        right.push_back(A[i].idx_set[j]);
        mid = A[i].idx_set[j];
      }else{
        left.push_back(A[i].idx_set[j]);
      }
    }
    // SVD
    svd(A[i],left,right,U,V,S,MoveFromLeft);
    mid.resize(S.size());
    bond_dims[i+1] = S.size();
    A[i] = U;
    A[i].idx_set.back() = mid;
    A[i+1] = std::move(V*A[i+1]);
    A[i+1].idx_set[0] = mid;
  }
  center = length-1;
}
template void dTensorTrain<double, 1>::lc();
template void dTensorTrain<double, 2>::lc();
template void dTensorTrain<std::complex<double>, 1>::lc();
template void dTensorTrain<std::complex<double>, 2>::lc();


template <typename T, unsigned N>
void dTensorTrain<T, N>::normalize(){
  assert(tensors_allocated);
  double nm = norm();
  if(center == -1)
    A[0] /= nm;
  else
    A[center] /= nm;
}
template void dTensorTrain<double, 1>::normalize();
template void dTensorTrain<double, 2>::normalize();
template void dTensorTrain<std::complex<double>, 1>::normalize();
template void dTensorTrain<std::complex<double>, 2>::normalize();


template <typename T, unsigned N>
double dTensorTrain<T, N>::norm(){
  assert(tensors_allocated);
  dTensorTrain<T, N> t;
  t = *this;
  t.rc();
  return t.A[0].norm();
}
template double dTensorTrain<double, 1>::norm();
template double dTensorTrain<double, 2>::norm();
template double dTensorTrain<std::complex<double>, 1>::norm();
template double dTensorTrain<std::complex<double>, 2>::norm();


//---------------------------------------------------------------------------
// Adjust canonical center
template <typename T, unsigned N>
double dTensorTrain<T, N>::position(int site){
  assert(tensors_allocated);
  assert(site>=0 && site<int(length));
  // Initialize center position if not in canonical form
  if(center == -1) rc();
  dtensor<T> U,V;
  vector<double> S;
  while( center!=site ){
    vector<dtensor_index> left;
    vector<dtensor_index> right;
    dtensor_index mid;
    if(center>site){
      // Move center to the left
      string LinkName = "ID"+to_string(_id)+"Link"+to_string(center);
      for (size_t j = 0; j < A[center].rank; j++) {
        string idx_name = A[center].idx_set[j].name();
        if (idx_name == LinkName) {
          left.push_back(A[center].idx_set[j]);
          mid = A[center].idx_set[j];
        }else{
          right.push_back(A[center].idx_set[j]);
        }
      }
      svd(A[center],left,right,U,V,S,MoveFromRight);
      mid.resize(S.size());
      bond_dims[center] = S.size();
      A[center] = V;
      A[center].idx_set[0] = mid;
      A[center-1] = std::move(A[center-1]*U);
      A[center-1].idx_set.back() = mid;
      --center;
    }else if(center<site){
      // Move center to the right
      string LinkName = "ID"+to_string(_id)+"Link"+to_string(center+1);
      for (size_t j = 0; j < A[center].rank; j++) {
        string idx_name = A[center].idx_set[j].name();
        if (idx_name == LinkName) {
          right.push_back(A[center].idx_set[j]);
          mid = A[center].idx_set[j];
        }else{
          left.push_back(A[center].idx_set[j]);
        }
      }
      svd(A[center],left,right,U,V,S,MoveFromLeft);
      mid.resize(S.size());
      bond_dims[center+1] = S.size();
      A[center] = U;
      A[center].idx_set.back() = mid;
      A[center+1] = std::move(V*A[center+1]);
      A[center+1].idx_set[0] = mid;
      ++center;
    }
  }
  double vNEE = 0.0;
  for(auto sg:S){
    if(sg>1e-24) vNEE -= sg*sg*std::log(sg*sg);
  }
  return vNEE;
}
template double dTensorTrain<double, 1>::position(int site);
template double dTensorTrain<double, 2>::position(int site);
template double dTensorTrain<std::complex<double>, 1>::position(int site);
template double dTensorTrain<std::complex<double>, 2>::position(int site);

#endif
