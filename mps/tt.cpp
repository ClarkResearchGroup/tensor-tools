#ifndef My_TENSORTRAIN_CLASS
#define My_TENSORTRAIN_CLASS

#include "tt.h"

template <typename T>
TensorTrain<T>::TensorTrain() {
  tensors_allocated = false;
  TTID = unsigned(1e8*drand48()+1e5*drand48());
}
template TensorTrain<double>::TensorTrain();
template TensorTrain< std::complex<double> >::TensorTrain();

template <typename T>
TensorTrain<T>::TensorTrain(const TensorTrain<T>& other){
  tensors_allocated = false;
  TTID = unsigned(1e8*drand48()+1e5*drand48());
  if(other.tensors_allocated){
    setLength(other.length);
    setIndexSize(other.index_size);
    setBondDims(other.bond_dims);
    allocateTensors();
    for (int i = 0; i < length; i++){
      for (int j = 0; j < index_size; j++) {
        M[i][j] = other.M[i][j];
      }
    }
    tensors_allocated = true;
  }
}
template TensorTrain<double>::TensorTrain(const TensorTrain<double>& other);
template TensorTrain< std::complex<double> >::TensorTrain(const TensorTrain< std::complex<double> >& other);

template <typename T>
TensorTrain<T>::TensorTrain(TensorTrain<T>&& other){
  tensors_allocated = false;
  TTID = unsigned(1e8*drand48()+1e5*drand48());
  if(other.tensors_allocated){
    setLength(other.length);
    setIndexSize(other.index_size);
    setBondDims(other.bond_dims);
    allocateTensors();
    M = other.M;
    other.M = nullptr;
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
}
template TensorTrain<double>::TensorTrain(TensorTrain<double>&& other);
template TensorTrain< std::complex<double> >::TensorTrain(TensorTrain< std::complex<double> >&& other);

template <typename T>
void TensorTrain<T>::setLength(int L) {
  if(!tensors_allocated) {
    length = L;
    bond_dims.resize(length+1);
  }
  assert(L>1);
}
template void TensorTrain<double>::setLength(int L);
template void TensorTrain< std::complex<double> >::setLength(int L);

template <typename T>
void TensorTrain<T>::setIndexSize(int s){
  if(!tensors_allocated) index_size = s;
  assert(s>0);
}
template void TensorTrain<double>::setIndexSize(int s);
template void TensorTrain< std::complex<double> >::setIndexSize(int s);

template <typename T>
void TensorTrain<T>::setBondDim(int bd){
  assert(bd>0);
  bond_dims.at(0) = 1;
  for (int i = 1; i < length; i++) bond_dims.at(i) = bd;
  bond_dims.at(length) = 1;
}
template void TensorTrain<double>::setBondDim(int bd);
template void TensorTrain< std::complex<double> >::setBondDim(int bd);

template <typename T>
void TensorTrain<T>::setBondDims(std::vector<int> bds){
  for (int i = 0; i < length+1; i++) {
    bond_dims[i] = bds.at(i);
    assert(bond_dims[i] > 0);
  }
  assert(bond_dims[0]==1 && bond_dims[length]==1);
}
template void TensorTrain<double>::setBondDims(std::vector<int> bds);
template void TensorTrain< std::complex<double> >::setBondDims(std::vector<int> bds);

template <typename T>
void TensorTrain<T>::allocateTensors(){
  if(!tensors_allocated) {
    M = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>* [length];
    for (int i = 0; i < length; i++) {
      M[i] = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> [index_size];
      for (int j = 0; j < index_size; j++) {
        M[i][j].setZero(bond_dims[i],bond_dims[i+1]);
      }
    }
    tensors_allocated = true;
  }
}
template void TensorTrain<double>::allocateTensors();
template void TensorTrain< std::complex<double> >::allocateTensors();

template <typename T>
void TensorTrain<T>::freeTensors(){
  if(tensors_allocated) {
    if(M!=nullptr){
      for (int i = 0; i < length; i++) delete [] M[i];
      delete [] M;
    }
    tensors_allocated = false;
  }
}
template void TensorTrain<double>::freeTensors();
template void TensorTrain< std::complex<double> >::freeTensors();

template <typename T>
TensorTrain<T>::~TensorTrain(){
  freeTensors();
}
template TensorTrain<double>::~TensorTrain();
template TensorTrain< std::complex<double> >::~TensorTrain();

template <typename T>
void TensorTrain<T>::setZero(){
  if(tensors_allocated) {
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < index_size; j++) {
        M[i][j].setZero(bond_dims[i],bond_dims[i+1]);
      }
    }
  }
}
template void TensorTrain<double>::setZero();
template void TensorTrain< std::complex<double> >::setZero();

template <typename T>
void TensorTrain<T>::setZero(int bd){
  setBondDim(bd);
  if(tensors_allocated) {
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < index_size; j++) {
        M[i][j].setZero(bond_dims[i],bond_dims[i+1]);
      }
    }
  }
}
template void TensorTrain<double>::setZero(int bd);
template void TensorTrain< std::complex<double> >::setZero(int bd);

template <typename T>
void TensorTrain<T>::setRandom(){
  if(tensors_allocated) {
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < index_size; j++) {
        M[i][j].setRandom(bond_dims[i],bond_dims[i+1]);
      }
    }
  }
}
template void TensorTrain<double>::setRandom();
template void TensorTrain< std::complex<double> >::setRandom();

template <typename T>
void TensorTrain<T>::setRandom(int bd){
  setBondDim(bd);
  if(tensors_allocated) {
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < index_size; j++) {
        M[i][j].setRandom(bond_dims[i],bond_dims[i+1]);
      }
    }
  }
}
template void TensorTrain<double>::setRandom(int bd);
template void TensorTrain< std::complex<double> >::setRandom(int bd);

template <typename T>
void TensorTrain<T>::print(int level){
  std::cout << length << " " << index_size << '\n';
  for(auto v : bond_dims) std::cout << v << " ";
  std::cout<<std::endl;
  if(level>0){
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < index_size; j++) {
        std::cout << M[i][j] << '\n';
      }
      std::cout << '\n';
    }
  }
}
template void TensorTrain<double>::print(int level);
template void TensorTrain< std::complex<double> >::print(int level);

template <typename T>
TensorTrain<T>& TensorTrain<T>::operator = (const TensorTrain<T>& other){
  assert(other.tensors_allocated);
  if(this!=&other){
    freeTensors();
    setLength(other.length);
    setIndexSize(other.index_size);
    setBondDims(other.bond_dims);
    allocateTensors();
    for (int i = 0; i < length; i++){
      for (int j = 0; j < index_size; j++) {
        M[i][j] = other.M[i][j];
      }
    }
  }
  return *this;
}
template TensorTrain<double>& TensorTrain<double>::operator = (const TensorTrain<double>& other);
template TensorTrain< std::complex<double> >& TensorTrain< std::complex<double> >::operator = (const TensorTrain< std::complex<double> >& other);

template <typename T>
TensorTrain<T>& TensorTrain<T>::operator = (TensorTrain<T>&& other){
  assert(other.tensors_allocated);
  if(this!=&other){
    freeTensors();
    setLength(other.length);
    setIndexSize(other.index_size);
    setBondDims(other.bond_dims);
    M = other.M;
    other.M = nullptr;
    tensors_allocated = true;
    other.tensors_allocated = false;
  }
  return *this;
}
template TensorTrain<double>& TensorTrain<double>::operator = (TensorTrain<double>&& other);
template TensorTrain< std::complex<double> >& TensorTrain< std::complex<double> >::operator = (TensorTrain< std::complex<double> >&& other);

template <typename T>
TensorTrain<T>& TensorTrain<T>::operator += (const TensorTrain<T>& A){
  assert(tensors_allocated && A.tensors_allocated);
  TensorTrain<T> t;
  t.setLength(length);
  t.setIndexSize(index_size);
  t.bond_dims[0] = 1;
  for (int i = 1; i < length; i++){
    t.bond_dims[i] = bond_dims[i] + A.bond_dims[i];
  }
  t.bond_dims[length] = 1;
  t.allocateTensors();
  for (int j = 0; j < index_size; j++) {
    // First site
    t.M[0][j].block(0,0,1,bond_dims[1]) = M[0][j];
    t.M[0][j].block(0,bond_dims[1],1,A.bond_dims[1]) = A.M[0][j];
    // Last site
    t.M[length-1][j].block(0,0,bond_dims[length-1],1) = M[length-1][j];
    t.M[length-1][j].block(bond_dims[length-1],0,A.bond_dims[length-1],1) = A.M[length-1][j];
  }
  for (int i = 1; i < length-1; i++){
    for (int j = 0; j < index_size; j++) {
      t.M[i][j].block(0,0,bond_dims[i],bond_dims[i+1]) = M[i][j];
      t.M[i][j].block(bond_dims[i],bond_dims[i+1],A.bond_dims[i],A.bond_dims[i+1]) = A.M[i][j];
    }
  }
  *this = t;
  return *this;
}
template TensorTrain<double>& TensorTrain<double>::operator += (const TensorTrain<double>& A);
template TensorTrain< std::complex<double> >& TensorTrain< std::complex<double> >::operator += (const TensorTrain< std::complex<double> >& A);

template <typename T>
TensorTrain<T>& TensorTrain<T>::operator -= (const TensorTrain<T>& A){
  assert(tensors_allocated && A.tensors_allocated);
  TensorTrain<T> t;
  t.setLength(length);
  t.setIndexSize(index_size);
  t.bond_dims[0] = 1;
  for (int i = 1; i < length; i++){
    t.bond_dims[i] = this->bond_dims[i] + A.bond_dims[i];
  }
  t.bond_dims[length] = 1;
  t.allocateTensors();
  for (int j = 0; j < index_size; j++) {
    // First site
    t.M[0][j].block(0,0,1,bond_dims[1]) = M[0][j];
    t.M[0][j].block(0,bond_dims[1],1,A.bond_dims[1]) = -A.M[0][j];
    // Last site
    t.M[length-1][j].block(0,0,bond_dims[length-1],1) = M[length-1][j];
    t.M[length-1][j].block(bond_dims[length-1],0,A.bond_dims[length-1],1) = A.M[length-1][j];
  }
  for (int i = 1; i < length-1; i++){
    for (int j = 0; j < index_size; j++) {
      t.M[i][j].block(0,0,bond_dims[i],bond_dims[i+1]) = M[i][j];
      t.M[i][j].block(bond_dims[i],bond_dims[i+1],A.bond_dims[i],A.bond_dims[i+1]) = A.M[i][j];
    }
  }
  *this = t;
  return *this;
}
template TensorTrain<double>& TensorTrain<double>::operator -= (const TensorTrain<double>& A);
template TensorTrain< std::complex<double> >& TensorTrain< std::complex<double> >::operator -= (const TensorTrain< std::complex<double> >& A);

template <typename T>
TensorTrain<T> TensorTrain<T>::operator + (const TensorTrain<T>& A){
  assert(tensors_allocated && A.tensors_allocated);
  TensorTrain<T> t;
  t.setLength(length);
  t.setIndexSize(index_size);
  t.bond_dims[0] = 1;
  for (int i = 1; i < length; i++){
    t.bond_dims[i] = this->bond_dims[i] + A.bond_dims[i];
  }
  t.bond_dims[length] = 1;
  t.allocateTensors();
  for (int j = 0; j < index_size; j++) {
    // First site
    t.M[0][j].block(0,0,1,bond_dims[1]) = M[0][j];
    t.M[0][j].block(0,bond_dims[1],1,A.bond_dims[1]) = A.M[0][j];
    // Last site
    t.M[length-1][j].block(0,0,bond_dims[length-1],1) = M[length-1][j];
    t.M[length-1][j].block(bond_dims[length-1],0,A.bond_dims[length-1],1) = A.M[length-1][j];
  }
  for (int i = 1; i < length-1; i++){
    for (int j = 0; j < index_size; j++) {
      t.M[i][j].block(0,0,bond_dims[i],bond_dims[i+1]) = M[i][j];
      t.M[i][j].block(bond_dims[i],bond_dims[i+1],A.bond_dims[i],A.bond_dims[i+1]) = A.M[i][j];
    }
  }
  return t;
}
template TensorTrain<double> TensorTrain<double>::operator + (const TensorTrain<double>& A);
template TensorTrain< std::complex<double> > TensorTrain< std::complex<double> >::operator + (const TensorTrain< std::complex<double> >& A);

template <typename T>
TensorTrain<T> TensorTrain<T>::operator - (const TensorTrain<T>& A){
  assert(tensors_allocated && A.tensors_allocated);
  TensorTrain<T> t;
  t.setLength(length);
  t.setIndexSize(index_size);
  t.bond_dims[0] = 1;
  for (int i = 1; i < length; i++){
    t.bond_dims[i] = this->bond_dims[i] + A.bond_dims[i];
  }
  t.bond_dims[length] = 1;
  t.allocateTensors();
  for (int j = 0; j < index_size; j++) {
    // First site
    t.M[0][j].block(0,0,1,bond_dims[1]) = M[0][j];
    t.M[0][j].block(0,bond_dims[1],1,A.bond_dims[1]) = -A.M[0][j];
    // Last site
    t.M[length-1][j].block(0,0,bond_dims[length-1],1) = M[length-1][j];
    t.M[length-1][j].block(bond_dims[length-1],0,A.bond_dims[length-1],1) = A.M[length-1][j];
  }
  for (int i = 1; i < length-1; i++){
    for (int j = 0; j < index_size; j++) {
      t.M[i][j].block(0,0,bond_dims[i],bond_dims[i+1]) = M[i][j];
      t.M[i][j].block(bond_dims[i],bond_dims[i+1],A.bond_dims[i],A.bond_dims[i+1]) = A.M[i][j];
    }
  }
  return t;
}
template TensorTrain<double> TensorTrain<double>::operator - (const TensorTrain<double>& A);
template TensorTrain< std::complex<double> > TensorTrain< std::complex<double> >::operator - (const TensorTrain< std::complex<double> >& A);

template <typename T>
TensorTrain<T>& TensorTrain<T>::operator *= (const T& c){
  assert(tensors_allocated);
  for (int j = 0; j < index_size; j++) {
    // First site
    M[0][j] *= c;
  }
  return *this;
}
template TensorTrain<double>& TensorTrain<double>::operator *= (const double& c);
template TensorTrain< std::complex<double> >& TensorTrain< std::complex<double> >::operator *= (const  std::complex<double> & c);

template <typename T>
TensorTrain<T>& TensorTrain<T>::operator /= (const T& c){
  assert(tensors_allocated);
  for (int j = 0; j < index_size; j++) {
    // First site
    M[0][j] /= c;
  }
  return *this;
}
template TensorTrain<double>& TensorTrain<double>::operator /= (const double& c);
template TensorTrain< std::complex<double> >& TensorTrain< std::complex<double> >::operator /= (const  std::complex<double> & c);

template <typename T>
TensorTrain<T> TensorTrain<T>::operator * (const T& c){
  assert(tensors_allocated);
  TensorTrain<T> t;
  t = *this;
  for (int j = 0; j < index_size; j++) {
    // First site
    t.M[0][j] *= c;
  }
  return t;
}
template TensorTrain<double> TensorTrain<double>::operator * (const double& c);
template TensorTrain< std::complex<double> > TensorTrain< std::complex<double> >::operator * (const  std::complex<double> & c);

template <typename T>
TensorTrain<T> TensorTrain<T>::operator / (const T& c){
  assert(tensors_allocated);

  TensorTrain<T> t;
  t = *this;
  for (int j = 0; j < index_size; j++) {
    // First site
    t.M[0][j] /= c;
  }
  return t;
}
template TensorTrain<double> TensorTrain<double>::operator / (const double& c);
template TensorTrain< std::complex<double> > TensorTrain< std::complex<double> >::operator / (const  std::complex<double> & c);

template <typename T>
void TensorTrain<T>::save(std::string fn){
  assert(tensors_allocated);
#ifdef USE_EZH5
  ezh5::File fh5 (fn, H5F_ACC_TRUNC);
  fh5["length"] = length;
  fh5["index_size"] = index_size;
  fh5["bond_dims"] = bond_dims;
  for (int i = 0; i < length; i++){
    for (int j = 0; j < index_size; j++) {
      std::string tensor_name = "M_site"+std::to_string(i)+"_"+std::to_string(j);
      fh5[tensor_name] = M[i][j];
    }
  }
#else
  std::ofstream fout;
  fout.precision(15);
  fout.open(fn.c_str());
  fout<<length<<" "<<index_size<<std::endl;
  for(auto v : bond_dims) fout<<v<<" ";
  fout<<std::endl;
  // MPS matrices
  for(int i = 0; i < length; ++i){
    for(int j = 0; j < index_size; ++j){
      for(int k = 0; k < M[i][j].size(); ++k) {
        fout<<M[i][j].data()[k]<<" ";
      }
      fout<<std::endl;
    }
  }
  fout.close();
#endif
}
template void TensorTrain<double>::save(std::string fn);
template void TensorTrain< std::complex<double> >::save(std::string fn);

template <typename T>
void TensorTrain<T>::load(std::string fn){
#ifdef USE_EZH5
  freeTensors();
  ezh5::File fh5 (fn, H5F_ACC_RDONLY);
  fh5["length"] >> length;
  fh5["index_size"] >> index_size;
  fh5["bond_dims"] >> bond_dims;
  allocateTensors();
  for (int i = 0; i < length; i++){
    for (int j = 0; j < index_size; j++) {
      std::string tensor_name = "M_site"+std::to_string(i)+"_"+std::to_string(j);
      fh5[tensor_name] >> M[i][j];
    }
  }
#else
  freeTensors();
  std::ifstream fin;
  fin.open(fn);
  fin>>length;
  fin>>index_size;
  bond_dims.resize(length+1);
  for (int i = 0; i < length+1; i++) {
    fin>>bond_dims[i];
  }
  allocateTensors();
  double temp;
  std::vector<double> vec;
  for(int i = 0; i < length; ++i){
    for(int j = 0; j < index_size; ++j){
      vec.clear();
      for(int k = 0; k < bond_dims[i]*bond_dims[i+1]; ++k) {
        fin>>temp;
        vec.push_back(temp);
      }
      M[i][j] = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(&vec[0],bond_dims[i],bond_dims[i+1]);
    }
  }
  fin.close();
#endif
}
template void TensorTrain<double>::load(std::string fn);
template void TensorTrain< std::complex<double> >::load(std::string fn);

template <typename T>
void TensorTrain<T>::rc(){
  assert(tensors_allocated);
  int row, col;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, Q, R;
  for(int i = length-1; i > 0; --i){
    row=M[i][0].rows();
    col=M[i][0].cols();
    TM.resize(row,index_size*col);
    for(int j = 0; j < index_size; j++) {
      TM.block(0,j*col,row,col) = M[i][j];
    }
    TM.adjointInPlace();
    QR(TM,Q);
    R.noalias() = TM.adjoint()*Q;
    Q.adjointInPlace();
    for(int j = 0; j < index_size; j++) {
      int min_row = std::min(row,int(Q.cols()));
      M[i][j].setZero();
      M[i][j].block(0,0,min_row,col)=Q.block(0,j*col,min_row,col);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
      tempM.noalias()=M[i-1][j]*R;
      M[i-1][j].setZero();
      M[i-1][j].block(0,0,tempM.rows(),tempM.cols())=tempM;
    }
  }
}
template void TensorTrain<double>::rc();
template void TensorTrain< std::complex<double> >::rc();

template <typename T>
void TensorTrain<T>::lc(){
	assert(tensors_allocated);
  int row, col;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, Q, R;
  for(int i = 0; i < length-1; i++){
    row=M[i][0].rows();
    col=M[i][0].cols();
    TM.resize(index_size*row,col);
    for(int j = 0; j < index_size; j++) {
      TM.block(j*row,0,row,col)=M[i][j];
    }
    QR(TM,Q);
    R.noalias() = Q.adjoint()*TM;
    for(int j = 0; j < index_size; j++) {
      int min_col = std::min(col,int(Q.rows()));
      M[i][j].setZero();
      M[i][j].block(0,0,row,min_col)=Q.block(j*row,0,row,min_col);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
      tempM.noalias()=R*M[i+1][j];
      M[i+1][j].setZero();
      M[i+1][j].block(0,0,tempM.rows(),tempM.cols())=tempM;
    }
  }
}
template void TensorTrain<double>::lc();
template void TensorTrain< std::complex<double> >::lc();

template <typename T>
void TensorTrain<T>::normalize(){
	assert(tensors_allocated);
  int row, col;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, Q, R;
  for(int i = length-1; i >= 0; --i){
    row=M[i][0].rows();
    col=M[i][0].cols();
    TM.resize(row,index_size*col);
    for(int j = 0; j < index_size; j++) {
      TM.block(0,j*col,row,col) = M[i][j];
    }
    TM.adjointInPlace();
    QR(TM,Q);
    R.noalias() = TM.adjoint()*Q;
    Q.adjointInPlace();
    for(int j = 0; j < index_size; j++) {
      int min_row = std::min(row,int(Q.cols()));
      M[i][j].setZero();
      M[i][j].block(0,0,min_row,col)=Q.block(0,j*col,min_row,col);
      if(i!=0) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
        tempM.noalias()=M[i-1][j]*R;
        M[i-1][j].setZero();
        M[i-1][j].block(0,0,tempM.rows(),tempM.cols())=tempM;
      }
    }
  }
}
template void TensorTrain<double>::normalize();
template void TensorTrain< std::complex<double> >::normalize();

template <typename T>
double TensorTrain<T>::norm(){
  assert(tensors_allocated);
  TensorTrain<T> t;
  t = *this;
  double nm=0.0;
  int row, col;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, Q, R;
  for(int i = length-1; i >= 0; --i){
    row=t.M[i][0].rows();
    col=t.M[i][0].cols();
    TM.resize(row,index_size*col);
    for(int j = 0; j < index_size; j++) {
      TM.block(0,j*col,row,col) = t.M[i][j];
    }
    TM.adjointInPlace();
    QR(TM,Q);
    R.noalias() = TM.adjoint()*Q;
    Q.adjointInPlace();
    for(int j = 0; j < index_size; j++) {
      int min_row = std::min(row,int(Q.cols()));
      t.M[i][j].setZero();
      t.M[i][j].block(0,0,min_row,col)=Q.block(0,j*col,min_row,col);
      if(i!=0) {
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
        tempM.noalias()=M[i-1][j]*R;
        t.M[i-1][j].setZero();
        t.M[i-1][j].block(0,0,tempM.rows(),tempM.cols())=tempM;
      }else{
        nm = std::abs(R(0,0));
      }
    }
  }
  return nm;
}
template double TensorTrain<double>::norm();
template double TensorTrain< std::complex<double> >::norm();

template <typename T>
void TensorTrain<T>::moveLeft(int site, bool print_EE){
  assert(tensors_allocated);
  if(site==length/2){
    // Calculate EE at mid bond
    int bD = *std::max_element(bond_dims.begin(), bond_dims.end());
    double * sv = new double [index_size*bD]();
    int row, col, tDim;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, U, V;
    row=M[site][0].rows();
    col=M[site][0].cols();
    TM.resize(row,index_size*col);
    tDim = std::min(TM.rows(),TM.cols());
    for(int tid = 0; tid < index_size; tid++){
      TM.block(0,tid*col,row,col) = M[site][tid].block(0,0,row,col);
    }
    SVD(TM,index_size*bD,sv,U,V,'l');
    double EE = 0;
    for(int j = 0; j < tDim; ++j){
      if(sv[j]>1e-15) EE -= sv[j]*sv[j]*log(sv[j]*sv[j]);
      // std::cout << sv[j] << " ";
    }
    // std::cout<<endl;
    if(print_EE) std::cout<<EE<<" ";
    // std::cout<<endl;
    for(int tid = 0; tid < index_size; tid++){
      M[site][tid].setZero();
      M[site][tid].block(0,0,std::min(row,tDim),col) = V.block(0,tid*col,std::min(row,tDim),col);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
      tempM.noalias() = M[site-1][tid] * U;
      M[site-1][tid].setZero();
      M[site-1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
    }
    delete [] sv;
  }else{
    int row, col;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, Q, R;
    row=M[site][0].rows();
    col=M[site][0].cols();
    TM.resize(row,index_size*col);
    for(int j = 0; j < index_size; j++) {
      TM.block(0,j*col,row,col) = M[site][j];
    }
    TM.adjointInPlace();
    QR(TM,Q);
    R.noalias() = TM.adjoint()*Q;
    Q.adjointInPlace();
    for(int j = 0; j < index_size; j++) {
      int min_row = std::min(row,int(Q.cols()));
      M[site][j].setZero();
      M[site][j].block(0,0,min_row,col)=Q.block(0,j*col,min_row,col);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
      tempM.noalias()=M[site-1][j]*R;
      M[site-1][j].setZero();
      M[site-1][j].block(0,0,tempM.rows(),tempM.cols())=tempM;
    }
  }
}
template void TensorTrain<double>::moveLeft(int site, bool print_EE);
template void TensorTrain< std::complex<double> >::moveLeft(int site, bool print_EE);

template <typename T>
void TensorTrain<T>::moveRight(int site, bool print_EE){
  assert(tensors_allocated);
  if(site==length/2-1){
    // Calculate EE at mid bond
    int bD = *std::max_element(bond_dims.begin(), bond_dims.end());
    double * sv = new double [index_size*bD]();
    int row, col, tDim;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, U, V;
    row=M[site][0].rows();
    col=M[site][0].cols();
    TM.resize(index_size*row,col);
    tDim = std::min(TM.rows(),TM.cols());
    for(int tid = 0; tid < index_size; tid++){
      TM.block(tid*row,0,row,col) = M[site][tid].block(0,0,row,col);
    }
    SVD(TM,index_size*bD,sv,U,V,'r');
    double EE = 0;
    for(int j = 0; j < tDim; ++j){
      if(sv[j]>1e-15) EE -= sv[j]*sv[j]*log(sv[j]*sv[j]);
      // std::cout << sv[j] << " ";
    }
    // std::cout<<endl;
    if(print_EE) std::cout<<EE<<" ";
    // std::cout<<endl;
    for(int tid = 0; tid < index_size; tid++){
      M[site][tid].setZero();
      M[site][tid].block(0,0,row,std::min(col,tDim)) = U.block(tid*row,0,row,std::min(col,tDim));
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
      tempM.noalias() = V * M[site+1][tid];
      M[site+1][tid].setZero();
      M[site+1][tid].block(0,0,tempM.rows(),tempM.cols())=tempM;
    }
    delete [] sv;
  }else{
    int row, col;
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, Q, R;
    row=M[site][0].rows();
    col=M[site][0].cols();
    TM.resize(index_size*row,col);
    for(int j = 0; j < index_size; j++) {
      TM.block(j*row,0,row,col)=M[site][j];
    }
    QR(TM,Q);
    R.noalias() = Q.adjoint()*TM;
    for(int j = 0; j < index_size; j++) {
      int min_col = std::min(col,int(Q.rows()));
      M[site][j].setZero();
      M[site][j].block(0,0,row,min_col)=Q.block(j*row,0,row,min_col);
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tempM;
      tempM.noalias()=R*M[site+1][j];
      M[site+1][j].setZero();
      M[site+1][j].block(0,0,tempM.rows(),tempM.cols())=tempM;
    }
  }
}
template void TensorTrain<double>::moveRight(int site, bool print_EE);
template void TensorTrain< std::complex<double> >::moveRight(int site, bool print_EE);

template <typename T>
void TensorTrain<T>::svd(double tol, bool dry_run){
  rc();
  int row, col;
  double svd_tol = tol;
  std::vector<double> sv;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> TM, U, V;
  // Left to Right
  for(int i = 0; i < length-1; i++) {
    row=M[i][0].rows();
    col=M[i][0].cols();
    TM.resize(index_size*row,col);
    for(int j = 0; j < index_size; j++) {
      TM.block(j*row,0,row,col)=M[i][j];
    }
    SVD(TM,sv,U,V,svd_tol,dry_run);
    for(int j = 0; j < index_size; j++){
      M[i][j]   = U.block(j*row,0,row,U.cols());
      M[i+1][j] = V * M[i+1][j];
    }
  }
  // Right to Left
  for(int i = length-1; i > 0; i--) {
    row=M[i][0].rows();
    col=M[i][0].cols();
    TM.resize(row,index_size*col);
    for(int j = 0; j < index_size; j++) {
      TM.block(0,j*col,row,col)=M[i][j];
    }
    TM.adjointInPlace();
    SVD(TM,sv,U,V,svd_tol,dry_run);
    U.adjointInPlace();
    V.adjointInPlace();
    for(int j = 0; j < index_size; j++){
      M[i][j]   = U.block(0,j*col,U.rows(),col);
      M[i-1][j] = M[i-1][j] * V;
    }
  }
  // calibrate bond_dims
  for(int i = 0; i < length; ++i){
    bond_dims[i] = M[i][0].rows();
  }
}
template void TensorTrain<double>::svd(double tol, bool dry_run);
template void TensorTrain< std::complex<double> >::svd(double tol, bool dry_run);

#endif
