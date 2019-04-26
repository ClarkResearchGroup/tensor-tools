#ifndef DENSE_TENSOR_OPERATIONS
#define DENSE_TENSOR_OPERATIONS

#include "dtensor_op.h"

template <typename T>
void svd(dtensor<T>& A,
         vector<dtensor_index>& left,
         vector<dtensor_index>& right,
         dtensor<T>& U, dtensor<T>& V, vector<double>& S,
         int direction)
{
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
  
  A.__T[indA.c_str()].svd(U.__T[indU.c_str()],_S[indS.c_str()],V.__T[indV.c_str()]); //,R+3); //Does this work independent of what's in U and V
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
  S.resize(np);
  //S = std::move(std::vector<double>(np,reinterpret_cast<double*>(data)));
  std::copy(S.begin(),S.end(),data);
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
    cerr<<endl;*/


    /*cerr<<"Length U: "<<i<<" is  "<<U.__T.lens[j]<<endl;
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

  exit(1);

  exit(1);
  

  // Permute dtensor
  /*unsigned r=1, c=1;
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
  for (auto i : perm)
    cerr<<i<<endl;
  exit(1);
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
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L');
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R');
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }*/
  assert(1==2);
}
template void svd(dtensor<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction);
template void svd(dtensor< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction);



template <typename T>
void svd_bond(dtensor<T>& combined, dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
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
  //  assert(1==2);
  //cerr<<"Pre SVD: "<<endl;
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
  S.resize(np);
  //S = std::move(std::vector<double>(np,reinterpret_cast<double*>(data)));
  std::copy(S.begin(),S.end(),data);
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

  exit(1);

  exit(1);
  

  // Permute dtensor
  /*unsigned r=1, c=1;
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
  for (auto i : perm)
    cerr<<i<<endl;
  exit(1);
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
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L');
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R');
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }*/
  assert(1==2);
}


template void svd(dtensor<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd(dtensor< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff, long unsigned K);



#endif
