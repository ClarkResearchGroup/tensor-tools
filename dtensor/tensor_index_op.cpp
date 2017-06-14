#ifndef INDEX_OP_CLASS_FOR_DENSE_TENSOR
#define INDEX_OP_CLASS_FOR_DENSE_TENSOR

#include "tensor_index_op.h"

void index_sets_union(const vector<tensor_index>& Ac, const vector<tensor_index>& Bc, vector<tensor_index>& res){
  res.resize(Ac.size()+Bc.size());
  vector<tensor_index>::iterator it;
  vector<tensor_index> A(Ac.begin(), Ac.end());
  vector<tensor_index> B(Bc.begin(), Bc.end());
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  // make sure indices are unique
  it = std::unique (A.begin(), A.end());
  A.resize(std::distance(A.begin(),it));
  it = std::unique (B.begin(), B.end());
  B.resize(std::distance(B.begin(),it));
  it = std::set_union(A.begin(), A.end(), B.begin(), B.end(), res.begin());
  res.resize(it-res.begin());
}

void index_sets_difference(const vector<tensor_index>& Ac, const vector<tensor_index>& Bc, vector<tensor_index>& res){
  vector<tensor_index>::iterator it;
  vector<tensor_index> A(Ac.begin(), Ac.end());
  vector<tensor_index> B(Bc.begin(), Bc.end());
  // make sure indices are unique
  it = std::unique (A.begin(), A.end());
  A.resize(std::distance(A.begin(),it));
  it = std::unique (B.begin(), B.end());
  B.resize(std::distance(B.begin(),it));
  res.clear();
  for (size_t i = 0; i < A.size(); i++) {
    bool duplicated = false;
    for (size_t j = 0; j < B.size(); j++) {
      if(A[i]==B[j]){
        duplicated = true;
        break;
      }
    }
    if(!duplicated){
       res.push_back(A[i]);
    }
  }
  for (size_t i = 0; i < B.size(); i++) {
    bool duplicated = false;
    for (size_t j = 0; j < A.size(); j++) {
      if(B[i]==A[j]){
        duplicated = true;
        break;
      }
    }
    if(!duplicated){
       res.push_back(B[i]);
    }
  }
}

void index_sets_intersection(const vector<tensor_index>& Ac, const vector<tensor_index>& Bc, vector<tensor_index>& res){
  res.resize(Ac.size()+Bc.size());
  vector<tensor_index>::iterator it;
  vector<tensor_index> A(Ac.begin(), Ac.end());
  vector<tensor_index> B(Bc.begin(), Bc.end());
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  // make sure indices are unique
  it = std::unique (A.begin(), A.end());
  A.resize(std::distance(A.begin(),it));
  it = std::unique (B.begin(), B.end());
  B.resize(std::distance(B.begin(),it));
  it = std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), res.begin());
  res.resize(it-res.begin());
}

#endif
