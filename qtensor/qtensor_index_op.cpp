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

#ifndef INDEX_OP_CLASS_FOR_QUANTUM_NUMBERED_TENSOR
#define INDEX_OP_CLASS_FOR_QUANTUM_NUMBERED_TENSOR

#include "qtensor_index_op.h"

void find_index_permutation(vector<qtensor_index>& from, vector<qtensor_index>& to, vector<unsigned>& perm){
  assert(from.size()==to.size());
  if(perm.size()>0) perm.clear();
  // brute force O(n^2)
  int n = to.size();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if(to[i] == from[j]){
        perm.push_back(j);
        break;
      }
    }
  }
  assert(int(perm.size())==n);
}

void index_sets_union(const vector<qtensor_index>& Ac, const vector<qtensor_index>& Bc, vector<qtensor_index>& res){
  res.resize(Ac.size()+Bc.size());
  vector<qtensor_index>::iterator it;
  vector<qtensor_index> A(Ac.begin(), Ac.end());
  vector<qtensor_index> B(Bc.begin(), Bc.end());
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

void index_sets_difference(const vector<qtensor_index>& Ac, const vector<qtensor_index>& Bc, vector<qtensor_index>& res){
  vector<qtensor_index>::iterator it;
  vector<qtensor_index> A(Ac.begin(), Ac.end());
  vector<qtensor_index> B(Bc.begin(), Bc.end());
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

void index_sets_intersection(const vector<qtensor_index>& Ac, const vector<qtensor_index>& Bc, vector<qtensor_index>& res){
  res.resize(Ac.size()+Bc.size());
  vector<qtensor_index>::iterator it;
  vector<qtensor_index> A(Ac.begin(), Ac.end());
  vector<qtensor_index> B(Bc.begin(), Bc.end());
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
