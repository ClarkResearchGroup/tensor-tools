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

#ifndef DENSE_TENSOR_OPERATIONS
#define DENSE_TENSOR_OPERATIONS

#include "dtensor_op.h"

template <typename T>
void qr(dtensor<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& Q, dtensor<T>& R)
{
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
  QR(r, c, A._T.data(), Q._T.data(), R._T.data());
}
template void qr(dtensor<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& Q, dtensor<double>& R);
template void qr(dtensor< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& Q, dtensor< std::complex<double> >& R);

template <typename T>
void qr(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& Q, dtensor<T>& R)
{
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
  QR(r, c, A._T.data(), Q._T.data(), R._T.data());
}
template void qr(dtensor_view<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& Q, dtensor<double>& R);
template void qr(dtensor_view< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& Q, dtensor< std::complex<double> >& R);


template <typename T>
void svd(dtensor<T>& A,
         vector<dtensor_index>& left,
         vector<dtensor_index>& right,
         dtensor<T>& U, dtensor<T>& V, vector<double>& S,
         int direction)
{
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
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L');
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R');
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }
}
template void svd(dtensor<double>& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction);
template void svd(dtensor< std::complex<double> >& A,vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction);


template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction)
{
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
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'L');
    V = U; V.dag(); V.conj(); V.idx_set.back().prime(); V = std::move(V*A); V.idx_set[0].prime(-1);
  }
  else if(direction==MoveFromRight){
    V.reset(V_idx_set);
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R');
    U = V; U.dag(); U.conj(); U.idx_set[0].prime(); U = std::move(A*U); U.idx_set.back().prime(-1);
  }
}
template void svd(dtensor_view<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction);
template void svd(dtensor_view< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction);


template <typename T>
void svd(dtensor<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff)
{
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
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R', cutoff);
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

template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff)
{
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
    SVD(r, c, A._T.data(), U._T.data(), S, V._T.data(), 'R', cutoff);
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
template void svd(dtensor_view< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff);

template <typename T>
void svd(dtensor<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
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
  }
}
template void svd(dtensor<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd(dtensor< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff, long unsigned K);

template <typename T>
void svd(dtensor_view<T>& A,
        vector<dtensor_index>& left,
        vector<dtensor_index>& right,
        dtensor<T>& U, dtensor<T>& V, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
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
  }
}
template void svd(dtensor_view<double>& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor<double>& U, dtensor<double>& V, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd(dtensor_view< std::complex<double> >& A, vector<dtensor_index>& left, vector<dtensor_index>& right, dtensor< std::complex<double> >& U, dtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff, long unsigned K);


template <typename T>
void svd_bond(dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction)
{
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
template void svd_bond(dtensor< std::complex<double> >& A_left, dtensor< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction);


template <typename T>
void svd_bond(dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction)
{
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
template void svd_bond(dtensor_view< std::complex<double> >& A_left, dtensor_view< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction);



template <typename T>
void svd_bond(dtensor<T>& A_left, dtensor<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
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


template <typename T>
void svd_bond(dtensor_view<T>& A_left, dtensor_view<T>& A_right,
        dtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
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
template void svd_bond(dtensor_view< std::complex<double> >& A_left, dtensor_view< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);


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
void svd_bond(dtensor<T>& combined, dtensor_view<T>& A_left, dtensor_view<T>& A_right,
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
template void svd_bond(dtensor<double>& combined, dtensor_view<double>& A_left, dtensor_view<double>& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(dtensor< std::complex<double> >& combined, dtensor_view< std::complex<double> >& A_left, dtensor_view< std::complex<double> >& A_right, dtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);


#endif
