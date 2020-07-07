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

#ifndef QUANTUM_NUMBERED_TENSOR_OPERATIONS
#define QUANTUM_NUMBERED_TENSOR_OPERATIONS

#include "qtensor_op.h"

template <typename T>
void qr(qtensor<T>& A,
        vector<qtensor_index>& left,
        vector<qtensor_index>& right,
        qtensor<T>& Q, qtensor<T>& R)
{
  // Permute qtensor, group left indices, and right indices
  vector<qtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
  }
  for(auto v : right){
    new_idx_set.push_back(v);
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Accumulate legal quantum numbers for the mid bond
  set<int> mid_Q_set;
  unordered_map<int, set<uint_vec> > l_index_qi;
  unordered_map<int, set<uint_vec> > r_index_qi;
  for (size_t i = 0; i < A.block.size(); i++) {
    int mid_QN = 0;
    uint_vec l_qi, r_qi;
    for (size_t j = 0; j < left.size(); j++) {
      if(left[j].arrow()==Inward){
        mid_QN += A.block_index_qn[i][j];
      }else{
        mid_QN -= A.block_index_qn[i][j];
      }
      l_qi.push_back(A.block_index_qi[i][j]);
    }
    for (size_t j = 0; j < right.size(); j++) {
      r_qi.push_back(A.block_index_qi[i][j+left.size()]);
    }
    mid_Q_set.insert(mid_QN);
    l_index_qi[mid_QN].insert(l_qi);
    r_index_qi[mid_QN].insert(r_qi);
  }
  // Set up the initial mid bond
  unordered_map<int, str_vec>  l_qn_str_map;
  unordered_map<int, str_vec>  r_qn_str_map;
  unordered_map<int, uint_vec> l_qn_sizes_map;
  unordered_map<int, uint_vec> r_qn_sizes_map;
  vector<int> mid_Q(mid_Q_set.begin(), mid_Q_set.end());
  uint_vec mid_QDim(mid_Q.size());
  for (size_t i = 0; i < mid_Q.size(); i++) {
    int q = mid_Q[i];
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    unsigned row = 0;
    unsigned col = 0;
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      unsigned s1 = 1;
      string l_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        s1 *= left[j].qdim(l_qi[j]);
        l_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
      }
      l_qn_str_map[q].push_back(l_qn_str);
      l_qn_sizes_map[q].push_back(s1);
      row += s1;
    }
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      unsigned s1 = 1;
      string r_qn_str;
      for (size_t j = 0; j < r_qi.size(); j++) {
        s1 *= right[j].qdim(r_qi[j]);
        r_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
      }
      r_qn_str_map[q].push_back(r_qn_str);
      r_qn_sizes_map[q].push_back(s1);
      col += s1;
    }
    mid_QDim[i] = std::min(row,col);
  }
  qtensor_index mid(Outward);
  mid.addQN(mid_Q, mid_QDim);
  // Set up Q and R
  vector<qtensor_index> Q_idx_set(left);
  Q_idx_set.push_back(mid);
  vector<qtensor_index> R_idx_set;
  mid.dag();
  R_idx_set.push_back(mid);
  mid.dag();
  for(auto v : right){
    R_idx_set.push_back(v);
  }
  Q.reset(Q_idx_set);
  Q._initted = true;
  R.reset(R_idx_set);
  R._initted = true;
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim = mid_QDim[ii];
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      uint_vec Q_block_index_qi;
      uint_vec Q_block_index_qd;
      int_vec  Q_block_index_qn;
      unsigned Q_block_size = 1;
      string   Q_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        Q_block_index_qi.push_back(l_qi[j]);
        Q_block_index_qd.push_back(left[j].qdim(l_qi[j]));
        Q_block_index_qn.push_back(left[j].qn(l_qi[j]));
        Q_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
        Q_block_size *= Q_block_index_qd.back();
      }
      // last bond
      Q_block_index_qi.push_back(ii);
      Q_block_index_qd.push_back(dim);
      Q_block_index_qn.push_back(q);
      Q_qn_str += (to_string(q)+" ");
      Q_block_size *= Q_block_index_qd.back();
      // build Q
      Q.block.emplace_back(Q_block_size);
      Q.block_index_qn.push_back(Q_block_index_qn);
      Q.block_index_qd.push_back(Q_block_index_qd);
      Q.block_index_qi.push_back(Q_block_index_qi);
      Q.block_id_by_qn_str[Q_qn_str] = Q.block.size()-1;
    }
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      uint_vec R_block_index_qi;
      uint_vec R_block_index_qd;
      int_vec  R_block_index_qn;
      unsigned R_block_size = 1;
      string   R_qn_str;
      // first bond
      R_block_index_qi.push_back(ii);
      R_block_index_qd.push_back(dim);
      R_block_index_qn.push_back(q);
      R_qn_str += (to_string(q)+" ");
      R_block_size *= R_block_index_qd.back();
      for (size_t j = 0; j < r_qi.size(); j++) {
        R_block_index_qi.push_back(r_qi[j]);
        R_block_index_qd.push_back(right[j].qdim(r_qi[j]));
        R_block_index_qn.push_back(right[j].qn(r_qi[j]));
        R_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
        R_block_size *= R_block_index_qd.back();
      }
      // build Q
      R.block.emplace_back(R_block_size);
      R.block_index_qn.push_back(R_block_index_qn);
      R.block_index_qd.push_back(R_block_index_qd);
      R.block_index_qi.push_back(R_block_index_qi);
      R.block_id_by_qn_str[R_qn_str] = R.block.size()-1;
    }
  }
  // QR block by block
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    // Get big block size
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    unsigned dim = mid_QDim[ii];
    vector<T> mA(row*col);
    vector<T> mQ(row*dim);
    vector<T> mR(dim*col);
    // copy data from qtensor blocks to vector mA
    unsigned c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      unsigned c_col = 0;
      for (size_t j = 0; j < r_qn_str_map[q].size(); j++) {
        string qn_str = l_qn_str_map[q][i] + r_qn_str_map[q][j];
        unsigned A_block = A.block_id_by_qn_str[qn_str];
        for (size_t k = 0; k < A.block[A_block].size(); k++) {
          unsigned b_row = k%l_qn_sizes_map[q][i];
          unsigned b_col = k/l_qn_sizes_map[q][i];
          mA[(c_row+b_row) + (c_col+b_col)*row] = A.block[A_block][k];
        }
        c_col += r_qn_sizes_map[q][j];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    // QR
    QR(row, col, mA.data(), mQ.data(), mR.data());
    // copy data from mQ to qtensor Q
    c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      string Q_qn_str = l_qn_str_map[q][i] + to_string(q) + " ";
      unsigned Q_block = Q.block_id_by_qn_str[Q_qn_str];
      for (size_t j = 0; j < Q.block[Q_block].size(); j++) {
        unsigned b_row = j%l_qn_sizes_map[q][i];
        unsigned b_col = j/l_qn_sizes_map[q][i];
        Q.block[Q_block][j] = mQ[(c_row+b_row) + b_col*row];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    // copy data from mR to qtensor R
    unsigned c_col = 0;
    for (size_t i = 0; i < r_qn_str_map[q].size(); i++) {
      string R_qn_str = to_string(q) + " " + r_qn_str_map[q][i];
      unsigned R_block = R.block_id_by_qn_str[R_qn_str];
      for (size_t j = 0; j < R.block[R_block].size(); j++) {
        unsigned b_row = j%dim;
        unsigned b_col = j/dim;
        R.block[R_block][j] = mR[b_row + (b_col+c_col)*dim];
      }
      c_col += r_qn_sizes_map[q][i];
    }
  }
}
template void qr(qtensor<double>& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor<double>& Q, qtensor<double>& R);
template void qr(qtensor< std::complex<double> >& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor< std::complex<double> >& Q, qtensor< std::complex<double> >& R);


template <typename T>
void svd(qtensor<T>& A,
         vector<qtensor_index>& left,
         vector<qtensor_index>& right,
         qtensor<T>& U, qtensor<T>& V, vector<double>& S,
         int direction)
{
  S.clear();
  // Permute qtensor, group left indices, and right indices
  vector<qtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
  }
  for(auto v : right){
    new_idx_set.push_back(v);
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Accumulate legal quantum numbers for the mid bond
  set<int> mid_Q_set;
  unordered_map<int, set<uint_vec> > l_index_qi;
  unordered_map<int, set<uint_vec> > r_index_qi;
  for (size_t i = 0; i < A.block.size(); i++) {
    int mid_QN = 0;
    uint_vec l_qi, r_qi;
    for (size_t j = 0; j < left.size(); j++) {
      if(left[j].arrow()==Inward){
        mid_QN += A.block_index_qn[i][j];
      }else{
        mid_QN -= A.block_index_qn[i][j];
      }
      l_qi.push_back(A.block_index_qi[i][j]);
    }
    for (size_t j = 0; j < right.size(); j++) {
      r_qi.push_back(A.block_index_qi[i][j+left.size()]);
    }
    mid_Q_set.insert(mid_QN);
    l_index_qi[mid_QN].insert(l_qi);
    r_index_qi[mid_QN].insert(r_qi);
  }
  // Set up the initial mid bond
  unordered_map<int, str_vec>  l_qn_str_map;
  unordered_map<int, str_vec>  r_qn_str_map;
  unordered_map<int, uint_vec> l_qn_sizes_map;
  unordered_map<int, uint_vec> r_qn_sizes_map;
  vector<int> mid_Q(mid_Q_set.begin(), mid_Q_set.end());
  uint_vec mid_QDim(mid_Q.size());
  for (size_t i = 0; i < mid_Q.size(); i++) {
    int q = mid_Q[i];
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    unsigned row = 0;
    unsigned col = 0;
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      unsigned s1 = 1;
      string l_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        s1 *= left[j].qdim(l_qi[j]);
        l_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
      }
      l_qn_str_map[q].push_back(l_qn_str);
      l_qn_sizes_map[q].push_back(s1);
      row += s1;
    }
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      unsigned s1 = 1;
      string r_qn_str;
      for (size_t j = 0; j < r_qi.size(); j++) {
        s1 *= right[j].qdim(r_qi[j]);
        r_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
      }
      r_qn_str_map[q].push_back(r_qn_str);
      r_qn_sizes_map[q].push_back(s1);
      col += s1;
    }
    mid_QDim[i] = std::min(row,col);
  }
  qtensor_index mid(Outward);
  mid.addQN(mid_Q, mid_QDim);
  // Set up U and V index
  vector<qtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<qtensor_index> V_idx_set;
  mid.dag();
  V_idx_set.push_back(mid);
  mid.dag();
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  U.reset(U_idx_set);
  U._initted = true;
  V.reset(V_idx_set);
  V._initted = true;
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim = mid_QDim[ii];
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      uint_vec U_block_index_qi;
      uint_vec U_block_index_qd;
      int_vec  U_block_index_qn;
      unsigned U_block_size = 1;
      string   U_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        U_block_index_qi.push_back(l_qi[j]);
        U_block_index_qd.push_back(left[j].qdim(l_qi[j]));
        U_block_index_qn.push_back(left[j].qn(l_qi[j]));
        U_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
        U_block_size *= U_block_index_qd.back();
      }
      // last bond
      U_block_index_qi.push_back(ii);
      U_block_index_qd.push_back(dim);
      U_block_index_qn.push_back(q);
      U_qn_str += (to_string(q)+" ");
      U_block_size *= U_block_index_qd.back();
      // build U
      U.block.emplace_back(U_block_size);
      U.block_index_qn.push_back(U_block_index_qn);
      U.block_index_qd.push_back(U_block_index_qd);
      U.block_index_qi.push_back(U_block_index_qi);
      U.block_id_by_qn_str[U_qn_str] = U.block.size()-1;
    }
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      uint_vec V_block_index_qi;
      uint_vec V_block_index_qd;
      int_vec  V_block_index_qn;
      unsigned V_block_size = 1;
      string   V_qn_str;
      // first bond
      V_block_index_qi.push_back(ii);
      V_block_index_qd.push_back(dim);
      V_block_index_qn.push_back(q);
      V_qn_str += (to_string(q)+" ");
      V_block_size *= V_block_index_qd.back();
      for (size_t j = 0; j < r_qi.size(); j++) {
        V_block_index_qi.push_back(r_qi[j]);
        V_block_index_qd.push_back(right[j].qdim(r_qi[j]));
        V_block_index_qn.push_back(right[j].qn(r_qi[j]));
        V_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
        V_block_size *= V_block_index_qd.back();
      }
      // build Q
      V.block.emplace_back(V_block_size);
      V.block_index_qn.push_back(V_block_index_qn);
      V.block_index_qd.push_back(V_block_index_qd);
      V.block_index_qi.push_back(V_block_index_qi);
      V.block_id_by_qn_str[V_qn_str] = V.block.size()-1;
    }
  }
  // SVD block by block
  vector< vector<double> > SQ(mid_Q.size());
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    // Get big block size
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    unsigned dim = mid_QDim[ii];
    SQ[ii].resize(dim);
    vector<T> mA(row*col);
    vector<T> mU(row*dim);
    vector<T> mV(dim*col);
    // copy data from qtensor blocks to vector mA
    unsigned c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      unsigned c_col = 0;
      for (size_t j = 0; j < r_qn_str_map[q].size(); j++) {
        string qn_str = l_qn_str_map[q][i] + r_qn_str_map[q][j];
        unsigned A_block = A.block_id_by_qn_str[qn_str];
        for (size_t k = 0; k < A.block[A_block].size(); k++) {
          unsigned b_row = k%l_qn_sizes_map[q][i];
          unsigned b_col = k/l_qn_sizes_map[q][i];
          mA[(c_row+b_row) + (c_col+b_col)*row] = A.block[A_block][k];
        }
        c_col += r_qn_sizes_map[q][j];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    // SVD
    SVD(row, col, mA.data(), mU.data(), SQ[ii], mV.data());
    if(direction==MoveFromLeft){
      for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < col; j++) {
          mV[i+j*dim] *= SQ[ii][i];
        }
      }
    }else if(direction==MoveFromRight){
      for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < dim; j++) {
          mU[i+j*row] *= SQ[ii][j];
        }
      }
    }
    // copy data from mU to qtensor U
    c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      string U_qn_str = l_qn_str_map[q][i] + to_string(q) + " ";
      unsigned U_block = U.block_id_by_qn_str[U_qn_str];
      for (size_t j = 0; j < U.block[U_block].size(); j++) {
        unsigned b_row = j%l_qn_sizes_map[q][i];
        unsigned b_col = j/l_qn_sizes_map[q][i];
        U.block[U_block][j] = mU[(c_row+b_row) + b_col*row];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    // copy data from mU to atensor V
    unsigned c_col = 0;
    for (size_t i = 0; i < r_qn_str_map[q].size(); i++) {
      string V_qn_str = to_string(q) + " " + r_qn_str_map[q][i];
      unsigned V_block = V.block_id_by_qn_str[V_qn_str];
      for (size_t j = 0; j < V.block[V_block].size(); j++) {
        unsigned b_row = j%dim;
        unsigned b_col = j/dim;
        V.block[V_block][j] = mV[b_row + (b_col+c_col)*dim];
      }
      c_col += r_qn_sizes_map[q][i];
    }
  }
  // Collect singular values
  for(auto Sv : SQ){
    for ( auto v : Sv){
      S.push_back(v);
    }
  }
}
template void svd(qtensor<double>& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor<double>& U, qtensor<double>& V, vector<double>& S, int direction);
template void svd(qtensor< std::complex<double> >& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor< std::complex<double> >& U, qtensor< std::complex<double> >& V, vector<double>& S, int direction);


template <typename T>
void svd(qtensor<T>& A,
        vector<qtensor_index>& left,
        vector<qtensor_index>& right,
        qtensor<T>& U, qtensor<T>& V, vector<double>& S,
        int direction, double cutoff)
{
  S.clear();
  // Permute qtensor, group left indices, and right indices
  vector<qtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
  }
  for(auto v : right){
    new_idx_set.push_back(v);
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Accumulate legal quantum numbers for the mid bond
  set<int> mid_Q_set;
  unordered_map<int, set<uint_vec> > l_index_qi;
  unordered_map<int, set<uint_vec> > r_index_qi;
  for (size_t i = 0; i < A.block.size(); i++) {
    int mid_QN = 0;
    uint_vec l_qi, r_qi;
    for (size_t j = 0; j < left.size(); j++) {
      if(left[j].arrow()==Inward){
        mid_QN += A.block_index_qn[i][j];
      }else{
        mid_QN -= A.block_index_qn[i][j];
      }
      l_qi.push_back(A.block_index_qi[i][j]);
    }
    for (size_t j = 0; j < right.size(); j++) {
      r_qi.push_back(A.block_index_qi[i][j+left.size()]);
    }
    mid_Q_set.insert(mid_QN);
    l_index_qi[mid_QN].insert(l_qi);
    r_index_qi[mid_QN].insert(r_qi);
  }
  // Set up the initial mid bond
  unordered_map<int, str_vec>  l_qn_str_map;
  unordered_map<int, str_vec>  r_qn_str_map;
  unordered_map<int, uint_vec> l_qn_sizes_map;
  unordered_map<int, uint_vec> r_qn_sizes_map;
  vector<int> mid_Q(mid_Q_set.begin(), mid_Q_set.end());
  uint_vec mid_QDim(mid_Q.size());
  for (size_t i = 0; i < mid_Q.size(); i++) {
    int q = mid_Q[i];
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    unsigned row = 0;
    unsigned col = 0;
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      unsigned s1 = 1;
      string l_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        s1 *= left[j].qdim(l_qi[j]);
        l_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
      }
      l_qn_str_map[q].push_back(l_qn_str);
      l_qn_sizes_map[q].push_back(s1);
      row += s1;
    }
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      unsigned s1 = 1;
      string r_qn_str;
      for (size_t j = 0; j < r_qi.size(); j++) {
        s1 *= right[j].qdim(r_qi[j]);
        r_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
      }
      r_qn_str_map[q].push_back(r_qn_str);
      r_qn_sizes_map[q].push_back(s1);
      col += s1;
    }
    mid_QDim[i] = std::min(row,col);
  }
  // SVD block by block
  vector< vector<double> > SQ(mid_Q.size());
  uint_vec new_QDim(mid_QDim.size());
  vector< vector<T> > mA(mid_Q.size());
  vector< vector<T> > mU(mid_Q.size());
  vector< vector<T> > mV(mid_Q.size());
  // SVD block by block
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim = mid_QDim[ii];
    // Get big block size
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    SQ[ii].resize(dim);
    mA[ii].resize(row*col);
    mU[ii].resize(row*dim);
    mV[ii].resize(dim*col);
    // copy data from qtensor blocks to vector mA
    unsigned c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      unsigned c_col = 0;
      for (size_t j = 0; j < r_qn_str_map[q].size(); j++) {
        string qn_str = l_qn_str_map[q][i] + r_qn_str_map[q][j];
        unsigned A_block = A.block_id_by_qn_str[qn_str];
        for (size_t k = 0; k < A.block[A_block].size(); k++) {
          unsigned b_row = k%l_qn_sizes_map[q][i];
          unsigned b_col = k/l_qn_sizes_map[q][i];
          mA[ii][(c_row+b_row) + (c_col+b_col)*row] = A.block[A_block][k];
        }
        c_col += r_qn_sizes_map[q][j];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    // SVD
    SVD(row, col, mA[ii].data(), mU[ii].data(), SQ[ii], mV[ii].data());
    if(direction==MoveFromLeft){
      for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < col; j++) {
          mV[ii][i+j*dim] *= SQ[ii][i];
        }
      }
    }else if(direction==MoveFromRight){
      for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < dim; j++) {
          mU[ii][i+j*row] *= SQ[ii][j];
        }
      }
    }
  }
  // Tune dimensions of the mid bond
  double nm = 0.0;
  vector< pair<double, unsigned> > S_and_block_idx;
  for (size_t i = 0; i < SQ.size(); i++) {
    for (size_t j = 0; j < SQ[i].size(); j++) {
      S_and_block_idx.push_back(std::make_pair(SQ[i][j], i));
      nm += SQ[i][j] * SQ[i][j];
    }
  }
  std::sort(S_and_block_idx.begin(), S_and_block_idx.end(), pairCompare);
  double cumsum = 0;
  long unsigned accBD = 0;
  for (size_t i = 0; i < S_and_block_idx.size(); i++) {
    cumsum += S_and_block_idx[i].first * S_and_block_idx[i].first / nm;
    new_QDim[S_and_block_idx[i].second] += 1;
    accBD += 1;
    if(1-cumsum<cutoff) break;
  }
  // set up mid bond
  qtensor_index mid(Outward);
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    if(new_QDim[ii]>0) mid.addQN(std::make_pair(mid_Q[ii], new_QDim[ii]));
  }
  // Set up U and V index
  vector<qtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<qtensor_index> V_idx_set;
  mid.dag();
  V_idx_set.push_back(mid);
  mid.dag();
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  U.reset(U_idx_set);
  U._initted = true;
  V.reset(V_idx_set);
  V._initted = true;
  unsigned mid_num = 0;
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim  = mid_QDim[ii];
    unsigned ndim = new_QDim[ii];
    if(ndim==0) continue;
    // initialize blocks
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      uint_vec U_block_index_qi;
      uint_vec U_block_index_qd;
      int_vec  U_block_index_qn;
      unsigned U_block_size = 1;
      string   U_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        U_block_index_qi.push_back(l_qi[j]);
        U_block_index_qd.push_back(left[j].qdim(l_qi[j]));
        U_block_index_qn.push_back(left[j].qn(l_qi[j]));
        U_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
        U_block_size *= U_block_index_qd.back();
      }
      // last bond
      U_block_index_qi.push_back(mid_num);
      U_block_index_qd.push_back(ndim);
      U_block_index_qn.push_back(q);
      U_qn_str += (to_string(q)+" ");
      U_block_size *= U_block_index_qd.back();
      // build U
      U.block.emplace_back(U_block_size);
      U.block_index_qn.push_back(U_block_index_qn);
      U.block_index_qd.push_back(U_block_index_qd);
      U.block_index_qi.push_back(U_block_index_qi);
      U.block_id_by_qn_str[U_qn_str] = U.block.size()-1;
    }
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      uint_vec V_block_index_qi;
      uint_vec V_block_index_qd;
      int_vec  V_block_index_qn;
      unsigned V_block_size = 1;
      string   V_qn_str;
      // first bond
      V_block_index_qi.push_back(mid_num);
      V_block_index_qd.push_back(ndim);
      V_block_index_qn.push_back(q);
      V_qn_str += (to_string(q)+" ");
      V_block_size *= V_block_index_qd.back();
      for (size_t j = 0; j < r_qi.size(); j++) {
        V_block_index_qi.push_back(r_qi[j]);
        V_block_index_qd.push_back(right[j].qdim(r_qi[j]));
        V_block_index_qn.push_back(right[j].qn(r_qi[j]));
        V_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
        V_block_size *= V_block_index_qd.back();
      }
      // build Q
      V.block.emplace_back(V_block_size);
      V.block_index_qn.push_back(V_block_index_qn);
      V.block_index_qd.push_back(V_block_index_qd);
      V.block_index_qi.push_back(V_block_index_qi);
      V.block_id_by_qn_str[V_qn_str] = V.block.size()-1;
    }
    ++mid_num;
    // copy data
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    // copy data from mU to qtensor U
    unsigned c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      string U_qn_str = l_qn_str_map[q][i] + to_string(q) + " ";
      unsigned U_block = U.block_id_by_qn_str.at(U_qn_str);
      for (size_t j = 0; j < U.block[U_block].size(); j++) {
        unsigned b_row = j%l_qn_sizes_map[q][i];
        unsigned b_col = j/l_qn_sizes_map[q][i];
        U.block[U_block][j] = mU[ii][(c_row+b_row) + b_col*row];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    unsigned c_col = 0;
    for (size_t i = 0; i < r_qn_str_map[q].size(); i++) {
      string V_qn_str = to_string(q) + " " + r_qn_str_map[q][i];
      unsigned V_block = V.block_id_by_qn_str.at(V_qn_str);
      for (size_t j = 0; j < V.block[V_block].size(); j++) {
        unsigned b_row = j%ndim;
        unsigned b_col = j/ndim;
        V.block[V_block][j] = mV[ii][b_row + (b_col+c_col)*dim];
      }
      c_col += r_qn_sizes_map[q][i];
    }
  }
  // Collect singular values
  for (size_t i = 0; i < SQ.size(); i++) {
    for (size_t j = 0; j < new_QDim[i]; j++) {
      S.push_back(SQ[i][j]);
    }
  }
}
template void svd(qtensor<double>& A, vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor<double>& U, qtensor<double>& V, vector<double>& S, int direction, double cutoff);
template void svd(qtensor< std::complex<double> >& A, vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor< std::complex<double> >& U, qtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff);


template <typename T>
void svd(qtensor<T>& A,
        vector<qtensor_index>& left,
        vector<qtensor_index>& right,
        qtensor<T>& U, qtensor<T>& V, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  S.clear();
  // Permute qtensor, group left indices, and right indices
  vector<qtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
  }
  for(auto v : right){
    new_idx_set.push_back(v);
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  A.permute(perm);
  // Accumulate legal quantum numbers for the mid bond
  set<int> mid_Q_set;
  unordered_map<int, set<uint_vec> > l_index_qi;
  unordered_map<int, set<uint_vec> > r_index_qi;
  for (size_t i = 0; i < A.block.size(); i++) {
    int mid_QN = 0;
    uint_vec l_qi, r_qi;
    for (size_t j = 0; j < left.size(); j++) {
      if(left[j].arrow()==Inward){
        mid_QN += A.block_index_qn[i][j];
      }else{
        mid_QN -= A.block_index_qn[i][j];
      }
      l_qi.push_back(A.block_index_qi[i][j]);
    }
    for (size_t j = 0; j < right.size(); j++) {
      r_qi.push_back(A.block_index_qi[i][j+left.size()]);
    }
    mid_Q_set.insert(mid_QN);
    l_index_qi[mid_QN].insert(l_qi);
    r_index_qi[mid_QN].insert(r_qi);
  }
  // Set up the initial mid bond
  unordered_map<int, str_vec>  l_qn_str_map;
  unordered_map<int, str_vec>  r_qn_str_map;
  unordered_map<int, uint_vec> l_qn_sizes_map;
  unordered_map<int, uint_vec> r_qn_sizes_map;
  vector<int> mid_Q(mid_Q_set.begin(), mid_Q_set.end());
  uint_vec mid_QDim(mid_Q.size());
  for (size_t i = 0; i < mid_Q.size(); i++) {
    int q = mid_Q[i];
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    unsigned row = 0;
    unsigned col = 0;
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      unsigned s1 = 1;
      string l_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        s1 *= left[j].qdim(l_qi[j]);
        l_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
      }
      l_qn_str_map[q].push_back(l_qn_str);
      l_qn_sizes_map[q].push_back(s1);
      row += s1;
    }
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      unsigned s1 = 1;
      string r_qn_str;
      for (size_t j = 0; j < r_qi.size(); j++) {
        s1 *= right[j].qdim(r_qi[j]);
        r_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
      }
      r_qn_str_map[q].push_back(r_qn_str);
      r_qn_sizes_map[q].push_back(s1);
      col += s1;
    }
    mid_QDim[i] = std::min(row,col);
  }
  // SVD block by block
  vector< vector<double> > SQ(mid_Q.size());
  uint_vec new_QDim(mid_QDim.size());
  vector< vector<T> > mA(mid_Q.size());
  vector< vector<T> > mU(mid_Q.size());
  vector< vector<T> > mV(mid_Q.size());
  // SVD block by block
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim = mid_QDim[ii];
    // Get big block size
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    SQ[ii].resize(dim);
    mA[ii].resize(row*col);
    mU[ii].resize(row*dim);
    mV[ii].resize(dim*col);
    // copy data from qtensor blocks to vector mA
    unsigned c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      unsigned c_col = 0;
      for (size_t j = 0; j < r_qn_str_map[q].size(); j++) {
        string qn_str = l_qn_str_map[q][i] + r_qn_str_map[q][j];
        unsigned A_block = A.block_id_by_qn_str[qn_str];
        for (size_t k = 0; k < A.block[A_block].size(); k++) {
          unsigned b_row = k%l_qn_sizes_map[q][i];
          unsigned b_col = k/l_qn_sizes_map[q][i];
          mA[ii][(c_row+b_row) + (c_col+b_col)*row] = A.block[A_block][k];
        }
        c_col += r_qn_sizes_map[q][j];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    // normalize A
    // double A_nm = 0;
    // for (size_t i = 0; i < mA[ii].size(); i++) {
    //   A_nm += std::real(cconj(mA[ii][i])*mA[ii][i]);
    // }
    // A_nm = std::sqrt(A_nm);
    // for (size_t i = 0; i < mA[ii].size(); i++) {
    //   mA[ii][i] *= A_nm;
    // }
    // SVD
    SVD(row, col, mA[ii].data(), mU[ii].data(), SQ[ii], mV[ii].data());
    // reset singular values norms
    // for (size_t i = 0; i < dim; i++) {
    //   SQ[ii][i] /= A_nm;
    // }
    // merge S into V or U
    if(direction==MoveFromLeft){
      for (size_t i = 0; i < dim; i++) {
        for (size_t j = 0; j < col; j++) {
          mV[ii][i+j*dim] *= SQ[ii][i];
        }
      }
    }else if(direction==MoveFromRight){
      for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < dim; j++) {
          mU[ii][i+j*row] *= SQ[ii][j];
        }
      }
    }
  }
  // Tune dimensions of the mid bond
  double nm = 0.0;
  vector< pair<double, unsigned> > S_and_block_idx;
  for (size_t i = 0; i < SQ.size(); i++) {
    for (size_t j = 0; j < SQ[i].size(); j++) {
      S_and_block_idx.push_back(std::make_pair(SQ[i][j], i));
      nm += SQ[i][j] * SQ[i][j];
    }
  }
  std::sort(S_and_block_idx.begin(), S_and_block_idx.end(), pairCompare);
  double cumsum = 0;
  long unsigned accBD = 0;
  for (size_t i = 0; i < S_and_block_idx.size(); i++) {
    cumsum += S_and_block_idx[i].first * S_and_block_idx[i].first / nm;
    new_QDim[S_and_block_idx[i].second] += 1;
    accBD += 1;
    if(1-cumsum<cutoff || accBD>=K) break;
  }
  // set up mid bond
  qtensor_index mid(Outward);
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    if(new_QDim[ii]>0) mid.addQN(std::make_pair(mid_Q[ii], new_QDim[ii]));
  }
  // Set up U and V index
  vector<qtensor_index> U_idx_set(left);
  U_idx_set.push_back(mid);
  vector<qtensor_index> V_idx_set;
  mid.dag();
  V_idx_set.push_back(mid);
  mid.dag();
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  U.reset(U_idx_set);
  U._initted = true;
  V.reset(V_idx_set);
  V._initted = true;
  unsigned mid_num = 0;
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim  = mid_QDim[ii];
    unsigned ndim = new_QDim[ii];
    if(ndim==0) continue;
    // initialize blocks
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      uint_vec U_block_index_qi;
      uint_vec U_block_index_qd;
      int_vec  U_block_index_qn;
      unsigned U_block_size = 1;
      string   U_qn_str;
      for (size_t j = 0; j < l_qi.size(); j++) {
        U_block_index_qi.push_back(l_qi[j]);
        U_block_index_qd.push_back(left[j].qdim(l_qi[j]));
        U_block_index_qn.push_back(left[j].qn(l_qi[j]));
        U_qn_str += (to_string(left[j].qn(l_qi[j]))+" ");
        U_block_size *= U_block_index_qd.back();
      }
      // last bond
      U_block_index_qi.push_back(mid_num);
      U_block_index_qd.push_back(ndim);
      U_block_index_qn.push_back(q);
      U_qn_str += (to_string(q)+" ");
      U_block_size *= U_block_index_qd.back();
      // build U
      U.block.emplace_back(U_block_size);
      U.block_index_qn.push_back(U_block_index_qn);
      U.block_index_qd.push_back(U_block_index_qd);
      U.block_index_qi.push_back(U_block_index_qi);
      U.block_id_by_qn_str[U_qn_str] = U.block.size()-1;
    }
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      uint_vec V_block_index_qi;
      uint_vec V_block_index_qd;
      int_vec  V_block_index_qn;
      unsigned V_block_size = 1;
      string   V_qn_str;
      // first bond
      V_block_index_qi.push_back(mid_num);
      V_block_index_qd.push_back(ndim);
      V_block_index_qn.push_back(q);
      V_qn_str += (to_string(q)+" ");
      V_block_size *= V_block_index_qd.back();
      for (size_t j = 0; j < r_qi.size(); j++) {
        V_block_index_qi.push_back(r_qi[j]);
        V_block_index_qd.push_back(right[j].qdim(r_qi[j]));
        V_block_index_qn.push_back(right[j].qn(r_qi[j]));
        V_qn_str += (to_string(right[j].qn(r_qi[j]))+" ");
        V_block_size *= V_block_index_qd.back();
      }
      // build Q
      V.block.emplace_back(V_block_size);
      V.block_index_qn.push_back(V_block_index_qn);
      V.block_index_qd.push_back(V_block_index_qd);
      V.block_index_qi.push_back(V_block_index_qi);
      V.block_id_by_qn_str[V_qn_str] = V.block.size()-1;
    }
    ++mid_num;
    // copy data
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    // copy data from mU to qtensor U
    unsigned c_row = 0;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      string U_qn_str = l_qn_str_map[q][i] + to_string(q) + " ";
      unsigned U_block = U.block_id_by_qn_str.at(U_qn_str);
      for (size_t j = 0; j < U.block[U_block].size(); j++) {
        unsigned b_row = j%l_qn_sizes_map[q][i];
        unsigned b_col = j/l_qn_sizes_map[q][i];
        U.block[U_block][j] = mU[ii][(c_row+b_row) + b_col*row];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    unsigned c_col = 0;
    for (size_t i = 0; i < r_qn_str_map[q].size(); i++) {
      string V_qn_str = to_string(q) + " " + r_qn_str_map[q][i];
      unsigned V_block = V.block_id_by_qn_str.at(V_qn_str);
      for (size_t j = 0; j < V.block[V_block].size(); j++) {
        unsigned b_row = j%ndim;
        unsigned b_col = j/ndim;
        V.block[V_block][j] = mV[ii][b_row + (b_col+c_col)*dim];
      }
      c_col += r_qn_sizes_map[q][i];
    }
  }
  // Collect singular values
  for (size_t i = 0; i < SQ.size(); i++) {
    for (size_t j = 0; j < new_QDim[i]; j++) {
      S.push_back(SQ[i][j]);
    }
  }
}
template void svd(qtensor<double>& A, vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor<double>& U, qtensor<double>& V, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd(qtensor< std::complex<double> >& A, vector<qtensor_index>& left, vector<qtensor_index>& right, qtensor< std::complex<double> >& U, qtensor< std::complex<double> >& V, vector<double>& S, int direction, double cutoff, long unsigned K);


template <typename T>
void svd_bond(qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction)
{
  qtensor<T> U,V;
  qtensor<T> combined = std::move(A_left * A_right);
  vector<qtensor_index> left;
  vector<qtensor_index> right;
  qtensor_index dag_mid(mid); dag_mid.dag();
  string tag  = mid.tag();
  string dtag = dag_mid.tag();
  // Separate qtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != dtag) {
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
  // Restore the correct name, type, and prime level
  U.idx_set.back() = mid;
  V.idx_set[0] = dag_mid;
  // Copy
  A_right = V;
  A_left = U;
}
template void svd_bond(qtensor<double>& A_left, qtensor<double>& A_right, qtensor_index& mid, vector<double>& S, int direction);
template void svd_bond(qtensor< std::complex<double> >& A_left, qtensor< std::complex<double> >& A_right, qtensor_index& mid, vector<double>& S, int direction);


template <typename T>
void svd_bond(qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  qtensor<T> U,V;
  qtensor<T> combined = std::move(A_left * A_right);
  vector<qtensor_index> left;
  vector<qtensor_index> right;
  qtensor_index dag_mid(mid); dag_mid.dag();
  string tag  = mid.tag();
  string dtag = dag_mid.tag();
  // Separate qtensor_index
  for (size_t j = 0; j < A_right.rank; j++) {
    string idx_tag = A_right.idx_set[j].tag();
    if (idx_tag != dtag) {
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
  // Restore the correct name, type, and prime level
  U.idx_set.back().rename(mid.name());
  U.idx_set.back().retype(mid.type());
  V.idx_set[0].rename(dag_mid.name());
  V.idx_set[0].retype(dag_mid.type());
  if(U.idx_set.back().level() != mid.level()){
    U.idx_set.back().prime(int(mid.level())-int(U.idx_set.back().level()));
  }
  if(V.idx_set[0].level() != dag_mid.level()){
    V.idx_set[0].prime(int(dag_mid.level())-int(V.idx_set[0].level()));
  }
  // Copy
  A_right = V;
  A_left = U;
}
template void svd_bond(qtensor<double>& A_left, qtensor<double>& A_right, qtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(qtensor< std::complex<double> >& A_left, qtensor< std::complex<double> >& A_right, qtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);


template <typename T>
void svd_bond(qtensor<T>& combined, qtensor<T>& A_left, qtensor<T>& A_right,
        qtensor_index& mid, vector<double>& S,
        int direction, double cutoff, long unsigned K)
{
  qtensor<T> U,V;
  vector<qtensor_index> left;
  vector<qtensor_index> right;
  qtensor_index dag_mid(mid); dag_mid.dag();
  // Separate qtensor_index
  for (size_t j = 0; j < A_left.rank; j++) {
    string qidx_tag = A_left.idx_set[j].tag();
    if (qidx_tag != mid.tag()) {
      left.push_back(A_left.idx_set[j]);
    }
  }
  for (size_t j = 0; j < A_right.rank; j++) {
    string qidx_tag = A_right.idx_set[j].tag();
    if (qidx_tag != dag_mid.tag()) {
      right.push_back(A_right.idx_set[j]);
    }
  }
  // SVD
  svd(combined,left,right,U,V,S,direction,cutoff,K);
  // Restore the correct name and prime level
  U.idx_set.back().rename(mid.name());
  U.idx_set.back().retype(mid.type());
  V.idx_set[0].rename(dag_mid.name());
  V.idx_set[0].retype(dag_mid.type());
  if(U.idx_set.back().level() != mid.level()){
    U.idx_set.back().prime(int(mid.level())-int(U.idx_set.back().level()));
  }
  if(V.idx_set[0].level() != dag_mid.level()){
    V.idx_set[0].prime(int(dag_mid.level())-int(V.idx_set[0].level()));
  }
  // Copy
  A_right = V;
  A_left = U;
}
template void svd_bond(qtensor<double>& combined, qtensor<double>& A_left, qtensor<double>& A_right, qtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(qtensor< std::complex<double> >& combined, qtensor< std::complex<double> >& A_left, qtensor< std::complex<double> >& A_right, qtensor_index& mid, vector<double>& S, int direction, double cutoff, long unsigned K);

#endif
