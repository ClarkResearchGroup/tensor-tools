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

#ifndef MPS_MPO_METHODS
#define MPS_MPO_METHODS

#include "mps_mpo_methods.h"

template <typename T, unsigned N>
void sum(vector< dTensorTrain<T, N> >& x, dTensorTrain<T, N>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose){
  assert(x.size()>=1);
  assert(max_iter>0);
  assert(max_bd>0);
  assert(cutoff>0 && tol>0);
  // Number of TensorTrain objects to sum up
  unsigned num = x.size();
  // Properties of the systems
  unsigned length = x[0].length;
  unsigned phy_dim = x[0].phy_dim;
  // Initialize res
  res.freeTensors();
  res.setLength(length);
  res.setPhysicalDim(phy_dim);
  res.setBondDim(2);
  res.allocateTensors();
  res.setRandom();
  res.rc();
  // Containers for holding environment tensors, num * length
  vector< vector< dtensor<T> > > TL(num);
  vector< vector< dtensor<T> > > TR(num);
  for (size_t i = 0; i < num; i++) {
    assert(x[i].tensors_allocated);
    TL[i].resize(length);
    TR[i].resize(length);
  }
  // Build environment from right (TR)
  dtensor<T> t1,t2,t3,t4;
  for (size_t i = 0; i < num; i++) {
    // Left ends
    {
      dtensor<T> left_ends({res.A[0].idx_set[0], x[i].A[0].idx_set[0]});
      left_ends.idx_set[0].primeLink();
      left_ends.setOne();
      left_ends.dag();
      TL[i][0] = left_ends;
    }
    // right ends
    {
      dtensor<T> right_ends({res.A[length-1].idx_set.back(), x[i].A[length-1].idx_set.back()});
      right_ends.idx_set[0].primeLink();
      right_ends.setOne();
      right_ends.dag();
      TR[i][length-1] = right_ends;
    }
    for (size_t site = length-1; site > 0; site--) {
      t1 = res.A[site];
      t1.dag(); t1.conj(); t1.primeLink();
      t2 = x[i].A[site];
      t3 = std::move(t1*TR[i][site]);
      TR[i][site-1] = std::move(t3*t2);
    }
  }
  int step = 0;
  double delta = 1.0;
  // Start the cycle of updating site and TL/TR, using two-site algorithm
  while(step<max_iter){
    vector<double> S;
    dtensor<T> U,V;
    // Left to Right
    for (size_t site = 0; site < length-2; site++) {
      dtensor<T> agg;
      dtensor<T> org = std::move(res.A[site]*res.A[site+1]); org.dag(); org.conj();
      for (size_t i = 0; i < num; i++) {
        t2 = std::move(x[i].A[site]*x[i].A[site+1]);
        t3 = std::move(t2*TR[i][site+1]);
        t1 = std::move(TL[i][site]*t3);
        if(i==0)
          agg = t1;
        else
          agg += t1;
      }
      agg.primeLink(-1); agg.dag();
      // svd
      vector<dtensor_index> left;
      vector<dtensor_index> right;
      dtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site+1);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          left.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site+1].rank; i++) {
        string idx_name = res.A[site+1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site+1].idx_set[i];
        }else{
          right.push_back(res.A[site+1].idx_set[i]);
        }
      }
      svd(agg,left,right,U,V,S,MoveFromLeft,cutoff,max_bd);
      mid.resize(S.size());
      res.bond_dims[site+1] = S.size();
      res.A[site] = U;
      res.A[site].idx_set.back() = mid;
      res.A[site+1] = V;
      res.A[site+1].idx_set[0] = mid;
      res.center = site+1;
      double overlap = std::abs(org.contract(agg));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink();
      for (size_t i = 0; i < num; i++) {
        t2 = x[i].A[site];
        t3 = std::move(t1*TL[i][site]);
        TL[i][site+1] = std::move(t3*t2);
      }
    }
    // Right to Left
    for (size_t site = length-1; site > 1; site--) {
      dtensor<T> agg;
      dtensor<T> org = std::move(res.A[site-1]*res.A[site]);  org.dag(); org.conj();
      for (size_t i = 0; i < num; i++) {
        t2 = std::move(x[i].A[site-1]*x[i].A[site]);
        t3 = std::move(t2*TR[i][site]);
        t1 = std::move(TL[i][site-1]*t3);
        if(i==0)
          agg = t1;
        else
          agg += t1;
      }
      agg.primeLink(-1); agg.dag();
      // svd
      vector<dtensor_index> left;
      vector<dtensor_index> right;
      dtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          right.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site-1].rank; i++) {
        string idx_name = res.A[site-1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site-1].idx_set[i];
        }else{
          left.push_back(res.A[site-1].idx_set[i]);
        }
      }
      svd(agg,left,right,U,V,S,MoveFromRight,cutoff,max_bd);
      mid.resize(S.size());
      res.bond_dims[site] = S.size();
      res.A[site-1] = U;
      res.A[site-1].idx_set.back() = mid;
      res.A[site] = V;
      res.A[site].idx_set[0] = mid;
      res.center = site-1;
      double overlap = std::abs(org.contract(agg));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.conj(); t1.dag(); t1.primeLink();
      for (size_t i = 0; i < num; i++) {
        t2 = x[i].A[site];
        t3 = std::move(t1*TR[i][site]);
        TR[i][site-1] = std::move(t3*t2);
      }
    }
    if(delta<tol) break;
    ++step;
  }
}
template void sum(vector< dTensorTrain<double, 1> >& x, dTensorTrain<double, 1>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void sum(vector< dTensorTrain<double, 2> >& x, dTensorTrain<double, 2>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void sum(vector< dTensorTrain<std::complex<double>, 1> >& x, dTensorTrain<std::complex<double>, 1>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void sum(vector< dTensorTrain<std::complex<double>, 2> >& x, dTensorTrain<std::complex<double>, 2>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);


template <typename T, unsigned N>
void sum(vector< qTensorTrain<T, N> >& x, qTensorTrain<T, N>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose){
  assert(x.size()>=1);
  assert(max_iter>0);
  assert(max_bd>0);
  assert(cutoff>0 && tol>0);
#ifndef DNDEBUG
  bool same_qn = true;
  int xQ = x[0].totalQ;
  for (size_t i = 0; i < x.size(); i++) {
    if(x[i].totalQ != xQ){
      same_qn = false;
      break;
    }
  }
  if(!same_qn) std::cout << "Input vector of qMPS/qMPO have different total quantum numbers!" << '\n';
  assert(same_qn);
#endif
  // Number of TensorTrain objects to sum up
  unsigned num = x.size();
  // Properties of the systems
  unsigned length = x[0].length;
  unsigned phy_dim = x[0].phy_dim;
  // Initialize res
  res.freeTensors();
  res.setLength(length);
  res.setPhysicalDim(phy_dim);
  for (size_t i = 0; i < length+1; i++) {
    res.bond_dims[i] = 1;
  }
  res.phy_qn = x[0].phy_qn;
  res.totalQ = x[0].totalQ;
  res.allocateTensors();
  res.setRandom();
  res.rc();
  // Containers for holding environment tensors, num * length
  vector< vector< qtensor<T> > > TL(num);
  vector< vector< qtensor<T> > > TR(num);
  for (size_t i = 0; i < num; i++) {
    assert(x[i].tensors_allocated);
    TL[i].resize(length);
    TR[i].resize(length);
  }
  // Build environment from right (TR)
  qtensor<T> t1,t2,t3,t4;
  for (size_t i = 0; i < num; i++) {
    // left ends
    {
      qtensor<T> left_ends({res.A[0].idx_set[0], x[i].A[0].idx_set[0]});
      left_ends.idx_set[0].primeLink(); left_ends.idx_set[1].dag();
      left_ends.initBlock();
      left_ends.setOne();
      TL[i][0] = left_ends;
    }
    //right ends
    {
      qtensor<T> right_ends({res.A[length-1].idx_set.back(), x[i].A[length-1].idx_set.back()});
      right_ends.idx_set[0].primeLink(); right_ends.idx_set[1].dag();
      right_ends.initBlock();
      right_ends.setOne();
      TR[i][length-1] = right_ends;
    }
    for (size_t site = length-1; site > 0 ; site--) {
      t1 = res.A[site];
      t1.dag(); t1.conj(); t1.primeLink();
      t2 = x[i].A[site];
      t3 = std::move(t1*TR[i][site]);
      TR[i][site-1] = std::move(t3*t2);
    }
  }
  int step = 0;
  double delta = 1.0;
  // Start the cycle of updating site and TL/TR, using two-site algorithm
  while(step<max_iter){
    vector<double> S;
    qtensor<T> U,V;
    // Left to Right
    for (size_t site = 0; site < length-2; site++) {
      qtensor<T> agg;
      qtensor<T> org = std::move(res.A[site]*res.A[site+1]); org.dag(); org.conj();
      for (size_t i = 0; i < num; i++) {
        t2 = std::move(x[i].A[site]*x[i].A[site+1]);
        t3 = std::move(t2*TR[i][site+1]);
        t1 = std::move(TL[i][site]*t3);
        if(i==0)
          agg = t1;
        else
          agg += t1;
      }
      agg.primeLink(-1);
      // svd
      vector<qtensor_index> left;
      vector<qtensor_index> right;
      qtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site+1);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          left.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site+1].rank; i++) {
        string idx_name = res.A[site+1].idx_set[i].name();
        if (idx_name != LinkName) {
          right.push_back(res.A[site+1].idx_set[i]);
        }
      }
      svd(agg,left,right,U,V,S,MoveFromLeft,cutoff,max_bd);
      res.bond_dims[site+1] = S.size();
      res.A[site] = U;
      res.A[site].idx_set.back().rename(mid.name());
      res.A[site].idx_set.back().prime(mid.level()-res.A[site].idx_set.back().level());
      res.A[site+1] = V;
      res.A[site+1].idx_set[0].rename(mid.name());
      res.A[site+1].idx_set[0].prime(mid.level()-res.A[site+1].idx_set[0].level());
      res.center = site+1;
      double overlap = std::abs(org.contract(agg));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink();
      for (size_t i = 0; i < num; i++) {
        t2 = x[i].A[site];
        t3 = std::move(t1*TL[i][site]);
        TL[i][site+1] = std::move(t3*t2);
      }
    }
    // Right to Left
    for (size_t site = length-1; site > 1; site--) {
      qtensor<T> agg;
      qtensor<T> org = std::move(res.A[site-1]*res.A[site]);  org.dag(); org.conj();
      for (size_t i = 0; i < num; i++) {
        t2 = std::move(x[i].A[site-1]*x[i].A[site]);
        t3 = std::move(t2*TR[i][site]);
        t1 = std::move(TL[i][site-1]*t3);
        if(i==0)
          agg = t1;
        else
          agg += t1;
      }
      agg.primeLink(-1);
      // svd
      vector<qtensor_index> left;
      vector<qtensor_index> right;
      qtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name != LinkName) {
          right.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site-1].rank; i++) {
        string idx_name = res.A[site-1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site-1].idx_set[i];
        }else{
          left.push_back(res.A[site-1].idx_set[i]);
        }
      }
      svd(agg,left,right,U,V,S,MoveFromRight,cutoff,max_bd);
      res.bond_dims[site] = S.size();
      res.A[site-1] = U;
      res.A[site-1].idx_set.back().rename(mid.name());
      res.A[site-1].idx_set.back().prime(mid.level()-res.A[site-1].idx_set.back().level());
      res.A[site] = V;
      res.A[site].idx_set[0].rename(mid.name());
      res.A[site].idx_set[0].prime(mid.level()-res.A[site].idx_set[0].level());
      res.center = site-1;
      double overlap = std::abs(org.contract(agg));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.conj(); t1.dag(); t1.primeLink();
      for (size_t i = 0; i < num; i++) {
        t2 = x[i].A[site];
        t3 = std::move(t1*TR[i][site]);
        TR[i][site-1] = std::move(t3*t2);
      }
    }
    if(delta<tol) break;
    ++step;
  }
}
template void sum(vector< qTensorTrain<double, 1> >& x, qTensorTrain<double, 1>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void sum(vector< qTensorTrain<double, 2> >& x, qTensorTrain<double, 2>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void sum(vector< qTensorTrain<std::complex<double>, 1> >& x, qTensorTrain<std::complex<double>, 1>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void sum(vector< qTensorTrain<std::complex<double>, 2> >& x, qTensorTrain<std::complex<double>, 2>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);


template <typename T>
void mult(MPO<T>& A, MPO<T>& B, MPO<T>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose){
  assert(A.tensors_allocated && B.tensors_allocated);
  assert(A.length == B.length);
  assert(A.phy_dim == B.phy_dim);
  assert(max_iter>0);
  assert(max_bd>0);
  assert(cutoff>0 && tol>0);
  unsigned length = A.length;
  unsigned phy_dim = B.phy_dim;
  // Initialize res
  res.freeTensors();
  res.setLength(length);
  res.setPhysicalDim(phy_dim);
  res.setBondDim(2);
  res.allocateTensors();
  res.setRandom();
  res.rc();
  // Containers for holding environment tensors
  vector< dtensor<T> > TL(length);
  vector< dtensor<T> > TR(length);
  // Build environment from right (TR)
  dtensor<T> t1,t2,t3,t4,t5;
  // Left ends
  {
    dtensor<T> left_ends({res.A[0].idx_set[0], A.A[0].idx_set[0], B.A[0].idx_set[0]});
    left_ends.idx_set[0].primeLink(2);
    left_ends.idx_set[1].primeLink(1);
    left_ends.setOne();
    left_ends.dag();
    TL[0] = left_ends;
  }
  // right ends
  {
    dtensor<T> right_ends({res.A[length-1].idx_set.back(), A.A[length-1].idx_set.back(), B.A[length-1].idx_set.back()});
    right_ends.idx_set[0].primeLink(2);
    right_ends.idx_set[1].primeLink(1);
    right_ends.setOne();
    right_ends.dag();
    TR[length-1] = right_ends;
  }
  for (size_t site = length-1; site > 0; site--) {
    t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink(2); t1.primeSite(1); t1.mapPrime(1,0,Site);
    t2 = A.A[site];   t2.prime(1);
    t3 = B.A[site];
    t4 = std::move(t1*TR[site]);
    t5 = std::move(t4*t2);
    TR[site-1] = std::move(t5*t3);
  }
  int step = 0;
  double delta = 1.0;
  // Start the cycle of updating site and TL/TR, using two-site algorithm
  while(step<max_iter){
    vector<double> S;
    dtensor<T> U,V;
    // Left to Right
    for (size_t site = 0; site < length-2; site++) {
      // save an original copy
      dtensor<T> org = std::move(res.A[site]*res.A[site+1]);
      org.dag(); org.conj();
      // get the updated tensor
      dtensor<T> upd;
      t1 = A.A[site+1]; t1.prime();
      t2 = B.A[site+1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site+1]);
      t1 = A.A[site]; t1.prime();
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site]*t4);
      upd.primeLink(-2); upd.mapPrime(2,1);
      // svd
      vector<dtensor_index> left;
      vector<dtensor_index> right;
      dtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site+1);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          left.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site+1].rank; i++) {
        string idx_name = res.A[site+1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site+1].idx_set[i];
        }else{
          right.push_back(res.A[site+1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromLeft,cutoff,max_bd);
      mid.resize(S.size());
      res.bond_dims[site+1] = S.size();
      res.A[site] = U;
      res.A[site].idx_set.back() = mid;
      res.A[site+1] = V;
      res.A[site+1].idx_set[0] = mid;
      res.center = site+1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink(2); t1.primeSite(1); t1.mapPrime(1,0,Site);
      t2 = A.A[site];   t2.prime(1);
      t3 = B.A[site];
      t4 = std::move(t1*TL[site]);
      t5 = std::move(t4*t2);
      TL[site+1] = std::move(t5*t3);
    }
    // Right to Left
    for (size_t site = length-1; site > 1; site--) {
      // save an original copy
      dtensor<T> org = std::move(res.A[site-1]*res.A[site]);
      org.dag(); org.conj();
      // update
      dtensor<T> upd;
      t1 = A.A[site]; t1.prime();
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site]);
      t1 = A.A[site-1]; t1.prime();
      t2 = B.A[site-1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site-1]*t4);
      upd.primeLink(-2); upd.mapPrime(2,1);
      // svd
      vector<dtensor_index> left;
      vector<dtensor_index> right;
      dtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          right.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site-1].rank; i++) {
        string idx_name = res.A[site-1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site-1].idx_set[i];
        }else{
          left.push_back(res.A[site-1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromRight,cutoff,max_bd);
      mid.resize(S.size());
      res.bond_dims[site] = S.size();
      res.A[site-1] = U;
      res.A[site-1].idx_set.back() = mid;
      res.A[site] = V;
      res.A[site].idx_set[0] = mid;
      res.center = site-1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink(2); t1.primeSite(1); t1.mapPrime(1,0,Site);
      t2 = A.A[site];   t2.prime(1);
      t3 = B.A[site];
      t4 = std::move(t1*TR[site]);
      t5 = std::move(t4*t2);
      TR[site-1] = std::move(t5*t3);
    }
    if(delta<tol) break;
    ++step;
  }
}
template void mult(MPO<double>& A, MPO<double>& B, MPO<double>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void mult(MPO< std::complex<double> >& A, MPO< std::complex<double> >& B, MPO< std::complex<double> >& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);



template <typename T>
void mult(qMPO<T>& A, qMPO<T>& B, qMPO<T>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose){
  assert(A.tensors_allocated && B.tensors_allocated);
  assert(A.length == B.length);
  assert(A.phy_dim == B.phy_dim);
  assert(max_iter>0);
  assert(max_bd>0);
  assert(cutoff>0 && tol>0);
  unsigned length = A.length;
  unsigned phy_dim = B.phy_dim;
  // Initialize res
  res.freeTensors();
  res.setLength(length);
  res.setPhysicalDim(phy_dim);
  for (size_t i = 0; i < length+1; i++) {
    res.bond_dims[i] = 1;
  }
  res.phy_qn = A.phy_qn;
  res.totalQ = A.totalQ;
  res.allocateTensors();
  res.setRandom();
  res.rc();
  // Containers for holding environment tensors, num * length
  vector< qtensor<T> > TL(length);
  vector< qtensor<T> > TR(length);
  // Build environment from right (TR)
  qtensor<T> t1,t2,t3,t4,t5;
  // Left ends
  {
    qtensor<T> left_ends({res.A[0].idx_set[0], A.A[0].idx_set[0], B.A[0].idx_set[0]});
    left_ends.idx_set[0].primeLink(2); left_ends.idx_set[0].dag();
    left_ends.idx_set[1].primeLink(1);
    left_ends.initBlock();
    left_ends.setOne();
    left_ends.dag();
    TL[0] = left_ends;
  }
  // right ends
  {
    qtensor<T> right_ends({res.A[length-1].idx_set.back(), A.A[length-1].idx_set.back(), B.A[length-1].idx_set.back()});
    right_ends.idx_set[0].primeLink(2); right_ends.idx_set[0].dag();
    right_ends.idx_set[1].primeLink(1);
    right_ends.initBlock();
    right_ends.setOne();
    right_ends.dag();
    TR[length-1] = right_ends;
  }
  for (size_t site = length-1; site > 0; site--) {
    t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink(2); t1.primeSite(1); t1.mapPrime(1,0,Site);
    t2 = A.A[site];   t2.prime(1);
    t3 = B.A[site];
    t4 = std::move(t1*TR[site]);
    t5 = std::move(t4*t2);
    TR[site-1] = std::move(t5*t3);
  }
  int step = 0;
  double delta = 1.0;
  // Start the cycle of updating site and TL/TR, using two-site algorithm
  while(step<max_iter){
    vector<double> S;
    qtensor<T> U,V;
    // Left to Right
    for (size_t site = 0; site < length-2; site++) {
      // save an original copy
      qtensor<T> org = std::move(res.A[site]*res.A[site+1]);
      org.dag(); org.conj();
      // get the updated tensor
      qtensor<T> upd;
      t1 = A.A[site+1]; t1.prime();
      t2 = B.A[site+1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site+1]);
      t1 = A.A[site]; t1.prime();
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site]*t4);
      upd.primeLink(-2); upd.mapPrime(2,1);
      // svd
      vector<qtensor_index> left;
      vector<qtensor_index> right;
      qtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site+1);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          left.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site+1].rank; i++) {
        string idx_name = res.A[site+1].idx_set[i].name();
        if (idx_name != LinkName) {
          right.push_back(res.A[site+1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromLeft,cutoff,max_bd);
      res.bond_dims[site+1] = S.size();
      res.A[site] = U;
      res.A[site].idx_set.back().rename(mid.name());
      res.A[site].idx_set.back().prime(mid.level()-res.A[site].idx_set.back().level());
      res.A[site+1] = V;
      res.A[site+1].idx_set[0].rename(mid.name());
      res.A[site+1].idx_set[0].prime(mid.level()-res.A[site+1].idx_set[0].level());
      res.center = site+1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink(2); t1.primeSite(1); t1.mapPrime(1,0,Site);
      t2 = A.A[site];   t2.prime(1);
      t3 = B.A[site];
      t4 = std::move(t1*TL[site]);
      t5 = std::move(t4*t2);
      TL[site+1] = std::move(t5*t3);
    }
    // Right to Left
    for (size_t site = length-1; site > 1; site--) {
      // save an original copy
      qtensor<T> org = std::move(res.A[site-1]*res.A[site]);
      org.dag(); org.conj();
      // update
      qtensor<T> upd;
      t1 = A.A[site]; t1.prime();
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site]);
      t1 = A.A[site-1]; t1.prime();
      t2 = B.A[site-1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site-1]*t4);
      upd.primeLink(-2); upd.mapPrime(2,1);
      // svd
      vector<qtensor_index> left;
      vector<qtensor_index> right;
      qtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name != LinkName) {
          right.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site-1].rank; i++) {
        string idx_name = res.A[site-1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site-1].idx_set[i];
        }else{
          left.push_back(res.A[site-1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromRight,cutoff,max_bd);
      res.bond_dims[site] = S.size();
      res.A[site-1] = U;
      res.A[site-1].idx_set.back().rename(mid.name());
      res.A[site-1].idx_set.back().prime(mid.level()-res.A[site-1].idx_set.back().level());
      res.A[site] = V;
      res.A[site].idx_set[0].rename(mid.name());
      res.A[site].idx_set[0].prime(mid.level()-res.A[site].idx_set[0].level());
      res.center = site-1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.primeLink(2); t1.primeSite(1); t1.mapPrime(1,0,Site);
      t2 = A.A[site];   t2.prime(1);
      t3 = B.A[site];
      t4 = std::move(t1*TR[site]);
      t5 = std::move(t4*t2);
      TR[site-1] = std::move(t5*t3);
    }
    if(delta<tol) break;
    ++step;
  }
}
template void mult(qMPO<double>& A, qMPO<double>& B, qMPO<double>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void mult(qMPO< std::complex<double> >& A, qMPO< std::complex<double> >& B, qMPO< std::complex<double> >& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);


template <typename T>
void mult(MPO<T>& A, MPS<T>& B, MPS<T>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose){
  assert(A.tensors_allocated && B.tensors_allocated);
  assert(A.length == B.length);
  assert(A.phy_dim == B.phy_dim);
  assert(max_iter>0);
  assert(max_bd>0);
  assert(cutoff>0 && tol>0);
  unsigned length = A.length;
  unsigned phy_dim = B.phy_dim;
  // Initialize res
  res.freeTensors();
  res.setLength(length);
  res.setPhysicalDim(phy_dim);
  res.setBondDim(2);
  res.allocateTensors();
  res.setRandom();
  res.rc();
  // Containers for holding environment tensors
  vector< dtensor<T> > TL(length);
  vector< dtensor<T> > TR(length);
  // Build environment from right (TR)
  dtensor<T> t1,t2,t3,t4,t5;
  // Left ends
  {
    dtensor<T> left_ends({res.A[0].idx_set[0], A.A[0].idx_set[0], B.A[0].idx_set[0]});
    left_ends.idx_set[0].primeLink();
    left_ends.setOne();
    left_ends.dag();
    TL[0] = left_ends;
  }
  // right ends
  {
    dtensor<T> right_ends({res.A[length-1].idx_set.back(), A.A[length-1].idx_set.back(), B.A[length-1].idx_set.back()});
    right_ends.idx_set[0].primeLink();
    right_ends.setOne();
    right_ends.dag();
    TR[length-1] = right_ends;
  }
  for (size_t site = length-1; site > 0; site--) {
    t1 = res.A[site]; t1.dag(); t1.conj(); t1.prime();
    t2 = A.A[site];
    t3 = B.A[site];
    t4 = std::move(t1*TR[site]);
    t5 = std::move(t4*t2);
    TR[site-1] = std::move(t5*t3);
  }
  int step = 0;
  double delta = 1.0;
  // Start the cycle of updating site and TL/TR, using two-site algorithm
  while(step<max_iter){
    vector<double> S;
    dtensor<T> U,V;
    // Left to Right
    for (size_t site = 0; site < length-2; site++) {
      // save an original copy
      dtensor<T> org = std::move(res.A[site]*res.A[site+1]);
      org.dag(); org.conj();
      // get the updated tensor
      dtensor<T> upd;
      t1 = A.A[site+1];
      t2 = B.A[site+1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site+1]);
      t1 = A.A[site];
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site]*t4);
      upd.prime(-1);
      // svd
      vector<dtensor_index> left;
      vector<dtensor_index> right;
      dtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site+1);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          left.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site+1].rank; i++) {
        string idx_name = res.A[site+1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site+1].idx_set[i];
        }else{
          right.push_back(res.A[site+1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromLeft,cutoff,max_bd);
      mid.resize(S.size());
      res.bond_dims[site+1] = S.size();
      res.A[site] = U;
      res.A[site].idx_set.back() = mid;
      res.A[site+1] = V;
      res.A[site+1].idx_set[0] = mid;
      res.center = site+1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.prime();
      t2 = A.A[site];
      t3 = B.A[site];
      t4 = std::move(t1*TL[site]);
      t5 = std::move(t4*t2);
      TL[site+1] = std::move(t5*t3);
    }
    // Right to Left
    for (size_t site = length-1; site > 1; site--) {
      // save an original copy
      dtensor<T> org = std::move(res.A[site-1]*res.A[site]);
      org.dag(); org.conj();
      // update
      dtensor<T> upd;
      t1 = A.A[site];
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site]);
      t1 = A.A[site-1];
      t2 = B.A[site-1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site-1]*t4);
      upd.prime(-1);
      // svd
      vector<dtensor_index> left;
      vector<dtensor_index> right;
      dtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          right.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site-1].rank; i++) {
        string idx_name = res.A[site-1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site-1].idx_set[i];
        }else{
          left.push_back(res.A[site-1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromRight,cutoff,max_bd);
      mid.resize(S.size());
      res.bond_dims[site] = S.size();
      res.A[site-1] = U;
      res.A[site-1].idx_set.back() = mid;
      res.A[site] = V;
      res.A[site].idx_set[0] = mid;
      res.center = site-1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.prime();
      t2 = A.A[site];
      t3 = B.A[site];
      t4 = std::move(t1*TR[site]);
      t5 = std::move(t4*t2);
      TR[site-1] = std::move(t5*t3);
    }
    if(delta<tol) break;
    ++step;
  }
}
template void mult(MPO<double>& A, MPS<double>& B, MPS<double>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void mult(MPO< std::complex<double> >& A, MPS< std::complex<double> >& B, MPS< std::complex<double> >& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);


template <typename T>
void mult(qMPO<T>& A, qMPS<T>& B, qMPS<T>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose){
  assert(A.tensors_allocated && B.tensors_allocated);
  assert(A.length == B.length);
  assert(A.phy_dim == B.phy_dim);
  assert(max_iter>0);
  assert(max_bd>0);
  assert(cutoff>0 && tol>0);
  unsigned length = A.length;
  unsigned phy_dim = B.phy_dim;
  // Initialize res
  res.freeTensors();
  res.setLength(length);
  res.setPhysicalDim(phy_dim);
  for (size_t i = 0; i < length+1; i++) {
    res.bond_dims[i] = 1;
  }
  res.phy_qn = A.phy_qn;
  res.totalQ = A.totalQ;
  res.allocateTensors();
  res.setRandom();
  res.rc();
  // Containers for holding environment tensors, num * length
  vector< qtensor<T> > TL(length);
  vector< qtensor<T> > TR(length);
  // Build environment from right (TR)
  qtensor<T> t1,t2,t3,t4,t5;
  // Left ends
  {
    qtensor<T> left_ends({res.A[0].idx_set[0], A.A[0].idx_set[0], B.A[0].idx_set[0]});
    left_ends.idx_set[0].primeLink(); left_ends.idx_set[0].dag();
    left_ends.initBlock();
    left_ends.setOne();
    left_ends.dag();
    TL[0] = left_ends;
  }
  // right ends
  {
    qtensor<T> right_ends({res.A[length-1].idx_set.back(), A.A[length-1].idx_set.back(), B.A[length-1].idx_set.back()});
    right_ends.idx_set[0].primeLink(); right_ends.idx_set[0].dag();
    right_ends.initBlock();
    right_ends.setOne();
    right_ends.dag();
    TR[length-1] = right_ends;
  }
  for (size_t site = length-1; site > 0; site--) {
    t1 = res.A[site]; t1.dag(); t1.conj(); t1.prime();
    t2 = A.A[site];
    t3 = B.A[site];
    t4 = std::move(t1*TR[site]);
    t5 = std::move(t4*t2);
    TR[site-1] = std::move(t5*t3);
  }
  int step = 0;
  double delta = 1.0;
  // Start the cycle of updating site and TL/TR, using two-site algorithm
  while(step<max_iter){
    vector<double> S;
    qtensor<T> U,V;
    // Left to Right
    for (size_t site = 0; site < length-2; site++) {
      // save an original copy
      qtensor<T> org = std::move(res.A[site]*res.A[site+1]);
      org.dag(); org.conj();
      // get the updated tensor
      qtensor<T> upd;
      t1 = A.A[site+1];
      t2 = B.A[site+1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site+1]);
      t1 = A.A[site];
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site]*t4);
      upd.prime(-1);
      // svd
      vector<qtensor_index> left;
      vector<qtensor_index> right;
      qtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site+1);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site].idx_set[i];
        }else{
          left.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site+1].rank; i++) {
        string idx_name = res.A[site+1].idx_set[i].name();
        if (idx_name != LinkName) {
          right.push_back(res.A[site+1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromLeft,cutoff,max_bd);
      res.bond_dims[site+1] = S.size();
      res.A[site] = U;
      res.A[site].idx_set.back().rename(mid.name());
      res.A[site].idx_set.back().prime(mid.level()-res.A[site].idx_set.back().level());
      res.A[site+1] = V;
      res.A[site+1].idx_set[0].rename(mid.name());
      res.A[site+1].idx_set[0].prime(mid.level()-res.A[site+1].idx_set[0].level());
      res.center = site+1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.prime();
      t2 = A.A[site];
      t3 = B.A[site];
      t4 = std::move(t1*TL[site]);
      t5 = std::move(t4*t2);
      TL[site+1] = std::move(t5*t3);
    }
    // Right to Left
    for (size_t site = length-1; site > 1; site--) {
      // save an original copy
      qtensor<T> org = std::move(res.A[site-1]*res.A[site]);
      org.dag(); org.conj();
      // update
      qtensor<T> upd;
      t1 = A.A[site];
      t2 = B.A[site];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*TR[site]);
      t1 = A.A[site-1];
      t2 = B.A[site-1];
      t3 = std::move(t2*t1);
      t4 = std::move(t3*t4);
      upd = std::move(TL[site-1]*t4);
      upd.prime(-1);
      // svd
      vector<qtensor_index> left;
      vector<qtensor_index> right;
      qtensor_index mid;
      string LinkName = "ID"+to_string(res._id)+"Link"+to_string(site);
      for (size_t i = 0; i < res.A[site].rank; i++) {
        string idx_name = res.A[site].idx_set[i].name();
        if (idx_name != LinkName) {
          right.push_back(res.A[site].idx_set[i]);
        }
      }
      for (size_t i = 0; i < res.A[site-1].rank; i++) {
        string idx_name = res.A[site-1].idx_set[i].name();
        if (idx_name == LinkName) {
          mid = res.A[site-1].idx_set[i];
        }else{
          left.push_back(res.A[site-1].idx_set[i]);
        }
      }
      svd(upd,left,right,U,V,S,MoveFromRight,cutoff,max_bd);
      res.bond_dims[site] = S.size();
      res.A[site-1] = U;
      res.A[site-1].idx_set.back().rename(mid.name());
      res.A[site-1].idx_set.back().prime(mid.level()-res.A[site-1].idx_set.back().level());
      res.A[site] = V;
      res.A[site].idx_set[0].rename(mid.name());
      res.A[site].idx_set[0].prime(mid.level()-res.A[site].idx_set[0].level());
      res.center = site-1;
      double overlap = std::abs(org.contract(upd));
      double nm = org.norm();
      delta = std::sqrt(std::abs(overlap - nm*nm))/nm;
      if(verbose) std::cout << "delta = " << delta << '\n';
      // update env
      t1 = res.A[site]; t1.dag(); t1.conj(); t1.prime();
      t2 = A.A[site];
      t3 = B.A[site];
      t4 = std::move(t1*TR[site]);
      t5 = std::move(t4*t2);
      TR[site-1] = std::move(t5*t3);
    }
    if(delta<tol) break;
    ++step;
  }
}
template void mult(qMPO<double>& A, qMPS<double>& B, qMPS<double>& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);
template void mult(qMPO< std::complex<double> >& A, qMPS< std::complex<double> >& B, qMPS< std::complex<double> >& res, int max_iter, double cutoff, int max_bd, double tol, bool verbose);

//Below functions are based off of ITensor code
//////////////////////////////////////////////////////////////////////
//returns  error,newStart
std::tuple<double,int> determineCutoff(double* evals,int size,int maxm,
                                          double cutoff, bool absCutoff=false, bool rescale=true){
  if(maxm<0) maxm = size;
  //Note: evals is sorted smallest->largest

  //if all zeros
  if(size==1)            return std::make_tuple(0.,0.);
  //for convenience make negatives zero
  for(int i=0;i<size;i++){
    if(evals[i]>=0.) break;
    evals[i] = 0.;
  }
  if(evals[size-1]==0.0) return std::make_tuple(0.,0.);
  double error = 0.0;
  int n=0;
  //determine error from just maxm truncation
  while ( size-n > maxm){ error += evals[n++]; }
  
  if(absCutoff)
    while(evals[n] < cutoff) { error += evals[n++]; }
  else{
    double scale = 1.0;
    if(rescale){ //normalize error
      scale = 0.0;
      for(int i=0;i<size;i++) scale+=evals[i];
    }
   while(error+evals[n] < cutoff*scale && n<size){
    error += evals[n++];
   }
  error /= scale; 
  }
  return std::make_tuple(error,n);

} 
template<typename T>
void addIndex(dtensor<T>& A, string name){
  /*A.idx_set.emplace_back(1,name,Link);
  ++A.rank;*/
  dtensor<T> contractor({ {1,name,Link}});
  contractor.setOne();
  A = std::move(contractor*A);
  /*auto tmp_set = A.idx_set;
  tmp_set.emplace_back(1,name,Link);
  A = std::move(dtensor<T>(*/

}

template void addIndex(dtensor<double>& A, string name);
template void addIndex(dtensor<complex<double>>& A, string name);

template<typename T>
void removeIndex(dtensor<T>& A, string name){
  /*for(int l=0;l<A.rank;l++){
    if(A.idx_set[l].name()==name){
      assert(A.idx_set[l].size()==1); //only for dangling
      A.idx_set.erase(A.idx_set.begin()+l);
      break;
    }
  }
  --A.rank;*/
  dtensor<T> contractor({ {1,name,Link}});
  contractor.setOne();
  A = std::move(A*contractor);
}
template void removeIndex(dtensor<double>& A, string name);
template void removeIndex(dtensor<complex<double>>& A, string name);

template<typename T>
MPS<T> exactApplyMPO(MPO<T> & K, MPS<T> & psi,double cutoff,int maxm, bool verbose){
  //TODO: allow const multiply
  assert(K.length==psi.length);
  //TODO: check that K and psi have the same sites
  int L = K.length;
  //HACK, remove dangling indices
  removeIndex(psi.A[0],"ID"+to_string(psi._id)+"Link0");
  removeIndex(psi.A[L-1],"ID"+to_string(psi._id)+"Link"+to_string(L));
  removeIndex(K.A[0],"ID"+to_string(K._id)+"Link0");
  removeIndex(K.A[L-1],"ID"+to_string(K._id)+"Link"+to_string(L));
  MPS<T> res=psi;
  //build environment tensors
  auto E = std::vector<dtensor<T> >(L);
  {
    MPS<T> psic = psi;
    MPO<T> Kc   = K;
    for(int j=0;j<L;j++){
      if(j==0){
        /*
         * auto ci = commonIndex(psi.A(1),psi.A(2),linkType);
         *psic.Aref(j) = dag(mapprime(psi.A(j),siteType,0,2,ci,0,plev));
         *ci = commonIndex(Kc.A(1),Kc.A(2),linkType);
         *Kc.Aref(j) = dag(mapprime(K.A(j),siteType,0,2,ci,0,plev));*/
        //find link index
        std::vector<dtensor_index> commonBonds;
        index_sets_intersection(psi.A[0].idx_set, psi.A[1].idx_set, commonBonds);
        assert(commonBonds.size()==1);
        psic.A[j].mapPrime(0,2,Site);
        psic.A[j].mapPrime(commonBonds,0,1717);

        index_sets_intersection(Kc.A[0].idx_set, Kc.A[1].idx_set, commonBonds);
        assert(commonBonds.size()==1);
        Kc.A[j].mapPrime(0,2,Site);
        Kc.A[j].mapPrime(commonBonds,0,1717);
      
      }
      else{
        psic.A[j].mapPrime(0,2,Site);
        psic.A[j].mapPrime(0,1717,Link);
        Kc.A[j].mapPrime(0,2,Site);
        Kc.A[j].mapPrime(0,1717,Link);
      }
    }
    E[0] = std::move(psi.A[0]*K.A[0]*Kc.A[0]*psic.A[0]);
    for(int j=1;j<L-1;j++){
      E[j] = std::move(E[j-1]*psi.A[j]*K.A[j]*Kc.A[j]*psic.A[j]);
    }
  }
  if(verbose) std::cerr<<"made Enviro"<<std::endl;
  //done making enviro
  auto O = std::move(psi.A[L-1]*K.A[L-1]);
  O.noPrime(Site);
  
  //auto Otemp = O;
  //Otemp.prime(1717);
  //auto rho = std::move(E[L-2]* O * Otemp);
  auto rho = std::move(E[L-2]*O);
       O.prime(1717);
       rho = std::move(rho*O);
       O.prime(-1717);
  //cerr<<rho.norm()<<" "<<rho.contract(rho)<<endl;
  //DIAG will destroy rho and replace with eigenvectors
  //but we need to init TBLIS so we keep rho around
  int matrixSize = 1;
  for(auto it = rho.idx_set.begin();it!=rho.idx_set.end();++it){
      if(it->level()==0) matrixSize *= it->size();
  }
  double* evals     = new double [matrixSize];
  double error; int newStart;
  DIAGD(matrixSize, rho._T.data(), evals);
  //for(int l=0;l<matrixSize;l++){ cerr<<evals[l]<<" ";} cerr<<endl;
  std::tie(error,newStart) = determineCutoff(evals,matrixSize,maxm,cutoff);
  //cerr<<"Err: "<<error<< " "<<newStart<<endl;
  delete[] evals;
  unsigned newm = 1;
  auto it = rho.idx_set.begin();
  while (it != rho.idx_set.end()){  
    if(it->level() !=0 ){
      newm *= it->size();
      it = rho.idx_set.erase(it);
    }
    else{ ++it; }
  }
  newm -= newStart;
  rho.idx_set.emplace_back(newm,"a"+to_string(L-1),Link);
  res.bond_dims[L-1] = newm;
  //res.A[L-1] = std::move(dtensor<double>(rho.idx_set,rho._T.data()+(matrixSize*newStart)));
  res.A[L-1].resize(rho.idx_set);
  std::copy(rho._T.data()+(matrixSize*newStart),rho._T.data()+(matrixSize*matrixSize),res.A[L-1]._T.data());
  //cerr<<L-1<< " "<<res.A[L-1].norm()<<" "<<res.A[L-1].contract(res.A[L-1])<<endl;
  assert(res.A[L-1].rank==2);
  O = std::move(O*res.A[L-1]*psi.A[L-2]*K.A[L-2]);
  /*O = std::move(O*res.A[L-1]);
  O = std::move(O*psi.A[L-2]);
  O = std::move(O*K.A[L-2]);*/
  O.noPrime(Site);
  
  for(int j = L-2; j > 0; --j){
    //Otemp = O; Otemp.prime(1717);
    //rho = std::move(E[j-1]*O*Otemp);
    rho = std::move(E[j-1]*O);
    O.prime(1717);
    rho = std::move(rho*O);
    O.prime(-1717);
    //cerr<<j<< " "<<rho.norm()<<" "<<rho.contract(rho)<<endl;
    //cerr<<j<<endl;
    //if(j>L/2) rho.print(0);
    
    matrixSize = 1;
    for(auto it = rho.idx_set.begin();it!=rho.idx_set.end();++it){
      if(it->level()==0) matrixSize *= it->size();
    }
    evals     = new double [matrixSize];
    DIAGD(matrixSize, rho._T.data(), evals);
    std::tie(error,newStart) = determineCutoff(evals,matrixSize,maxm,cutoff);
    if(verbose) std::cerr<<j<<" Err: "<<error<< " "<<newStart<<std::endl;
    delete[] evals;
    //convert indices from primed to new link
    newm = 1;
    auto it = rho.idx_set.begin();
    while (it != rho.idx_set.end()){  
      if(it->level() !=0){
        newm *= it->size();
        it = rho.idx_set.erase(it);
      }
      else{ ++it; }
    }
    newm -=newStart;
    rho.idx_set.emplace_back(newm,"a"+to_string(j),Link);
    rho.rank = rho.idx_set.size();
    res.bond_dims[j] = newm;
    //res.A[j] = std::move(dtensor<double>(rho.idx_set,
    //                                     rho._T.data()+(matrixSize*newStart) ));

    res.A[j].resize(rho.idx_set); 
    std::copy(rho._T.data()+(matrixSize*newStart),rho._T.data()+(matrixSize*matrixSize),res.A[j]._T.data());
    //cerr<<j<< " "<<res.A[j].norm()<<" "<<res.A[j].contract(res.A[j])<<endl;
    //if(j>L/2) res.A[j].print();

    O = std::move(O*res.A[j]*psi.A[j-1]*K.A[j-1]);
    /*O = std::move(O*res.A[j]);
    O = std::move(O*psi.A[j-1]);
    O = std::move(O*K.A[j-1]);*/
    O.noPrime(Site);
  }

  //O /= O.norm();
  res.A[0] = std::move(O);
  res.center = 0;
  //rename things because this is annoying
  for(int j=0;j<L;j++){
    auto& this_idx = res.A[j].idx_set;
    auto it = this_idx.begin();
    while (it!=this_idx.end()){
      if(it->name()=="a"+to_string(j)){
        it->rename("ID"+to_string(res._id)+"Link"+to_string(j));
      }
      if(it->name()=="a"+to_string(j+1)){
        it->rename("ID"+to_string(res._id)+"Link"+to_string(j+1));
      }
      ++it;
    }
  }
  //HACK
  addIndex(psi.A[0],"ID"+to_string(psi._id)+"Link0");
  addIndex(psi.A[L-1],"ID"+to_string(psi._id)+"Link"+to_string(L));
  addIndex(res.A[0],"ID"+to_string(res._id)+"Link0");
  addIndex(res.A[L-1],"ID"+to_string(res._id)+"Link"+to_string(L));
  addIndex(K.A[0],"ID"+to_string(K._id)+"Link0");
  addIndex(K.A[L-1],"ID"+to_string(K._id)+"Link"+to_string(L));
  /*for(int j=0;j<L;j++){
    uint_vec perm;
    vector<dtensor_index> new_idx_set = res.A[j].idx_set;
    std::sort(new_idx_set.begin(),new_idx_set.end());
    find_index_permutation(res.A[j].idx_set, new_idx_set, perm);
    res.A[j].permute(perm);
  }*/
  if(verbose){
    std::cerr<<"Res:"<<std::endl;
    res.print(1);
    std::cerr<<"!!"<<std::endl;
  }
  return res;
}
MPS<double> exactApplyMPO(MPO<double> & K, MPS<double> & psi,double cutoff,int maxm, bool verbose);
//MPS<complex<double>> exactApplyMPO(MPO<complex<double>> & K, MPS<complex<double>> & psi,double cutoff,int maxm, bool verbose){

  
template <typename T>  
T errorMPOProd(MPS<T> & psi2, MPO<T> & K, MPS<T> & psi1){
  T err = overlap(psi2,psi2);
  err += -2.*std::real(overlap(psi2,K,psi1));
  //create K\dagger
  MPO<T> Kd = K;
  for(unsigned j=0;j<K.length;j++){
    //swap 0,1, on sites
    Kd.A[j].mapPrime(0,3,Site);
    Kd.A[j].mapPrime(1,0,Site);
    Kd.A[j].mapPrime(3,1,Site);
  }
  err /= overlap(psi1,Kd,K,psi1);
  err = std::sqrt(abs(1.0+err));//abs needed for underflow*/
  return err;
}
template double errorMPOProd(MPS<double> & psi2, MPO<double> & K, MPS<double> & psi1);


#endif
