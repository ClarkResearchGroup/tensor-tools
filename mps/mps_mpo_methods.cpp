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


#endif
