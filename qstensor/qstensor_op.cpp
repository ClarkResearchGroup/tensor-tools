#ifndef QUANTUM_NUMBERED_TENSOR_SPARSE_OPERATIONS
#define QUANTUM_NUMBERED_TENSOR_SPARSE_OPERATIONS

#include "qstensor_op.h"

template <typename T>
void qr(qstensor<T>& A,
        vector<qtensor_index>& left,
        vector<qtensor_index>& right,
        qstensor<T>& Q, qstensor<T>& R)
{
  assert(1==2);
  // Permute qstensor, group left indices, and right indices
  /*vector<qtensor_index> new_idx_set;
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
    // copy data from qstensor blocks to vector mA
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
    // copy data from mQ to qstensor Q
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
    // copy data from mR to qstensor R
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
  }*/
}
template void qr(qstensor<double>& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qstensor<double>& Q, qstensor<double>& R);
template void qr(qstensor< std::complex<double> >& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qstensor< std::complex<double> >& Q, qstensor< std::complex<double> >& R);



template <typename T>
void svd(qstensor<T>& A,
         vector<qtensor_index>& left,
         vector<qtensor_index>& right,
         qstensor<T>& U, qstensor<T>& V, qtensor<T>& S,
         int direction,double cutoff, unsigned K)
{
  S.clearBlock();
  // Permute qstensor, group left indices, and right indices
  vector<qtensor_index> new_idx_set;
  for(auto v : left){
    new_idx_set.push_back(v);
  }
  for(auto v : right){
    new_idx_set.push_back(v);
  }
  uint_vec perm;
  find_index_permutation(A.idx_set, new_idx_set, perm);
  //A.permute(perm);
  // Accumulate legal quantum numbers for the mid bond
  unordered_map< string, unsigned > A_block_id_by_qn_str;
  for(size_t i=0;i<A.block_index_qn.size();i++){
    string A_qn_str;
    for (size_t j = 0; j < A.rank; j++) {
      A_qn_str += (to_string(A.block_index_qn[i][perm[j]])+" ");
    }
    A_block_id_by_qn_str[A_qn_str] = i;
  }
  set<int> mid_Q_set;
  unordered_map<int, set<uint_vec> > l_index_qi;
  unordered_map<int, set<uint_vec> > r_index_qi;
  for (size_t i = 0; i < A.block_index_qn.size(); i++) {
    int mid_QN = 0;
    uint_vec l_qi, r_qi;
    for (size_t j = 0; j < left.size(); j++) {
      if(left[j].arrow()==Inward){
        mid_QN += A.block_index_qn[i][perm[j]];
      }else{
        mid_QN -= A.block_index_qn[i][perm[j]];
      }
      l_qi.push_back(A.block_index_qi[i][perm[j]]);
    }
    for (size_t j = 0; j < right.size(); j++) {
      r_qi.push_back(A.block_index_qi[i][perm[j+left.size()]]);
    }
    mid_Q_set.insert(mid_QN);
    l_index_qi[mid_QN].insert(l_qi);
    r_index_qi[mid_QN].insert(r_qi);
  }
  // Set up the initial mid bond
  unordered_map<int, str_vec>  l_qn_str_map;
  unordered_map<int, str_vec>  r_qn_str_map;
  unordered_map<int, int_vec> l_qn_sizes_map;
  unordered_map<int, int_vec> r_qn_sizes_map;
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
  vector<qtensor_index> S_idx_set = {mid};
  vector< pair<double, unsigned> > S_and_block_idx;
  double nm = 0.0; //sum of all singular values

  unordered_map<string,char> charMap;
  auto indA = A.getIndices(charMap);
  auto indU = indicesToChar(U_idx_set,charMap);
  auto indV = indicesToChar(V_idx_set,charMap);
  auto indS = string(1,charMap[mid.tagNoArrow()]);

  uint_vec new_QDim(mid_QDim.size(),0);
  vector<CTF::Matrix<T>> mU(mid_Q.size());
  vector<CTF::Matrix<T>> mV(mid_Q.size());
  vector<CTF::Vector<T>> mS(mid_Q.size());

  vector<unordered_map<int,int64_t> > blockOffsets; //find corner of block in dense
  getOffsets(A,blockOffsets);
 
  // SVD block by block
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim = mid_QDim[ii];
    unsigned row = 0;
    unsigned col = 0;
    for (auto sz : l_qn_sizes_map[q]){
      row += sz;
    }
    for (auto sz : r_qn_sizes_map[q]){
      col += sz;
    }
    CTF::Matrix<T> mA(row,col);
    CTF::Matrix<T>& _U = mU[ii];
    CTF::Matrix<T>& _V = mV[ii];
    CTF::Vector<T>& _S = mS[ii];
    /*if(mA.get_tot_size(false)==1 && false){ //handle single element tensor manually so we don't do a bunch of reshapes etc
      assert(l_qn_str_map[q].size()==1 && r_qn_str_map[q].size()==1);
      unsigned A_block = A_block_id_by_qn_str[l_qn_str_map[q][0] + r_qn_str_map[q][0]];
      vector<int64_t> starts(A.rank);
      vector<int64_t> ends(A.rank);
      for(int j=0;j<A.rank;j++){
        starts[j] = blockOffsets[j][A.block_index_qn[A_block][j]];
        ends[j]   = starts[j]+A.block_index_qd[A_block][j];
        assert(ends[j]-starts[j] == 1);
      }
      vector<int64_t> lens = {1}; 
      _S = std::move(A._T.slice(starts.data(),ends.data()).reshape(1,lens.data()));
      //new_QDim[ii] = _S.len;
			int64_t ns;
			T* s_data;
			_S.read_all(&ns,&s_data);
			for (size_t j = 0; j < ns; j++) {
				S_and_block_idx.push_back(std::make_pair(s_data[j], ii));
				nm += s_data[j]*s_data[j];
			}
			delete[] s_data;
      continue;
    }*/

    std::vector<int64_t> offsetA(A.rank,0);
    vector<int64_t> ends(A.rank);
    int c_row = 0; int c_col = 0; 
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      c_col = 0;
      for (size_t j = 0; j < r_qn_str_map[q].size(); j++) {
        string qn_str = l_qn_str_map[q][i] + r_qn_str_map[q][j];
        unsigned A_block = A_block_id_by_qn_str[qn_str];

        std::vector<int64_t> offsetmA = {c_row,c_col};
        std::vector<int64_t> endsmA   = {c_row+l_qn_sizes_map[q][i],c_col+r_qn_sizes_map[q][j]};
        
        for(int k=0;k<A.rank;k++){
          offsetA[k] = blockOffsets[k][A.block_index_qn[A_block][k]];
          ends[k]    = offsetA[k]+A.block_index_qd[A_block][k];
        }
        mA.slice(offsetmA.data(),endsmA.data(),0.,A._T,offsetA.data(),ends.data(),1.);
        c_col += r_qn_sizes_map[q][j];
      }
      c_row += l_qn_sizes_map[q][i];
    }
    assert(mA.get_tot_size(false)==row*col);
    if(row==1 && col==2){ //brute force workaround for scalapack bug
      _U = CTF::Matrix<T>(1,1);
      _U["ij"] = 1.;
      //_V = CTF::Matrix<T>(1,2);
      _S = CTF::Vector<T>(1);
      _V = mA;
      double d=0.;
      CTF::Scalar<double> tot;
      using cmplx = std::complex<double>;
      if(std::is_same<T,cmplx>::value){
        tot[""] += CTF::Function<double,cmplx>([](cmplx r){ return std::real(cconj(r)*r);})(mA["ij"]);
        d=tot;
        d=sqrt(d);
      } else{
        tot[""]=std::real(mA.norm2());
        d=tot;
      }
      _S["i"] = d;
      if(d!=0.) _V["ij"] *= 1./d;
      else    { _V["ij"] = 1./sqrt(2);} //this is what lapack seems to do

    } else{
        mA.svd(_U,_S,_V,K,0);
    }
    //new_QDim[ii] = _S.len;
    if(cutoff==0 and K==0) assert(_S.get_tot_size(false)==dim);

    if(direction==MoveFromLeft){
      _V["ba"] = _S["b"]*_V["ba"];
    } else if(direction==MoveFromRight){
      _U["ab"] = _S["b"]*_U["ab"];
    }
    //TODO: make this scale better
    int64_t ns;
    T* s_data;
    _S.read_all(&ns,&s_data);
    for (size_t j = 0; j < ns; j++) {
      S_and_block_idx.push_back(std::make_pair(s_data[j], ii));
      nm += s_data[j]*s_data[j];
    }
    delete[] s_data;
  }
  // Tune dimensions of the mid bond
  std::sort(S_and_block_idx.begin(), S_and_block_idx.end(), pairCompare);
  double cumsum = 0;
  long unsigned accBD = 0;
  //TODO: have a cutoff free version
  for (size_t i = 0; i < S_and_block_idx.size(); i++) {
    cumsum += S_and_block_idx[i].first * S_and_block_idx[i].first / nm;
		if( S_and_block_idx[i].first==0. &&cutoff!=0.) break;
    new_QDim[S_and_block_idx[i].second] += 1;
    accBD += 1;
    if( (1-cumsum<cutoff &&cutoff!=0.) || (accBD>=K && K!=0) ) break;
  }
  //perr<<"Trunc err="<< std::max(1-cumsum,0.) <<" "<<endl;
  // set up mid bond
  qtensor_index midCut(Outward);
  for(size_t ii=0;ii<mid_Q.size();ii++){
    if(new_QDim[ii]>0) midCut.addQN(std::make_pair(mid_Q[ii], new_QDim[ii]));
  }
  // Set up U and V index
  U_idx_set = vector<qtensor_index>(left);
  U_idx_set.push_back(midCut);
  V_idx_set.clear();
  midCut.dag();
  V_idx_set.push_back(midCut);
  midCut.dag();
  for(auto v : right){
    V_idx_set.push_back(v);
  }
  U.reset(U_idx_set);
  U._initted = true;
  V.reset(V_idx_set);
  V._initted = true;
  S_idx_set[0] = midCut;
  S.reset(S_idx_set);
  S._block.reserve(new_QDim.size());
  S._initted = true;
  unsigned mid_num = 0;
  for (size_t ii = 0; ii < mid_Q.size(); ii++) {
    int q = mid_Q[ii];
    unsigned dim  = mid_QDim[ii];
    int ndim      = new_QDim[ii]; //cast to int to get rid of warnings
    if(ndim==0) continue;
    const set<uint_vec>& l_qi_set = l_index_qi[q];
    for (auto i1 = l_qi_set.begin(); i1 != l_qi_set.end(); ++i1){
      const uint_vec& l_qi = *i1;
      uint_vec U_block_index_qi;
      int_vec U_block_index_qd;
      int_vec  U_block_index_qn;
      int_vec  U_idx_sizes;
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
      U._block.emplace_back(U_idx_set.size(),U_block_index_qd.data());
      U.block_index_qn.push_back(U_block_index_qn);
      U.block_index_qd.push_back(U_block_index_qd);
      U.block_index_qi.push_back(U_block_index_qi);
      U.block_id_by_qn_str[U_qn_str] = U._block.size()-1;
    }
    const set<uint_vec>& r_qi_set = r_index_qi[q];
    for (auto i1 = r_qi_set.begin(); i1 != r_qi_set.end(); ++i1){
      const uint_vec& r_qi = *i1;
      uint_vec V_block_index_qi;
      int_vec V_block_index_qd;
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
      V._block.emplace_back(V_idx_set.size(),V_block_index_qd.data());
      V.block_index_qn.push_back(V_block_index_qn);
      V.block_index_qd.push_back(V_block_index_qd);
      V.block_index_qi.push_back(V_block_index_qi);
      V.block_id_by_qn_str[V_qn_str] = V._block.size()-1;
    }
    mid_num++;
    CTF::Matrix<T>& _U = mU[ii];
    CTF::Matrix<T>& _V = mV[ii];
    CTF::Vector<T>& _S = mS[ii];
    if(_S.len < ndim) assert(1==2); //nonsense?
    if(_S.len != ndim){ _S = _S.slice(0,ndim-1);} 
    S._block.emplace_back(_S);
    //for slicing
    vector<int> offU(U.rank,0);
    vector<int> offV(V.rank,0);
    /*if(l_qn_str_map[q].size()==1 && r_qn_str_map[q].size()==1 && l_qn_sizes_map[q][0]==1 && r_qn_sizes_map[q][0]==1 && false){ //handle single element tensor manually so we don't do a bunch of reshapes etc
      assert(l_qn_str_map[q].size()==1 && r_qn_str_map[q].size()==1);
      //unsigned A_block = A.block_id_by_qn_str[l_qn_str_map[q][0] + r_qn_str_map[q][0]];
      unsigned U_block = U.block_id_by_qn_str[l_qn_str_map[q][0]+to_string(q)+" "];
      unsigned V_block = V.block_id_by_qn_str[to_string(q)+" "+r_qn_str_map[q][0]];
      if(direction==MoveFromLeft){
          U._block[U_block][indU.c_str()] = 1.;
          V._block[V_block][indV.c_str()] = _S[indS.c_str()];
      } else if(direction==MoveFromRight){
          U._block[U_block][indU.c_str()] = _S[indS.c_str()];
          V._block[V_block][indV.c_str()] = 1.;
      }
      continue;
    }*/
    int c_row = 0;
    int c_col = 0;
    //perr<<"U slice"<<endl;
    for (size_t i = 0; i < l_qn_str_map[q].size(); i++) {
      string U_qn_str = l_qn_str_map[q][i] + to_string(q) + " ";
      unsigned U_block = U.block_id_by_qn_str[U_qn_str];
      int    tot_size = U._block[U_block].get_tot_size(false);
      vector<int> start = {c_row,c_col};
      vector<int> end   = {c_row+(tot_size-1)%l_qn_sizes_map[q][i]+1, c_col+tot_size/l_qn_sizes_map[q][i]};
      //perr<< "  $("<<row<<","<<col<<") "<<i<<" "<<start[0]<<","<<start[1]<<" "<<end[0]<<","<<end[1]<<" "<<tot_size<<endl;
      //U._block[U_block] = (_U.slice(start.data(),end.data()));
      //U._block[U_block] = (U._block[U_block].reshape(U.rank,U.block_index_qd[U_block].data()));
      U._block[U_block].slice(offU.data(),U.block_index_qd[U_block].data(),0.,_U,start.data(),end.data(),1.);
      c_row += l_qn_sizes_map[q][i];
    }
    //perr<<endl;
    //perr<<"V slice"<<endl;
    // copy data from mU to atensor V
    c_col=0;
    c_row=0;
    for (size_t i = 0; i < r_qn_str_map[q].size(); i++) {
      string V_qn_str = to_string(q) + " " + r_qn_str_map[q][i];
      unsigned V_block = V.block_id_by_qn_str[V_qn_str];
      int    tot_size = V._block[V_block].get_tot_size(false);
      vector<int> start = {c_row,c_col};
      vector<int> end   = {c_row+(tot_size-1)%ndim+1, c_col+(tot_size)/ndim};
      //perr<< " %("<<row<<","<<col<<") "<<i<<" "<<start[0]<<","<<start[1]<<" "<<end[0]<<","<<end[1]<<" "<<tot_size<<endl;
      //V._block[V_block] = (_V.slice(start.data(),end.data()));
      //V._block[V_block] = (V._block[V_block].reshape(V.rank,V.block_index_qd[V_block].data()));
      V._block[V_block].slice(offV.data(),V.block_index_qd[V_block].data(),0.,_V,start.data(),end.data(),1.);
      c_col += r_qn_sizes_map[q][i];
    }
    //perr<<endl;

  }//end midQ loop*/

#ifndef NDEBUG
  for(int ii=0;ii<mid_num;ii++){
    assert(U.rank==U._block[ii].order);
    for(int l=0;l<U.rank;l++){
      assert(U._block[ii].lens[l]==U.block_index_qd[ii][l]);
    }
    assert(V.rank==V._block[ii].order);
    for(int l=0;l<V.rank;l++){
      assert(V._block[ii].lens[l]==V.block_index_qd[ii][l]);
    }
  }
  //perr<<"Passed SVD"<<endl;
#endif
  //convert blocks in U,V back to sparse
  blockToSparse(U,U._T);
  blockToSparse(V,V._T);
  U._block.clear();
  V._block.clear();

}
template void svd(qstensor<double>& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qstensor<double>& U, qstensor<double>& V, qtensor<double>& S, int direction,double cutoff, unsigned K);

template<> void svd(qstensor< std::complex<double> >& A,vector<qtensor_index>& left, vector<qtensor_index>& right, qstensor< std::complex<double> >& U, qstensor< std::complex<double> >& V, qtensor<std::complex<double> >& S, int direction,double cutoff,unsigned K) { 
  assert(1==2); //CTF has broken things for complex<double> so for the meantime get it working for double
};

template <typename T>
void svd_bond(qstensor<T>& A_left, qstensor<T>& A_right,
        qtensor_index& mid, qtensor<T>& S,
        int direction, double cutoff, long unsigned K)
{
  qstensor<T> U,V;
  qstensor<T> combined = std::move(A_left * A_right);
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
template void svd_bond(qstensor<double>& A_left, qstensor<double>& A_right, qtensor_index& mid, qtensor<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(qstensor< std::complex<double> >& A_left, qstensor< std::complex<double> >& A_right, qtensor_index& mid, qtensor<std::complex<double> >& S, int direction, double cutoff, long unsigned K);


template <typename T>
void svd_bond(qstensor<T>& combined, qstensor<T>& A_left, qstensor<T>& A_right,
        qtensor_index& mid, qtensor<T>& S,
        int direction, double cutoff, long unsigned K)
{
  qstensor<T> U,V;
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
template void svd_bond(qstensor<double>& combined, qstensor<double>& A_left, qstensor<double>& A_right, qtensor_index& mid, qtensor<double>& S, int direction, double cutoff, long unsigned K);
template void svd_bond(qstensor< std::complex<double> >& combined, qstensor< std::complex<double> >& A_left, qstensor< std::complex<double> >& A_right, qtensor_index& mid, qtensor<std::complex<double> >& S, int direction, double cutoff, long unsigned K);

#endif
