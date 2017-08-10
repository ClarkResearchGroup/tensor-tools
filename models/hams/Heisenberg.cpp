#ifndef SPIN_HALF_HEISENBERG_MODEL
#define SPIN_HALF_HEISENBERG_MODEL

#include "Heisenberg.h"

template <typename T>
void Heisenberg<T>::addOperators(MPO<T>& H, unsigned site, unsigned r, unsigned c, string op, double val){
  // if(site==0 && r==1){
  //   for (size_t k = 0; k < 2; k++) {
  //     for (size_t b = 0; b < 2; b++) {
  //       H.A[site]._T.data()[k + 2*b + 2*2*c] = val * (*_s).d_bra_op_ket(b, op, k);
  //     }
  //   }
  // }else if(site==H.length-1 && c==0){
  //   for (size_t k = 0; k < 2; k++) {
  //     for (size_t b = 0; b < 2; b++) {
  //       H.A[site]._T.data()[r + 5*k + 5*2*b] = val * (*_s).d_bra_op_ket(b, op, k);
  //     }
  //   }
  // }else if(site>0 and site<H.length-1){
  //   for (size_t k = 0; k < 2; k++) {
  //     for (size_t b = 0; b < 2; b++) {
  //       H.A[site]._T.data()[r + 5*k + 5*2*b + 5*2*2*c] = val * (*_s).d_bra_op_ket(b, op, k);
  //     }
  //   }
  // }
    for (size_t k = 0; k < 2; k++) {
      for (size_t b = 0; b < 2; b++) {
        unsigned s0 = H.A[site]._T.stride(0);
        unsigned s1 = H.A[site]._T.stride(1);
        unsigned s2 = H.A[site]._T.stride(2);
        unsigned s3 = H.A[site]._T.stride(3);
        if(site==0){
          assert(k*s1 + b*s2 + s3*c < H.A[site].size);
          H.A[site]._T.data()[k*s1 + b*s2 + s3*c] = val * (*_s).d_bra_op_ket(b, op, k);
        }else{
          assert(s0*r + k*s1 + b*s2 + s3*c < H.A[site].size);
          H.A[site]._T.data()[s0*r + k*s1 + b*s2 + s3*c] = val * (*_s).d_bra_op_ket(b, op, k);
        }
      }
    }
}
template void Heisenberg<double>::addOperators(MPO<double>& H, unsigned site, unsigned r, unsigned c, string op, double val);
template void Heisenberg< std::complex<double> >::addOperators(MPO< std::complex<double> >& H, unsigned site, unsigned r, unsigned c, string op, double val);


template <typename T>
void Heisenberg<T>::addOperators(qMPO<T>& H, unsigned site, unsigned r, unsigned c, string op, double val, int Qi, int Qo){
  for (size_t k = 0; k < 2; k++) {
    for (size_t b = 0; b < 2; b++) {
      string qn_str;
      qn_str += (to_string(Qi)+" ");
      qn_str += (to_string(H.phy_qn[k])+" ");
      qn_str += (to_string(H.phy_qn[b])+" ");
      qn_str += (to_string(Qo)+" ");
      if(H.A[site].block_id_by_qn_str.count(qn_str) > 0){
        unsigned block_pos = H.A[site].block_id_by_qn_str[qn_str];
        unsigned s0 = 1;
        unsigned s1 = s0 * H.A[site].idx_set[0].qdim(H.A[site].block_index_qi[block_pos][0]);
        unsigned s2 = s1 * H.A[site].idx_set[1].qdim(H.A[site].block_index_qi[block_pos][1]);
        unsigned s3 = s2 * H.A[site].idx_set[2].qdim(H.A[site].block_index_qi[block_pos][2]);
        if(site==0){
          assert(s3*c < H.A[site].block[block_pos].size());
          H.A[site].block[block_pos][s3*c] = val * (*_s).d_bra_op_ket(b, op, k);
        }else{
          assert(s0*r + s3*c < H.A[site].block[block_pos].size());
          H.A[site].block[block_pos][s0*r + s3*c] = val * (*_s).d_bra_op_ket(b, op, k);
        }
      }
    }
  }
}
template void Heisenberg<double>::addOperators(qMPO<double>& H, unsigned site, unsigned r, unsigned c, string op, double val, int Qi, int Qo);
template void Heisenberg< std::complex<double> >::addOperators(qMPO< std::complex<double> >& H, unsigned site, unsigned r, unsigned c, string op, double val, int Qi, int Qo);


template <typename T>
void Heisenberg<T>::buildHam(MPO<T>&  H){
  // Resize the MPO
  MPO<T> A((*_s).N(), (*_s).phy_dim(), 5);
  A.setZero();
  H = A;
  // Setup MPO
  for (size_t site = 0; site < H.length; site++) {
    // Identity block
    if(site>0)          addOperators(H, site, 0, 0, "Id", 1.0);
    if(site<H.length-1) addOperators(H, site, 1, 1, "Id", 1.0);
    if(_dh!=nullptr)    addOperators(H, site, 1, 0, "Sz", _dh[site]);
    if(_tE!=0)          addOperators(H, site, 1, 0, "Id", -_tE/H.length);
    // Sz
    if(site>0)          addOperators(H, site, 2, 0, "Sz", 1.0);
    if(site<H.length-1) addOperators(H, site, 1, 2, "Sz", 1.0);
    // Sp Sm
    if(site>0)          addOperators(H, site, 3, 0, "S+", 1.0);
    if(site<H.length-1) addOperators(H, site, 1, 3, "S-", 0.5);
    // Sm Sp
    if(site>0)          addOperators(H, site, 4, 0, "S-", 1.0);
    if(site<H.length-1) addOperators(H, site, 1, 4, "S+", 0.5);
  }
}
template void Heisenberg<double>::buildHam(MPO<double>&  H);
template void Heisenberg< std::complex<double> >::buildHam(MPO< std::complex<double> >&  H);


template <typename T>
void Heisenberg<T>::buildHam(qMPO<T>& H){
  qMPO<T> A(_s);
  A.setZero();
  H = A;
  // Resize qMPO
  H.A.clear();
  string Link_name_pref = "ID"+to_string(H._id)+"Link";
  string Site_name_pref = "Site";
  for (size_t i = 0; i < H.length; i++) {
    string left_link_name  = Link_name_pref+to_string(i);
    string right_link_name = Link_name_pref+to_string(i+1);
    string site_name       = Site_name_pref+to_string(i);
    H.A.push_back(
      std::move(
        qtensor<T>(
          {Inward, Outward, Inward, Outward},
          {left_link_name, site_name, site_name, right_link_name},
          {Link, Site, Site, Link},
          {0,0,1,0})
      )
    );
    // Set QNs
    // Left Link
    if(i==0){
      H.A[i].addQNtoIndex(0, std::make_pair(0, 1));
    }else{
      H.A[i].addQNtoIndex(0, std::make_pair( 0, 3));
      H.A[i].addQNtoIndex(0, std::make_pair(-2, 1));
      H.A[i].addQNtoIndex(0, std::make_pair(+2, 1));
    }
    // Site
    for (size_t j = 0; j < H.phy_qn.size(); j++) {
      H.A[i].addQNtoIndex(1, std::make_pair(H.phy_qn[j], 1));
    }
    // Site
    for (size_t j = 0; j < H.phy_qn.size(); j++) {
      H.A[i].addQNtoIndex(2, std::make_pair(H.phy_qn[j], 1));
    }
    // Right Link
    if(i==H.length-1){
      H.A[i].addQNtoIndex(3, std::make_pair(0, 1));
    }else{
      H.A[i].addQNtoIndex(3, std::make_pair( 0, 3));
      H.A[i].addQNtoIndex(3, std::make_pair(-2, 1));
      H.A[i].addQNtoIndex(3, std::make_pair(+2, 1));
    }
    // Set up blocks
    H.A[i].initBlock();
    H.A[i].setZero();
  }
  // Setup MPO
  for (size_t site = 0; site < H.length; site++) {
    // Identity block
    if(site>0)          addOperators(H, site, 0, 0, "Id", 1.0, 0, 0);
    if(site<H.length-1) addOperators(H, site, 1, 1, "Id", 1.0, 0, 0);
    if(_dh!=nullptr)    addOperators(H, site, 1, 0, "Sz", _dh[site], 0, 0);
    if(_tE!=0)          addOperators(H, site, 1, 0, "Id", -_tE/H.length, 0, 0);
    if(site>0)          addOperators(H, site, 2, 0, "Sz", 1.0, 0, 0);
    if(site<H.length-1) addOperators(H, site, 1, 2, "Sz", 1.0, 0, 0);
    // S+ in
    if(site>0)          addOperators(H, site, 0, 0, "S+", 1.0, -2, 0);
    // S- in
    if(site>0)          addOperators(H, site, 0, 0, "S-", 1.0, +2, 0);
    // S- out
    if(site<H.length-1) addOperators(H, site, 1, 0, "S-", 0.5, 0, -2);
    // S+ out
    if(site<H.length-1) addOperators(H, site, 1, 0, "S+", 0.5, 0, +2);
  }
}
template void Heisenberg<double>::buildHam(qMPO<double>& H);
template void Heisenberg< std::complex<double> >::buildHam(qMPO< std::complex<double> >& H);

#endif