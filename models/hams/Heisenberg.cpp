#ifndef SPIN_HALF_HEISENBERG_MODEL
#define SPIN_HALF_HEISENBERG_MODEL

#include "Heisenberg.h"

struct SiteQN{
    SiteTerm st;
    QN_t q;

    SiteQN() {}
    SiteQN(SiteTerm const& st_, QN_t const& q_): st(st_),q(q_){}
};

string
startTerm(const string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cdagup","Adagup*F"),
           make_pair("Cup","Aup*F"),
           make_pair("Cdagdn","Adagdn"),
           make_pair("Cdn","Adn"),
           make_pair("C","A*F"), //A*F is -A, so essentially a trick for putting in a -1
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }

string
endTerm(const string& op)
    {
    static array<pair<string,string>,6>
           rewrites =
           {{
           make_pair("Cup","Aup"),
           make_pair("Cdagup","Adagup"),
           make_pair("Cdn","F*Adn"),
           make_pair("Cdagdn","F*Adagdn"),
           make_pair("C","A"),
           make_pair("Cdag","Adag")
           }};
    for(auto& p : rewrites)
        {
        if(p.first == op) return p.second;
        }
    return op;
    }
template <typename T>
void Heisenberg<T>::addOperators(MPO<T>& H, unsigned site, unsigned r, unsigned c, string op, double val){
  //perr<<"Adding @"<<site<<" r="<<r<<" c="<<c<<" op="<<op<<" val="<<val<<endl;
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
    auto s0 = 1;
    auto s1 = H.A[site].idx_set[0].size()*s0;
    auto s2 = H.A[site].idx_set[1].size()*s1;
    auto s3 = H.A[site].idx_set[2].size()*s2;
    for (size_t k = 0; k < H.phy_dim; k++) {
      for (size_t b = 0; b < H.phy_dim; b++) {
        if(site==0){
          // assert(k*s1 + b*s2 + s3*c < H.A[site].size);

          //vector<long int> inds = {0,k,b,c};
          //TODO: throw if unsigned long int doesnt fit in long int
          vector<long int> inds = {static_cast<long int>(k*s1 + b*s2 + s3*c)};

          const vector<T> myData = {val * (*_s).d_bra_op_ket(b, op, k)};

          //cerr<<site<<" "<<val * (*_s).d_bra_op_ket(b, op, k)<<endl;
          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          if(rank==0)
            H.A[site].__T.write(1, inds.data(), myData.data());
          else
            H.A[site].__T.write(0, inds.data(), myData.data());
        }else{
          //  assert(s0*r + k*s1 + b*s2 + s3*c < H.A[site].size);

          //vector<long int> inds = {r,k,b,c};
          vector<long int> inds = {static_cast<long int>(s0*r + k*s1 + b*s2 + s3*c)};
          const vector<T> myData = {val * (*_s).d_bra_op_ket(b, op, k)};
          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          if(rank==0)
            H.A[site].__T.write(1, inds.data(), myData.data());
          else
            H.A[site].__T.write(0, inds.data(), myData.data());
          //	  H.A[site].__T(r ,k,b,c) = val * (*_s).d_bra_op_ket(b, op, k);
        }
      }
    }
}
template void Heisenberg<double>::addOperators(MPO<double>& H, unsigned site, unsigned r, unsigned c, string op, double val);
template void Heisenberg< std::complex<double> >::addOperators(MPO< std::complex<double> >& H, unsigned site, unsigned r, unsigned c, string op, double val);


template <typename T>
void Heisenberg<T>::addOperators(qMPO<T>& H, unsigned site, unsigned r, unsigned c, string op, double val, QN_t Qi, QN_t Qo){
  bool added = false;
  for (size_t k = 0; k < H.phy_dim; k++) {
    for (size_t b = 0; b < H.phy_dim; b++) {
      string qn_str;
      qn_str += (to_string(Qi)+" ");
      qn_str += (to_string(H.phy_qn[k])+" ");
      qn_str += (to_string(H.phy_qn[b])+" ");
      qn_str += (to_string(Qo)+" ");
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if(H.A[site].block_id_by_qn_str.count(qn_str) > 0){
        added=true;
        unsigned block_pos = H.A[site].block_id_by_qn_str[qn_str];
        unsigned s0 = 1;
        unsigned s1 = s0 * H.A[site].idx_set[0].qdim(H.A[site].block_index_qi[block_pos][0]);
        unsigned s2 = s1 * H.A[site].idx_set[1].qdim(H.A[site].block_index_qi[block_pos][1]);
        unsigned s3 = s2 * H.A[site].idx_set[2].qdim(H.A[site].block_index_qi[block_pos][2]);

        const vector<T> myData = {val * (*_s).d_bra_op_ket(b, op, k)};
        if(site==0){
          assert(s3*c < H.A[site]._block[block_pos].get_tot_size(false));
          vector<long int> inds = {static_cast<long int>(s3*c)};
          //H.A[site].block[block_pos][s3*c] = val * (*_s).d_bra_op_ket(b, op, k);
          if(rank==0)
            H.A[site]._block[block_pos].write(1, inds.data(), myData.data());
          else
            H.A[site]._block[block_pos].write(0, inds.data(), myData.data());
        }else{
          assert(s0*r + s3*c < H.A[site]._block[block_pos].get_tot_size(false));
          //H.A[site].block[block_pos][s0*r + s3*c] = val * (*_s).d_bra_op_ket(b, op, k);
          vector<long int> inds = {static_cast<long int>(s0*r+s3*c)};
          if(rank==0)
            H.A[site]._block[block_pos].write(1, inds.data(), myData.data());
          else
            H.A[site]._block[block_pos].write(0, inds.data(), myData.data());
        }
      }
    }
  }
  //if(added) perr<<"Adding @"<<site<<" r="<<r<<" c="<<c<<" op="<<op<<" val="<<val<<" Qi="<<to_string(Qi)<<" Qo="<<to_string(Qo)<<endl;
  //else perr<<"NOT Adding @"<<site<<" r="<<r<<" c="<<c<<" op="<<op<<" val="<<val<<" Qi="<<to_string(Qi)<<" Qo="<<to_string(Qo)<<endl;
  assert(added);
}
template void Heisenberg<double>::addOperators(qMPO<double>& H, unsigned site, unsigned r, unsigned c, string op, double val, QN_t Qi, QN_t Qo);
template void Heisenberg< std::complex<double> >::addOperators(qMPO< std::complex<double> >& H, unsigned site, unsigned r, unsigned c, string op, double val, QN_t Qi, QN_t Qo);


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
void Heisenberg<T>::buildHam(AutoMPO& ampo, MPO<T>&  H){
  auto N = (*_s).N();
  //partially taken from ITensor(v2) AutoMPO.cc
  auto IL = SiteTerm("IL",0);
  auto HL = SiteTerm("HL",0);
  vector<vector<SiteQN>> basis(N+1);

  for(unsigned n=0;n<N;++n){
    basis.at(n).emplace_back(IL,0);
  }
  for(unsigned n=1;n<=N;++n){
    basis.at(n).emplace_back(HL,0);
  }
  const QN_t Zero = 0; //QN type
  
  //Fill up the basis array at each site with the unique operator types occuring on the site
  //unique including their coefficient
  //and starting a string of operators (i.e. first op of an HTerm)
  for(auto& ht : ampo.terms()){
    for(auto n=ht.first().i+1; n<= ht.last().i+1;++n){
      auto& bn = basis.at(n);
      auto test_has_first = [&ht](SiteQN const& sq){return sq.st == ht.first();};
      bool has_first = (std::find_if(bn.begin(),bn.end(),test_has_first)!=bn.end() );
      if(!has_first){
        if(true){
          bn.emplace_back(ht.first(),- _s->div(ht.first().op));
        }else { bn.emplace_back(ht.first(),Zero); }
      }
    }
  }
  if(true){
    auto qn_comp = [&Zero](const SiteQN& sq1,const SiteQN& sq2){
                    //first two if statements are to artificially make
                    //Zero QN come first in the sort
                    if(sq1.q==Zero      && sq2.q != Zero) return true;
                    else if(sq2.q==Zero && sq1.q !=Zero)  return false;
                    return sq1.q < sq2.q;
                  };
    for(auto& bn : basis) std::sort(bn.begin(),bn.end(),qn_comp);
  }

  auto links = vector<vector<quantum_number>>(N+1);
  // first: qn;    second: dimension
  auto inqn = vector<quantum_number>{}; //IndexQN
  for(unsigned n=0;n<=N;++n){
    auto& bn = basis.at(n);
    inqn.clear();
    QN_t currq = bn.front().q;
    unsigned currm = 0;
    int count = 0;
    for(auto& sq : bn){
      if(sq.q == currq){ ++currm; }
      else{
        inqn.emplace_back(currq,currm);
        currq = sq.q;
        currm = 1;
      }
    }
    //TODO: make more effecient
    inqn.emplace_back(currq,currm);
    links.at(n) = std::move(inqn);
  }
  //create arrays indexed by lattice sites
  //for lattice sites "j", ht_by_n[j] contains all HTerms, operator strings
  //which begin on, end on, or cross site "j"
  auto ht_by_n = vector<vector<HTerm>>(N+1);
  for(auto& ht : ampo.terms()){
    for(auto& st: ht.ops)
      ht_by_n.at(st.i+1).push_back(ht);
  }
  vector<unsigned> bdList(N+1);
  for(unsigned i=0;i<N+1;i++) bdList[i] = basis.at(i).size();
  bdList[0] = 1; bdList.back()=1;
  unsigned maxbd = *std::max_element(bdList.begin(),bdList.end());
  perr<<"Warning! MPO max bond dim is "<<maxbd<<endl;
  MPO<T> A(_s, bdList);
  A.setZero();
  H = A;
  assert(ht_by_n[0].size()==0);
  for (size_t i = 0; i < H.length; i++) {
    auto& bn1 = basis.at(i);
    auto& bn = basis.at(i+1);
    unsigned c_max = (i==N-1)?1:bn.size(); //itensor bug/feature where the last tensor isn't a vector
    for(unsigned r_=0;r_<bn1.size();r_++){
      for(unsigned c_=0;c_<c_max;c_++){
        auto& rst = bn1.at(r_).st;
        auto& cst = bn.at(c_).st;
        unsigned r = (i==0)? 1 : r_; //technically row zero but we define it as 1
        unsigned c = c_;//(i==N-1)? c_: c_;
        //start a new operator string
        if(cst.i==(i) && rst == IL){
          auto op = startTerm(cst.op);
          if(op!="HL" && op!="IL"){ addOperators(H, i, r, c, op, 1.0); }
        }
        if(cst == rst){ 
          if(isFermionic(cst))  addOperators(H, i, r, c, "F", 1.0);
          else                  addOperators(H, i, r, c, "Id", 1.0);
        }
        if(cst==HL){
          for(const auto& ht:  ht_by_n.at(i+1)){
            if(rst==ht.first() && ht.last().i==(i)){
              auto op = endTerm(ht.last().op);
              addOperators(H,i,r,c,op,ht.coef.real());
            }
          }
        }
        if(rst ==IL && cst ==HL){
          for(const auto& ht : ht_by_n.at(i+1)){
            if(ht.first().i==ht.last().i){
              addOperators(H,i,r,c,ht.first().op,ht.coef.real());
            }
          }
        }
      }
    }//end rc loop
  }
}
template void Heisenberg<double>::buildHam(AutoMPO& ampo,MPO<double>&  H);
template void Heisenberg< std::complex<double> >::buildHam(AutoMPO& ampo,MPO< std::complex<double> >&  H);
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

  //vector<unordered_map<int,unsigned>> qnToMaxIdx(N+1);
  //vector<unordered_map<unsigned,unsigned>> rcToIdx(N+1);
unsigned basisToIndex(QN_t qni,unordered_map<QN_t,unsigned>& qnToMaxIdx,
                      unordered_map<unsigned,unsigned>& rcToIdx, unsigned i){
  if(qnToMaxIdx.count(qni)==0) qnToMaxIdx[qni] = 1;
  if(rcToIdx.count(i)==0){
    rcToIdx[i] = qnToMaxIdx[qni]++; //always make something 1 at minimum
  }
  return rcToIdx[i]-1;
}
template <typename T>
void Heisenberg<T>::buildHam(AutoMPO& ampo, qMPO<T>& H){
  auto N = (*_s).N();
  qMPO<T> A(_s);
  A.setZero();
  H = A;
  // Resize qMPO
  H.A.clear();
  string Link_name_pref = "ID"+to_string(H._id)+"Link";
  string Site_name_pref = "Site";
  /*for (size_t i = 0; i < H.length; i++) {
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
  }*/
  // Setup MPO
  /*for (size_t site = 0; site < H.length; site++) {
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
  }*/
  auto IL = SiteTerm("IL",0);
  auto HL = SiteTerm("HL",0);
  vector<vector<SiteQN>> basis(N+1);

  for(unsigned n=0;n<N;++n){
    basis.at(n).emplace_back(IL,0);
  }
  for(unsigned n=1;n<=N;++n){
    basis.at(n).emplace_back(HL,0);
  }
  const QN_t Zero = 0; //QN type
  
  //Fill up the basis array at each site with the unique operator types occuring on the site
  //unique including their coefficient
  //and starting a string of operators (i.e. first op of an HTerm)
  for(auto& ht : ampo.terms()){
    for(auto n=ht.first().i+1; n<= ht.last().i+1;++n){
      auto& bn = basis.at(n);
      auto test_has_first = [&ht](SiteQN const& sq){return sq.st == ht.first();};
      bool has_first = (std::find_if(bn.begin(),bn.end(),test_has_first)!=bn.end() );
      if(!has_first){
        if(true){
          bn.emplace_back(ht.first(),- _s->div(ht.first().op));
        }else { bn.emplace_back(ht.first(),Zero); }
      }
    }
  }
  if(true){
    auto qn_comp = [&Zero](const SiteQN& sq1,const SiteQN& sq2){
                    //first two if statements are to artificially make
                    //Zero QN come first in the sort
                    if(sq1.q==Zero && sq2.q != Zero) return true;
                    else if(sq2.q==Zero && sq1.q !=Zero) return false;
                    return sq1.q < sq2.q;
                  };
    for(auto& bn : basis) std::sort(bn.begin(),bn.end(),qn_comp);
  }

  auto links = vector<vector<quantum_number>>(N+1);
  // first: qn;    second: dimension
  auto inqn = vector<quantum_number>{}; //IndexQN
  for(unsigned n=0;n<=N;++n){
    auto& bn = basis.at(n);
    inqn.clear();
    QN_t currq = bn.front().q;
    unsigned currm = 0;
    int count = 0;
    for(auto& sq : bn){
      if(sq.q == currq){ ++currm; }
      else{
        inqn.emplace_back(currq,currm);
        currq = sq.q;
        currm = 1;
      }
    }
    //TODO: make more effecient
    inqn.emplace_back(currq,currm);
    links.at(n) = std::move(inqn);
  }
  //create arrays indexed by lattice sites
  //for lattice sites "j", ht_by_n[j] contains all HTerms, operator strings
  //which begin on, end on, or cross site "j"
  auto ht_by_n = vector<vector<HTerm>>(N+1);
  for(auto& ht : ampo.terms()){
    for(auto& st: ht.ops)
      ht_by_n.at(st.i+1).push_back(ht);
  }
  vector<unordered_map<QN_t,unsigned>> qnToMaxIdx(N+1);
  vector<unordered_map<unsigned,unsigned>> rcToIdx(N+1);
  assert(ht_by_n[0].size()==0);
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
    for(auto& row_qn : links.at(i)){
      H.A[i].addQNtoIndex(0,row_qn);
    }
    if(i!=N-1)
      for(auto& col_qn : links.at(i+1)){
        H.A[i].addQNtoIndex(3,col_qn);
      }
    else
      H.A[i].addQNtoIndex(3,std::make_pair(0,1));
    // Site
    for (size_t j = 0; j < H.phy_qn.size(); j++) {
      H.A[i].addQNtoIndex(1, std::make_pair(H.phy_qn[j], 1));
      H.A[i].addQNtoIndex(2, std::make_pair(H.phy_qn[j], 1));
    }
    H.A[i].initBlock();
    H.A[i].setZero();

    auto& bn1 = basis.at(i);
    auto& bn = basis.at(i+1);
    unsigned c_max = (i==N-1)?1:bn.size(); //itensor bug/feature where the last tensor isn't a vector
    for(unsigned r_ =0;r_<bn1.size();r_++){
      for(unsigned c_=0;c_<c_max;c_++){
        auto& rst = bn1.at(r_).st;
        auto& cst = bn.at(c_).st;
        auto qr  = -bn1.at(r_).q;
        auto qc  = -bn.at(c_).q;
        unsigned r = (i==0)? 1 : r_; //technically row zero but we define it as 1
        unsigned c = c_;//(i==N-1)? c_: c_;
        //start a new operator string
        if(cst.i==(i) && rst == IL){
          auto op = startTerm(cst.op);
          if(op!="HL" && op!="IL"){
          unsigned thisR = basisToIndex(-qr,qnToMaxIdx.at(i),rcToIdx.at(i),r);
          unsigned thisC = basisToIndex(qc,qnToMaxIdx.at(i+1),rcToIdx.at(i+1),c);
          if(i==0){ thisR = 1;  }
          addOperators(H, i, thisR, thisC, op, 1.0, qr,qc);
          }
        }
        if(cst == rst){
          unsigned thisR = basisToIndex(-qr,qnToMaxIdx.at(i),rcToIdx.at(i),r);
          unsigned thisC = basisToIndex(qc,qnToMaxIdx.at(i+1),rcToIdx.at(i+1),c);
            if(isFermionic(cst)) addOperators(H, i, thisR, thisC, "F", 1.0, qr, qc);
            else                 addOperators(H, i, thisR, thisC, "Id", 1.0, qr, qc);
        }
        if(cst==HL){
          for(const auto& ht:  ht_by_n.at(i+1)){
            if(rst==ht.first() && ht.last().i==(i)){
              auto op = endTerm(ht.last().op);
              unsigned thisR = basisToIndex(-qr,qnToMaxIdx.at(i),rcToIdx.at(i),r);
              unsigned thisC = basisToIndex(qc,qnToMaxIdx.at(i+1),rcToIdx.at(i+1),c);
              addOperators(H,i,thisR,thisC,op,ht.coef.real(),qr,qc);
            }
          }
        }
        if(rst ==IL && cst ==HL){
          for(const auto& ht : ht_by_n.at(i+1)){
            if(ht.first().i==ht.last().i){
              unsigned thisR = basisToIndex(-qr,qnToMaxIdx.at(i),rcToIdx.at(i),r);
              unsigned thisC = basisToIndex(qc,qnToMaxIdx.at(i+1),rcToIdx.at(i+1),c);
              addOperators(H,i,thisR,thisC,ht.first().op,ht.coef.real(),qr,qc); 
            }
          }
        }
      }
    }//end rc loop
  }
}
template void Heisenberg<double>::buildHam(AutoMPO& ampo, qMPO<double>& H);
template void Heisenberg< std::complex<double> >::buildHam(AutoMPO& ampo, qMPO< std::complex<double> >& H);
#endif
