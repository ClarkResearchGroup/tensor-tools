#ifndef ELECTRON_SITES_TYPE
#define ELECTRON_SITES_TYPE

#include "electron.h"

double electron::d_bra_op_ket(unsigned bra, string op, unsigned ket){
  assert(bra<4 && ket<4);
  if(op=="Id" || op=="I"){
    if(ket==bra) return 1;
    return 0;
  }else if(op=="Sz"){
    if(ket==Dn && bra==Dn) return -0.5;
    if(ket==Up && bra==Up) return +0.5;
    return 0;
  }else if(op=="Sx"){
    if(ket==Dn && bra==Up) return +0.5;
    if(ket==Up && bra==Dn) return +0.5;
    return 0;
  }else if(op=="S+" || op=="Sp"){
    if(ket==Dn && bra==Up) return +1;
    return 0;
  }else if(op=="S-" || op=="Sm"){
    if(ket==Up && bra==Dn) return +1;
    return 0;
  }else if(op=="Nup"){
    if(ket==Up   && bra==Up)   return +1;
    if(ket==UpDn && bra==UpDn) return +1;
    return 0;
  }else if(op=="Ndn"){
    if(ket==Dn   && bra==Dn)   return +1;
    if(ket==UpDn && bra==UpDn) return +1;
    return 0;
  }else if(op=="Nupdn"){
    if(ket==UpDn && bra==UpDn) return +1;
    return 0;
  }else if(op=="Ntot"){
    if(ket==Up   && bra==Up)   return +1;
    if(ket==Dn   && bra==Dn)   return +1;
    if(ket==UpDn && bra==UpDn) return +2;
    return 0;
  }else if(op=="Cup"){
    if(ket==Up   && bra==Emp) return +1;
    if(ket==UpDn && bra==Dn)  return +1;
    return 0;
  }else if(op=="Cdagup"){
    if(ket==Emp && bra==Up)   return +1;
    if(ket==Dn  && bra==UpDn) return +1;
    return 0;
  }else if(op=="Cdn"){
    if(ket==Dn   && bra==Emp) return +1;
    if(ket==UpDn && bra==Up)  return -1;
    return 0;
  }else if(op=="Cdagdn"){
    if(ket==Emp && bra==Dn)    return +1;
    if(ket==Up  && bra==UpDn)  return -1;
    return 0;
  }else if(op=="Aup"){
    if(ket==Up   && bra==Emp) return +1;
    if(ket==UpDn && bra==Dn)  return +1;
    return 0;
  }else if(op=="Adagup"){
    if(ket==Emp && bra==Up)   return +1;
    if(ket==Dn  && bra==UpDn) return +1;
    return 0;
  }else if(op=="Adn"){
    if(ket==Dn   && bra==Emp) return +1;
    if(ket==UpDn && bra==Up)  return +1;
    return 0;
  }else if(op=="Adagdn"){
    if(ket==Emp && bra==Dn)    return +1;
    if(ket==Up  && bra==UpDn)  return +1;
    return 0;
  }else if(op=="F"){
    if(ket==Emp  && bra==Emp)  return +1;
    if(ket==Up   && bra==Up)   return -1;
    if(ket==Dn   && bra==Dn)   return -1;
    if(ket==UpDn && bra==UpDn) return +1;
    return 0;
  }else if(op=="Fup"){
    if(ket==Emp  && bra==Emp)  return +1;
    if(ket==Up   && bra==Up)   return -1;
    if(ket==Dn   && bra==Dn)   return +1;
    if(ket==UpDn && bra==UpDn) return -1;
    return 0;
  }else if(op=="Fdn"){
    if(ket==Emp  && bra==Emp)  return +1;
    if(ket==Up   && bra==Up)   return +1;
    if(ket==Dn   && bra==Dn)   return -1;
    if(ket==UpDn && bra==UpDn) return -1;
    return 0;
  }else if(op[0]=='F' && op[1]=='*'){
    string op2 = op.substr(2);
    double tot = 0.;
    for (unsigned k=0;k<phy_dim();++k)
      tot+= d_bra_op_ket(bra,"F",k)*d_bra_op_ket(k,op2,ket);
    return tot;
  }else if(op.size()>3 && op.substr( op.length() - 2 )=="*F"){
    string op2 = op.substr(0,op.length() - 2);
    double tot = 0.;
    for (unsigned k=0;k<phy_dim();++k)
      tot+= d_bra_op_ket(bra,op2,k)*d_bra_op_ket(k,"F",ket);
    return tot;
  }else{
    perr << "(electron class) Operator " << op << " is not supported!" << '\n';
    assert(1==2);
    return 0;
  }
}

std::complex<double> electron::c_bra_op_ket(unsigned bra, string op, unsigned ket){
  assert(bra<4 && ket<4);
  if(op=="Sy"){
    if(ket==Dn && bra==Up) return std::complex<double>(0, -0.5);
    if(ket==Up && bra==Dn) return std::complex<double>(0, +0.5);
    return 0;
  }else{
    return std::complex<double>(d_bra_op_ket(bra,op,ket),0.);
  }
}

uint_vec electron::product_state(str_vec& st){
  unsigned m = 0;
  uint_vec ps;
  for(auto s : st){
    if(s=="Up"){
      ps.push_back(1);
    }else if(s=="Dn"){
      ps.push_back(0);
    }else if(s=="Emp"){
      ps.push_back(2);
    }else if(s=="UpDn"){
      ps.push_back(3);
    }else{
      perr << "Input produst state string is incompatible with spinhalf sites!" << '\n';
      abort();
    }
    ++m;
  }
  assert(m == _N);
  return ps;
}
/*
 * Sz (Sz=0,Nf=0)
 * S+ (Sz=2,Nf=0)
 * S- (Sz=-2,Nf=0)
 * Id (Sz=0,Nf=0)
 * Nup (Sz=0,Nf=0)
 * Ndn (Sz=0,Nf=0)
 * Ntot (Sz=0,Nf=0)
 * Cup (Sz=-1,Nf=-1)
 * Cdagup (Sz=1,Nf=1)
 * Cdn (Sz=1,Nf=-1)
 * Cdagdn (Sz=-1,Nf=1)
 * Aup (Sz=-1,Nf=-1)
 * Adagup (Sz=1,Nf=1)
 * Adn (Sz=1,Nf=-1)
 * Adagdn (Sz=-1,Nf=1)
 * F (Sz=0,Nf=0)
 * Fup (Sz=0,Nf=0)
 * Fdn (Sz=0,Nf=0)
 * */
QN_t electron::div(string op){
  if(op=="Sz")
    return {0,0};
  if(op=="S+")
    return {2,0};
  if(op=="S-")
    return {-2,0};
  if(op=="Id")
    return 0;
  if(op=="Nup" || op=="Ndn" || op=="Ntot" || op=="Nupdn")
    return 0;
  if(op=="Cup")
    return {-1,-1};
  if(op=="Cdagup")
    return {1,1};
  if(op=="Cdn")
    return {1,-1};
  if(op=="Cdagdn")
    return {-1,1};
  
  perr<<"Bad string:"<<op<<std::endl;
  assert(1==2); //Not found op string!
  return 0;
}

#endif
