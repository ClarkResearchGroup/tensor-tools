#ifndef ELECTRON_SITES_TYPE_HEADER
#define ELECTRON_SITES_TYPE_HEADER

#include "sites.h"

class electron : public abstract_sites
{
public:
  electron() : abstract_sites(0,4),Dn(0),Up(1),Emp(2),UpDn(3) {} 
  electron(unsigned L) : abstract_sites(L,4), Dn(0), Up(1),Emp(2),UpDn(3) {
    abstract_sites::_phyQN.push_back({-1,1});
    abstract_sites::_phyQN.push_back({+1,1});
    abstract_sites::_phyQN.push_back({0,0});
    abstract_sites::_phyQN.push_back({0,2});
  }
  ~electron(){}

  double d_bra_op_ket(unsigned bra, string op, unsigned ket);
  std::complex<double> c_bra_op_ket(unsigned bra, string op, unsigned ket);
  uint_vec product_state(str_vec& st);
  QN_t div(string op);

  const unsigned Dn;
  const unsigned Up;
  const unsigned Emp;
  const unsigned UpDn;
};

#endif
