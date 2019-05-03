#ifndef SPIN_HALF_SITES_TYPE_HEADER
#define SPIN_HALF_SITES_TYPE_HEADER

#include "sites.h"

class spinhalf : public abstract_sites
{
public:
  spinhalf() : abstract_sites(0,2),Dn(0),Up(1) {} 
  spinhalf(unsigned L) : abstract_sites(L,2), Dn(0), Up(1) {
    abstract_sites::_phyQN.push_back(-1);
    abstract_sites::_phyQN.push_back(+1);
  }
  ~spinhalf(){}

  double d_bra_op_ket(unsigned bra, string op, unsigned ket);
  std::complex<double> c_bra_op_ket(unsigned bra, string op, unsigned ket);
  uint_vec product_state(str_vec& st);

  const unsigned Dn;
  const unsigned Up;
};

#endif
