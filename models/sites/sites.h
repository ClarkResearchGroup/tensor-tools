#ifndef ABSTRACT_SITES_HEADER
#define ABSTRACT_SITES_HEADER

#include "../../util/types_and_headers.h"

class abstract_sites{
public:

  abstract_sites(unsigned L, unsigned pD) : _N(L) , _phyDim(pD) {}
  ~abstract_sites(){}

  virtual double d_bra_op_ket(unsigned bra, string op, unsigned ket) = 0;
  virtual std::complex<double> c_bra_op_ket(unsigned bra, string op, unsigned ket) = 0;
  virtual uint_vec product_state(str_vec& st) = 0;
  virtual QN_t div(string op)=0;

  unsigned N() {return _N;}
  unsigned phy_dim() {return _phyDim;}
  vector<QN_t> phy_qn() {return _phyQN;}

protected:

  unsigned _N;
  unsigned _phyDim;
  vector<QN_t> _phyQN;
};

#endif
