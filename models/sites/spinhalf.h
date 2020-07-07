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

#ifndef SPIN_HALF_SITES_TYPE_HEADER
#define SPIN_HALF_SITES_TYPE_HEADER

#include "sites.h"

class spinhalf : public abstract_sites
{
public:
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
