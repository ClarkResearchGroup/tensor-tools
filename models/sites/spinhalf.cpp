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

#ifndef SPIN_HALF_SITES_TYPE
#define SPIN_HALF_SITES_TYPE

#include "spinhalf.h"

double spinhalf::d_bra_op_ket(unsigned bra, string op, unsigned ket){
  assert(bra<2 && ket<2);
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
  }else{
    std::cout << "(spinhalf class) Operator " << op << " is not supported!" << '\n';
    return 0;
  }
}

std::complex<double> spinhalf::c_bra_op_ket(unsigned bra, string op, unsigned ket){
  assert(bra<2 && ket<2);
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
  }else if(op=="Sy"){
    if(ket==Dn && bra==Up) return std::complex<double>(0, -0.5);
    if(ket==Up && bra==Dn) return std::complex<double>(0, +0.5);
    return 0;
  }else if(op=="S+" || op=="Sp"){
    if(ket==Dn && bra==Up) return +1;
    return 0;
  }else if(op=="S-" || op=="Sm"){
    if(ket==Up && bra==Dn) return +1;
    return 0;
  }else{
    std::cout << "(spinhalf class) Operator " << op << " is not supported!" << '\n';
    return 0;
  }
}

uint_vec spinhalf::product_state(str_vec& st){
  unsigned m = 0;
  uint_vec ps;
  for(auto s : st){
    if(s=="Up"){
      ps.push_back(1);
    }else if(s=="Dn"){
      ps.push_back(0);
    }else{
      std::cout << "Input produst state string is incompatible with spinhalf sites!" << '\n';
      abort();
    }
    ++m;
  }
  assert(m == _N);
  return ps;
}

#endif
