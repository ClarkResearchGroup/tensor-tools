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

#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "../dtensor/dtensor_all.h"
#include "../qtensor/qtensor_all.h"
#include "../mps/mps_all.h"

#include "../models/sites/spinhalf.h"
#include "../models/hams/Heisenberg.h"

#include "../algos/dmrg/dmrg.h"

#include "../util/timer.h"
using namespace std;

int main(int argc, char const *argv[]) {
  unsigned N = 30;
  spinhalf sites(N);
  str_vec ps;
  for (size_t i = 0; i < sites.N(); i++) {
    if(i%2==0)
      ps.push_back("Dn");
    else
      ps.push_back("Up");
  }

  std::cout.precision(8);

  {
    std::cout << "\n" << "Test MPS DMRG." << '\n';
    MPS< double > psi(&sites,1);
    psi.setRandom();
    psi.normalize();

    MPO< double > H;
    Heisenberg< double > HB(&sites);
    HB.buildHam(H);

    int nsweeps = 20;
    int maxm = 60;
    double cutoff = 1e-8;
    dmrg(psi, H, nsweeps, maxm, cutoff);

    psi.print();

    std::cout << psiHphi(psi,H,psi) << '\n';
  }

  {
    std::cout << "\n" << "Test qMPS DMRG." << '\n';
    qMPS< double > psi(&sites,ps);
    psi.normalize();

    qMPO< double > H;
    Heisenberg< double > HB(&sites);
    HB.buildHam(H);

    int nsweeps = 20;
    int maxm = 60;
    double cutoff = 1e-8;
    dmrg(psi, H, nsweeps, maxm, cutoff);

    psi.print();

    std::cout << psiHphi(psi,H,psi) << '\n';
  }
  //------------------------------------
  return 0;
}
