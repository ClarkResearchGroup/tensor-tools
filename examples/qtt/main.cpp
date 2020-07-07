/*
 * Copyright 2020 Ryan Levy, Xiongjie Yu, and Bryan K. Clark
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


#include "../../util/types_and_headers.h"

#include "../../linalg/lapack_wrapper.cpp"

#include "../../dtensor/dtensor_index.cpp"
#include "../../dtensor/dtensor_index_op.cpp"
#include "../../dtensor/dtensor.cpp"
#include "../../dtensor/dtensor_view.cpp"
#include "../../dtensor/dtensor_op.cpp"
#include "../../dtensor/big_dtensor.cpp"

#include "../../qtensor/qtensor_index.cpp"
#include "../../qtensor/qtensor_index_op.cpp"
#include "../../qtensor/qtensor.cpp"
#include "../../qtensor/qtensor_op.cpp"
#include "../../qtensor/big_qtensor.cpp"

#include "../../mps/tt.cpp"
#include "../../mps/qtt.cpp"
#include "../../mps/observables.cpp"
#include "../../mps/mps_mpo_methods.cpp"

#include "../../models/sites/sites.h"
#include "../../models/sites/spinhalf.cpp"
#include "../../models/hams/Heisenberg.cpp"

#include "../../linalg/tensor_davidson.cpp"
#include "../../algos/dmrg/dmrg.cpp"

#include "../../util/ezh5.cpp"
#include "../../util/timer.h"
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

  // {
  //   std::cout << "\n" << "Test qMPS product state initialization and save/load." << '\n';
  //
  //   qMPS<double> psi(&sites, ps);
  //   psi.print();
  //   psi.save("test_qmps.h5");
  //
  //   qMPS<double> phi;
  //   phi.load("test_qmps.h5");
  //   phi.print();
  // }
  //
  // {
  //   std::cout << "\n" << "Test qMPS norm() method." << '\n';
  //   qMPS<double> psi(&sites, 2);
  //   psi.setRandom();
  //   psi.print();
  //   std::cout << "Norm = " << psi.norm() << '\n';
  //   psi.rc();
  //   std::cout << "Norm = " << psi.norm() << '\n';
  //   psi.normalize();
  //   std::cout << "Norm = " << psi.norm() << '\n';
  //   psi.position(1);
  //   std::cout << "Norm = " << psi.norm() << '\n';
  //   psi.position(7);
  //   std::cout << "Norm = " << psi.norm() << '\n';
  //   psi.position(4);
  //   std::cout << "Norm = " << psi.norm() << '\n';
  //   psi.print();
  // }
  //
  // {
  //   std::cout << "\n" << "Test basic MPS observables." << '\n';
  //
  //   MPS<double> psi(&sites,20);
  //   psi.setRandom();
  //   psi.lc();
  //   psi.normalize();
  //   psi.print();
  //   std::cout << psi.norm() << '\n';
  //
  //   MPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   double t1 = get_clock();
  //   std::cout << psiHphi(psi,H,psi) << '\n';
  //
  //   std::cout << psiHKphi(psi,H,H,psi) << '\n';
  //   double t2 = get_clock();
  //
  //   std::cout << "MPS time = " << t2-t1 << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test basic qMPS observables." << '\n';
  //
  //   qMPS<double> psi(&sites,0);
  //   psi.setRandom();
  //   psi.lc();
  //   psi.normalize();
  //   std::cout << psi.norm() << '\n';
  //   psi.print();
  //   std::cout << psiphi(psi,psi) << '\n';
  //
  //   qMPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   double t1 = get_clock();
  //   std::cout << psiHphi(psi,H,psi) << '\n';
  //
  //   std::cout << psiHKphi(psi,H,H,psi) << '\n';
  //   double t2 = get_clock();
  //   std::cout << "qMPS time = " << t2-t1 << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test MPS summation." << '\n';
  //   MPS<double> p1(&sites,10);
  //   p1.setRandom();
  //   p1.normalize();
  //   p1.rc();
  //
  //   MPS<double> p2(&sites,10);
  //   p2.setRandom();
  //   p2.normalize();
  //   p2.rc();
  //
  //   MPS<double> p3(&sites,10);
  //   p3.setRandom();
  //   p3.normalize();
  //   p3.rc();
  //
  //   MPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   MPS<double> res;
  //
  //   vector< MPS<double> > x(3);
  //   x.at(0) = p1;
  //   x.at(1) = p2;
  //   x.at(2) = p3;
  //
  //   sum(x, res, 6, 1e-12, 24);
  //   res.print();
  //   std::cout << psiphi(p1,p1)+psiphi(p1,p2)+psiphi(p1,p3) << '\n';
  //   std::cout << psiphi(p1,res) << '\n';
  //
  //   std::cout << psiHphi(p1,H,p1)+psiHphi(p1,H,p2)+psiHphi(p1,H,p3) << '\n';
  //   std::cout << psiHphi(p1,H,res) << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test qMPS summation." << '\n';
  //   qMPS<double> p1(&sites,0);
  //   p1.setRandom();
  //   p1.normalize();
  //   p1.rc();
  //
  //   qMPS<double> p2(&sites,0);
  //   p2.setRandom();
  //   p2.normalize();
  //   p2.rc();
  //
  //   qMPS<double> p3(&sites,0);
  //   p3.setRandom();
  //   p3.normalize();
  //   p3.rc();
  //
  //   qMPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   qMPS<double> res;
  //
  //   vector< qMPS<double> > x(3);
  //   x.at(0) = p1;
  //   x.at(1) = p2;
  //   x.at(2) = p3;
  //
  //   sum(x, res, 6, 1e-12, 50);
  //   res.print();
  //   std::cout << psiphi(p1,p1)+psiphi(p1,p2)+psiphi(p1,p3) << '\n';
  //   std::cout << psiphi(p1,res) << '\n';
  //
  //   std::cout << psiHphi(p1,H,p1)+psiHphi(p1,H,p2)+psiHphi(p1,H,p3) << '\n';
  //   std::cout << psiHphi(p1,H,res) << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test MPO mult." << '\n';
  //   MPS<double> psi(&sites,10);
  //   psi.setRandom();
  //   psi.normalize();
  //   psi.rc();
  //   psi.print();
  //   std::cout << psi.norm() << '\n';
  //   std::cout << psiphi(psi,psi) << '\n';
  //
  //   MPS<double> phi(&sites,10);
  //   phi.setRandom();
  //   phi.normalize();
  //   phi.lc();
  //   phi.print();
  //   std::cout << phi.norm() << '\n';
  //   std::cout << psiphi(phi,phi) << '\n';
  //
  //   MPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   MPO<double> HS;
  //   mult(H, H, HS, 6, 1e-12, 15);
  //   HS.print();
  //
  //   std::cout << psiHKphi(psi,H,H,psi) << '\n';
  //   std::cout << psiHKphi(phi,H,H,phi) << '\n';
  //   std::cout << psiHKphi(psi,H,H,phi) << '\n';
  //
  //   std::cout << psiHphi(psi,HS,psi) << '\n';
  //   std::cout << psiHphi(phi,HS,phi) << '\n';
  //   std::cout << psiHphi(psi,HS,phi) << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test qMPO mult." << '\n';
  //   qMPS<double> psi(&sites,0);
  //   psi.setRandom();
  //   psi.normalize();
  //   psi.rc();
  //   psi.print();
  //   std::cout << psi.norm() << '\n';
  //   std::cout << psiphi(psi,psi) << '\n';
  //
  //   qMPS<double> phi(&sites,0);
  //   phi.setRandom();
  //   phi.normalize();
  //   phi.lc();
  //   phi.print();
  //   std::cout << phi.norm() << '\n';
  //   std::cout << psiphi(phi,phi) << '\n';
  //
  //   qMPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   qMPO<double> HS;
  //   mult(H, H, HS, 6, 1e-12, 10);
  //   HS.print();
  //
  //   std::cout << psiHKphi(psi,H,H,psi) << '\n';
  //   std::cout << psiHKphi(phi,H,H,phi) << '\n';
  //   std::cout << psiHKphi(psi,H,H,phi) << '\n';
  //
  //   std::cout << psiHphi(psi,HS,psi) << '\n';
  //   std::cout << psiHphi(phi,HS,phi) << '\n';
  //   std::cout << psiHphi(psi,HS,phi) << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test MPO-MPS mult." << '\n';
  //   MPS<double> psi(&sites,10);
  //   psi.setRandom();
  //   psi.normalize();
  //   psi.rc();
  //   psi.print();
  //   std::cout << psi.norm() << '\n';
  //   std::cout << psiphi(psi,psi) << '\n';
  //
  //   MPS<double> phi(&sites,10);
  //   phi.setRandom();
  //   phi.normalize();
  //   phi.lc();
  //   phi.print();
  //   std::cout << phi.norm() << '\n';
  //   std::cout << psiphi(phi,phi) << '\n';
  //
  //   MPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   MPS<double> res;
  //   mult(H, phi, res, 6, 1e-12, 30);
  //   res.print();
  //
  //   std::cout << psiHphi(psi,H,phi) << '\n';
  //
  //   std::cout << psiphi(psi,res) << '\n';
  // }
  //
  // {
  //   std::cout << "\n" << "Test qMPO-qMPS mult." << '\n';
  //   qMPS<double> psi(&sites,0);
  //   psi.setRandom();
  //   psi.normalize();
  //   psi.rc();
  //   psi.print();
  //   std::cout << psi.norm() << '\n';
  //   std::cout << psiphi(psi,psi) << '\n';
  //
  //   qMPS<double> phi(&sites,0);
  //   phi.setRandom();
  //   phi.normalize();
  //   phi.lc();
  //   phi.print();
  //   std::cout << phi.norm() << '\n';
  //   std::cout << psiphi(phi,phi) << '\n';
  //
  //   qMPO<double> H;
  //   Heisenberg<double> HB(&sites);
  //   HB.buildHam(H);
  //
  //   qMPS<double> res;
  //   mult(H, phi, res, 8, 1e-30, 24);
  //   res.print();
  //
  //   std::cout << psiHphi(psi,H,phi) << '\n';
  //
  //   std::cout << psiphi(psi,res) << '\n';
  // }

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
    double cutoff = 2e-11;
    dmrg(psi, H, nsweeps, maxm, cutoff);

    psi.print();

    std::cout << psiHphi(psi,H,psi) << '\n';
  }
  //------------------------------------
  return 0;
}
