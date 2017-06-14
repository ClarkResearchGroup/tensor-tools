#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>

#include "../../mps/tt.cpp"
// #include "../../mps/mps.cpp"
// #include "../../mps/mpo.cpp"
// #include "../../mps/observables.cpp"
// #include "../../models/hams.cpp"
#include "../../util/ezh5.cpp"
#include "../../linalg/lapack_wrapper.cpp"

using namespace std;

int main(int argc, char const *argv[]) {
  cout<<"//------------------------------------"<<endl;
  cout<<"This program tests the TensorTrain base class,"<<endl;
  cout<<"and the derived MPS and MPO classes."<<endl;
  cout<<"//------------------------------------"<<endl;
  //------------------------------------
  int L = 10, bd = 12, xs = 2;
  //------------------------------------
  TensorTrain< std::complex<double> > t1;
  t1.setLength(L);
  t1.setBondDim(bd);
  t1.setIndexSize(xs);
  t1.allocateTensors();
  t1.setRandom();
  std::cout << t1.TTID << '\n';
  std::cout << t1.norm() << '\n';
  std::cout << t1.norm() << '\n';
  std::cout << t1.norm() << '\n';
  std::cout << t1.norm() << '\n';
  t1.normalize();
  t1.lc();
  std::cout << t1.norm() << '\n' << '\n';
  t1.normalize();
  t1.moveRight(0);
  t1.moveRight(1);
  t1.moveRight(2);
  t1.moveRight(3);
  t1.moveRight(4);
  t1.moveRight(5);
  t1.moveRight(6);
  t1.moveRight(7);
  t1.moveRight(8);
  std::cout << t1.norm() << '\n';
  t1.moveLeft(9);
  t1.moveLeft(8);
  t1.moveLeft(7);
  t1.moveLeft(6);
  t1.moveLeft(5);
  t1.moveLeft(4);
  t1.moveLeft(3);
  t1.moveLeft(2);
  t1.moveLeft(1);
  std::cout << t1.norm() << '\n';
  t1.print();
  t1.save("test.h5");
  TensorTrain< std::complex<double> > t2(std::move(t1));
  std::cout << t2.TTID << '\n';
  t2.print();
  t1.load("test.h5");
  t1.print();

  // MPO H(L,xs,5);
  // buildHeisenberg(H);
  // psi.print(0);
  // std::cout << psiHphi(psi,H,psi) << '\n';
  // std::cout << psiphi(psi,psi) << '\n';
  // std::cout << trace(H) << '\n';
  // std::cout << var(H) << '\n';
  // std::vector<double> v;
  // EE(H,v);
  // for(auto val : v) std::cout << val << " ";
  // std::cout << '\n';
  //------------------------------------
  return 0;
}
