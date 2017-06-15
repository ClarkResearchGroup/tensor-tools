#include "../../dtensor/tensor_index.cpp"
#include "../../dtensor/tensor_index_op.cpp"
#include "../../dtensor/dtensor.cpp"
#include "../../mps/tt.cpp"
#include "../../mps/mps.cpp"
#include "../../mps/mpo.cpp"
#include "../../mps/observables.cpp"
#include "../../models/hams.cpp"
#include "../../util/ezh5.cpp"
#include "../../linalg/lapack_wrapper.cpp"

using namespace std;

int main(int argc, char const *argv[]) {
  cout<<"//------------------------------------"<<endl;
  cout<<"This program tests the TensorTrain base class,"<<endl;
  cout<<"and the derived MPS and MPO classes."<<endl;
  cout<<"//------------------------------------"<<endl;
  //------------------------------------
  int L = 10, bd = 30, xs = 2;
  //------------------------------------
  MPS<double> psi(L,xs,bd);
  MPS<double> phi(L,xs,bd);
  MPO<double> H(L,xs,5);
  buildHeisenberg(H);

  std::cout << psi.norm() << " " << phi.norm() << '\n';
  std::cout << psiphi(psi,psi) << " " << psiphi(phi,phi) << " " << psiphi(psi,phi) << '\n';
  std::cout << psiHphi(psi,H,psi) << '\n';
  std::cout << psiHphi(phi,H,phi) << '\n';
  std::cout << psiHphi(psi,H,phi) << '\n';

  psi.fit(phi,4,1e-6);
  std::cout << psiphi(psi,psi) << " " << psiphi(phi,phi) << " " << psiphi(psi,phi) << '\n';

  MPO<double> W(L,xs,5);
  W.fit(H,4,1e-6);
  MPO<double> V;
  V = H - W;
  std::cout << H.norm() << " " << W.norm() << " " << V.norm() << '\n';
  //------------------------------------
  return 0;
}
