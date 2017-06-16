#include "../../dtensor/tensor_index.cpp"
#include "../../dtensor/tensor_index_op.cpp"
#include "../../dtensor/dtensor.cpp"
#include "../../mps/tt.cpp"
#include "../../mps/mps.cpp"
#include "../../mps/mpo.cpp"
#include "../../mps/observables.cpp"
#include "../../mps/mps_mpo_ops.cpp"
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

  std::cout << psiphi(psi,psi) << " " << psiphi(phi,phi) << " " << psiphi(psi,phi) << '\n';
  psi.fit(phi,4,1e-6);
  std::cout << psiphi(psi,psi) << " " << psiphi(phi,phi) << " " << psiphi(psi,phi) << '\n';

  MPO<double> W(L,xs,5);
  W.fit(H,4,1e-6);
  MPO<double> V;
  V = H - W;
  std::cout << H.norm() << " " << W.norm() << " " << V.norm() << '\n';

  fitApplyMPO(H,psi,phi);
  std::cout << psiHphi(psi,H,psi) << " " << psiphi(psi,phi) << '\n';

  H.print();
  fitApplyMPO(H, H, W);
  W.print();
  V.setMPO(L,xs,1); V.setIdentity();
  V = V+W;
  V.print();
  std::cout << "trace(V) = " << trace(V) << '\n';
  //------------------------------------
  return 0;
}
