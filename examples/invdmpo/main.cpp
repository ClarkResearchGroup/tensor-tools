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
#include "../../algos/invdmpo/invdmpo.cpp"
#include "../../linalg/tensor_cg.cpp"

using namespace std;

void print_diag_elems(MPO<double>& H){
  int L = H.length;
  int phys[L];
  for (int i = 0; i < std::pow(2,L); i++) {
    for (int j = 0; j < L; j++) {
      phys[j] = (i>>j)%2;
    }
    Mxd tp;
    for (int j = 0; j < L; j++) {
      if(j==0){
        tp = H.M[j][phys[j]*2+phys[j]];
      }else{
        tp = tp * H.M[j][phys[j]*2+phys[j]];
      }
    }
    std::cout << tp(0,0) << " ";
  }
  std::cout << '\n';
  std::cout << '\n';
}

int main(int argc, char const *argv[]) {
  cout<<"//------------------------------------"<<endl;
  cout<<"This program tests the TensorTrain base class,"<<endl;
  cout<<"and the derived MPS and MPO classes."<<endl;
  cout<<"//------------------------------------"<<endl;
  //------------------------------------
  int L = 10, bd = 20, xs = 2;
  //------------------------------------
  MPS<double> psi(L,xs,bd);
  MPS<double> phi(L,xs,bd);
  MPO<double> H(L,xs,5), Hd(L,xs,5);
  buildHeisenberg(H); Hd = diagonal(H);
  MPO<double> U(L,xs,bd), W, V(L,xs,1);

  std::cout << "l2norm(Hd) = " << l2norm(Hd) << '\n';
  std::cout << "trace(Hd)  = " << trace(Hd) << '\n';
  fitApplyMPO(Hd, Hd, W);
  W *= 1.5;
  std::cout << "l2norm(0.1Hd^2) = " << l2norm(W) << '\n';
  std::cout << "trace(0.1Hd^2)  = " << trace(W) << '\n';
  V = V + W;
  std::cout << "l2norm(I+0.1Hd^2) = " << l2norm(V) << '\n';
  std::cout << "trace(I+0.1Hd^2)  = " << trace(V) << '\n';
  U.print();
  std::cout << "l2norm(U) = " << l2norm(U) << '\n';
  std::cout << "trace(U)  = " << trace(U) << '\n';
  std::cout << "Starting to invert 1 + 0.1 H^2 ..." << '\n';
  invdmpo(V, U, 20, 1e-10, 1e-14);

  print_diag_elems(V);
  print_diag_elems(U);
  //------------------------------------
  return 0;
}
