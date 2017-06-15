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
#include "../../algos/ddmrg/dmrg.cpp"
#include "../../algos/ddmrg/dmrg_diag.cpp"

using namespace std;

int main(int argc, char const *argv[]) {
  cout.precision(10);
  cout<<"//------------------------------------"<<endl;
  cout<<"This program tests the single-site DMRG method."<<endl;
  cout<<"It does not have conserved quantum numbers."<<endl;
  cout<<"//------------------------------------"<<endl;
  //------------------------------------
  int L = 10, bd = 30, xs = 2;
  //------------------------------------
  MPS< double > psi(L,xs,bd);
  MPO< double > H(L,xs,5);
  buildHeisenberg(H);
  //------------------------------------
  dmrg(psi, H, 100);
  //------------------------------------
  return 0;
}
