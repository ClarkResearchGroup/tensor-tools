#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "../dtensor/dtensor_all.h"
#include "../qtensor/qtensor_all.h"
#include "../mps/mps_all.h"

#include "../models/sites/spinhalf.h"
#include "../models/hams/Heisenberg.h"

#include "../algos/dmrg/dmrg.h"

#include "../util/timer.h"
#include <ctf.hpp>
using namespace std;

#include "tblis.h"
typedef vector<tblis::len_type>    len_vec;
typedef vector<tblis::label_type>  lab_vec;
#include <Eigen/Dense>
#include "hdf5.h"
#ifdef USE_HPTT
#include <hptt.h>
#endif

filebuf fbuf;
ostream perr(&fbuf);
filebuf fbufo;
ostream pout(&fbufo);

typedef pair<int, unsigned> quantum_number;
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  unsigned N = 10;
  spinhalf sites(N);
  str_vec ps;
  for (size_t i = 0; i < sites.N(); i++) {
    if(i%2==0)
      ps.push_back("Dn");
    else
      ps.push_back("Up");
  }

  //std::cout.precision(8);
  //std::cerr.precision(8);

  {
    CTF::World world(argc,argv);
    if (world.rank == 0){
      perr.rdbuf(cerr.rdbuf());
      pout.rdbuf(cout.rdbuf());
      perr.precision(8);
      pout.precision(8);
    }
    perr << "\n" << "Test MPS DMRG." << '\n';
    MPS< double > psi(&sites,2);
    ///    psi.setRandom();
    for (int i=0;i<psi.length;i++){
      psi.A[i].setOne();
    }
    /*auto psi2 = psi;
    cerr<<"NORM:"<<psiphi(psi,psi)<<endl;
    psi.normalize();
    cerr<<"NORM:"<<psiphi(psi,psi)<<endl;
    cerr<<"NORM:"<<psiphi(psi,psi2)<<endl;*/

    //exit(1);
    
    MPO< double > H;
    Heisenberg< double > HB(&sites);
    HB.buildHam(H);
    //H.print(2);
    auto num   = psiHphi(psi,H,psi);
    auto denom = psiphi(psi,psi);
    perr<<"Initial overlap "<<num/denom<<" "<<num<<" "<<denom <<endl;
    int nsweeps = 20;
    int maxm = 60;
    double cutoff = 1e-8;
    auto finalEnergy = dmrg(psi, H, nsweeps, maxm, cutoff);
    if(world.rank==0) psi.print();

    pout << "Final Energy = "<<finalEnergy << '\n';
  }
  MPI_Finalize();
  return 0;
  exit(1);
  /*{
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
  }*/
  //------------------------------------
  return 0;
}
