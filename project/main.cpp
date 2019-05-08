#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "../dtensor/dtensor_all.h"
#include "../qtensor/qtensor_all.h"
#include "../mps/mps_all.h"

#include "../models/sites/spinhalf.h"
#include "../models/hams/Heisenberg.h"

#include "../algos/dmrg/dmrg.h"

#include "../util/timer.h"

#include "../models/hams/AutoMPO.h"

#include "../models/lattice/square.h"

#include <ctf.hpp>
using namespace std;

//#include "tblis.h"
//typedef vector<tblis::len_type>    len_vec;
//typedef vector<tblis::label_type>  lab_vec;
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
  auto Nx=6;
  auto Ny=4;
  auto N =Nx*Ny;
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
    psi.setRandom();
    psi.normalize();
    perr<<"norm"<<endl;
    //    for (int i=0;i<psi.length;i++){
    //      psi.A[i].setOne();
    //    }
    /*auto psi2 = psi;
    cerr<<"NORM:"<<psiphi(psi,psi)<<endl;
    psi.normalize();
    cerr<<"NORM:"<<psiphi(psi,psi)<<endl;
    cerr<<"NORM:"<<psiphi(psi,psi2)<<endl;*/

    //exit(1);

    /*AutoMPO ampo(sites);
    for(int i=0;i<(int)N-1;i++){
      ampo+=0.5,"S+",i,"S-",i+1;
      ampo+=0.5,"S-",i,"S+",i+1;
      ampo+=    "Sz",i,"Sz",i+1;
    }
    for(int i=0;i<(int)N-3;i++){
      ampo+=0.5,"S+",i,"S-",i+3;
      ampo+=0.5,"S-",i,"S+",i+3;
      ampo+=    "Sz",i,"Sz",i+3;
    }
    for(int i=0;i<(int)N-2;i++){
      ampo+=-0.5,"S+",i,"S-",i+2;
      ampo+=-0.5,"S-",i,"S+",i+2;
      ampo+=-1.0,"Sz",i,"Sz",i+2;
    }
    for(int i=0;i<(int)N;i++)
      ampo+= -0.1,"Sz",i;
    ampo+=0.2,"Sz",1;*/

    //perr<<ampo<<endl;
    bool yperiodic=true;
    //auto lattice = squareLattice(Nx,Ny,yperiodic);
    auto lattice = squareNextNeighbor(Nx,Ny,yperiodic);

    auto J1 = 1.0;
    auto J2 = 0.5;

    AutoMPO ampo(sites);
    for(auto bnd: lattice){
      if(bnd.type=="1"){
        ampo+=J1/2,"S+",bnd.s1-1,"S-",bnd.s2-1;
        ampo+=J1/2,"S-",bnd.s1-1,"S+",bnd.s2-1;
        ampo+=J1  ,"Sz",bnd.s1-1,"Sz",bnd.s2-1;
      }
      if(bnd.type=="2"){
        ampo+=J2/2,"S+",bnd.s1-1,"S-",bnd.s2-1;
        ampo+=J2/2,"S-",bnd.s1-1,"S+",bnd.s2-1;
        ampo+=J2  ,"Sz",bnd.s1-1,"Sz",bnd.s2-1;
      }
    }
    
    MPO< double > H;
    //vector<double> h = {0,0.5,0,0.5,0,0.7,0.1,0.5,0.2,0.6,0.2,0.1,0.6,0.8,0.9,0.1};
    Heisenberg< double > HB(&sites);
    HB.buildHam(ampo,H);
    //HB.buildHam(H);
    //H.print(2);
    //Ha.print(2);
    auto num   = psiHphi(psi,H,psi);
    auto denom = psiphi(psi,psi);
    perr<<"Initial overlap "<<num/denom<<" "<<num<<" "<<denom <<endl;
    bool RUN = 1;
    int nsweeps = 2;
    /*int maxm = 60;
    double cutoff = 1e-8;
    int nsweeps = 5;*/
        std::vector<int> maxm = {50,80,100,150,150,200,200,200};
    //std::vector<int> maxm = {1000,1000,1000,1000,1000,1000,1000,1000};
    //    std::vector<int> maxm = {20000};
    std::vector<double> cutoff = {1E-6,1E-8,1E-10,1E-12,1E-12,1E-12,1E-12,1E-12};
    //    std::vector<double> cutoff = {0}; 
    std::vector<int> max_restart = {2};
    //TODO: noise
    if(RUN){
      auto finalEnergy = dmrg(psi, H, nsweeps, maxm, cutoff,max_restart);
      if(world.rank==0) psi.print();
      pout << "Final Energy = "<<finalEnergy << '\n';
      psi.save("output");
      pout <<"Saved!"<<endl;
    }
    else{
      MPS<double> psi2;
      psi2.load("output.h5");
      auto fe = psiHphi(psi2,H,psi2);
      perr<<"Loaded E ="<<fe<<endl;
      //do measurements
    }
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
