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


#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "../dtensor/dtensor_all.h"
#include "../qtensor/qtensor_all.h"
#include "../mps/mps_all.h"

#include "../models/sites/spinhalf.h"
#include "../models/sites/electron.h"
#include "../models/hams/Heisenberg.h"

#include "../algos/dmrg/dmrg.h"

#include "../util/timer.h"

#include "../models/hams/AutoMPO.h"

#include "../models/lattice/square.h"
#include "../models/lattice/triangular.h"

#include <ctf.hpp>
#include <sys/stat.h>
using namespace std;
inline bool check_existence(const string& name);

//#include <Eigen/Dense>
#include "hdf5.h"

filebuf fbuf;
ostream perr(&fbuf);
filebuf fbufo;
ostream pout(&fbufo);

template <typename T, template <typename,unsigned> class TensorType>
void fixTensors(TensorType<T,1>& psi,TensorType<T,2>& H){
  uint_vec perm;
  string Link_name_pref = "ID"+to_string(psi._id)+"Link";
  string Site_name_pref = "Site";
  for(unsigned l=0;l<psi.length;l++){
    string left_link_name  = Link_name_pref+to_string(l);
    string right_link_name = Link_name_pref+to_string(l+1);
    string site_name       = Site_name_pref+to_string(l);
    vector<qtensor_index> inds(3);
    for(auto& ix : psi.A[l].idx_set){
      if(ix.name()==left_link_name) inds[0] = ix;
      else if(ix.name()==site_name) inds[1] = ix;
      else if(ix.name()==right_link_name) inds[2] = ix;
    }
    find_index_permutation(psi.A[l].idx_set,inds, perm);
    psi.A[l].permute(perm);
  }
  //check if we need to dagger H
  string site_name = Site_name_pref+to_string(0);
  auto& psiSet = psi.A[0].idx_set;
  auto& HSet   = H.A[0].idx_set;
  auto it  = find_if(psiSet.begin(), psiSet.end(), [&](qtensor_index& qi){ return qi.name()==site_name;});
  auto itH = find_if(HSet.begin(),   HSet.end(),   [&](qtensor_index& qi){ return (qi.name()==site_name) && (qi.level()==0);});
  assert(it!=psiSet.end() && itH != HSet.end());
  if( it->arrow() == itH->arrow()){
    for(auto& Ai: H.A) Ai.dag();
  }

}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  {
    CTF::World world(argc,argv);
    if (world.rank == 0){
      perr.rdbuf(cerr.rdbuf());
      pout.rdbuf(cout.rdbuf());
      perr.precision(12);
      pout.precision(12);
    }
    if(argc==1) {
      perr<<"Error: Need filename!"<<endl;
      assert(1==2);
    }
    string inputName = string(argv[1]);
    ifstream in;
    in.open(inputName);
    assert(in);
    string line; stringstream linest;

    /*
     *              input file format
     * type
     * Nx Ny
     * J1 J2
     * fname
     * prefix
     * nsweeps
     * maxm maxm maxm
     * cutoff cutoff cutoff
     * restart restart restart
     *
     */

    string type; getline(in,type);
    getline(in,line); linest.clear(); linest.str(line);
    int Nx,Ny;
    linest >> Nx >> Ny;
    auto N = Nx*Ny;
    double J1,J2;
    getline(in,line); linest.clear(); linest.str(line); 
    perr << "J1 line=" << line << endl;
    linest >> J1 >> J2;
    string fname="";  
    getline(in,line); if(!line.empty()) fname=line;
    string pref=""; 
    getline(in,line); if(!line.empty()) pref=line; 
    int nsweeps = 0; 
    getline(in,line); linest.clear(); linest.str(line);
    linest >> nsweeps; 

    vector<int> maxm;
    vector<double> cutoff;
    vector<int> max_restart;
    getline(in,line); linest.clear(); linest.str(line);
    while(getline(linest,line,' ')){
      maxm.push_back(std::stoi(line));
    }
    getline(in,line); linest.clear(); linest.str(line);
    while(getline(linest,line,' ')){
      cutoff.push_back(std::stof(line));
    }
    getline(in,line); linest.clear(); linest.str(line);
    while(getline(linest,line,' ')){
      max_restart.push_back(std::stoi(line));
    }
    in.close();

    if(world.rank==0)
      printf("type=%s N=(%i,%i) J1=%f J2=%f\nfile=%s prefix=%s\nnsweeps=%i \n",
          type.c_str(),Nx,Ny,J1,J2,fname.c_str(),pref.c_str(),nsweeps);
    perr<<"maxm:";        for(auto m: maxm)        perr<<m<<" "; perr<<"\n";
    perr<<"cutoff:";      for(auto c: cutoff)      perr<<c<<" "; perr<<"\n";
    perr<<"max_restart:"; for(auto r: max_restart) perr<<r<<" "; perr<<"\n";
    //spinhalf sites(N);
    electron sites(N);
    str_vec ps;
    for (size_t i = 0; i < sites.N(); i++) {
      if(i%2==0)
        ps.push_back("Dn");
      else
        ps.push_back("Up");
    }

    AutoMPO ampo(sites);
    bool yperiodic = true;
    auto lattice = triangularLattice(Nx,Ny,yperiodic);
    auto t = 1.0;
    auto U = J2;

    for(auto bnd : lattice){
      //hopping terms in the 1DEG
      ampo+=-t,"Cdagup",bnd.s1-1,"Cup",bnd.s2-1;
      ampo+=-t,"Cdagdn",bnd.s1-1,"Cdn",bnd.s2-1;
      ampo+=-t,"Cdagup",bnd.s2-1,"Cup",bnd.s1-1;
      ampo+=-t,"Cdagdn",bnd.s2-1,"Cdn",bnd.s1-1;
    }
    //-------------------------
    for(int i=0;i<N;i++)
      ampo+=U,"Nupdn",i;

    if(type=="d"){
      pout << "\n" << "Dense MPS DMRG (Linear Heisenberg)" << '\n';
      
      MPS<double> psi(&sites,ps);
      //psi.load(fname); //TODO: fix pref
      psi.print();
      MPO< double > H;
      Heisenberg< double > HB(&sites);
      HB.buildHam(H);
      dmrg(psi, H, nsweeps, maxm, cutoff, max_restart);
      psi.print();

    }
    if(type=="q"){
      pout << "\n" << "qMPS Fermionic DMRG" << '\n';
      
      qMPS< double > psi(&sites,ps);
      bool yperiodic=true;
      auto lattice = triangularLattice(Nx,Ny,yperiodic);

      auto t = 1.0;
      auto U = J2;

      AutoMPO ampo(sites);
      for(auto bnd : lattice){
        //hopping terms in the 1DEG
        ampo+=-t,"Cdagup",bnd.s1-1,"Cup",bnd.s2-1;
        ampo+=-t,"Cdagdn",bnd.s1-1,"Cdn",bnd.s2-1;
        ampo+=-t,"Cdagup",bnd.s2-1,"Cup",bnd.s1-1;
        ampo+=-t,"Cdagdn",bnd.s2-1,"Cdn",bnd.s1-1;
      }
      //-------------------------
      for(int i=0;i<N;i++)
        ampo+=U,"Nupdn",i;

      string postfix = to_string(Nx)+"x"+to_string(Ny)+"_U"+to_string(U);
      auto sp = "psi_"+postfix;
      qMPO< double > H;
      Heisenberg< double > HB(&sites);
      perr<<"making H"<<endl;
      HB.buildHam(ampo,H);

      psi.load(fname,pref);
      psi.print();

      assert(H._id != psi._id);
      fixTensors(psi,H);

      dmrg(psi, H, nsweeps, maxm, cutoff, max_restart);

    }
    if(type=="qs"){
      pout << "\n" << "qsMPS Fermionic DMRG" << '\n';
      
      qsMPS<double> psi(&sites,ps);
      psi.load(fname,pref);
      psi.print();
      qMPO< double > H;
      Heisenberg< double > HB(&sites);
      HB.buildHam(ampo,H);
      qsMPO< double > Hq = H;

      assert(H._id != psi._id);
      fixTensors(psi,Hq);

      dmrg(psi, Hq, nsweeps, maxm, cutoff, max_restart);

      psi.print();
      for(int l=0;l<N;l++){
        auto& Al = psi.A[l];
        perr<<Al._T.nnz_tot<<","<<(double)Al._T.nnz_tot/(Al._T.get_tot_size(false))<<endl;
        for(size_t i=0;i<Al.block_index_qd.size();i++)
          perr<<"   ("<<Al.block_index_qd[i][0]<<","
              << Al.block_index_qd[i][1]<<","<<Al.block_index_qd[i][2]<<")\n";
      }
    }
    if(type=="qToqs"){
      pout << "\n" << "q to qsMPS Fermionic DMRG" << '\n';
      
      qMPO< double > H;
      Heisenberg< double > HB(&sites);
      HB.buildHam(ampo,H);
      H.load("H_36_1.h5");
      H.print();
      qsMPO< double > Hq = H;

      qMPS<double> psi(&sites);
      psi.load(fname,pref);
      psi.print();
      qsMPS<double> psiq = psi;
      psiq.print();

     
      assert(H._id != psi._id);
      fixTensors(psiq,Hq);
      
      dmrg(psiq, Hq, nsweeps, maxm, cutoff, max_restart);

      psiq.print();
      /*for(int l=0;l<N;l++){
        auto& Al = psiq.A[l];
        perr<<Al._T.nnz_tot<<","<<(double)Al._T.nnz_tot/(Al._T.get_tot_size(false))<<endl;
        perr<<"   ";
        for(size_t i=0;i<Al.block_index_qd.size();i++)
          perr<<"("<<Al.block_index_qd[i][0]<<","
              << Al.block_index_qd[i][1]<<","<<Al.block_index_qd[i][2]<<") ";
        perr<<'\n';
      }*/
    }
  }
  MPI_Finalize();
  return 0;
  exit(1);
  //------------------------------------
  return 0;
}

inline bool check_existence(const std::string& name){
  struct stat buffer;
  return (stat (name.c_str(),&buffer) == 0);
}
