#include "../util/types_and_headers.h"
#include "../linalg/lapack_wrapper.h"
#include "../dtensor/dtensor_all.h"
#include "../qtensor/qtensor_all.h"
#include "../mps/mps_all.h"

#include "../models/sites/spinhalf.h"
#include "../models/hams/Heisenberg.h"

#include "../algos/dmrg/dmrg.h"

#include "../util/timer.h"

#include "randomTools.h"

#include <ctf.hpp>

using namespace std;

template<typename T>
void addIndex(dtensor<T>& A, string name){
  /*A.idx_set.emplace_back(1,name,Link);
  ++A.rank;*/
  dtensor<T> contractor({ {1,name,Link}});
  contractor.setOne();
  A = std::move(A*contractor);
}

template<typename T>
void removeIndex(dtensor<T>& A, string name){
  /*for(int l=0;l<A.rank;l++){
    if(A.idx_set[l].name()==name){
      assert(A.idx_set[l].size()==1); //only for dangling
      A.idx_set.erase(A.idx_set.begin()+l);
      break;
    }
  }
  --A.rank;*/
  dtensor<T> contractor({ {1,name,Link}});
  contractor.setOne();
  A = std::move(A*contractor);
}
template<typename T>
MPS<T> exactApplyMPO(MPO<T> & K, MPS<T> & psi,double cutoff=1E-13,int maxm=-1){
  //TODO: allow const multiply
  MPS<T> res=psi;
  assert(K.length==psi.length);
  //TODO: check that K and psi have the same sites
  int L = K.length;
  cerr<< overlap(psi,K,psi)<<endl;
  //HACK, remove dangling indices
  removeIndex(psi.A[0],"ID"+to_string(psi._id)+"Link0");
  removeIndex(psi.A[L-1],"ID"+to_string(psi._id)+"Link"+to_string(L));
  removeIndex(K.A[0],"ID"+to_string(K._id)+"Link0");
  removeIndex(K.A[L-1],"ID"+to_string(K._id)+"Link"+to_string(L));
  cerr<<"N:"<<K.A[L-1].norm()<<endl;
  //build environment tensors
  auto E = std::vector<dtensor<T> >(L);
  {
    MPS<T> psic = psi;
    MPO<T> Kc   = K;
    for(int j=0;j<L;j++){
      if(j==0){
        /*
         * auto ci = commonIndex(psi.A(1),psi.A(2),linkType);
         *psic.Aref(j) = dag(mapprime(psi.A(j),siteType,0,2,ci,0,plev));
         *ci = commonIndex(Kc.A(1),Kc.A(2),linkType);
         *Kc.Aref(j) = dag(mapprime(K.A(j),siteType,0,2,ci,0,plev));*/
        //find link index
        std::vector<dtensor_index> commonBonds;
        index_sets_intersection(psi.A[0].idx_set, psi.A[1].idx_set, commonBonds);
        assert(commonBonds.size()==1);
        psic.A[j].mapPrime(0,2,Site);
        psic.A[j].mapPrime(commonBonds,0,1717);

        index_sets_intersection(Kc.A[0].idx_set, Kc.A[1].idx_set, commonBonds);
        assert(commonBonds.size()==1);
        Kc.A[j].mapPrime(0,2,Site);
        Kc.A[j].mapPrime(commonBonds,0,1717);
      
      }
      else{
        psic.A[j].mapPrime(0,2,Site);
        psic.A[j].mapPrime(0,1717,Link);
        Kc.A[j].mapPrime(0,2,Site);
        Kc.A[j].mapPrime(0,1717,Link);
      }
    }
    E[0] = (psi.A[0]*K.A[0]*Kc.A[0]*psic.A[0]);
    for(int j=1;j<L-1;j++){
      E[j] = (E[j-1]*psi.A[j]*K.A[j]*Kc.A[j]*psic.A[j]);
    }
  }
  std::cerr<<"made Enviro"<<endl;
  //done making enviro
  auto O = std::move(psi.A[L-1]*K.A[L-1]);
  O.noPrime(Site);
  //O.print();
  
  auto Otemp = O;
  Otemp.prime(1717);
  //the ends are pesky rank 3, just make sure they contract
  /*if(Otemp.rank==4){
    auto di = dtensor_index({1},{"ID"+to_string(psi._id)+"Link"+to_string(L)},{Link},1717);
    Otemp.mapPrime(di,1717,0);
  }*/
  //Otemp.print();
  auto rho = E[L-2]* O * Otemp;
  //rho.print();
  //DIAG will destroy rho and replace with eigenvectors
  //so let's just put it where we want it
  //res.A[L-1] = rho;
  int matrixSize = 1;
  for(auto it = rho.idx_set.begin();it!=rho.idx_set.end();++it){
      if(it->level()==0) matrixSize *= it->size();
  }
  double* evals     = new double [matrixSize];
  //DIAG(matrixSize, res.A[L-1]._T.data(), evals);
  DIAG(matrixSize, rho._T.data(), evals);
  //cerr<<"$"<<rho.contract(rho)<<endl;
  //rho.print(2);
  //for(int l=0;l<matrixSize;l++){ cerr<<evals[l]<<" ";} cerr<<endl;
  //TODO determine truncation
  delete[] evals;
  //we always map the primes into links
  /*for(int l=0;l<res.A[L-1].rank;l++){
    if(res.A[L-1].idx_set[l].level()!=0){
      res.A[L-1].idx_set[l] = dtensor_index(res.A[L-1].idx_set[l].size(),"a"+to_string(L-1),Link);
      res.bond_dims[L-1] = res.A[L-1].idx_set[l].size();
      break;
    }
  }*/
  //res.A[L-1].idx_set = rho.idx_set;
  unsigned newm = 1;
  auto it = rho.idx_set.begin();
  while (it != rho.idx_set.end()){  
    if(it->level() !=0 ){
      newm *= it->size();
      it = rho.idx_set.erase(it);
    }
    else{ ++it; }
  }
  rho.idx_set.push_back(dtensor_index(newm,"a"+to_string(L-1),Link));
  //rho.rank = res.A[L-1].idx_set.size();
  res.bond_dims[L-1] = newm;
  //std::reverse(rho.idx_set.begin(),rho.idx_set.end());
  res.A[L-1] = std::move(dtensor<double>(rho.idx_set,rho._T.data()));
  res.A[L-1].print(1);
  assert(res.A[L-1].rank==2);
  O = O*res.A[L-1]*psi.A[L-2]*K.A[L-2];
  O.noPrime(Site);
  
  for(int j = L-2; j > 0; --j){
    Otemp = O; Otemp.prime(1717);
    rho = E[j-1]*O*Otemp;
    //res.A[j] = rho;
    //cerr<<j<<endl;
    //if(j>L/2) rho.print(0);
    
    matrixSize = 1;
    for(auto it = rho.idx_set.begin();it!=rho.idx_set.end();++it){
      if(it->level()==0) matrixSize *= it->size();
    }
    evals     = new double [matrixSize];
    //TODO: check that permutation is correct...
    //DIAG(matrixSize, res.A[j]._T.data(), evals);
    DIAG(matrixSize, rho._T.data(), evals);
    //for(int l=0;l<matrixSize;l++){ cerr<<evals[l]<<" ";} cerr<<endl;
    delete[] evals;
    //convert indices from primed to new link
    //res.A[j].idx_set = rho.idx_set;
    newm = 1;
    auto it = rho.idx_set.begin();
    while (it != rho.idx_set.end()){  
      if(it->level() !=0){
        newm *= it->size();
        it = rho.idx_set.erase(it);
      }
      else{ ++it; }
    }
    rho.idx_set.emplace_back(newm,"a"+to_string(j),Link);
    //std::reverse(rho.idx_set.begin(),rho.idx_set.end());
    rho.rank = rho.idx_set.size();
    //res.bond_dims[j] = newm;
    res.A[j] = std::move(dtensor<double>(rho.idx_set,rho._T.data()));
    //std::reverse(res.A[j].idx_set.begin(),res.A[j].idx_set.end());

    //if(j>L/2) res.A[j].print();

    O = O*res.A[j]*psi.A[j-1]*K.A[j-1];
    O.noPrime(Site);
  }

  O.print();
  //O /= O.norm();
  //std::reverse(O.idx_set.begin(),O.idx_set.end());
  res.A[0] = std::move(O);
  res.center = 0;
  //rename things because this is annoying
  for(int j=0;j<L;j++){
    auto& this_idx = res.A[j].idx_set;
    auto it = this_idx.begin();
    while (it!=this_idx.end()){
      if(it->name()=="a"+to_string(j)){
        it->rename("ID"+to_string(res._id)+"Link"+to_string(j));
      }
      if(it->name()=="a"+to_string(j+1)){
        it->rename("ID"+to_string(res._id)+"Link"+to_string(j+1));
      }
      ++it;
    }
  }
  //HACK
  addIndex(psi.A[0],"ID"+to_string(psi._id)+"Link0");
  addIndex(psi.A[L-1],"ID"+to_string(psi._id)+"Link"+to_string(L));
  addIndex(res.A[0],"ID"+to_string(res._id)+"Link0");
  addIndex(res.A[L-1],"ID"+to_string(res._id)+"Link"+to_string(L));
  addIndex(K.A[0],"ID"+to_string(K._id)+"Link0");
  addIndex(K.A[L-1],"ID"+to_string(K._id)+"Link"+to_string(L));
  /*for(int j=0;j<L;j++){
    uint_vec perm;
    vector<dtensor_index> new_idx_set = res.A[j].idx_set;
    std::sort(new_idx_set.begin(),new_idx_set.end());
    find_index_permutation(res.A[j].idx_set, new_idx_set, perm);
    res.A[j].permute(perm);
  }*/
  cerr<<"Res:"<<endl;
  //res.print(1);
  cerr<<"!!"<<endl;
  //res= psi;
  return res;
}

template <typename T>  
T errorMPOProd(MPS<T> & psi2, MPO<T> & K, MPS<T> & psi1){
  T err = overlap(psi2,psi2);
  //cerr<<"psi2 err "<<err<<endl;
   err += -2.*real(overlap(psi2,K,psi1));
  //cerr<<"err"<<err<<endl;
   //create K\dagger
   MPO<T> Kd = K;
   for(int j=0;j<K.length;j++){
     //swap 0,1, on sites
    Kd.A[j].mapPrime(0,3,Site);
    Kd.A[j].mapPrime(1,0,Site);
    Kd.A[j].mapPrime(3,1,Site);
   }
   err /= overlap(psi1,Kd,K,psi1);
  //cerr<<"err"<<err<<endl;
   err = std::sqrt(1.0+err);
   return err;

}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    CTF::World world(argc,argv);

  int N=10;
  int bonddim=2;
  int localdim=2;
  int mps_maxm = 2200;
  double r = 0.5; //chooses 1 vs d
  int seed = 0;
  int seedr= 0;
  int steps=2;
  if(argc == 1){
    std::cerr<<"Using all defaults"<<std::endl;
  }
  else if(argc >17 || (argc-1)%2){
    std::cout<<"Wrong input format!"<<std::endl;
    std::cout<<"-N [sites] -BD [bond dimension of random MPO] -Neel [0/1]"<<endl;
    return 1;
  }
  else{
    for(int i=1;i<argc;i+=2){
      if(std::string(argv[i]) == "-N")          N = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-BD")   bonddim = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-LD")   localdim = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-r") r = std::stof(argv[i+1]);
      if(std::string(argv[i]) == "-S")       seed = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-maxm") mps_maxm = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-steps") steps = std::stoi(argv[i+1]);
    }
  }  

  if(seed!=0){
    //readFromFile("sites_"+to_string(N)+"_"+to_string(seed),sites);
    //readFromFile("mps_"+to_string(N)+"_"+to_string(seed),psi);
  }
  else{
    seed  = make_seed();
    seedr = make_new_seed();
    //writeToFile("sites_"+to_string(N)+"_"+to_string(seed),sites);
    //writeToFile("sites_"+to_string(N)+"_"+to_string(seed)+"_"+to_string(seedr),sites);
  }
	cerr<<"N="<<N<<" bond_dim="<<bonddim<<" r="<<r<<" maxm="<<mps_maxm<<endl;
	cerr<<"Seed="<<seed<<" seedr="<<seedr<<endl;
	//std::default_random_engine generator(seed);
	std::mt19937 generator(seed);
	std::mt19937 generatorR(seedr);
	//ITensor seed
	//seedRNG(seed);
	std::normal_distribution<double> distribution(0,1.0);
	std::uniform_real_distribution<double> rdist(0.0, 1.); //for coin flips
	auto returnGauss = [&distribution,&generator]() { return distribution(generator); };


  int gaussNum = 0;
  auto testGauss = [&gaussNum](){ return gaussNum++;};

  spinhalf sites(N);
  std::cout.precision(8);

	MPS< double > psi(&sites,1);
  if(true){
    vector<string> ps;
    for (size_t i = 0; i < sites.N(); i++) {
      if(i%2==0)
        ps.push_back("Up");
      else
        ps.push_back("Dn");
    }
    psi = MPS<double>(&sites,ps);

  }
  else{
    psi.setRandom();
  }
	psi.normalize();
  cerr<<"norm:"<<overlap(psi,psi)<<endl;
  //psi.print(1);

	string mpoBDName = "mpoSize_"+to_string(N)+"_"+to_string(seed)+"_"+to_string(seedr);
	//ofstream fmpo(mpoBDName.c_str());

	string mpoLDName = "mpoLocalSize_"+to_string(N)+"_"+to_string(seed)+"_"+to_string(seedr);
	//ofstream fmpoL(mpoLDName.c_str());
	cout<<"Generating"<<endl;
	auto currentState = psi;
	int spot = N/2;
	double SvN=0.0;
  vector<double> S;
	std::vector<dtensor_index> links(N+1);
	std::vector<dtensor_index> newSites(N);
	std::vector<dtensor_index> oldSites(N);
	//convert from last time
	for(int i=0;i<N;i++) oldSites.at(i)=dtensor_index(2,"Site"+to_string(i),Site);
	
	MPO<double> O;
  O.setLength(N);
  O.setPhysicalDim(2);
  O.A.resize(N);
  O.bond_dims.resize(N+1);

  string Link_name_pref = "ID"+to_string(O._id)+"Link";
  for(int i=0;i<steps;i++){

    for(int l=0;l<N;l++){
      int thisDim = std::round(rdist(generatorR)); 
      thisDim = l==0 ? 1:2;//max(1,thisDim);
      links.at(l) = dtensor_index(thisDim,Link_name_pref+to_string(l),Link);
      O.bond_dims[l]=thisDim;
      //if (l!=0) fmpo << thisDim<<" ";

      thisDim = std::round(rdist(generatorR)); 
      thisDim = 2;//max(1,thisDim);
      newSites.at(l) = dtensor_index(thisDim,"Site"+to_string(l),Site);
      //fmpoL << thisDim<<" ";
    }
    links.at(N) = dtensor_index(1,Link_name_pref+to_string(N),Link);
    O.bond_dims[N]=1;
      //fmpo << thisDim<<" ";
   //fmpo<<endl;
   //fmpoL<<endl;
    
    //for testing
    newSites=oldSites;

    for(int l=0;l<N;l++){
      O.A[l]=std::move( dtensor<double>({links[l],oldSites[l],prime(newSites[l]),links[l+1]}) );
    }
    //O.A[0] = std::move( dtensor<double>({links[0],oldSites[0],prime(newSites[0])}) );
    //O.A[N-1] = std::move( dtensor<double>({links[N-2],oldSites[N-1],prime(newSites[N-1])}) );
    O.tensors_allocated=true;

    //O.A[0].print(1); //as seen here, by creating the tensor this way we're filling with 0s for no reason...
                    //TODO see if faster to pass array to be copied instead of filling with 0s
    for(int l=0;l<N;l++){
      //O.A[l].setRandom();
      //O.A[l].generate(testGauss);
      O.A[l].generate(returnGauss);
      cerr<<l<<" "<<O.A[l].norm()<<endl;
    }
    currentState.position(spot);
    //dtensor<double> x = std::move(currentState.A[spot]*currentState.A[spot+1]);

    //x.print();
    auto left = currentState.A[spot];
    auto right = currentState.A[spot+1];
    //left.print();
    //right.print();
    vector<dtensor_index> commonBonds;
    index_sets_intersection(right.idx_set, left.idx_set, commonBonds);
    assert(commonBonds.size()==1);
    /*for(int l=0;l<commonBonds.size();l++){
      cerr<<commonBonds[l].tag()<<"!"<<endl;
    }*/
    //cerr<<overlap(currentState,currentState)<<endl;

    //svd_bond(left,right,commonBonds[0],S,MoveFromLeft);
    //cerr<<overlap(currentState,currentState)<<endl;
    //svd(x,left,right,U,V,S,MoveFromLeft);
    SvN = 0.0;
    for(auto sg:S){
      if(sg>1e-60) SvN -= sg*sg*std::log(sg*sg);
    }
    cerr<<"SvN="<<SvN<<endl;

    MPS<double> newState = exactApplyMPO(O,currentState,1e-50,mps_maxm);
    newState.normalize();
    cerr<<"Err:"<< errorMPOProd(newState,O,currentState)<<" ov:"<<overlap(newState,newState)<<endl;
    currentState=newState;
    currentState.normalize();
    auto& T = currentState.A[0];
    cerr<<"NN:"<<T.contract(T)<<endl;
    vector<int> lengths; for(auto x: T.idx_set) lengths.push_back(x.size());
    CTF::Tensor<> TT(T.rank,lengths.data(),world);
    long int npair = T.size;
    vector<long int> idxs(npair);
    std::iota (std::begin(idxs), std::end(idxs), 0);
    TT.write(npair,idxs.data(),T._T.data());
    //T.print(1);
    //TT.print();
    cerr<<"TT:"<<TT.norm2()*TT.norm2()<<endl;
    //currentState.rc();
    cerr<<currentState.position(1)<<endl;

    //removeIndex(currentState.A[0],"ID"+to_string(currentState._id)+"Link0");
    //addIndex(currentState.A[0],"ID"+to_string(currentState._id)+"Link0");
    /*auto com = currentState.A[0]*currentState.A[1];
    com.print();
    {
    auto left = currentState.A[0];
    auto right = currentState.A[1];
    vector<dtensor_index> commonBonds;
    index_sets_intersection(right.idx_set, left.idx_set, commonBonds);
    assert(commonBonds.size()==1);
    svd_bond(left,right,commonBonds[0],S,MoveFromLeft);
    SvN = 0.0;
    for(auto sg:S){
      if(sg>1e-60) SvN -= sg*sg*std::log(sg*sg);
    }
    cerr<<"SvN="<<SvN<<endl;
    }*/
    //currentState.position(1);
    //currentState.print(1);
    cerr<<overlap(currentState,currentState)<< " "<<overlap(currentState,newState)<<endl;
    //cerr<<overlap(currentState,currentState)<<" " <<overlap(newState,newState)<<endl;
    

    dtensor<double> R({2,2,2,2},{"Site","Site","left","left"},{Site,Site,Link,Link},{0,1,0,1});
    R.generate(testGauss);
    R.print(1);
    auto dR = R.diagonal();
    dR.print(1);

    return 1;

  }

  

  //------------------------------------
  return 0;
}
