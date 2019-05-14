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

//#include <ctf.hpp>
#include <thread>

using namespace std;

//returns  error,newStart
std::tuple<double,int> determineCutoff(double* evals,int size,int maxm,
                                          double cutoff, bool absCutoff=false, bool rescale=true){
  if(maxm<0) maxm = size;
  //Note: evals is sorted smallest->largest

  //if all zeros
  if(size==1)            return std::make_tuple(0.,0.);
  //for convenience make negatives zero
  for(int i=0;i<size;i++){
    if(evals[i]>=0.) break;
    evals[i] = 0.;
  }
  if(evals[size-1]==0.0) return std::make_tuple(0.,0.);
  double error = 0.0;
  int n=0;
  //determine error from just maxm truncation
  while ( size-n > maxm){ error += evals[n++]; }
  
  if(absCutoff)
    while(evals[n] < cutoff) { error += evals[n++]; }
  else{
    double scale = 1.0;
    if(rescale){ //normalize error
      scale = 0.0;
      for(int i=0;i<size;i++) scale+=evals[i];
    }
   while(error+evals[n] < cutoff*scale && n<size){
    error += evals[n++];
   }
  error /= scale; 
  }
  return std::make_tuple(error,n);

} 
template<typename T>
void addIndex(dtensor<T>& A, string name){
  /*A.idx_set.emplace_back(1,name,Link);
  ++A.rank;*/
  dtensor<T> contractor({ {1,name,Link}});
  contractor.setOne();
  A = std::move(contractor*A);
  /*auto tmp_set = A.idx_set;
  tmp_set.emplace_back(1,name,Link);
  A = std::move(dtensor<T>(*/

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
MPS<T> exactApplyMPO(MPO<T> & K, MPS<T> & psi,double cutoff=1E-13,int maxm=-1, bool verbose=false){
  //TODO: allow const multiply
  MPS<T> res=psi;
  assert(K.length==psi.length);
  //TODO: check that K and psi have the same sites
  int L = K.length;
  //HACK, remove dangling indices
  removeIndex(psi.A[0],"ID"+to_string(psi._id)+"Link0");
  removeIndex(psi.A[L-1],"ID"+to_string(psi._id)+"Link"+to_string(L));
  removeIndex(K.A[0],"ID"+to_string(K._id)+"Link0");
  removeIndex(K.A[L-1],"ID"+to_string(K._id)+"Link"+to_string(L));
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
    E[0] = std::move(psi.A[0]*K.A[0]*Kc.A[0]*psic.A[0]);
    for(int j=1;j<L-1;j++){
      E[j] = std::move(E[j-1]*psi.A[j]*K.A[j]*Kc.A[j]*psic.A[j]);
    }
  }
  if(verbose) std::cerr<<"made Enviro"<<endl;
  //done making enviro
  auto O = std::move(psi.A[L-1]*K.A[L-1]);
  O.noPrime(Site);
  
  //auto Otemp = O;
  //Otemp.prime(1717);
  //auto rho = std::move(E[L-2]* O * Otemp);
  auto rho = std::move(E[L-2]*O);
       O.prime(1717);
       rho = std::move(rho*O);
       O.prime(-1717);
  //cerr<<rho.norm()<<" "<<rho.contract(rho)<<endl;
  //DIAG will destroy rho and replace with eigenvectors
  //but we need to init TBLIS so we keep rho around
  int matrixSize = 1;
  for(auto it = rho.idx_set.begin();it!=rho.idx_set.end();++it){
      if(it->level()==0) matrixSize *= it->size();
  }
  double* evals     = new double [matrixSize];
  double error; int newStart;
  DIAGD(matrixSize, rho._T.data(), evals);
  //for(int l=0;l<matrixSize;l++){ cerr<<evals[l]<<" ";} cerr<<endl;
  std::tie(error,newStart) = determineCutoff(evals,matrixSize,maxm,cutoff);
  //cerr<<"Err: "<<error<< " "<<newStart<<endl;
  delete[] evals;
  unsigned newm = 1;
  auto it = rho.idx_set.begin();
  while (it != rho.idx_set.end()){  
    if(it->level() !=0 ){
      newm *= it->size();
      it = rho.idx_set.erase(it);
    }
    else{ ++it; }
  }
  newm -= newStart;
  rho.idx_set.push_back(dtensor_index(newm,"a"+to_string(L-1),Link));
  res.bond_dims[L-1] = newm;
  res.A[L-1] = std::move(dtensor<double>(rho.idx_set,rho._T.data()+(matrixSize*newStart)));
  //res.A[L-1].print(1);
  assert(res.A[L-1].rank==2);
  O = std::move(O*res.A[L-1]*psi.A[L-2]*K.A[L-2]);
  /*O = std::move(O*res.A[L-1]);
  O = std::move(O*psi.A[L-2]);
  O = std::move(O*K.A[L-2]);*/
  O.noPrime(Site);
  
  for(int j = L-2; j > 0; --j){
    //Otemp = O; Otemp.prime(1717);
    //rho = std::move(E[j-1]*O*Otemp);
    rho = std::move(E[j-1]*O);
    O.prime(1717);
    rho = std::move(rho*O);
    O.prime(-1717);
    //cerr<<j<< " "<<rho.norm()<<" "<<rho.contract(rho)<<endl;
    //cerr<<j<<endl;
    //if(j>L/2) rho.print(0);
    
    matrixSize = 1;
    for(auto it = rho.idx_set.begin();it!=rho.idx_set.end();++it){
      if(it->level()==0) matrixSize *= it->size();
    }
    evals     = new double [matrixSize];
    DIAGD(matrixSize, rho._T.data(), evals);
    std::tie(error,newStart) = determineCutoff(evals,matrixSize,maxm,cutoff);
    if(verbose) cerr<<j<<" Err: "<<error<< " "<<newStart<<endl;
    delete[] evals;
    //convert indices from primed to new link
    newm = 1;
    auto it = rho.idx_set.begin();
    while (it != rho.idx_set.end()){  
      if(it->level() !=0){
        newm *= it->size();
        it = rho.idx_set.erase(it);
      }
      else{ ++it; }
    }
    newm -=newStart;
    rho.idx_set.emplace_back(newm,"a"+to_string(j),Link);
    rho.rank = rho.idx_set.size();
    res.bond_dims[j] = newm;
    res.A[j] = std::move(dtensor<double>(rho.idx_set,
                                         rho._T.data()+(matrixSize*newStart) ));
    //cerr<<j<< " "<<res.A[j].norm()<<" "<<res.A[j].contract(res.A[j])<<endl;
    //if(j>L/2) res.A[j].print();

    O = std::move(O*res.A[j]*psi.A[j-1]*K.A[j-1]);
    /*O = std::move(O*res.A[j]);
    O = std::move(O*psi.A[j-1]);
    O = std::move(O*K.A[j-1]);*/
    O.noPrime(Site);
  }

  //O /= O.norm();
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
  if(verbose){
    cerr<<"Res:"<<endl;
    res.print(1);
    cerr<<"!!"<<endl;
  }
  return res;
}

template <typename T>  
T errorMPOProd(MPS<T> & psi2, MPO<T> & K, MPS<T> & psi1){
  T err = overlap(psi2,psi2);
  err += -2.*real(overlap(psi2,K,psi1));
  //create K\dagger
  MPO<T> Kd = K;
  for(int j=0;j<K.length;j++){
    //swap 0,1, on sites
    Kd.A[j].mapPrime(0,3,Site);
    Kd.A[j].mapPrime(1,0,Site);
    Kd.A[j].mapPrime(3,1,Site);
  }
  err /= overlap(psi1,Kd,K,psi1);
  err = std::sqrt(abs(1.0+err));//abs needed for underflow
  return err;

}

int main(int argc, char *argv[]) {
    //MPI_Init(&argc, &argv);
    //CTF::World world(argc,argv);

  int N=10;
  int bonddim=2;
  int localdim=2;
  int mps_maxm = 2200;
  double r = 0.1; //chooses 1 vs d
  double mu = 0.5;
  int seed = 0;
  int seedr= 0;
  int steps=10;
  int trial = 0;
  if(argc == 1){
    std::cerr<<"Using all defaults"<<std::endl;
  }
  else if(argc >19 || (argc-1)%2){
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
      if(std::string(argv[i]) == "-mu") mu = std::stof(argv[i+1]);
      if(std::string(argv[i]) == "-S")       seed = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-maxm") mps_maxm = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-steps") steps = std::stoi(argv[i+1]);
      if(std::string(argv[i]) == "-trial") trial = std::stoi(argv[i+1]);
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
  cerr<<"N="<<N<<" bond_dim="<<bonddim<<" r="<<r<<" mu="<<mu<<" maxm="<<mps_maxm<<endl;
  cerr<<"Seed="<<seed<<" seedr="<<seedr<<endl;
	//std::default_random_engine generator(seed);
	std::mt19937 generator(seed); _gen().seed(seed);
	std::mt19937 generatorR(seedr);
  //index IDs
  std::hash<std::thread::id> hasher;
  static thread_local std::mt19937 indexID(std::clock() + hasher(std::this_thread::get_id()));
	//ITensor seed
	//seedRNG(seed);
	std::normal_distribution<double> distribution(0,1.0);
	//std::uniform_real_distribution<double> rdist(0.0, 1.); //for coin flips
  std::lognormal_distribution<> rdist(mu,r); //hardcode lognormal
	auto returnGauss = [&distribution,&generator]() { return distribution(generator); };


  int gaussNum = 0;
  auto testGauss = [&gaussNum](){ return gaussNum++;};

  spinhalf sites(N);
  std::cout.precision(8);

	MPS< double > psi(&sites,1);
  if(false){
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

	string mpoBDName = "mpoSize_"+to_string(N)+"_"+to_string(seed)+"_"+to_string(seedr);
	ofstream fmpo(mpoBDName.c_str());

	string mpoLDName = "mpoLocalSize_"+to_string(N)+"_"+to_string(seed)+"_"+to_string(seedr);
	ofstream fmpoL(mpoLDName.c_str());
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
      thisDim = l==0 ? 1:max(1,thisDim);
      links.at(l) = dtensor_index(thisDim,Link_name_pref+to_string(l),Link);
      O.bond_dims[l]=thisDim;
      if (l!=0) fmpo << thisDim<<" ";

      thisDim = std::round(rdist(generatorR)); 
      thisDim = max(1,thisDim);
      newSites.at(l) = dtensor_index(thisDim,to_string(indexID())+"_Site"+to_string(l),Site);
      fmpoL << thisDim<<" ";
    }
    links.at(N) = dtensor_index(1,Link_name_pref+to_string(N),Link);
    O.bond_dims[N]=1;
    fmpo<<endl;
    fmpoL<<endl;

    for(int l=0;l<N;l++){
      O.A[l]=std::move( dtensor<double>({links[l],oldSites[l],prime(newSites[l]),links[l+1]}) );
    }
    O.tensors_allocated=true;

    //O.A[0].print(1); //as seen here, by creating the tensor this way we're filling with 0s for no reason...
                    //TODO see if faster to pass array to be copied instead of filling with 0s
    for(int l=0;l<N;l++){
      //O.A[l].generate(testGauss);
      O.A[l].generate(returnGauss);
    }
    //Heisenberg< double > HB(&sites);
    //HB.buildHam(O);
    //std::cerr<<"1Svd="<<currentState.position(spot)<<std::endl;
    currentState.position(spot);

    int maxm = *std::max_element(currentState.bond_dims.begin(),currentState.bond_dims.end());
    auto& A_left = currentState.A[spot-1];
    auto& A_right = currentState.A[spot];
    vector<dtensor_index> commonbonds;
    index_sets_intersection(A_right.idx_set, A_left.idx_set, commonbonds);
    assert(commonbonds.size()==1);
    //svd_bond(left,right,commonbonds[0],S,MoveFromRight);


    dtensor<double> U,V;
    dtensor<double> combined = std::move(A_left * A_right);
    vector<dtensor_index> left;
    vector<dtensor_index> right;
    string tag = commonbonds[0].tag();
    // Separate dtensor_index
    for (size_t j = 0; j < A_right.rank; j++) {
      string idx_tag = A_right.idx_set[j].tag();
      if (idx_tag != tag) {
        right.push_back(A_right.idx_set[j]);
      }
    }
    for (size_t j = 0; j < A_left.rank; j++) {
      string idx_tag = A_left.idx_set[j].tag();
      if (idx_tag != tag) {
        left.push_back(A_left.idx_set[j]);
      }
    }
    // SVD
    svd(combined,left,right,U,V,S,MoveFromRight);
    SvN = 0.0;
    for(auto sg:S){
      if(sg>1e-60) SvN -= sg*sg*std::log(sg*sg);
    }
    cerr<<maxm<<" "<<SvN<<endl;

    MPS<double> newState = std::move(exactApplyMPO(O,currentState,1e-50,mps_maxm));
    //cerr<<"made state"<<endl;
    for(int l=0;l<N;l++){
      assert(newState.A[l].rank == newState.A[l]._T.dimension());
      for(int k=0;k<newState.A[l].rank;k++)
        assert(newState.A[l].idx_set[k].size()==newState.A[l]._T.length(k));
    }
    //cerr<<"Err:"<< errorMPOProd(newState,O,currentState)<<" ov:"<<overlap(newState,newState)<<" E:"<<overlap(newState,O,newState)<<endl;
    newState.normalize();
    currentState=newState;
    //cerr<<overlap(currentState,currentState)<< " "<<overlap(currentState,newState)<<endl;
    oldSites.swap(newSites);
  } //end steps
  fmpo.close();
  fmpoL.close();

  //------------------------------------
  return 0;
}
