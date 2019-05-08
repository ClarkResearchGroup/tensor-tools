#ifndef DMRG_CLASS
#define DMRG_CLASS

#include "dmrg.h"
#include "../../util/timer.h"

template <typename T>
void buildEnv(MPS<T>& psi, MPO<T>& H, std::vector< dtensor<T> >& TR, std::vector< dtensor<T> >& TL){
	// Left ends
	{
		dtensor<T> left_ends({psi.A[0].idx_set[0], H.A[0].idx_set[0], psi.A[0].idx_set[0]});
		left_ends.idx_set[0].primeLink();
		left_ends.setOne();
		left_ends.dag();
		TL[0] = left_ends;
	}
	// right ends
	{
		dtensor<T> right_ends({psi.A[psi.length-1].idx_set.back(), H.A[psi.length-1].idx_set.back(), psi.A[psi.length-1].idx_set.back()});
		right_ends.idx_set[0].primeLink();
		right_ends.setOne();
		right_ends.dag();
		TR[psi.length-1] = right_ends;
	}
	dtensor<T> t1, t2, t3, t4, t5;
	for (size_t i = psi.length-1; i > 0; i--) {
		t1 = psi.A[i]; t1.conj(); t1.dag(); t1.prime();
		t2 = H.A[i];
		t3 = psi.A[i];
		t4 = std::move(t1*TR[i]);
		t5 = std::move(t4*t2);
		TR[i-1] = std::move(t5*t3);
	}
}
template void buildEnv(MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL);
template void buildEnv(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL);


template <typename T>
void buildEnv(qMPS<T>& psi, qMPO<T>& H, std::vector< qtensor<T> >& TR, std::vector< qtensor<T> >& TL){
	// Left ends
	{
		qtensor<T> left_ends({psi.A[0].idx_set[0], H.A[0].idx_set[0], psi.A[0].idx_set[0]});
		left_ends.idx_set[0].primeLink(); left_ends.idx_set[0].dag();
		left_ends.initBlock();
		left_ends.setOne();
		left_ends.dag();
		TL[0] = left_ends;
	}
	// right ends
	{
		qtensor<T> right_ends({psi.A[psi.length-1].idx_set.back(), H.A[psi.length-1].idx_set.back(), psi.A[psi.length-1].idx_set.back()});
		right_ends.idx_set[0].primeLink(); right_ends.idx_set[0].dag();
		right_ends.initBlock();
		right_ends.setOne();
		right_ends.dag();
		TR[psi.length-1] = right_ends;
	}
	qtensor<T> t1, t2, t3, t4, t5;
	for (size_t i = psi.length-1; i > 0; i--) {
		t1 = psi.A[i]; t1.conj(); t1.dag(); t1.prime();
		t2 = H.A[i];
		t3 = psi.A[i];
		t4 = std::move(t1*TR[i]);
		t5 = std::move(t4*t2);
		TR[i-1] = std::move(t5*t3);
	}
}
template void buildEnv(qMPS<double>& psi, qMPO<double>& H, std::vector< qtensor<double> >& TR, std::vector< qtensor<double> >& TL);
template void buildEnv(qMPS< std::complex<double> >& psi, qMPO< std::complex<double> >& H, std::vector< qtensor< std::complex<double> > >& TR, std::vector< qtensor< std::complex<double> > >& TL);


template <typename T>
void updateSite(MPS<T>& psi, MPO<T>& H, std::vector< dtensor<T> >& TR, std::vector< dtensor<T> >& TL, const unsigned& site, T& energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart, Timer &davidsonTimer){
	if(direction==MoveFromLeft){
		// Set up big_dtensor for two site optimization
		big_dtensor<T> A;
		A.setLeft (&TL[site]);
		A.addMid(&H.A[site]);
		A.addMid(&H.A[site+1]);
		A.setRight(&TR[site+1]);
		// Davidson Eigen solver (two-site)
		dtensor<T> x = std::move(psi.A[site]*psi.A[site+1]);
		davidsonTimer.Start();
		energy = tensor_davidson(A, x, search_space_size, max_restart, 1e-12, mode);
		davidsonTimer.Stop();
		// SVD and move canonical center
		//vector<double> S;
		dtensor<T> S;
		dtensor_index mid;
		string LinkName = "ID"+to_string(psi._id)+"Link"+to_string(site+1);
		for (size_t i = 0; i < psi.A[site].rank; i++) {
			string idx_name = psi.A[site].idx_set[i].name();
			if (idx_name == LinkName) {
				mid = psi.A[site].idx_set[i];
				break;
			}
		}
		svd_bond(x,psi.A[site],psi.A[site+1],mid,S,MoveFromLeft,cutoff,max_bd);
		psi.bond_dims[site+1] = S.size;
		psi.center = site+1;
		// EE
		if(site==psi.length/2 - 1){
      double vNEE = calcEntropy(S);
			perr << vNEE << " ";
		}
	}else{
		// Set up big_dtensor for two site optimization
		big_dtensor<T> A;
		A.setLeft (&TL[site-1]);
		A.addMid(&H.A[site-1]);
		A.addMid(&H.A[site]);
		A.setRight(&TR[site]);
		// Davidson Eigen solver (two-site)
		dtensor<T> x = std::move(psi.A[site-1]*psi.A[site]);
		davidsonTimer.Start();
		energy = tensor_davidson(A, x, search_space_size, max_restart, 1e-12, mode);
		davidsonTimer.Stop();
		// SVD and move canonical center
		//vector<double> S;
    dtensor<T> S;
		dtensor_index mid;
		string LinkName = "ID"+to_string(psi._id)+"Link"+to_string(site);
		for (size_t i = 0; i < psi.A[site-1].rank; i++) {
			string idx_name = psi.A[site-1].idx_set[i].name();
			if (idx_name == LinkName) {
				mid = psi.A[site-1].idx_set[i];
				break;
			}
		}
		svd_bond(x,psi.A[site-1],psi.A[site],mid,S,MoveFromRight,cutoff,max_bd);
		psi.bond_dims[site] = S.size;
		psi.center = site - 1;
		// EE
		if(site==psi.length/2){
      double vNEE = calcEntropy(S);
      perr << vNEE << " ";
		}
	}
}
template void updateSite(MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, const unsigned& site, double& energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart,Timer &t);
template void updateSite(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, const unsigned& site,  std::complex<double> & energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart, Timer &t);


template <typename T>
void updateSite(qMPS<T>& psi, qMPO<T>& H, std::vector< qtensor<T> >& TR, std::vector< qtensor<T> >& TL, const unsigned& site, T& energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart, Timer &diag){
	if(direction==MoveFromLeft){
		// Set up big_dtensor for two site optimization
		big_qtensor<T> A;
		A.setLeft (&TL[site]);
		A.addMid(&H.A[site]);
		A.addMid(&H.A[site+1]);
		A.setRight(&TR[site+1]);
		// Davidson Eigen solver (two-site)
		qtensor<T> x = std::move(psi.A[site]*psi.A[site+1]);
		diag.Start();
		energy = tensor_davidson(A, x, search_space_size, max_restart, 1e-12, mode);
		diag.Stop();
		// SVD and move canonical center
		vector<double> S;
		qtensor_index mid;
		string LinkName = "ID"+to_string(psi._id)+"Link"+to_string(site+1);
		for (size_t i = 0; i < psi.A[site].rank; i++) {
			string idx_name = psi.A[site].idx_set[i].name();
			if (idx_name == LinkName) {
				mid = psi.A[site].idx_set[i];
				break;
			}
		}
		svd_bond(x,psi.A[site],psi.A[site+1],mid,S,MoveFromLeft,cutoff,max_bd);
		psi.bond_dims[site+1] = S.size();
		psi.center = site+1;
		// EE
		if(site==psi.length/2 - 1){
			double vNEE = 0.0;
			std::cout << vNEE << " ";
		}
	}else{
		// Set up big_dtensor for two site optimization
		big_qtensor<T> A;
		A.setLeft (&TL[site-1]);
		A.addMid(&H.A[site-1]);
		A.addMid(&H.A[site]);
		A.setRight(&TR[site]);
		// Davidson Eigen solver (two-site)
		qtensor<T> x = std::move(psi.A[site-1]*psi.A[site]);
		diag.Start();
		energy = tensor_davidson(A, x, search_space_size, max_restart, 1e-12, mode);
		diag.Stop();
		// SVD and move canonical center
		vector<double> S;
		qtensor_index mid;
		string LinkName = "ID"+to_string(psi._id)+"Link"+to_string(site);
		for (size_t i = 0; i < psi.A[site-1].rank; i++) {
			string idx_name = psi.A[site-1].idx_set[i].name();
			if (idx_name == LinkName) {
				mid = psi.A[site-1].idx_set[i];
				break;
			}
		}
		svd_bond(x,psi.A[site-1],psi.A[site],mid,S,MoveFromRight,cutoff,max_bd);
		psi.bond_dims[site] = S.size();
		psi.center = site - 1;
		// EE
		if(site==psi.length/2){
			double vNEE = 0.0;
			std::cout << vNEE << " ";
		}
	}
}
template void updateSite(qMPS<double>& psi, qMPO<double>& H, std::vector< qtensor<double> >& TR, std::vector< qtensor<double> >& TL, const unsigned& site, double& energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart,Timer &t);
template void updateSite(qMPS< std::complex<double> >& psi, qMPO< std::complex<double> >& H, std::vector< qtensor< std::complex<double> > >& TR, std::vector< qtensor< std::complex<double> > >& TL, const unsigned& site,  std::complex<double> & energy, int& direction, int max_bd, double cutoff, char mode, int search_space_size, int max_restart,Timer &t);


template <typename T>
void updateEnv(MPS<T>& psi, MPO<T>& H, std::vector< dtensor<T> >& TR, std::vector< dtensor<T> >& TL, const unsigned& site, int& direction){
	int phy = psi.phy_dim;
	unsigned L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5;
	t1 = psi.A[site]; t1.conj(); t1.dag(); t1.prime();
	t2 = H.A[site];
	t3 = psi.A[site];
	////////////////////////////////////////////////
	if(direction == MoveFromLeft){

		t4 = std::move(t1*TL[site]);
		t5 = std::move(t4*t2);
		TL[site+1] = std::move(t5*t3);

	}else if(direction == MoveFromRight){
		t4 = std::move(t1*TR[site]);
		t5 = std::move(t4*t2);
		TR[site-1] = std::move(t5*t3);
	}
}
template void updateEnv(MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, const unsigned& site, int& direction);
template void updateEnv(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, const unsigned& site, int& direction);


template <typename T>
void updateEnv(qMPS<T>& psi, qMPO<T>& H, std::vector< qtensor<T> >& TR, std::vector< qtensor<T> >& TL, const unsigned& site, int& direction){
	int phy = psi.phy_dim;
	unsigned L = psi.length;
	qtensor<T> t1, t2, t3, t4, t5;
	t1 = psi.A[site]; t1.conj(); t1.dag(); t1.prime();
	t2 = H.A[site];
	t3 = psi.A[site];
	////////////////////////////////////////////////
	if(direction == MoveFromLeft){
		t4 = std::move(t1*TL[site]);
		t5 = std::move(t4*t2);
		TL[site+1] = std::move(t5*t3);
	}else if(direction == MoveFromRight){
		t4 = std::move(t1*TR[site]);
		t5 = std::move(t4*t2);
		TR[site-1] = std::move(t5*t3);
	}
}
template void updateEnv(qMPS<double>& psi, qMPO<double>& H, std::vector< qtensor<double> >& TR, std::vector< qtensor<double> >& TL, const unsigned& site, int& direction);
template void updateEnv(qMPS< std::complex<double> >& psi, qMPO< std::complex<double> >& H, std::vector< qtensor< std::complex<double> > >& TR, std::vector< qtensor< std::complex<double> > >& TL, const unsigned& site, int& direction);


template <typename T>
T dmrg(MPS<T>& psi, MPO<T>& H, int num_sweeps, int max_bd, double cutoff, char mode, int search_space_size, int max_restart,int start_sweep){
  Timer davidsonTimer("Davidson");
	int L = H.length;
	int direction, site=0;
	if(start_sweep%2==0) psi.position(0); else psi.position(L-1);


	psi.normalize();

	//exit(1);
	T Energy = psiHphi(psi, H, psi);
	////////////////////////////////////////////
	// Environment tensors
	std::vector< dtensor<T> > TR(L);
	std::vector< dtensor<T> > TL(L);
	buildEnv(psi, H, TR, TL);
	////////////////////////////////////////////////
	// Repeat Nsweep
	if(start_sweep==0) perr<<"# Sweep # Mid bond EE # Energy # Time(s) # Davidson(s)"<<std::endl;
  Timer t("time");
	for(int l = start_sweep; l < num_sweeps; l++) {
    t.Start();
		perr<<l/2<<((l%2==0)? "_L" : "_R")<<"\t ";
		direction = ((l%2==0)? MoveFromLeft : MoveFromRight); // determine the direction
		for(int i = 0; i < L-2; i++) // direction change happens at the last site of any sweep
		{
			if(direction==MoveFromLeft)  site = i;
			if(direction==MoveFromRight) site = L-1-i;
			updateSite(psi, H, TR, TL, site, Energy, direction, max_bd, cutoff, mode, search_space_size, max_restart,davidsonTimer);
			updateEnv(psi, H, TR, TL, site, direction);

		}
    t.Stop();
		perr<<Energy<<"\t"<<t.Time()<<"\t"<<davidsonTimer.Time()<<std::endl;
		davidsonTimer.Clear();
    t.Clear();
	}
	//psi.position(0);
	psi.normalize();
	return Energy;
}
template double dmrg(MPS<double>& psi, MPO<double>& H, int num_sweeps, int max_bd, double cutoff, char mode, int search_space_size, int max_restart, int start_sweep);
template std::complex<double> dmrg(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, int num_sweeps, int max_bd, double cutoff, char mode, int search_space_size, int max_restart,int start_sweep);


template <typename T>
T dmrg(MPS<T>& psi, MPO<T>& H, int num_sweeps, const std::vector<int>& max_bd, const std::vector<double>& cutoff, const std::vector<int>& max_restart){
  //assert(max_bd.size() == cutoff.size());
  //assert(max_restart.size() == max_bd.size());

  T Energy;
  for(int sweep =0;sweep<num_sweeps;sweep++){
    //allow sweeps to loop through
    int lbd = sweep%(max_bd.size());
    int lc  = sweep%(cutoff.size());
    int lmr = sweep%(max_restart.size());
    //always do two sweeps so that we go left to right
    Energy = dmrg(psi,H,2*(sweep+1),
                  max_bd[lbd],cutoff[lc],'S',3,max_restart[lmr],2*sweep);
  }
  if(num_sweeps%2==0) psi.position(0); else psi.position(psi.length-1);
 return Energy; 
}
template double dmrg(MPS<double>& psi, MPO<double>& H, int num_sweeps, const std::vector<int>& max_bd, const std::vector<double>& cutoff, const std::vector<int>& max_restart);

template std::complex<double> dmrg(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, int num_sweeps, 
                                   const std::vector<int>& max_bd, const std::vector<double>& cutoff_vec,const std::vector<int>& max_restart);

template <typename T>
T dmrg(qMPS<T>& psi, qMPO<T>& H, int num_sweeps, int max_bd, double cutoff, char mode, int search_space_size, int max_restart,int start_sweep){
  Timer davidsonTimer("Davidon");

	int L = H.length;
	int direction, site=0;
	psi.position(0);
	psi.normalize();
	T Energy = psiHphi(psi, H, psi);
	/*std::cout<<"Initial energy of the MPS: "<<Energy<<std::endl;*/
	////////////////////////////////////////////
	// Environment tensors
	std::vector< qtensor<T> > TR(L);
	std::vector< qtensor<T> > TL(L);
	buildEnv(psi, H, TR, TL);
	////////////////////////////////////////////////
	// Repeat Nsweep
	if (start_sweep==0) std::cout<<"# Sweep # Mid bond EE # Energy #"<<std::endl;
	for(int l = start_sweep; l < num_sweeps; l++) {
		std::cout<<l<<"\t ";
		direction = ((l%2==0)? MoveFromLeft : MoveFromRight); // determine the direction
		for(int i = 0; i < L-2; i++) // direction change happens at the last site of any sweep
		{
			if(direction==MoveFromLeft)  site = i;
			if(direction==MoveFromRight) site = L-1-i;
			updateSite(psi, H, TR, TL, site, Energy, direction, max_bd, cutoff, mode, search_space_size, max_restart,davidsonTimer);
			updateEnv(psi, H, TR, TL, site, direction);

		}
		std::cout<<Energy<<std::endl;
		perr<<"Davidson Time: "<<davidsonTimer.Time()<<endl;
		davidsonTimer.Clear();
	}
	psi.position(0);
	psi.normalize();
	return Energy;
}
template double dmrg(qMPS<double>& psi, qMPO<double>& H, int num_sweeps, int max_bd, double cutoff, char mode, int search_space_size, int max_restart, int start_sweep);
template std::complex<double> dmrg(qMPS< std::complex<double> >& psi, qMPO< std::complex<double> >& H, int num_sweeps, int max_bd, double cutoff, char mode, int search_space_size, int max_restart, int start_sweep);

template <typename T>
T dmrg(qMPS<T>& psi, qMPO<T>& H, int num_sweeps, const std::vector<int>& max_bd, const std::vector<double>& cutoff, const std::vector<int>& max_restart){
  assert(max_bd.size() == cutoff.size());
  assert(max_restart.size() == max_bd.size());

  T Energy;
  for(int sweep =0;sweep<num_sweeps;sweep++){
    //allow sweeps to loop through
    int l = sweep%(max_bd.size());
    //always do two sweeps so that we go left to right
    Energy = dmrg(psi,H,2*(sweep+1),max_bd[l],cutoff[l],'S',3,max_restart[l],2*sweep);
  }
 return Energy; 
}
template double dmrg(qMPS<double>& psi, qMPO<double>& H, int num_sweeps, const std::vector<int>& max_bd, const std::vector<double>& cutoff, const std::vector<int>& max_restart);

template std::complex<double> dmrg(qMPS< std::complex<double> >& psi, qMPO< std::complex<double> >& H, int num_sweeps, 
                                   const std::vector<int>& max_bd, const std::vector<double>& cutoff_vec,const std::vector<int>& max_restart);
#endif
