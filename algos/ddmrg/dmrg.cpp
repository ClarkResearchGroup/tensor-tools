#ifndef DMRG_CLASS
#define DMRG_CLASS

#include "dmrg.h"

template <typename T>
void buildEnv(MPS<T>& psi, MPO<T>& H, std::vector< dtensor<T> >& TR, std::vector< dtensor<T> >& TL){
	int phy = psi.index_size;
	int L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5;
	for (int i = psi.length-1; i > 0; i--) {
		t1 = std::move( psi.tensorize(i) ); t1.dag();
		t2 = std::move( H.tensorize(i) );
		t3 = std::move( psi.tensorize(i) );
		if(i==psi.length-1){
			t4 = std::move(t1*t2);
			TR[i] = std::move(t4*t3);
		}else{
			t4 = std::move(TR[i+1]*t1);
			t5 = std::move(t4*t2);
			TR[i] = std::move(t5*t3);
		}
	}
}
template void buildEnv(MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL);
template void buildEnv(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL);

template <typename T>
void updateSite(const char& direction, MPS<T>& psi, MPO<T>& H, std::vector< dtensor<T> >& TR, std::vector< dtensor<T> >& TL, const int& site, T& EN, char whichEn){
	int phy = psi.index_size;
	int r = psi.bond_dims[site];
	int c = psi.bond_dims[site+1];
	int n = phy * r * c;
	int nev = 1;
	////////////////////////////////////////////////
	// solve sparse eigen problem
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M_in(r,c*phy);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> M_out(r,c*phy);
	for(int j = 0; j < phy; j++) {
		M_in.block(0,j*c,r,c) = psi.M[site][j];
	}
	T **Evecs, *Evals;
	Evals = new T[nev]();
	Evecs = new T*[nev];
	Evecs[0] = M_out.data();
	if(whichEn=='H')
		arpack_L(n, nev, Evals, Evecs, M_in.data(), psi, H, TR, TL, site);
	else
		arpack_S(n, nev, Evals, Evecs, M_in.data(), psi, H, TR, TL, site);
	EN = Evals[0];
	for(int tid=0; tid<phy; tid++){
		psi.M[site][tid] = M_out.block(0,tid*c,r,c);
	}
	delete [] Evals;
	delete [] Evecs;
	////////////////////////////////////////////////
	if(direction == 'r'){
		psi.moveRight(site,true);
	}else if(direction == 'l'){
		psi.moveLeft(site,true);
	}
}
template void updateSite(const char& direction, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, const int& site, double& EN, char whichEn);
template void updateSite(const char& direction, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, const int& site,  std::complex<double> & EN, char whichEn);

template <typename T>
void updateEnv(const char& direction, MPS<T>& psi, MPO<T>& H, std::vector< dtensor<T> >& TR, std::vector< dtensor<T> >& TL, const int& site){
	int phy = psi.index_size;
	int L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5;
	t1 = std::move( psi.tensorize(site) ); t1.dag();
	t2 = std::move( H.tensorize(site) );
	t3 = std::move( psi.tensorize(site) );
	////////////////////////////////////////////////
	if(direction == 'r'){
		if(site==0){
			t4 = std::move(t1*t2);
			TL[site] = std::move(t4*t3);
		}else{
			t4 = std::move(TL[site-1]*t1);
			t5 = std::move(t4*t2);
			TL[site] = std::move(t5*t3);
		}
	}else if(direction == 'l'){
		if(site==L-1){
			t4 = std::move(t1*t2);
			TR[site] = std::move(t4*t3);
		}else{
			t4 = std::move(TR[site+1]*t1);
			t5 = std::move(t4*t2);
			TR[site] = std::move(t5*t3);
		}
	}
}
template void updateEnv(const char& direction, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, const int& site);
template void updateEnv(const char& direction, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, const int& site);

template <typename T>
T dmrg(MPS<T>& psi, MPO<T>& H, int Nsweep, char whichEn){
	int L = H.length;
	////////////////////////////////////////////
	psi.normalize();
	MPO<T> HS = H;
	HS.square();
	////////////////////////////////////////////
	char direc;
	int site;
	T Energy = psiHphi(psi, H, psi);
	std::cout<<"Initial energy of the MPS: "<<Energy<<std::endl;
	////////////////////////////////////////////
	// Environment tensors
	std::vector< dtensor<T> > TR(L);
	std::vector< dtensor<T> > TL(L);
	////////////////////////////////////////////////
	psi.rc();
	buildEnv(psi, H, TR, TL);
	////////////////////////////////////////////////
	// Repeat Nsweep
	std::cout<<"# DMRG Sweep # Mid bond EE # Energy # Std #"<<std::endl;
	for(int l = 0; l < Nsweep; l++) {
		std::cout<<l<<" ";
		direc = ((l%2==0)?'r':'l'); // determine the direction
		for(int i = 0; i < L-1; i++) // direction change happens at the last site of any sweep
		{
			if(direc=='r') site = i;
			if(direc=='l') site = L-1-i;
			updateSite(direc, psi, H, TR, TL, site, Energy, whichEn);
			updateEnv(direc, psi, H, TR, TL, site);
		}
		T SE = psiHphi(psi,HS,psi);
		T E  = psiHphi(psi,H,psi);
		std::cout<<E<<" "<<std::sqrt(std::abs(SE-E*E))<<std::endl;
	}
	psi.rc();
	return Energy;
}
template double dmrg(MPS<double>& psi, MPO<double>& H, int Nsweep, char whichEn);
template std::complex<double> dmrg(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, int Nsweep, char whichEn);

#endif
