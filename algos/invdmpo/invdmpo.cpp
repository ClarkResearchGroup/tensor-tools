#ifndef DIAGONAL_MPO_INVERSION
#define DIAGONAL_MPO_INVERSION

#include "invdmpo.h"

template <typename T>
void invdmpo (MPO<T>& H, MPO<T>& invH, int Nsweep, double cutoff, double tol){
	int L = H.length;
	char direc;
	int site;
	////////////////////////////////////////////
	T delta = 1;
	MPO<T> HS = H;
	HS.square();
	HS.print();
	std::cout << "l2norm((0.1Hd^2)^2) = " << l2norm(HS) << '\n';
	std::cout << "trace((0.1Hd^2)^2)  = " << trace(HS) << '\n';
	MPO<T> HS_approx;
	HS_approx.fit(HS,10,cutoff);
	HS_approx.print();
	std::cout << "l2norm((0.1Hd^2)^2) = " << l2norm(HS_approx) << '\n';
	std::cout << "trace((0.1Hd^2)^2)  = " << trace(HS_approx) << '\n';
	MPS<T> phi = DiagonalMPOAsMPS(H);
	MPS<T> psi = DiagonalMPOAsMPS(invH);
	MPO<T> IdM(L,H.phy_size,1); IdM.setIdentity();
	////////////////////////////////////////////
	// Matrix side environment
	std::vector< dtensor<T> > MR(L);
	std::vector< dtensor<T> > ML(L);
	// Vector side environment
	std::vector< dtensor<T> > VR(L);
	std::vector< dtensor<T> > VL(L);
	////////////////////////////////////////////////
	buildEnv(psi, phi, HS, MR, ML, VR, VL);
	////////////////////////////////////////////////
	// Repeat Nsweep
	std::cout<<"# Sweep # |invD*D-I| #"<<std::endl;
	for(int l = 0; l < Nsweep; l++) {
		std::cout<<l<<" ";
		direc = ((l%2==0)?'r':'l'); // determine the direction
		for(int i = 0; i < L-1; i++) // direction change happens at the last site of any sweep
		{
			if(direc=='r') site = i;
			if(direc=='l') site = L-1-i;
			updateSite(direc, psi, phi, HS, MR, ML, VR, VL, site, tol);
			updateEnv(direc, psi, phi, HS, MR, ML, VR, VL, site);
		}
		invH = MPSAsDiagonalMPO(psi);
		MPO<T> res;
		fitApplyMPO(invH, H, res, 'I', true, cutoff);
		res = res - IdM;
		std::cout<<res.norm()<<std::endl;
	}
	invH = MPSAsDiagonalMPO(psi);
}
template void invdmpo (MPO<double>& H, MPO<double>& invH, int Nsweep, double cutoff, double tol);
template void invdmpo (MPO< std::complex<double> >& H, MPO< std::complex<double> >& invH, int Nsweep, double cutoff, double tol);


template <typename T>
void buildEnv(MPS<T>& psi, MPS<T>& phi, MPO<T>& H, std::vector<dtensor<T> >& MR, std::vector<dtensor<T> >& ML, std::vector<dtensor<T> >& VR, std::vector<dtensor<T> >& VL){
	int phy = psi.index_size;
	int L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5, t6;
	for (int i = psi.length-1; i > 0; i--) {
		t1 = std::move( psi.tensorize(i) ); t1.dag();
		t2 = std::move( H.tensorize(i) );
		t3 = std::move( psi.tensorize(i) );
		t4 = std::move( phi.tensorize(i) ); t4.primeSite();
		// Matrix side
		if(i==psi.length-1){
			t5 = std::move(t1*t2);
			MR[i] = std::move(t5*t3);
		}else{
			t5 = std::move(t1*MR[i+1]);
			t6 = std::move(t5*t2);
			MR[i] = std::move(t6*t3);
		}
		// Vector side
		if(i==psi.length-1){
			VR[i] = std::move(t1*t4);
		}else{
			t5 = std::move(t1*VR[i+1]);
			VR[i] = std::move(t5*t4);
		}
		// std::cout<<"Site "<<i<<std::endl;
		// MR[i].print();
		// VR[i].print();
	}
}
template void buildEnv(MPS<double>& psi, MPS<double>& phi, MPO<double>& H, std::vector<dtensor<double> >& MR, std::vector<dtensor<double> >& ML, std::vector<dtensor<double> >& VR, std::vector<dtensor<double> >& VL);
template void buildEnv(MPS< std::complex<double> >& psi, MPS< std::complex<double> >& phi, MPO< std::complex<double> >& H, std::vector<dtensor< std::complex<double> > >& MR, std::vector<dtensor< std::complex<double> > >& ML, std::vector<dtensor< std::complex<double> > >& VR, std::vector<dtensor< std::complex<double> > >& VL);

template <typename T>
void updateSite(const char& direction, MPS<T>& psi, MPS<T>& phi, MPO<T>& H, std::vector<dtensor<T> >& MR, std::vector<dtensor<T> >& ML, std::vector<dtensor<T> >& VR, std::vector<dtensor<T> >& VL, const int& site, double tol){
	// std::cout << "Hello, update!" << '\n';
	int max_iter = 400;
	int L = psi.length;
	////////////////////////////////////////////////
	// Linear equation problem Ax=b
	dtensor<T> x, b, W, t1, t2, t3, null_tensor;
	x  = std::move( psi.tensorize(site) );
	t1 = std::move( phi.tensorize(site) ); t1.primeSite();
	W  = std::move( H.tensorize(site) );
	if(site==0){
		b  = std::move( VR[site+1]*t1 );
	}else if(site==psi.length-1){
		b  = std::move( VL[site-1]*t1 );
	}else{
		t2 = std::move( VR[site+1]*t1 );
		b  = std::move( VL[site-1]*t2 );
	}
	b.prime(-1);
	if(site==0){
		tensor_CG(null_tensor, MR[site+1], W, x, b, max_iter, tol, false);
	}else if(site==L-1){
		tensor_CG(ML[site-1], null_tensor, W, x, b, max_iter, tol, false);
	}else{
		tensor_CG(ML[site-1], MR[site+1], W, x, b, max_iter, tol, false);
	}
	// x.normalize();
	////////////////////////////////////////////////
	// copy data from dtensor to MPS
	int idx = 0;
	for (int j = 0; j < psi.index_size; j++) {
		for (int k = 0; k < psi.M[site][j].size(); k++) {
			psi.M[site][j].data()[k] = x._T.data()[idx];
			++idx;
		}
	}
	// move to neighbor
	if(direction == 'r'){
		psi.moveRight(site,false);
	}else if(direction == 'l'){
		psi.moveLeft(site,false);
	}
}
template void updateSite(const char& direction, MPS<double>& psi, MPS<double>& phi, MPO<double>& H, std::vector<dtensor<double> >& MR, std::vector<dtensor<double> >& ML, std::vector<dtensor<double> >& VR, std::vector<dtensor<double> >& VL, const int& site, double tol);
template void updateSite(const char& direction, MPS< std::complex<double> >& psi, MPS< std::complex<double> >& phi, MPO< std::complex<double> >& H, std::vector<dtensor< std::complex<double> > >& MR, std::vector<dtensor< std::complex<double> > >& ML, std::vector<dtensor< std::complex<double> > >& VR, std::vector<dtensor< std::complex<double> > >& VL, const int& site, double tol);

template <typename T>
void updateEnv(const char& direction, MPS<T>& psi, MPS<T>& phi, MPO<T>& H, std::vector<dtensor<T> >& MR, std::vector<dtensor<T> >& ML, std::vector<dtensor<T> >& VR, std::vector<dtensor<T> >& VL, const int& site){
	int phy = psi.index_size;
	int L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5, t6;
	t1 = std::move( psi.tensorize(site) ); t1.dag();
	t2 = std::move( H.tensorize(site) );
	t3 = std::move( psi.tensorize(site) );
	t4 = std::move( phi.tensorize(site) ); t4.primeSite();
	////////////////////////////////////////////////
	if(direction == 'r'){
		// Matrix side
		if(site==0){
			t5 = std::move(t1*t2);
			ML[site] = std::move(t5*t3);
		}else{
			t5 = std::move(ML[site-1]*t1);
			t6 = std::move(t5*t2);
			ML[site] = std::move(t6*t3);
		}
		// Vector side
		if(site==0){
			VL[site] = std::move(t1*t4);
		}else{
			t5 = std::move(VL[site-1]*t1);
			VL[site] = std::move(t5*t4);
		}
	}else if(direction == 'l'){
		// Matrix side
		if(site==L-1){
			t5 = std::move(t1*t2);
			MR[site] = std::move(t5*t3);
		}else{
			t5 = std::move(t1*MR[site+1]);
			t6 = std::move(t5*t2);
			MR[site] = std::move(t6*t3);
		}
		// Vector side
		if(site==L-1){
			VR[site] = std::move(t1*t4);
		}else{
			t5 = std::move(t1*VR[site+1]);
			VR[site] = std::move(t5*t4);
		}
	}
}
template void updateEnv(const char& direction, MPS<double>& psi, MPS<double>& phi, MPO<double>& H, std::vector<dtensor<double> >& MR, std::vector<dtensor<double> >& ML, std::vector<dtensor<double> >& VR, std::vector<dtensor<double> >& VL, const int& site);
template void updateEnv(const char& direction, MPS< std::complex<double> >& psi, MPS< std::complex<double> >& phi, MPO< std::complex<double> >& H, std::vector<dtensor< std::complex<double> > >& MR, std::vector<dtensor< std::complex<double> > >& ML, std::vector<dtensor< std::complex<double> > >& VR, std::vector<dtensor< std::complex<double> > >& VL, const int& site);

#endif
