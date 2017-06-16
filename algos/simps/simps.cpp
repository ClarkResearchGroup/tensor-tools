#ifndef SHIFT_AND_INVERSE_MPS
#define SHIFT_AND_INVERSE_MPS

#include "simps.h"

template <typename T>
void simps (double tE, MPO<T>& H, MPS<T>& psi, int Nsweep){
	int L = H.length;
	char direc;
	int site;
	////////////////////////////////////////////
	T En, EnSq;
	MPO<T> HS = H;
	HS.square();
	////////////////////////////////////////////
	// Matrix side environment
	std::vector< dtensor<T> > MR(L);
	std::vector< dtensor<T> > ML(L);
	// Vector side environment
	std::vector< dtensor<T> > VR(L);
	std::vector< dtensor<T> > VL(L);
	////////////////////////////////////////////////
	psi.rc();
	buildEnv(psi, H, HS, MR, ML, VR, VL);
	////////////////////////////////////////////////
	// Repeat Nsweep
	std::cout<<"# Iteration # Mid-bond EE # Energy # Energy std-dev #"<<std::endl;
	for(int l = 0; l < Nsweep; l++) {
		std::cout<<l<<" ";
		direc = ((l%2==0)?'r':'l'); // determine the direction
		for(int i = 0; i < L-1; i++) // direction change happens at the last site of any sweep
		{
			if(direc=='r') site = i;
			if(direc=='l') site = L-1-i;
			updateSite(direc, psi, H, HS, MR, ML, VR, VL, site, En, EnSq);
			updateEnv(direc, psi, H, HS, MR, ML, VR, VL, site);
		}
		std::cout<<En+tE<<" "<<std::sqrt(std::abs(EnSq-En*En))<<std::endl;
	}
	psi.rc();
}
template void simps (double tE, MPO<double>& H, MPS<double>& psi, int Nsweep);
template void simps (double tE, MPO< std::complex<double> >& H, MPS< std::complex<double> >& psi, int Nsweep);


template <typename T>
void buildEnv(MPS<T>& psi, MPO<T>& H, MPO<T>& HSq, std::vector<dtensor<T> >& MR, std::vector<dtensor<T> >& ML, std::vector<dtensor<T> >& VR, std::vector<dtensor<T> >& VL){
	int phy = psi.index_size;
	int L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5, t6;
	for (int i = psi.length-1; i > 0; i--) {
		t1 = std::move( psi.tensorize(i) ); t3 = t1; t1.dag();
		t2 = std::move( HSq.tensorize(i) );
		t4 = std::move( H.tensorize(i) );
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
			t5 = std::move(t1*t4);
			VR[i] = std::move(t5*t3);
		}else{
			t5 = std::move(t1*VR[i+1]);
			t6 = std::move(t5*t4);
			VR[i] = std::move(t6*t3);
		}
		// std::cout<<"Site "<<i<<std::endl;
		// MR[i].print();
		// VR[i].print();
	}
}
template void buildEnv(MPS<double>& psi, MPO<double>& H, MPO<double>& HSq, std::vector<dtensor<double> >& MR, std::vector<dtensor<double> >& ML, std::vector<dtensor<double> >& VR, std::vector<dtensor<double> >& VL);
template void buildEnv(MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, MPO< std::complex<double> >& HSq, std::vector<dtensor< std::complex<double> > >& MR, std::vector<dtensor< std::complex<double> > >& ML, std::vector<dtensor< std::complex<double> > >& VR, std::vector<dtensor< std::complex<double> > >& VL);

template <typename T>
void updateSite(const char& direction, MPS<T>& psi, MPO<T>& H, MPO<T>& HSq, std::vector<dtensor<T> >& MR, std::vector<dtensor<T> >& ML, std::vector<dtensor<T> >& VR, std::vector<dtensor<T> >& VL, const int& site, T& En,  T& EnSq){
	// std::cout << "Hello, update!" << '\n';
	int phy = psi.index_size;
	int L = psi.length;
	int r = psi.bond_dims[site];
	int c = psi.bond_dims[site+1];
	int n = phy * r * c;
	int max_iter = 24;
	double tol = 1e-12;
	////////////////////////////////////////////////
	// Linear equation problem Ax=b
	dtensor<T> x, b1, b2, W1, W2, t1, t2, t3, t4, null_tensor;
	t1 = std::move( psi.tensorize(site) ); t2 = t1; x = t1; t1.dag();
	W1 = std::move( HSq.tensorize(site) );
	W2 = std::move( H.tensorize(site) );
	if(site==0){
		t3 = std::move( VR[site+1]*t2 );
		b2 = std::move( t3*W2 );

		t3 = std::move( MR[site+1]*t2 );
		b1 = std::move( t3*W1 );
	}else if(site==psi.length-1){
		t3 = std::move( VL[site-1]*t2 );
		b2 = std::move( t3*W2 );

		t3 = std::move( ML[site-1]*t2 );
		b1 = std::move( t3*W1 );
	}else{
		t3 = std::move( VR[site+1]*t2 );
		t4 = std::move( t3*W2 );
		b2 = std::move( VL[site-1]*t4 );

		t3 = std::move( MR[site+1]*t2 );
		t4 = std::move( t3*W1 );
		b1 = std::move( ML[site-1]*t4 );
	}
	b1.prime(-1);
	b2.prime(-1);
	// x.print();
	// b1.print();
	// b2.print();
	if(site==0){
		tensor_CG(null_tensor, MR[site+1], W1, x, b2, max_iter, tol);
	}else if(site==L-1){
		tensor_CG(ML[site-1], null_tensor, W1, x, b2, max_iter, tol);
	}else{
		tensor_CG(ML[site-1], MR[site+1], W1, x, b2, max_iter, tol);
	}
	x.normalize();
	EnSq = x.contract(b1);
	En   = x.contract(b2);
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
		psi.moveRight(site,true);
	}else if(direction == 'l'){
		psi.moveLeft(site,true);
	}
}
template void updateSite(const char& direction, MPS<double>& psi, MPO<double>& H, MPO<double>& HSq, std::vector<dtensor<double> >& MR, std::vector<dtensor<double> >& ML, std::vector<dtensor<double> >& VR, std::vector<dtensor<double> >& VL, const int& site, double& En, double& EnSq);
template void updateSite(const char& direction, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, MPO< std::complex<double> >& HSq, std::vector<dtensor< std::complex<double> > >& MR, std::vector<dtensor< std::complex<double> > >& ML, std::vector<dtensor< std::complex<double> > >& VR, std::vector<dtensor< std::complex<double> > >& VL, const int& site, std::complex<double>& En, std::complex<double>& EnSq);

template <typename T>
void updateEnv(const char& direction, MPS<T>& psi, MPO<T>& H, MPO<T>& HSq, std::vector<dtensor<T> >& MR, std::vector<dtensor<T> >& ML, std::vector<dtensor<T> >& VR, std::vector<dtensor<T> >& VL, const int& site){
	int phy = psi.index_size;
	int L = psi.length;
	dtensor<T> t1, t2, t3, t4, t5, t6;
	t1 = std::move( psi.tensorize(site) ); t3 = t1; t1.dag();
	t2 = std::move( HSq.tensorize(site) );
	t4 = std::move( H.tensorize(site) );
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
			t5 = std::move(t1*t4);
			VL[site] = std::move(t5*t3);
		}else{
			t5 = std::move(VL[site-1]*t1);
			t6 = std::move(t5*t4);
			VL[site] = std::move(t6*t3);
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
			t5 = std::move(t1*t4);
			VR[site] = std::move(t5*t3);
		}else{
			t5 = std::move(t1*VR[site+1]);
			t6 = std::move(t5*t4);
			VR[site] = std::move(t6*t3);
		}
	}
}
template void updateEnv(const char& direction, MPS<double>& psi, MPO<double>& H, MPO<double>& HSq, std::vector<dtensor<double> >& MR, std::vector<dtensor<double> >& ML, std::vector<dtensor<double> >& VR, std::vector<dtensor<double> >& VL, const int& site);
template void updateEnv(const char& direction, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, MPO< std::complex<double> >& HSq, std::vector<dtensor< std::complex<double> > >& MR, std::vector<dtensor< std::complex<double> > >& ML, std::vector<dtensor< std::complex<double> > >& VR, std::vector<dtensor< std::complex<double> > >& VL, const int& site);

#endif
