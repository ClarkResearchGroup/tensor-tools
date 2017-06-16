#ifndef My_HAMS
#define My_HAMS

#include "hams.h"

template <typename T>
void buildHeisenberg(MPO<T>& H, double tE) {
  assert(H.tensors_allocated);
	double J1=1, J2=1;
  double enps = tE/double(H.length);
	H.setZero(5);

	H.M[0][0]<<-enps,0,0,(J2)/2,1;
	H.M[0][1]<<0,0,(J1)/2,0,0;
	H.M[0][2]<<0,(J1)/2,0,0,0;
	H.M[0][3]<<-enps,0,0,-(J2)/2,1;
	for(int i = 1; i < H.length-1; i++) {
		H.M[i][0]<<1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0.5,0,0,0,0,  -enps,0,0,(J2)/2,1;
		H.M[i][1]<<0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,    0,0,(J1)/2,0,0;
		H.M[i][2]<<0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,    0,(J1)/2,0,0,0;
		H.M[i][3]<<1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  -0.5,0,0,0,0, -enps,0,0,-(J2)/2,1;
	}
	H.M[H.length-1][0]<<1,0,0,0.5,-enps;
	H.M[H.length-1][1]<<0,1,0,0,0;
	H.M[H.length-1][2]<<0,0,1,0,0;
	H.M[H.length-1][3]<<1,0,0,-0.5,-enps;
}
template void buildHeisenberg(MPO<double>& H, double tE);
template void buildHeisenberg(MPO< std::complex<double> >& H, double tE);

template <typename T>
void buildHeisenberg(MPO<T>& H, double tE, double* dJ, double* dh) {
  assert(H.tensors_allocated);
	double enps = tE/double(H.length);
	double J1=1, J2=1;
	H.setZero(5);

	H.M[0][0]<<-enps+0.5*dh[0],0,0,(J2+dJ[0])/2,1;
	H.M[0][1]<<0,0,(J1+dJ[0])/2,0,0;
	H.M[0][2]<<0,(J1+dJ[0])/2,0,0,0;
	H.M[0][3]<<-enps-0.5*dh[0],0,0,-(J2+dJ[0])/2,1;
	for(int i = 1; i < H.length-1; i++) {
		H.M[i][0]<<1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  0.5,0,0,0,0,  -enps+0.5*dh[i],0,0,(J2+dJ[i])/2,1;
		H.M[i][1]<<0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,    0,0,(J1+dJ[i])/2,0,0;
		H.M[i][2]<<0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,    0,(J1+dJ[i])/2,0,0,0;
		H.M[i][3]<<1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,  -0.5,0,0,0,0, -enps-0.5*dh[i],0,0,-(J2+dJ[i])/2,1;
	}
	H.M[H.length-1][0]<<1,0,0,0.5,-enps+0.5*dh[H.length-1];
	H.M[H.length-1][1]<<0,1,0,0,0;
	H.M[H.length-1][2]<<0,0,1,0,0;
	H.M[H.length-1][3]<<1,0,0,-0.5,-enps-0.5*dh[H.length-1];
}
template void buildHeisenberg(MPO<double>& H, double En, double* dJ, double* dh);
template void buildHeisenberg(MPO< std::complex<double> >& H, double En, double* dJ, double* dh);

template <typename T>
void buildSz(MPO<T>& H, int site) {
  assert(H.tensors_allocated);
	H.setZero(1);
	for(int i = 0; i < H.length; ++i){
		if(i!=site)
		{
			H.M[i][0].setIdentity();
			H.M[i][3].setIdentity();
		}else{
			H.M[i][0].setIdentity();
			H.M[i][0] *= 1;
			H.M[i][3].setIdentity();
			H.M[i][3] *= -1;
		}
	}
}
template void buildSz(MPO<double>& H, int site);
template void buildSz(MPO< std::complex<double> >& H, int site);

#endif
