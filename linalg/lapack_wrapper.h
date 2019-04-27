#ifndef My_LAPACK_WRAPPERS_H
#define My_LAPACK_WRAPPERS_H

#include "../util/types_and_headers.h"

// extern "C" void dgemm_(char*,cmakehar*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
extern "C" void dgeqrf_(int*, int*, double*, int*, double*, double*, int*, int*);
extern "C" void dorgqr_(int*, int*, int*, double*, int*, double*, double*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
extern "C" void dsytrs_(char*,int*,int*,double*,int*,int*,double*,int*,int*); // broken?
extern "C" void dsysv_(char*,int*,int*,double*,int*,int*,double*,int*,double*,int*,int*); // the only safe one to use
extern "C" void dposv_(char*,int*,int*,double*,int*,double*,int*,int*);
extern "C" void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*); // the only safe one to use
extern "C" void dsposv_(char*,int*,int*,double*,int*,double*,int*,double*,int*,double*,float*,int*,int*);
extern "C" void dsysvx_(char*,char*,int*,int*,double*,int*,double*,int*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*,int*);
extern "C" void dposvx_(char*,char*,int*,int*,double*,int*,double*,int*,char*,double*,double*,int*,double*,int*,double*,double*,double*,double*,int*,int*);
extern "C" void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
extern "C" void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
extern "C" void zgemm_(char*,char*,int*,int*,int*,std::complex<double>*,std::complex<double>*,int*,std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
extern "C" void zgeqrf_(int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);
extern "C" void zungqr_(int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);
extern "C" void zgesvd_(char*, char*, int*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, int*, double*, int*);
extern "C" void zheev_(char*,char*,int*,std::complex<double>*,int*,double*,std::complex<double>*,int*,double*,int*);
extern "C" void zhesv_(char*,int*,int*,std::complex<double>*,int*,int*,std::complex<double>*,int*,std::complex<double>*,int*,int*);
extern "C" void zgesv_(int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);
extern "C" void zhegv_(int*,char*,char*,int*,std::complex<double>*,int*,std::complex<double>*,int*,double*,std::complex<double>*,int*,double*,int*);




void MAT_VEC(int& m, int& k, int& n, double* A, double* B, double* C);
void MAT_VEC(int& m, int& k, int& n, std::complex<double>* A, std::complex<double>* B, std::complex<double>* C);

void DIAG(int m, double* A, double* evals);
void DIAG(int m, std::complex<double>* A, double* evals);

void ORTHO(int r, int c, double* A);
void ORTHO(int r, int c, std::complex<double>* A);

void QR(int r, int c, double* A, double* Q, double* R);
void QR(int r, int c, std::complex<double>* A, std::complex<double>* Q, std::complex<double>* R);

void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V);
void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V);
void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V, char direction);
void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V, char direction);
void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V, char direction, double cutoff);
void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V, char direction, double cutoff);
void SVD(int r, int c, double* A, double* U, vector<double>& S, double* V, char direction, int max_size);
void SVD(int r, int c, std::complex<double>* A, std::complex<double>* U, vector<double>& S, std::complex<double>* V, char direction, int max_size);

#endif
