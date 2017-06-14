#ifndef My_LAPACK_WRAPPERS_H
#define My_LAPACK_WRAPPERS_H

#include "../util/types_and_headers.h"

#ifdef USE_MKL
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#endif

#if !defined(USE_MKL)
// #if !defined(USE_MKL) && !defined(USE_TBLIS)
// double
// extern "C" void dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
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
// std::complex<double>
// extern "C" void zgemm_(char*,char*,int*,int*,int*,std::complex<double>*,std::complex<double>*,int*,std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
extern "C" void zgeqrf_(int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);
extern "C" void zungqr_(int*, int*, int*, std::complex<double>*, int*, std::complex<double>*, std::complex<double>*, int*, int*);
extern "C" void zgesvd_(char*, char*, int*, int*, std::complex<double>*, int*, double*, std::complex<double>*, int*, std::complex<double>*, int*, std::complex<double>*, int*, double*, int*);
extern "C" void zheev_(char*,char*,int*,std::complex<double>*,int*,double*,std::complex<double>*,int*,double*,int*);
extern "C" void zgemm_(char*,char*,int*,int*,int*,std::complex<double>*,std::complex<double>*,int*,std::complex<double>*,int*,std::complex<double>*,std::complex<double>*,int*);
extern "C" void zhesv_(char*,int*,int*,std::complex<double>*,int*,int*,std::complex<double>*,int*,std::complex<double>*,int*,int*);
extern "C" void zgesv_(int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);
extern "C" void zhegv_(int*,char*,char*,int*,std::complex<double>*,int*,std::complex<double>*,int*,double*,std::complex<double>*,int*,double*,int*);
#endif

void diag(double* M, double* evals, int nn);

void QR(const Mxd& MT, Mxd& Q);
void QR(const Mxc& MT, Mxc& Q);

void SVD(const Mxd& M, int ds, double * sv, Mxd& UM, Mxd& VTM, char direction);
void SVD(const Mxd& M, int ds, double * sv, Mxd& UM, Mxd& VTM, int truncD, double& svd_error);
void SVD(const Mxd& M, std::vector<double>& sv, Mxd& UM, Mxd& VTM, double& svd_cutoff, bool dry_run=true);
void SVD(const Mxc& M, int ds, double * sv, Mxc& UM, Mxc& VTM, char direction);
void SVD(const Mxc& M, int ds, double * sv, Mxc& UM, Mxc& VTM, int truncD, double& svd_error);
void SVD(const Mxc& M, std::vector<double>& sv, Mxc& UM, Mxc& VTM, double& svd_cutoff, bool dry_run=true);

void linearSolver(int ord, Mxd& A, Mxd& vec);

void indef_linearSolver(int ord, Mxd& A, Mxd& vec);

int def_linearSolver(int ord, Mxd& A, Mxd& vec);

void geindef_linearSolver(int ord, Mxd& A, Mxd& vec);

int sdef_linearSolver(int ord, Mxd& A, Mxd& vec);

void sindef_linearSolver(int ord, Mxd& A, Mxd& vec);

void sxdef_linearSolver(int ord, Mxd& A, Mxd& vec);

void rCSD(Mxd& A, Mxd& theta, Mxd& U1, Mxd& U2, Mxd& V1, Mxd& V2);

#endif
