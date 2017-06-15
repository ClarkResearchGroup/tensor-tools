#ifndef My_ARPACK_WRAPPERS_H
#define My_ARPACK_WRAPPERS_H


//-----------------------------------------------------------------------------
// Arpack -- REAL SYMMETRIC eigen pairs
extern "C" void dsaupd_(int *ido, char *bmat, int *n, char *which,
			int *nev, double *tol, double *resid, int *ncv,
			double *v1, int *ldv, int *iparam, int *ipntr,
			double *workd, double *workl, int *lworkl,
			int *info);

extern "C" void dseupd_(int *rvec, char *All, int *select, double *d,
			double *v2, int *ldv, double *sigma,
			char *bmat, int *n, char *which, int *nev,
			double *tol, double *resid, int *ncv, double *v3,
			int *ldv1, int *iparam, int *ipntr, double *workd,
      double *workl, int *lworkl, int *ierr);

//-----------------------------------------------------------------------------
// Arpack -- Complex SYMMETRIC eigen pairs
extern "C" void znaupd_(int *ido, char *bmat, int *n, char *which,
			int *nev, double *tol, std::complex<double> *resid, int *ncv,
			std::complex<double> *v1, int *ldv, int *iparam, int *ipntr,
			std::complex<double> *workd, std::complex<double> *workl, int *lworkl,
			double *rwork, int *info);

extern "C" void zneupd_(int *rvec, char *All, int *select, std::complex<double> *d,
			std::complex<double> *v2, int *ldv, std::complex<double> *sigma, std::complex<double> *workev,
			char *bmat, int *n, char *which, int *nev,
			double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *v3,
			int *ldv1, int *iparam, int *ipntr, std::complex<double> *workd,
			std::complex<double> *workl, int *lworkl, double *rwrok, int *ierr);

//-----------------------------------------------------------------------------
// Arpack -- NON-SYMMETRIC eigen pairs
extern "C" void dnaupd_(int *ido, char *bmat, int *n, char *which,
			int *nev, double *tol, double *resid, int *ncv,
			double *v1, int *ldv, int *iparam, int *ipntr,
			double *workd, double *workl, int *lworkl,
			int *info);
extern "C" void dneupd_(int* rvec, char* All, int *select, double *dr,
      double* di, double* z, int* ldz, double* sigmar, double* sigmai, double* workev,
      char* bmat, int* n, char* which, int* nev,
			double* tol, double* resid, int* ncv, double* v,
			int* ldv, int* iparam, int* ipntr, double* workd,
			double* workl, int* lworkl, int* ierr);
//-----------------------------------------------------------------------------

#endif
