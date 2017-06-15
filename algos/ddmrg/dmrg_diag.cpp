#ifndef My_DMRG_DIAG_CLASS
#define My_DMRG_DIAG_CLASS

#include "dmrg_diag.h"

void arpack_mv(int n, double *in, double *out, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, int site){
	dtensor<double> t1, t2, t3, t4, t5;
	t1 = std::move( psi.tensorize(site) ); // copy the indices from psi
	std::copy(in,in+n,t1._T.data()); // copy the in array into the tensor
	t2 = std::move( H.tensorize(site) );
	////////////////////////////////////////////////
	// Order of the multiplication is important below!!!
	// dtensor operator * preserves the order of the not-contracted indices
	// in order to get the out array to be in the right order without doing
	// additional transposition, the order of the multiplication needs to be
	// chosen carefully.
	if(site==0){
		t3 = std::move( TR[site+1]*t1 );
		t5 = std::move( t3*t2 );
	}else if(site==psi.length-1){
		t3 = std::move( TL[site-1]*t1 );
		t5 = std::move( t3*t2 );
	}else{
		t3 = std::move( TR[site+1]*t1 );
		t4 = std::move( t3*t2 );
		t5 = std::move( TL[site-1]*t4 );
	}
	////////////////////////////////////////////////
	std::copy(t5._T.data(),t5._T.data()+n,out);
}

void arpack_mv(int n, std::complex<double> *in, std::complex<double> *out, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, int site){
	dtensor< std::complex<double> > t1, t2, t3, t4, t5;
	t1 = std::move( psi.tensorize(site) ); // copy the indices from psi
	std::copy(in,in+n,t1._T.data()); // copy the in array into the tensor
	t2 = std::move( H.tensorize(site) );
	////////////////////////////////////////////////
	// Order of the multiplication is important below!!!
	// dtensor operator * preserves the order of the not-contracted indices
	// in order to get the out array to be in the right order without doing
	// additional transposition, the order of the multiplication needs to be
	// chosen carefully.
	if(site==0){
		t3 = std::move( TR[site+1]*t1 );
		t5 = std::move( t3*t2 );
	}else if(site==psi.length-1){
		t3 = std::move( TL[site-1]*t1 );
		t5 = std::move( t3*t2 );
	}else{
		t3 = std::move( TR[site+1]*t1 );
		t4 = std::move( t3*t2 );
		t5 = std::move( TL[site-1]*t4 );
	}
	////////////////////////////////////////////////
	// t5.dag();
	std::copy(t5._T.data(),t5._T.data()+n,out);
}

void arpack_S(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, int site){
	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "SA";
	double tol = 1e-8;
	double *resid;
	resid = initVec;
	// resid = new double[n];
	int ncv = 10*nev;
	if (ncv>n) ncv = n;
	double *v;
	int ldv = n;
	v = new double[ldv*ncv];
	int *iparam;
	iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = 100;
	iparam[6] = 1;
	int *ipntr;
	ipntr = new int[11];
	double *workd;
	workd = new double[3*n];
	double *workl;
	workl = new double[ncv*(ncv+8)+16];
	int lworkl = ncv*(ncv+8)+16;
	int info = 1;
	int rvec = 1;  // Changed from above
	int *select;
	select = new int[ncv];
	double *d;
	d = new double[2*ncv];
	double sigma;
	int ierr;

	do {
	dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);

	if ((ido==1)||(ido==-1)) arpack_mv(n, workd+ipntr[0]-1, workd+ipntr[1]-1, psi, H, TR, TL, site);
	} while ((ido==1)||(ido==-1));

	if (info<0) {
	   //  cout << "Error with dsaupd, info = " << info << "\n";
	//	 cout << "Check documentation in dsaupd\n\n";
	} else {
	dseupd_(&rvec, (char *) "A", select, d, v, &ldv, &sigma, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);

	if (ierr!=0) {
	//  cout << "Error with dseupd, info = " << ierr << "\n";
	 // cout << "Check the documentation of dseupd.\n\n";
	} else if (info==1) {
	//  cout << "Maximum number of iterations reached.\n\n";
	} else if (info==3) {
	//  cout << "No shifts could be applied during implicit\n";
	//  cout << "Arnoldi update, try increasing NCV.\n\n";
	}
	int i, j;
	for (i=0; i<nev; i++) Evals[i] = d[i];
	for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];

	// delete resid;
	delete [] v;
	delete [] iparam;
	delete [] ipntr;
	delete [] workd;
	delete [] workl;
	delete [] select;
	delete [] d;
	}
}

void arpack_L(int n, int nev, double *Evals, double **Evecs, double *initVec, MPS<double>& psi, MPO<double>& H, std::vector< dtensor<double> >& TR, std::vector< dtensor<double> >& TL, int site){
	int ido = 0;
	char bmat[2] = "I";
	char which[3] = "LA";
	double tol = 1e-8;
	double *resid;
	resid = initVec;
	// resid = new double[n];
	int ncv = 10*nev;
	if (ncv>n) ncv = n;
	double *v;
	int ldv = n;
	v = new double[ldv*ncv];
	int *iparam;
	iparam = new int[11];
	iparam[0] = 1;
	iparam[2] = 100;
	iparam[6] = 1;
	int *ipntr;
	ipntr = new int[11];
	double *workd;
	workd = new double[3*n];
	double *workl;
	workl = new double[ncv*(ncv+8)+16];
	int lworkl = ncv*(ncv+8)+16;
	int info = 1;
	int rvec = 1;  // Changed from above
	int *select;
	select = new int[ncv];
	double *d;
	d = new double[2*ncv];
	double sigma;
	int ierr;

	do {
	dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info);

	if ((ido==1)||(ido==-1)) arpack_mv(n, workd+ipntr[0]-1, workd+ipntr[1]-1, psi, H, TR, TL, site);
	} while ((ido==1)||(ido==-1));

	if (info<0) {
	   //  cout << "Error with dsaupd, info = " << info << "\n";
	//	 cout << "Check documentation in dsaupd\n\n";
	} else {
	dseupd_(&rvec, (char *) "A", select, d, v, &ldv, &sigma, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, &ierr);

	if (ierr!=0) {
	//  cout << "Error with dseupd, info = " << ierr << "\n";
	 // cout << "Check the documentation of dseupd.\n\n";
	} else if (info==1) {
	//  cout << "Maximum number of iterations reached.\n\n";
	} else if (info==3) {
	//  cout << "No shifts could be applied during implicit\n";
	//  cout << "Arnoldi update, try increasing NCV.\n\n";
	}
	int i, j;
	for (i=0; i<nev; i++) Evals[i] = d[i];
	for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];

	// delete resid;
	delete [] v;
	delete [] iparam;
	delete [] ipntr;
	delete [] workd;
	delete [] workl;
	delete [] select;
	delete [] d;
	}
}

void arpack_S(int n, int nev,  std::complex<double>  *Evals,  std::complex<double>  **Evecs,  std::complex<double>  *initVec, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, int site){
	int ido = 0;
  char bmat[2] = "I";
  char which[3] = "SR";
  double tol = 1e-8;
  std::complex<double> *resid;
	resid = initVec;
  // resid = new std::complex<double> [n];
  int ncv = 12*nev;
  if (ncv>n) ncv = n;
  std::complex<double> *v;
  int ldv = n;
  v = new std::complex<double> [ldv*ncv];
  int *iparam;
  iparam = new int[11];
  iparam[0] = 1;
  iparam[2] = 400;
  iparam[6] = 1;
  int *ipntr;
  ipntr = new int[14];
  std::complex<double> *workd;
  workd = new std::complex<double> [3*n];
  std::complex<double> *workl;
  workl = new std::complex<double> [3*ncv*(ncv+2)];
  int lworkl = 3*ncv*(ncv+2);
  int info = 0;
  int rvec = 1;  // Changed from above
  int *select;
  select = new int[ncv];
  std::complex<double> *d;
  d = new std::complex<double> [2*ncv];
  double *rwork = new double [ncv];
  std::complex<double> *workev = new std::complex<double> [2*ncv];
  std::complex<double> sigma;
  int ierr;

  do {
    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, rwork, &info);
	if ((ido==1)||(ido==-1)) arpack_mv(n, workd+ipntr[0]-1, workd+ipntr[1]-1, psi, H, TR, TL, site);
  } while ((ido==1)||(ido==-1));

  if (info<0) {
       //  std::cout << "Error with znaupd, info = " << info << "\n";
	//	 std::cout << "Check documentation in znaupd\n\n";
  } else {
    zneupd_(&rvec, (char *) "A", select, d, v, &ldv, &sigma, workev, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

    if (ierr!=0) {
    //  std::cout << "Error with zneupd, info = " << ierr << "\n";
     // std::cout << "Check the documentation of zneupd.\n\n";
    } else if (info==1) {
    //  std::cout << "Maximum number of iterations reached.\n\n";
    } else if (info==3) {
    //  std::cout << "No shifts could be applied during implicit\n";
    //  std::cout << "Arnoldi update, try increasing NCV.\n\n";
    }
		int i, j;
		for (i=0; i<nev; i++) Evals[i] = d[i];
		for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];

    // delete [] resid;
    delete [] v;
    delete [] iparam;
    delete [] ipntr;
    delete [] workd;
    delete [] workl;
    delete [] select;
    delete [] d;
		delete [] rwork;
		delete [] workev;
  }
}

void arpack_L(int n, int nev,  std::complex<double>  *Evals,  std::complex<double>  **Evecs,  std::complex<double>  *initVec, MPS< std::complex<double> >& psi, MPO< std::complex<double> >& H, std::vector< dtensor< std::complex<double> > >& TR, std::vector< dtensor< std::complex<double> > >& TL, int site){
	int ido = 0;
  char bmat[2] = "I";
  char which[3] = "LR";
  double tol = 1e-8;
  std::complex<double> *resid;
	resid = initVec;
  // resid = new std::complex<double> [n];
  int ncv = 12*nev;
  if (ncv>n) ncv = n;
  std::complex<double> *v;
  int ldv = n;
  v = new std::complex<double> [ldv*ncv];
  int *iparam;
  iparam = new int[11];
  iparam[0] = 1;
  iparam[2] = 400;
  iparam[6] = 1;
  int *ipntr;
  ipntr = new int[14];
  std::complex<double> *workd;
  workd = new std::complex<double> [3*n];
  std::complex<double> *workl;
  workl = new std::complex<double> [3*ncv*(ncv+2)];
  int lworkl = 3*ncv*(ncv+2);
  int info = 0;
  int rvec = 1;  // Changed from above
  int *select;
  select = new int[ncv];
  std::complex<double> *d;
  d = new std::complex<double> [2*ncv];
  double *rwork = new double [ncv];
  std::complex<double> *workev = new std::complex<double> [2*ncv];
  std::complex<double> sigma;
  int ierr;

  do {
    znaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, rwork, &info);
	if ((ido==1)||(ido==-1)) arpack_mv(n, workd+ipntr[0]-1, workd+ipntr[1]-1, psi, H, TR, TL, site);
  } while ((ido==1)||(ido==-1));

  if (info<0) {
       //  std::cout << "Error with znaupd, info = " << info << "\n";
	//	 std::cout << "Check documentation in znaupd\n\n";
  } else {
    zneupd_(&rvec, (char *) "A", select, d, v, &ldv, &sigma, workev, bmat,
	    &n, which, &nev, &tol, resid, &ncv, v, &ldv,
	    iparam, ipntr, workd, workl, &lworkl, rwork, &ierr);

    if (ierr!=0) {
    //  std::cout << "Error with zneupd, info = " << ierr << "\n";
     // std::cout << "Check the documentation of zneupd.\n\n";
    } else if (info==1) {
    //  std::cout << "Maximum number of iterations reached.\n\n";
    } else if (info==3) {
    //  std::cout << "No shifts could be applied during implicit\n";
    //  std::cout << "Arnoldi update, try increasing NCV.\n\n";
    }
    int i, j;
		for (i=0; i<nev; i++) Evals[i] = d[i];
		for (i=0; i<nev; i++) for (j=0; j<n; j++) Evecs[i][j] = v[i*n+j];

    // delete [] resid;
    delete [] v;
    delete [] iparam;
    delete [] ipntr;
    delete [] workd;
    delete [] workl;
    delete [] select;
    delete [] d;
		delete [] rwork;
		delete [] workev;
  }
}

#endif
